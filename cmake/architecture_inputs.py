"""基于编译产物与递归文件清单分析依赖。"""

from collections import defaultdict
from pathlib import Path
import json
import os
import re
import shlex
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "LMMC" / "cmake"))
from check_integers import load_lexer, source_tokens
from check_quality import IMPLEMENTATIONS, inventory, read_source, violation


def canonical(value, base=None):
    path = Path(value)
    if not path.is_absolute() and base is not None:
        path = Path(base) / path
    return path.resolve()


def cached_canonical(cache, value, base=None):
    key = (str(value), str(base) if base is not None else None)
    if key not in cache:
        cache[key] = canonical(value, base)
    return cache[key]


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def run_command(arguments, directory):
    result = subprocess.run(arguments, cwd=directory, capture_output=True,
                            text=True, encoding="utf-8", errors="replace", check=False)
    if result.returncode:
        raise ValueError(f"Command failed ({result.returncode}): {arguments!r}\n"
                         f"{result.stdout}\n{result.stderr}")
    return result.stdout


def production_inventory(root, report):
    # 复用品质检查的扫描与别名规则，仅收集 C++ 组件。
    return {source.path.resolve() for source in inventory(root, report, include_generated=True)
            if source.group == "production" and not source.relative.startswith("LMMC/")}


def command_arguments(entry):
    if "arguments" in entry:
        return entry["arguments"]
    arguments = shlex.split(entry["command"], posix=os.name != "nt")
    if os.name == "nt":
        arguments = [argument[1:-1] if argument.startswith('"') and argument.endswith('"')
                     else argument for argument in arguments]
    return arguments


def compile_architectures(arguments):
    architectures = []
    for index, argument in enumerate(arguments):
        if argument == "-arch":
            architectures.append(arguments[index + 1] if index + 1 < len(arguments) else None)
        elif argument.startswith("-arch="):
            architectures.append(argument.removeprefix("-arch=") or None)
    return architectures


def check_compile_architecture(manifest, entry, source, report):
    if manifest["system_name"] != "Darwin":
        return
    expected = manifest["osx_architectures"]
    actual = compile_architectures(command_arguments(entry))
    if actual != expected:
        violation(report, "compile.architecture", source,
                  "Each Darwin compiler record must select exactly the manifest architecture",
                  expected=expected, actual=actual)


def object_format_symbol(manifest, symbol):
    if manifest["system_name"] == "Darwin" and symbol.startswith("_"):
        return symbol[1:]
    return symbol


def check_direct_process_output(records, report):
    lexer, comment, string, _number = load_lexer()
    call_names = frozenset(("printf", "fprintf", "puts"))
    stream_names = frozenset(("cout", "cerr", "clog"))
    for record in records:
        tokens, _normalized = source_tokens(
            read_source(record["source"]), lexer, comment, string)
        for index, token in enumerate(tokens):
            name = token["value"]
            member_access = index > 0 and tokens[index - 1]["value"] in (".", "->")
            direct_call = (not member_access and name in call_names
                           and index + 1 < len(tokens)
                           and tokens[index + 1]["value"] == "(")
            if (not member_access and name in stream_names) or direct_call:
                violation(report, "source.direct_process_output", record["source"],
                          "Production code must not write directly to process output",
                          line=token["line"], symbol=name)


def object_output(entry):
    if entry.get("output"):
        return canonical(entry["output"], entry["directory"])
    arguments = command_arguments(entry)
    for index, argument in enumerate(arguments):
        if argument == "-o" and index + 1 < len(arguments):
            return canonical(arguments[index + 1], entry["directory"])
    raise ValueError(f"Compiler record has no object output: {entry.get('file')}")


def compile_objects(manifest, owners, database, report):
    root = canonical(manifest["root"])
    compiled = defaultdict(list)
    for entry in database:
        source = canonical(entry["file"], entry["directory"])
        if not source.is_relative_to(root / "src"):
            continue
        check_compile_architecture(manifest, entry, source, report)
        if source not in owners:
            violation(report, "ownership.unregistered_compile", source,
                      "Compiled first-party translation unit has no component owner")
            continue
        record = {"source": source, "path": object_output(entry),
                  "owner": owners[source], "directory": canonical(entry["directory"])}
        compiled[source].append(record)
    records = []
    for source, owner in owners.items():
        if source.suffix.lower() not in IMPLEMENTATIONS:
            continue
        candidates = compiled.get(source, [])
        if len(candidates) != 1:
            violation(report, "compile.membership", source,
                      "Each source must have exactly one compiler record", count=len(candidates))
            continue
        record = candidates[0]
        if f"{owner}.dir" not in record["path"].parts:
            violation(report, "compile.owner", source,
                      "Object output does not belong to its declared CMake target",
                      owner=owner, object=str(record["path"]))
        if not record["path"].is_file():
            violation(report, "compile.missing_object", source,
                      "Build the component before checking architecture", object=str(record["path"]))
            continue
        records.append(record)
    return records


def ninja_dependencies(manifest, records):
    build = canonical(manifest["build"])
    dependencies = {}
    paths = {}
    for offset in range(0, len(records), 24):
        objects = [record["path"].relative_to(build).as_posix()
                   for record in records[offset:offset + 24]]
        text = run_command([manifest["make_program"], "-t", "deps", *objects], build)
        current = None
        for line in text.splitlines():
            if not line.strip():
                current = None
                continue
            if ": #deps " in line:
                name, metadata = line.split(": #deps ", 1)
                if not metadata.endswith("(VALID)"):
                    raise ValueError(f"Stale compiler dependency record: {line}")
                current = cached_canonical(paths, name, build)
                dependencies[current] = []
            elif line.startswith("    ") and current is not None:
                dependencies[current].append(cached_canonical(paths, line[4:], build))
            else:
                raise ValueError(f"Invalid or missing Ninja dependency record: {line}")
    return dependencies


def make_dependency_words(text):
    """解析 Make 文件名转义与 Windows 盘符。"""
    text = text.replace("\\\r\n", "").replace("\\\n", "")
    separator = re.search(r":\s", text)
    if separator is None:
        raise ValueError("Compiler depfile has no target separator")
    text = text[separator.end():]
    words = []
    word = []
    index = 0
    while index < len(text):
        character = text[index]
        if character == "\\" and index + 1 < len(text) and text[index + 1] in " \t#\\:":
            index += 1
            word.append(text[index])
        elif character.isspace():
            if word:
                words.append("".join(word))
                word = []
        elif character == "$" and index + 1 < len(text) and text[index + 1] == "$":
            word.append("$")
            index += 1
        else:
            word.append(character)
        index += 1
    if word:
        words.append("".join(word))
    return words


def dependency_files(manifest, records):
    if manifest["generator"].startswith("Ninja"):
        return ninja_dependencies(manifest, records)
    dependencies = {}
    paths = {}
    for record in records:
        depfile = Path(str(record["path"]) + ".d")
        if not depfile.is_file():
            raise ValueError(f"Compiler depfile is missing: {depfile}")
        words = make_dependency_words(depfile.read_text(encoding="utf-8-sig"))
        dependencies[record["path"]] = [
            cached_canonical(paths, word, record["directory"]) for word in words]
    return dependencies


def check_dependency_freshness(record, dependencies, report, dependency_metadata):
    if record["source"] not in dependencies:
        violation(report, "dependency.source_missing", record["source"],
                  "Compiler dependencies do not contain the translation unit itself")
    object_time = record["path"].stat().st_mtime_ns
    for dependency in dependencies:
        if dependency not in dependency_metadata:
            exists = dependency.is_file()
            dependency_metadata[dependency] = (
                exists, dependency.stat().st_mtime_ns if exists else None)
        exists, dependency_time = dependency_metadata[dependency]
        if not exists:
            violation(report, "dependency.missing", dependency,
                      "A recorded compiler dependency no longer exists")
        elif dependency_time > object_time:
            violation(report, "dependency.stale", record["source"],
                      "Rebuild the object after its dependency changed", dependency=str(dependency))


def symbol_location(location, directory, paths):
    if not location or location.startswith("??"):
        return None
    path = re.sub(r":\d+(?: \(discriminator \d+\))?$", "", location)
    return cached_canonical(paths, path, directory)


def nm_records(manifest, records):
    build = canonical(manifest["build"])
    paths = {}
    for offset in range(0, len(records), 16):
        objects = [record["path"].relative_to(build).as_posix()
                   for record in records[offset:offset + 16]]
        text = run_command([manifest["nm"], "-A", "-l", "--extern-only",
                            "--format=posix", *objects], build)
        for line in text.splitlines():
            if not line.strip():
                continue
            if ": " not in line:
                raise ValueError(f"Unrecognized nm record: {line}")
            object_name, payload = line.split(": ", 1)
            fields, _, location = payload.partition("\t")
            match = re.fullmatch(r"(.*?) ([A-Za-z?])(?: [0-9a-fA-F]+){0,2}", fields.strip())
            if match is None:
                raise ValueError(f"Unrecognized nm symbol: {line}")
            yield (cached_canonical(paths, object_name, build),
                   object_format_symbol(manifest, match[1]), match[2],
                   symbol_location(location, build, paths))
