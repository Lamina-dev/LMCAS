#!/usr/bin/env python3
"""编译合法与非法证据图，验证架构、源码策略和产物检查。"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from architecture_inputs import (check_compile_architecture, compile_objects, dependency_files,
                                 make_dependency_words, nm_records, object_format_symbol,
                                 read_json, run_command, symbol_location)
from check_architecture import evaluate, is_first_party_symbol, load_owners


def create_sources(root, case):
    (root / "include").mkdir(parents=True)
    (root / "src/internal").mkdir(parents=True)
    header = "#pragma once\nnamespace LMCAS { int low(); int high(); }\n"
    header += "namespace LMCAS { template<class T> T fixture_identity(T x) { return x; } }\n"
    if case == "transitive_private_include":
        header += '#include "internal/upper.hpp"\n'
    (root / "include/allowed.hpp").write_text(header, encoding="utf-8")
    (root / "src/internal/upper.hpp").write_text(
        "#pragma once\nnamespace LMCAS { struct Upper { int value; }; }\n", encoding="utf-8")
    body = "return high();" if case == "upward_symbol" else "return fixture_identity(7);"
    prefix = '// std::cout and puts("inside comments are ignored")\n'
    extra = ""
    if case == "same_owner_textual_include":
        (root / "src/low_part.cpp").write_text(
            "namespace LMCAS { int low_part() { return 7; } }\n", encoding="utf-8")
        prefix += '#include "low_part.cpp"\n'
        body = "return low_part();"
        extra = " src/low_part.cpp"
    elif case == "direct_process_output":
        prefix += "#include <cstdio>\n"
        body = 'std::puts("compiled fixture output"); return fixture_identity(7);'
    (root / "src/low.cpp").write_text(
        prefix + '#include "allowed.hpp"\nnamespace LMCAS { int low() { '
        'const char* ignored = "printf("; (void)ignored; ' + body + " } }\n",
        encoding="utf-8")
    (root / "src/high.cpp").write_text(
        '#include "allowed.hpp"\nnamespace LMCAS { int high() { return fixture_identity(low()); } }\n',
        encoding="utf-8")
    if case == "unregistered_nested":
        (root / "src/nested").mkdir()
        (root / "src/nested/extra.cpp").write_text(
            "namespace LMCAS { int extra() { return 9; } }\n", encoding="utf-8")
        extra = " src/nested/extra.cpp"
    (root / "CMakeLists.txt").write_text(
        "cmake_minimum_required(VERSION 3.20)\nproject(ArchitectureFixture LANGUAGES CXX)\n"
        "set(CMAKE_EXPORT_COMPILE_COMMANDS ON)\n"
        f"add_library(low OBJECT src/low.cpp{extra})\nadd_library(high OBJECT src/high.cpp)\n"
        "foreach(target IN ITEMS low high)\n"
        "  target_compile_features(${target} PRIVATE cxx_std_17)\n"
        "  target_compile_options(${target} PRIVATE -g)\n"
        "  target_include_directories(${target} PRIVATE include src)\nendforeach()\n",
        encoding="utf-8")


def fixture_manifest(toolchain, root, build):
    files = [
        ("src/low.cpp", "low"), ("include/allowed.hpp", "low"),
        ("src/high.cpp", "high"), ("src/internal/upper.hpp", "high"),
    ]
    if (root / "src/low_part.cpp").is_file():
        files.append(("src/low_part.cpp", "low"))
    return {**{key: toolchain[key] for key in ("nm", "make_program", "generator")},
            "root": str(root), "build": str(build),
            "compile_commands": str(build / "compile_commands.json"),
            "system_name": toolchain["system_name"],
            "system_processor": toolchain["system_processor"],
            "osx_architectures": toolchain["osx_architectures"],
            "components": [{"name": "low", "rank": 0, "dependencies": []},
                           {"name": "high", "rank": 1, "dependencies": ["low"]}],
            "files": [{"path": str(root / path), "owner": owner}
                      for path, owner in files]}


class ArchitectureFixtureTests(unittest.TestCase):
    toolchain: dict

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="lmcas_architecture_")
        self.addCleanup(temporary.cleanup)
        # macOS resolves /var to /private/var in compiler evidence.
        self.directory = Path(temporary.name).resolve()

    def compile_fixture(self, case="valid"):
        root = self.directory / case
        build = self.directory / (case + "_build")
        create_sources(root, case)
        command = [self.toolchain["cmake"], "-S", str(root), "-B", str(build),
                   "-G", self.toolchain["generator"], "-DCMAKE_BUILD_TYPE=Debug",
                   "-DCMAKE_CXX_COMPILER=" + self.toolchain["compiler"],
                   "-DCMAKE_MAKE_PROGRAM=" + self.toolchain["make_program"]]
        if self.toolchain["system_name"] == "Darwin":
            command.append("-DCMAKE_OSX_ARCHITECTURES="
                           + ";".join(self.toolchain["osx_architectures"]))
        run_command(command, self.directory)
        run_command([self.toolchain["cmake"], "--build", str(build), "--parallel", "2"],
                    self.directory)
        return fixture_manifest(self.toolchain, root, build)

    def object_records(self, manifest):
        report = {"violations": []}
        owners, _ = load_owners(manifest, report)
        records = compile_objects(manifest, owners, read_json(manifest["compile_commands"]), report)
        self.assertEqual(report["violations"], [])
        return records

    def assert_darwin_architecture_rejected(self, arguments, actual):
        report = {"violations": []}
        manifest = {"system_name": "Darwin", "osx_architectures": ["arm64"]}
        check_compile_architecture(
            manifest, {"arguments": arguments}, "fixture.cpp", report)
        self.assertEqual([item["code"] for item in report["violations"]],
                         ["compile.architecture"])
        self.assertEqual(report["violations"][0]["expected"], ["arm64"])
        self.assertEqual(report["violations"][0]["actual"], actual)

    def test_darwin_compile_architecture_missing_rejected(self):
        self.assert_darwin_architecture_rejected(
            ["clang++", "-c", "fixture.cpp"], [])

    def test_darwin_compile_architecture_multiple_rejected(self):
        self.assert_darwin_architecture_rejected(
            ["clang++", "-arch", "arm64", "-arch=x86_64", "-c", "fixture.cpp"],
            ["arm64", "x86_64"])

    def test_darwin_compile_architecture_mismatch_rejected(self):
        self.assert_darwin_architecture_rejected(
            ["clang++", "-arch", "x86_64", "-c", "fixture.cpp"], ["x86_64"])

    def test_darwin_symbol_prefix_removes_exactly_one_underscore(self):
        manifest = {"system_name": "Darwin"}
        symbol = object_format_symbol(manifest, "__ZN5LMCAS3fooEv")
        self.assertEqual(symbol, "_ZN5LMCAS3fooEv")
        self.assertTrue(is_first_party_symbol(symbol))
        self.assertEqual(object_format_symbol(manifest, "___ZN5LMCAS3fooEv"),
                         "__ZN5LMCAS3fooEv")
        self.assertEqual(
            object_format_symbol({"system_name": "Linux"}, "__ZN5LMCAS3fooEv"),
            "__ZN5LMCAS3fooEv")

    def check_dependency_rule(self, case, expected):
        manifest = self.compile_fixture(case)
        build = Path(manifest["build"])
        manifest_path = build / "manifest.json"
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        report_path = build / "report.json"
        result = subprocess.run(
            [sys.executable, str(Path(__file__).with_name("check_architecture.py")),
             "--manifest", str(manifest_path), "--report", str(report_path)],
            capture_output=True, text=True, check=False)
        diagnostic = result.stdout + "\n" + result.stderr
        self.assertTrue(report_path.is_file(), diagnostic)
        report = read_json(report_path)
        self.assertEqual(result.returncode, 1 if expected else 0, diagnostic)
        self.assertEqual({item["code"] for item in report["violations"]}, expected, diagnostic)

    def test_valid_dependency_graph(self):
        self.check_dependency_rule("valid", set())

    def test_upward_symbol_rejected(self):
        self.check_dependency_rule("upward_symbol", {"symbol.direction"})

    def test_transitive_private_include_rejected(self):
        self.check_dependency_rule("transitive_private_include", {"include.direction"})

    def test_unregistered_nested_source_rejected(self):
        self.check_dependency_rule(
            "unregistered_nested", {"ownership.unregistered", "ownership.unregistered_compile"})

    def test_same_owner_textual_implementation_include_rejected(self):
        self.check_dependency_rule(
            "same_owner_textual_include", {"source.textual_implementation_include"})

    def test_direct_process_output_rejected(self):
        self.check_dependency_rule("direct_process_output", {"source.direct_process_output"})

    def test_shared_header_changes_are_observed(self):
        manifest = self.compile_fixture()
        records = self.object_records(manifest)
        root = Path(manifest["root"])
        header = root / "include/allowed.hpp"
        original = header.read_text(encoding="utf-8")
        original_time = header.stat().st_mtime_ns
        (root / "src/high.cpp").touch()
        run_command([self.toolchain["cmake"], "--build", manifest["build"], "--target", "high"], root)
        times = {record["owner"]: record["path"].stat().st_mtime_ns for record in records}
        self.assertLess(times["low"], times["high"])
        middle = times["low"] + (times["high"] - times["low"]) // 2
        os.utime(header, ns=(middle, middle))
        stale = evaluate(manifest)["violations"]
        self.assertEqual([(item["code"], item["path"]) for item in stale],
                         [("dependency.stale", str(root / "src/low.cpp"))])
        os.utime(header, ns=(original_time, original_time))
        self.assertTrue(evaluate(manifest)["passed"])
        header.unlink()
        missing = evaluate(manifest)["violations"]
        self.assertEqual(sum(item["code"] == "dependency.missing" and item["path"] == str(header)
                             for item in missing), 2, missing)
        header.write_text(original, encoding="utf-8")
        os.utime(header, ns=(original_time, original_time))
        self.assertTrue(evaluate(manifest)["passed"])

    def test_make_dependency_escaping(self):
        manifest = self.compile_fixture()
        records = self.object_records(manifest)
        header = Path(manifest["build"]) / "header with spaces.hpp"
        header.write_text("#pragma once\n", encoding="utf-8")
        for record in records:
            escaped = [str(path).replace("\\", "\\\\").replace(" ", "\\ ")
                       for path in (record["source"], header)]
            Path(str(record["path"]) + ".d").write_text(
                "target: " + " ".join(escaped) + "\n", encoding="utf-8")
        dependencies = dependency_files({**manifest, "generator": "Unix Makefiles"}, records)
        for owner, paths in dependencies.items():
            with self.subTest(owner=owner):
                self.assertIn(header, paths)
        self.assertEqual(
            make_dependency_words(r"C:\obj.o: C:\src\file.cpp path\ with\ spaces.hpp C\:/escaped.hpp"),
            [r"C:\src\file.cpp", "path with spaces.hpp", "C:/escaped.hpp"])

    def test_malformed_compiler_evidence_rejected(self):
        manifest = self.compile_fixture()
        build = Path(manifest["build"])
        source = build / "malformed_evidence.cpp"
        executable = build / "malformed_evidence.exe"
        source.write_text(
            '#include <cstdio>\nint main() { std::puts("malformed compiler evidence"); }\n',
            encoding="utf-8")
        run_command([self.toolchain["compiler"], str(source), "-o", str(executable)], build)
        checker = str(Path(__file__).with_name("check_architecture.py"))
        for key in ("make_program", "nm"):
            with self.subTest(tool=key):
                path = build / f"malformed_{key}.json"
                bad_manifest = {**manifest, key: str(executable)}
                if key == "make_program":
                    bad_manifest["generator"] = "Ninja"
                path.write_text(json.dumps(bad_manifest), encoding="utf-8")
                result = subprocess.run([sys.executable, "-B", checker, "--manifest", str(path)],
                                        capture_output=True, text=True, check=False)
                self.assertEqual(result.returncode, 2, result.stdout + "\n" + result.stderr)
                self.assertIn("malformed compiler evidence", result.stderr)

    def test_debug_locations_preserve_header_owner(self):
        manifest = self.compile_fixture()
        records = self.object_records(manifest)
        header = Path(manifest["root"]) / "include/allowed.hpp"
        symbols = [item for item in nm_records(manifest, records) if "fixture_identity" in item[1]]
        self.assertTrue(symbols, "Compiled fixture_identity definitions must be present")
        for symbol in symbols:
            with self.subTest(symbol=symbol[1]):
                self.assertEqual(symbol[3], header)
        self.assertEqual(
            symbol_location(f"{header}:2 (discriminator 1)", manifest["build"], {}), header)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    ArchitectureFixtureTests.toolchain = read_json(args.manifest)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ArchitectureFixtureTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
