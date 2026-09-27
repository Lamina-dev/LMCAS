#!/usr/bin/env python3
"""检查组件归属、编译包含关系与外部符号依赖。"""

import argparse
from collections import defaultdict
import json
from pathlib import Path
import re
import sys

from architecture_inputs import (IMPLEMENTATIONS, canonical, check_dependency_freshness,
                                 check_direct_process_output, compile_objects, dependency_files,
                                 nm_records, production_inventory, read_json, violation)
from architecture_public_tests import check_public_tests


def load_owners(manifest, report):
    root = canonical(manifest["root"])
    layers = {entry["name"]: entry for entry in manifest["components"]}
    memberships = defaultdict(list)
    for entry in manifest["files"]:
        path = canonical(entry["path"], root)
        memberships[path].append(entry["owner"])
    owners = {}
    for path, members in memberships.items():
        if len(members) != 1 or members[0] not in layers:
            violation(report, "ownership.unique", path,
                      "Every source/header needs exactly one known component", owners=members)
        else:
            owners[path] = members[0]
        if not path.is_file():
            violation(report, "ownership.missing", path, "Owned source/header does not exist")
    expected = production_inventory(root, report)
    for path in sorted(expected - owners.keys()):
        violation(report, "ownership.unregistered", path,
                  "Recursive first-party inventory has no component owner")
    for path in sorted(owners.keys() - expected):
        if path.is_relative_to(root):
            violation(report, "ownership.out_of_scope", path,
                      "Owned file is outside the shared production inventory")
    return owners, layers


def check_edge(caller, callee, layers, path, kind, report, **evidence):
    if caller == callee:
        return
    report["edges"].add((caller, callee, kind))
    if layers[callee]["rank"] >= layers[caller]["rank"]:
        violation(report, f"{kind}.direction", path,
                  "Dependency must point to a strictly lower component",
                  caller=caller, callee=callee, **evidence)
    elif callee not in layers[caller]["dependencies"]:
        violation(report, f"{kind}.undeclared", path,
                  "Observed component dependency is not declared in CMake",
                  caller=caller, callee=callee, **evidence)


def check_includes(manifest, records, owners, layers, report, dependency_metadata):
    dependencies = dependency_files(manifest, records)
    root = canonical(manifest["root"])
    for record in records:
        included = dependencies.get(record["path"], [])
        check_dependency_freshness(record, included, report, dependency_metadata)
        for path in included:
            if (path != record["source"] and path.suffix.lower() in IMPLEMENTATIONS
                    and path.is_relative_to(root / "src")):
                violation(report, "source.textual_implementation_include", record["source"],
                          "Production translation units must not include implementation files",
                          dependency=str(path))
            owner = owners.get(path)
            if owner is not None:
                check_edge(record["owner"], owner, layers, record["source"],
                           "include", report, dependency=str(path))
            elif path.is_relative_to(root / "src") or path.is_relative_to(root / "include"):
                violation(report, "include.unowned", record["source"],
                          "Compiler dependency has no first-party owner", dependency=str(path))


def is_first_party_symbol(symbol):
    # 匹配 Itanium ABI 顶层 LMCAS 命名空间及限定符、跳板、RTTI、虚表。
    return re.match(r"^_Z(?:N[KVROr]*|T[VISWH]N|T[hv][0-9_n]*N|GVZ?N|ZN)5LMCAS", symbol) is not None


def check_symbols(manifest, records, owners, layers, report):
    objects = {record["path"]: record for record in records}
    definitions = defaultdict(set)
    undefined = []
    for path, symbol, kind, location in nm_records(manifest, records):
        if not is_first_party_symbol(symbol):
            continue
        record = objects[path]
        if kind in ("U", "w", "v"):
            undefined.append((record, symbol))
        else:
            # 按调试位置将内联与模板定义归属到源头文件。
            definitions[symbol].add(owners.get(location, record["owner"]))
    for record, symbol in undefined:
        candidates = definitions[symbol]
        if len(candidates) != 1:
            violation(report, "symbol.unique_definition", record["source"],
                      "Undefined first-party symbol needs one definition owner",
                      symbol=symbol, owners=sorted(candidates))
            continue
        check_edge(record["owner"], next(iter(candidates)), layers, record["source"],
                   "symbol", report, symbol=symbol)
    report["undefined_first_party_symbols"] = len(undefined)
    report["defined_first_party_symbols"] = len(definitions)


def check_test_linkage(manifest, report):
    for target in manifest.get("tests", []):
        required, forbidden = (("lmcas", "lmcas_test_internal") if target["linkage"] == "public"
                               else ("lmcas_test_internal", "lmcas"))
        libraries = target["libraries"]
        if required not in libraries or forbidden in libraries:
            violation(report, "test.state_ownership", target["name"],
                      "Public tests use shared state; internal tests use only the object archive",
                      linkage=target["linkage"], libraries=libraries)


def evaluate(manifest):
    report = {"schema_version": 1, "root": manifest["root"], "violations": [], "edges": set()}
    owners, layers = load_owners(manifest, report)
    database = read_json(manifest["compile_commands"])
    records = compile_objects(manifest, owners, database, report)
    dependency_metadata = {}
    check_includes(manifest, records, owners, layers, report, dependency_metadata)
    check_direct_process_output(records, report)
    check_symbols(manifest, records, owners, layers, report)
    check_test_linkage(manifest, report)
    check_public_tests(manifest, database, report, dependency_metadata)
    report["compiled_objects"] = len(records)
    report["owned_files"] = len(owners)
    report["edges"] = [dict(zip(("caller", "callee", "kind"), edge))
                       for edge in sorted(report["edges"])]
    report["passed"] = not report["violations"]
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args(argv)
    try:
        report = evaluate(read_json(args.manifest))
    except (OSError, ValueError, KeyError) as error:
        print(f"Architecture evidence unavailable: {error}", file=sys.stderr)
        return 2
    if args.report:
        args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    for item in report["violations"]:
        print(json.dumps(item, ensure_ascii=False))
    print(f"Architecture: {report['compiled_objects']} objects, {report['owned_files']} files, "
          f"{len(report['edges'])} component edges, {len(report['violations'])} violations")
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
