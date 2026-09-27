"""检查公开接口测试仅使用公开头文件。"""

from collections import defaultdict

from architecture_inputs import (IMPLEMENTATIONS, canonical, check_compile_architecture,
                                 check_dependency_freshness, command_arguments, dependency_files,
                                 object_output, violation)


def public_test_objects(manifest, database, report):
    compiled = defaultdict(list)
    for entry in database:
        compiled[canonical(entry["file"], entry["directory"])].append(entry)
    records = []
    for target in manifest.get("tests", []):
        if target["linkage"] != "public":
            continue
        for value in target["sources"]:
            source = canonical(value)
            if source.suffix.lower() not in IMPLEMENTATIONS:
                continue
            candidates = compiled[source]
            if len(candidates) != 1:
                violation(report, "test.compile_membership", source,
                          "Public test needs exactly one compiler record", count=len(candidates))
                continue
            entry = candidates[0]
            check_compile_architecture(manifest, entry, source, report)
            path = object_output(entry)
            if not path.is_file():
                violation(report, "test.missing_object", source, "Build public tests before architecture validation")
                continue
            if f"{target['name']}.dir" not in path.parts:
                violation(report, "test.compile_owner", source,
                          "Public test object must belong to its registered executable")
            if any("LMCAS_STATIC_DEFINE" in argument for argument in command_arguments(entry)):
                violation(report, "test.static_definition", source,
                          "Public test must consume the shared-library declarations")
            records.append({"source": source, "path": path,
                            "directory": canonical(entry["directory"])})
    return records


def check_public_tests(manifest, database, report, dependency_metadata):
    records = public_test_objects(manifest, database, report)
    dependencies = dependency_files(manifest, records)
    private = canonical(manifest["root"]) / "src"
    for record in records:
        included = dependencies.get(record["path"], [])
        check_dependency_freshness(record, included, report, dependency_metadata)
        for path in included:
            if path.is_relative_to(private):
                violation(report, "test.private_dependency", record["source"],
                          "Public consumer includes a private implementation header", dependency=str(path))
    report["public_test_objects"] = len(records)
