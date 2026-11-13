#!/usr/bin/env python3
"""Stage configured code without executing analysis or overwriting changed files."""
import argparse
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, help="Absolute analysis working directory")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    root = Path(args.root).expanduser()
    if not root.is_absolute():
        parser.error("--root must be an absolute path")
    root = root.resolve()
    if root == repo or repo in root.parents or root in repo.parents:
        parser.error("Choose a separate working directory outside the source repository")
    # Paths are embedded in shell/R/Python source; reject characters requiring extra quoting.
    root_text = root.as_posix()
    if any(c in root_text for c in "\n\r\t'\"`$\\") or " " in root_text:
        parser.error("Use a working path without spaces, quotes or shell metacharacters")
    layout = json.loads((repo / "config/runtime_layout.json").read_text(encoding="utf-8"))
    plan = []
    for item in layout:
        source = (repo / item["file"]).resolve()
        if repo not in source.parents or not source.is_file():
            raise ValueError("Invalid source: " + item["file"])
        content = source.read_text(encoding="utf-8").replace("/path/to/chicken-aging", root_text)
        for relative in [item["runtime"]] + item.get("additional_runtime", []):
            target = (root / relative).resolve()
            if root not in target.parents:
                raise ValueError("Runtime target escapes working directory: " + relative)
            if target.exists() and target.read_text(encoding="utf-8") != content:
                raise FileExistsError("Preserving changed existing file: " + str(target))
            plan.append((target, content))
    targets = [str(p) for p, _ in plan]
    if len(targets) != len(set(targets)):
        raise ValueError("Duplicate runtime targets")
    print(f"{len(plan)} code/config files -> {root_text}")
    if args.dry_run:
        print("Dry run: no directories or files created; no analysis executed.")
        return
    for target, content in plan:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(content, encoding="utf-8", newline="\n")
    print(
        "Workspace code prepared. Supply inputs and activate the relevant environment before running a stage."
    )


if __name__ == "__main__":
    main()
