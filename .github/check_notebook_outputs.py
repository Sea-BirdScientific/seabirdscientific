"""Fails if any notebook has cell outputs or execution counts.

The example notebooks are committed with outputs cleared (see the PR template). Uses only
the standard library so it can run without installing the project.

Usage: python .github/check_notebook_outputs.py notebooks/*.ipynb
"""

import json
import sys
from pathlib import Path


def find_uncleared_cells(path: Path) -> list[int]:
    """Lists the code cells of a notebook that have outputs or an
    execution count

    :param path: path to the .ipynb file
    :return: 1-based indexes of the uncleared cells
    """
    notebook = json.loads(path.read_text(encoding="utf-8"))
    return [
        number
        for number, cell in enumerate(notebook.get("cells", []), start=1)
        if cell.get("cell_type") == "code"
        and (cell.get("outputs") or cell.get("execution_count") is not None)
    ]


def main(paths: list[str]) -> int:
    if not paths:
        print("usage: check_notebook_outputs.py NOTEBOOK [NOTEBOOK ...]", file=sys.stderr)
        return 2

    failed = False
    for path in map(Path, paths):
        cells = find_uncleared_cells(path)
        if cells:
            failed = True
            # GitHub Actions annotation, shown on the file in the PR
            print(
                f"::error file={path}::{path} has outputs or execution counts in cells "
                f"{', '.join(map(str, cells))}. Clear all outputs before committing."
            )

    if not failed:
        print(f"{len(paths)} notebooks have cleared outputs")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
