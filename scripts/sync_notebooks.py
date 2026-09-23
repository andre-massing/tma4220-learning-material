"""Sync jupytext-paired notebooks and fold ``hide-input`` cells in JupyterLab.

Code cells are folded in two independent places, each with its own marker:

* the MyST build (website, PDF) folds cells tagged ``hide-input``, and that tag
  lives in the ``.md`` file::

      ```{code-cell} ipython3
      :tags: [hide-input]
      ```

* JupyterLab folds cells whose metadata contains
  ``"jupyter": {"source_hidden": true}``, which lives in the ``.ipynb`` only.
  Writing it into the ``.md`` would make mystmd warn "unexpected option" for
  every such cell.

This script treats the ``hide-input`` tags in the ``.md`` as the single place to
edit.  For every Markdown file paired with a notebook (``formats`` containing
``ipynb`` in its jupytext frontmatter) it

1. runs ``jupytext --sync``,
2. sets ``source_hidden: true`` on every code cell tagged ``hide-input`` in the
   paired ``.ipynb``, and
3. touches the ``.md`` if the notebook changed, so that the ``.md`` is not
   older than the ``.ipynb``.  JupyterLab objects to a pair whose notebook is
   newer than its text file.

Usage::

    python scripts/sync_notebooks.py                   # all paired files
    python scripts/sync_notebooks.py chapter_02        # one directory
    python scripts/sync_notebooks.py --dry-run         # report only
    python scripts/sync_notebooks.py --prune           # also unfold untagged cells

For step 2 to be stable, the frontmatter must exclude the ``jupyter`` key and
keep ``tags``::

    jupytext:
      cell_metadata_filter: -jupyter

Files with tagged cells whose filter would drop the tags (e.g. ``-all``) or
write ``source_hidden`` back into the ``.md`` are skipped entirely, since
syncing them would already destroy the tags, and a warning explains what to
change.  The script never edits the filter itself, and exits with status 1 if
any file was skipped.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import jupytext
from jupytext.cli import jupytext as jupytext_cli
from jupytext.metadata_filter import metadata_filter_as_dict

ROOT = Path(__file__).resolve().parents[1]

#: The cell tag MyST uses to fold code, and the one this script mirrors.
FOLD_TAG = "hide-input"

#: Directories that never contain source files worth syncing.
SKIP_DIRS = {"_build", "node_modules", "templates", "scripts"}


def _skipped(path: Path) -> bool:
    parts = path.relative_to(ROOT).parts[:-1]
    return any(p.startswith(".") or p in SKIP_DIRS for p in parts)


def find_paired_markdown(targets: list[Path]) -> list[Path]:
    """Return all Markdown files under ``targets`` that are paired with a notebook."""
    candidates: list[Path] = []
    for target in targets:
        if target.is_file():
            candidates.append(target)
        else:
            candidates.extend(sorted(target.rglob("*.md")))

    paired = []
    for md in candidates:
        if md.suffix != ".md" or _skipped(md):
            continue
        formats = jupytext.read(md).metadata.get("jupytext", {}).get("formats", "")
        if "ipynb" in formats.split(","):
            paired.append(md)
    return paired


def _excludes(metadata_filter: dict, key: str) -> bool:
    """Whether ``metadata_filter`` keeps ``key`` out of the text file."""
    excluded = metadata_filter.get("excluded", [])
    additional = metadata_filter.get("additional", [])
    if excluded == "all":
        return not (additional == "all" or key in additional)
    return key in excluded


def filter_problems(md: Path) -> list[str]:
    """Explain why folding would be unsafe for ``md``, or return an empty list."""
    notebook = jupytext.read(md)
    raw = notebook.metadata.get("jupytext", {}).get("cell_metadata_filter")
    metadata_filter = metadata_filter_as_dict(raw)
    uses_tag = any(FOLD_TAG in cell.metadata.get("tags", []) for cell in notebook.cells)

    problems = []
    if uses_tag and _excludes(metadata_filter, "tags"):
        problems.append(
            f"cell_metadata_filter '{raw}' drops cell tags, so syncing from the "
            f".ipynb would remove '{FOLD_TAG}' from the .md"
        )
    if uses_tag and not _excludes(metadata_filter, "jupyter"):
        problems.append(
            f"cell_metadata_filter '{raw or '(none)'}' does not exclude 'jupyter', so "
            "source_hidden would be written into the .md and mystmd would warn"
        )
    if problems:
        problems.append("set 'cell_metadata_filter: -jupyter' under 'jupytext:' in the frontmatter")
    return problems


def fold_cells(ipynb: Path, prune: bool, dry_run: bool) -> tuple[int, int]:
    """Mirror the ``hide-input`` tags as ``source_hidden`` in ``ipynb``.

    Returns the number of cells folded and unfolded.
    """
    notebook = jupytext.read(ipynb)
    folded = unfolded = 0
    for cell in notebook.cells:
        if cell.cell_type != "code":
            continue
        tagged = FOLD_TAG in cell.metadata.get("tags", [])
        hidden = bool(cell.metadata.get("jupyter", {}).get("source_hidden"))
        if tagged and not hidden:
            cell.metadata.setdefault("jupyter", {})["source_hidden"] = True
            folded += 1
        elif prune and hidden and not tagged:
            cell.metadata["jupyter"].pop("source_hidden")
            if not cell.metadata["jupyter"]:
                cell.metadata.pop("jupyter")
            unfolded += 1
    if (folded or unfolded) and not dry_run:
        jupytext.write(notebook, ipynb)
    return folded, unfolded


def sync(md: Path) -> None:
    exit_code = jupytext_cli(["--sync", "--quiet", str(md)])
    if exit_code:
        raise RuntimeError(f"jupytext --sync failed for {md} (exit code {exit_code})")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        help="Markdown files or directories (default: the whole book)",
    )
    parser.add_argument("--dry-run", action="store_true", help="report without writing anything")
    parser.add_argument(
        "--prune",
        action="store_true",
        help=f"also unfold cells that are folded in JupyterLab but not tagged '{FOLD_TAG}'",
    )
    args = parser.parse_args(argv)

    targets = [p.resolve() for p in args.paths] or [ROOT]
    pairs = find_paired_markdown(targets)
    if not pairs:
        print("No Markdown files paired with a notebook found.")
        return 0

    warnings = 0
    for md in pairs:
        rel = md.relative_to(ROOT)
        ipynb = md.with_suffix(".ipynb")
        problems = filter_problems(md)

        if problems:
            # Do not sync: with such a filter the sync itself rewrites the .md
            # and destroys the very tags this script relies on.
            warnings += 1
            print(f"  {rel}: SKIPPED, neither synced nor folded")
            for problem in problems:
                print(f"      - {problem}")
            continue

        if not args.dry_run:
            sync(md)

        if not ipynb.exists():
            print(f"  {rel}: would create {ipynb.name} on sync" if args.dry_run else f"  {rel}: synced")
            continue

        folded, unfolded = fold_cells(ipynb, prune=args.prune, dry_run=args.dry_run)
        if (folded or unfolded) and not args.dry_run:
            # Keep the .md at least as recent as the patched .ipynb; JupyterLab
            # objects to a pair whose notebook is newer than its text file.
            md.touch()

        verb = "would fold" if args.dry_run else "folded"
        summary = f"{verb} {folded} cell(s)"
        if args.prune:
            summary += f", {'would unfold' if args.dry_run else 'unfolded'} {unfolded}"
        print(f"  {rel}: {'checked' if args.dry_run else 'synced'}, {summary}")

    print(f"\n{len(pairs)} paired file(s), {warnings} with warnings" + (" (dry run)" if args.dry_run else ""))
    return 1 if warnings else 0


if __name__ == "__main__":
    sys.exit(main())
