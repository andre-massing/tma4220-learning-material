"""Rewrite ```{math} directives as ``$$ ... $$`` blocks.

The ``{math}`` directive is not rendered by the VS Code MyST preview (its
directive plugin only knows admonitions, code, raw, image and figure), while
``$$ ... $$`` is.  Both are equivalent MyST for display math, so the source can
be converted once and for all::

    ```{math}
    :label: eq:poincare
    \\|v\\|_{L^2} \\leqslant C \\|\\nabla v\\|_{L^2}
    ```

becomes::

    $$
    \\|v\\|_{L^2} \\leqslant C \\|\\nabla v\\|_{L^2}
    $$ (eq:poincare)

Usage::

    python scripts/math_directive_to_dollars.py                 # dry run
    python scripts/math_directive_to_dollars.py --write
    python scripts/math_directive_to_dollars.py chapter_01 --diff

Blocks are always surrounded by blank lines, which a labelled ``$$`` requires:
without a preceding blank line MyST drops the label and leaks the ``(label)``
text into the page (MyST-Parser issue #897).  Indentation is preserved, so a
directive nested in a list item stays in that list item.

Directives carrying options other than ``:label:`` are left alone and reported,
since ``$$`` has nowhere to put them.
"""

from __future__ import annotations

import argparse
import difflib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

FENCE_RE = re.compile(r"^(?P<indent>[ \t]*)(?P<marker>`{3,}|~{3,})(?P<info>.*)$")
OPTION_RE = re.compile(r"^[ \t]*:(?P<name>[A-Za-z][A-Za-z0-9_-]*):[ \t]*(?P<value>.*?)[ \t]*$")

#: Fence info strings whose content is code, so nothing inside them is markdown.
CODE_DIRECTIVES = {"code", "code-block", "code-cell", "literalinclude", "raw"}


def _is_code_fence(info: str) -> bool:
    info = info.strip()
    if not info:
        return True
    if info.startswith("{"):
        return info[1:].split("}", 1)[0].strip() in CODE_DIRECTIVES
    return True


def _closes(line: str, char: str, length: int) -> bool:
    match = FENCE_RE.match(line.rstrip("\r\n"))
    if not match:
        return False
    marker = match.group("marker")
    return marker[0] == char and len(marker) >= length and not match.group("info").strip()


def convert(text: str, on_skip=None) -> tuple[str, int]:
    """Return the text with ``{math}`` directives rewritten, and how many."""
    lines = text.splitlines(keepends=True)
    newline = "\r\n" if text.endswith("\r\n") or "\r\n" in text[:200] else "\n"
    out: list[str] = []
    stack: list[tuple[str, int, bool]] = []  # marker char, length, is_code
    converted = 0
    i = 0

    while i < len(lines):
        line = lines[i]
        stripped = line.rstrip("\r\n")
        match = FENCE_RE.match(stripped)
        in_code = bool(stack) and stack[-1][2]

        if match:
            marker = match.group("marker")
            char, length = marker[0], len(marker)
            info = match.group("info").strip()
            top = stack[-1] if stack else None

            # The closing check must come first: a code fence has to be closed
            # even though we ignore everything else while inside it.
            if top is not None and char == top[0] and length >= top[1] and not info:
                stack.pop()
                out.append(line)
                i += 1
                continue

            if in_code:
                out.append(line)
                i += 1
                continue

            if info == "{math}":
                end = next(
                    (j for j in range(i + 1, len(lines)) if _closes(lines[j], char, length)),
                    None,
                )
                if end is not None:
                    block = _convert_block(
                        indent=match.group("indent"),
                        body=lines[i + 1 : end],
                        newline=newline,
                        on_skip=on_skip,
                        line_no=i + 1,
                    )
                    if block is not None:
                        if out and out[-1].strip():  # blank line before
                            out.append(newline)
                        out.extend(block)
                        if end + 1 < len(lines) and lines[end + 1].strip():
                            out.append(newline)  # blank line after
                        converted += 1
                        i = end + 1
                        continue

            stack.append((char, length, _is_code_fence(info)))

        out.append(line)
        i += 1

    return "".join(out), converted


def _convert_block(indent: str, body: list[str], newline: str, on_skip, line_no: int):
    """Build the ``$$`` replacement for one directive body, or None to skip it."""
    label = None
    while body:
        option = OPTION_RE.match(body[0].rstrip("\r\n"))
        if not option:
            break
        name, value = option.group("name"), option.group("value")
        if name not in {"label", "name"}:
            if on_skip:
                on_skip(line_no, f"unsupported option :{name}:")
            return None
        label = value or label
        body = body[1:]

    while body and not body[0].strip():  # drop padding inside the directive
        body = body[1:]
    while body and not body[-1].strip():
        body = body[:-1]
    if not body:
        if on_skip:
            on_skip(line_no, "empty directive")
        return None

    closing = f"{indent}$${f' ({label})' if label else ''}{newline}"
    return [f"{indent}$${newline}", *body, closing]


def collect_files(paths: list[Path]) -> list[Path]:
    files: list[Path] = []
    for path in paths:
        if path.is_dir():
            files.extend(sorted(p for p in path.rglob("*.md") if "_build" not in p.parts))
        elif path.suffix == ".md":
            files.append(path)
        else:
            raise SystemExit(f"not a markdown file or directory: {path}")
    return files


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("paths", nargs="*", type=Path, default=[ROOT])
    parser.add_argument("--write", action="store_true", help="rewrite the files in place")
    parser.add_argument("--diff", action="store_true", help="show a unified diff per file")
    args = parser.parse_args(argv)

    total = 0
    changed_files = 0
    for path in collect_files(args.paths):
        text = path.read_text()
        skipped: list[str] = []
        new_text, count = convert(
            text, on_skip=lambda line, why: skipped.append(f"  line {line}: {why}")
        )
        for note in skipped:
            print(f"{path}: skipped{note}", file=sys.stderr)
        if not count:
            continue
        changed_files += 1
        total += count
        print(f"{path}: {count} directive(s)")
        if args.diff:
            print(
                "".join(
                    difflib.unified_diff(
                        text.splitlines(keepends=True),
                        new_text.splitlines(keepends=True),
                        fromfile=str(path),
                        tofile=str(path),
                    )
                )
            )
        if args.write:
            path.write_text(new_text)

    verb = "converted" if args.write else "would convert"
    print(f"\n{verb} {total} directive(s) in {changed_files} file(s)")
    if total and not args.write:
        print("re-run with --write to apply")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
