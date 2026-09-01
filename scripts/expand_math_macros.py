"""Expand the MyST math macros of ``myst.yml`` into plain LaTeX, in place.

Macros are convenient to write but only MyST knows about them: the VS Code
preview, GitHub's markdown view and anything else that renders ``$...$`` with a
plain KaTeX/MathJax instance shows them unresolved.  This script lets you keep
writing ``\\mcT_h`` and later replace it by ``\\mathcal{T}_h`` in the source.

The macro definitions stay in ``myst.yml`` -- the script only reads them -- so
files that have already been expanded keep rendering, and new material can go on
using macros.

Usage::

    python scripts/expand_math_macros.py                 # dry run over the project
    python scripts/expand_math_macros.py --write         # rewrite the files
    python scripts/expand_math_macros.py chapter_02 --write
    python scripts/expand_math_macros.py chapter_01/lax_milgram.md --diff

What is left untouched:

* the YAML frontmatter,
* code fences (```` ```python ````, ```` ```{code-cell} ````, ...) and inline
  code spans, so that a ``\\PP`` inside a code cell survives,
* the macro definitions in ``myst.yml``.

Everything else is expanded, including the body of directives such as
``{prf:theorem}`` and ``{math}``.  Note that macros are expanded wherever they
appear in that text, not only between math delimiters: a macro outside math does
not render anyway, so there is nothing to lose and no need to guess where the
math regions are.

If a markdown file is paired with a notebook (jupytext), run ``jupytext --sync``
on it afterwards so the ``.ipynb`` picks up the expansion.
"""

from __future__ import annotations

import argparse
import difflib
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MYST_YML = ROOT / "myst.yml"

#: Fence info strings whose content is code rather than markdown.
CODE_DIRECTIVES = {
    "code",
    "code-block",
    "code-cell",
    "literalinclude",
    "raw",
}

MAX_EXPANSION_ROUNDS = 20


class MacroError(RuntimeError):
    pass


@dataclass(frozen=True)
class Macro:
    name: str  # without the leading backslash
    body: str
    n_args: int

    @classmethod
    def parse(cls, key: str, body: str) -> "Macro":
        name = key[1:] if key.startswith("\\") else key
        if not name.isalpha():
            raise MacroError(f"unsupported macro name {key!r}: expected letters only")
        used = [int(m) for m in re.findall(r"(?<!\\)#(\d)", body)]
        return cls(name=name, body=body, n_args=max(used, default=0))


def load_macros(myst_yml: Path, page_math: dict[str, str] | None = None) -> dict[str, Macro]:
    """Macros from ``myst.yml``, overridden by a page's own frontmatter."""
    cfg = yaml.safe_load(myst_yml.read_text()) or {}
    raw: dict[str, str] = dict(cfg.get("project", {}).get("math") or {})
    raw.update(page_math or {})
    return {m.name: m for m in (Macro.parse(k, v) for k, v in raw.items())}


# --------------------------------------------------------------------------
# splitting the document into expandable prose and protected code
# --------------------------------------------------------------------------

FENCE_RE = re.compile(r"^(?P<indent>\s*)(?P<marker>`{3,}|~{3,})(?P<info>.*)$")
FRONTMATTER_RE = re.compile(r"\A---\r?\n.*?^---[ \t]*\r?$\r?\n?", re.DOTALL | re.MULTILINE)


def _is_code_fence(info: str) -> bool:
    """True for fences holding code, False for MyST directives and math."""
    info = info.strip()
    if not info:
        return True
    if info.startswith("{"):
        directive = info[1:].split("}", 1)[0].strip()
        return directive in CODE_DIRECTIVES
    # ```python, ```text, ```{eval-rst} ... a bare language name means code.
    return True


CODE_SPAN_RE = re.compile(r"(?<!`)(`+)(?!`)(.+?)(?<!`)\1(?!`)")


def protected_spans(text: str) -> list[tuple[int, int]]:
    """Character ranges that must not be touched: frontmatter, code, code spans.

    Fences are tracked with a stack, so a directive fence (``{prf:theorem}``,
    ``{math}``, ...) may contain code fences and only the latter are protected.
    """
    spans: list[tuple[int, int]] = []

    frontmatter = FRONTMATTER_RE.match(text)
    if frontmatter:
        spans.append((0, frontmatter.end()))

    stack: list[tuple[str, int, bool, int]] = []  # char, length, is_code, content start
    offset = 0
    for line in text.splitlines(keepends=True):
        stripped = line.rstrip("\r\n")
        in_code = bool(stack) and stack[-1][2]
        match = FENCE_RE.match(stripped)
        is_fence_line = False
        if match:
            marker = match.group("marker")
            char, length = marker[0], len(marker)
            info = match.group("info").strip()
            top = stack[-1] if stack else None
            if top is not None and char == top[0] and length >= top[1] and not info:
                stack.pop()
                if top[2]:
                    spans.append((top[3], offset))
                is_fence_line = True
            elif not in_code:  # a fence inside a code block is just code
                stack.append((char, length, _is_code_fence(info), offset + len(line)))
                is_fence_line = True

        if not is_fence_line and not in_code:
            for span in CODE_SPAN_RE.finditer(stripped):
                # ``{math}`\mcT`` is a math role -- its content is math, not code.
                if stripped[: span.start()].rstrip().endswith("{math}"):
                    continue
                spans.append((offset + span.start(), offset + span.end()))
        offset += len(line)

    for _, _, is_code, start in stack:  # unterminated fences: protect to the end
        if is_code:
            spans.append((start, len(text)))

    return _merge(spans)


def _merge(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(spans):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


# --------------------------------------------------------------------------
# macro expansion
# --------------------------------------------------------------------------


def _read_argument(text: str, pos: int) -> tuple[str, int]:
    """Read one TeX argument starting at ``pos``; return it and the position after."""
    while pos < len(text) and text[pos] in " \t\n":
        pos += 1
    if pos >= len(text):
        raise MacroError("missing argument at end of input")
    if text[pos] == "{":
        depth, start = 0, pos
        while pos < len(text):
            char = text[pos]
            if char == "\\":
                pos += 2
                continue
            if char == "{":
                depth += 1
            elif char == "}":
                depth -= 1
                if depth == 0:
                    return text[start + 1 : pos], pos + 1
            pos += 1
        raise MacroError("unbalanced braces in macro argument")
    if text[pos] == "\\":  # a control sequence is a single argument token
        match = re.compile(r"\\[A-Za-z]+|\\.").match(text, pos)
        assert match is not None
        return match.group(), match.end()
    return text[pos], pos + 1


def _substitute(macro: Macro, args: list[str]) -> str:
    def replace(match: re.Match[str]) -> str:
        return args[int(match.group(1)) - 1]

    return re.sub(r"(?<!\\)#(\d)", replace, macro.body)


def expand(text: str, macros: dict[str, Macro], counts: Counter[str] | None = None) -> str:
    """Expand macros in ``text`` repeatedly, until nothing changes."""
    if not macros:
        return text
    # Longest name first so that \PPdc is not read as \PP followed by 'dc'.
    names = sorted(macros, key=len, reverse=True)
    pattern = re.compile(r"\\(" + "|".join(re.escape(n) for n in names) + r")(?![A-Za-z])")

    for _ in range(MAX_EXPANSION_ROUNDS):
        out, pos, changed = [], 0, False
        for match in pattern.finditer(text):
            if match.start() < pos:  # already consumed as an argument
                continue
            macro = macros[match.group(1)]
            cursor = match.end()
            try:
                args = []
                for _ in range(macro.n_args):
                    arg, cursor = _read_argument(text, cursor)
                    args.append(arg)
            except MacroError as exc:
                raise MacroError(f"\\{macro.name}: {exc}") from exc
            out.append(text[pos : match.start()])
            out.append(_substitute(macro, args))
            pos = cursor
            changed = True
            if counts is not None:
                counts[macro.name] += 1
        out.append(text[pos:])
        text = "".join(out)
        if not changed:
            return text
    raise MacroError(f"macro expansion did not terminate after {MAX_EXPANSION_ROUNDS} rounds")


def expand_document(text: str, macros: dict[str, Macro]) -> tuple[str, Counter[str]]:
    """Expand macros everywhere except in protected (code) regions."""
    counts: Counter[str] = Counter()
    pieces, pos = [], 0
    for start, end in protected_spans(text):
        pieces.append(expand(text[pos:start], macros, counts))
        pieces.append(text[start:end])
        pos = end
    pieces.append(expand(text[pos:], macros, counts))
    return "".join(pieces), counts


def page_macros(text: str) -> dict[str, str]:
    """The ``math:`` mapping of a page's own frontmatter, if any."""
    match = FRONTMATTER_RE.match(text)
    if not match:
        return {}
    body = match.group()[4:].rsplit("---", 1)[0]
    return (yaml.safe_load(body) or {}).get("math") or {}


# --------------------------------------------------------------------------
# command line
# --------------------------------------------------------------------------


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
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        default=[ROOT],
        help="markdown files or directories to process (default: the whole project)",
    )
    parser.add_argument("--write", action="store_true", help="rewrite the files in place")
    parser.add_argument("--diff", action="store_true", help="show a unified diff per file")
    parser.add_argument(
        "--myst-yml",
        type=Path,
        default=DEFAULT_MYST_YML,
        help=f"config holding the macro definitions (default: {DEFAULT_MYST_YML})",
    )
    args = parser.parse_args(argv)

    project_macros = load_macros(args.myst_yml)
    total: Counter[str] = Counter()
    changed_files = 0

    for path in collect_files(args.paths):
        text = path.read_text()
        macros = project_macros | {
            m.name: m for m in (Macro.parse(k, v) for k, v in page_macros(text).items())
        }
        try:
            new_text, counts = expand_document(text, macros)
        except MacroError as exc:
            print(f"{path}: {exc}", file=sys.stderr)
            return 1
        if new_text == text:
            continue

        changed_files += 1
        total.update(counts)
        summary = ", ".join(f"\\{name}x{n}" for name, n in counts.most_common())
        print(f"{path.relative_to(Path.cwd()) if path.is_absolute() else path}: {summary}")
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

    verb = "expanded" if args.write else "would expand"
    print(f"\n{verb} {sum(total.values())} macro uses in {changed_files} file(s)")
    if changed_files and not args.write:
        print("re-run with --write to apply")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
