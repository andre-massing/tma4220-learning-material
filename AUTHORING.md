# Authoring conventions

Notes for writing the lecture material in this directory. Not part of the book —
it is deliberately absent from the `toc` in `myst.yml`.

## Files and structure

* Chapters live in `chapter_NN/`, one MyST Markdown file per topic.
* A new file has to be registered **twice** in `myst.yml`: under `project.toc`
  (for the website) and under `project.exports[].articles` (for the PDF/LaTeX
  export, which share a YAML anchor `&book_articles`).
* Figures belong in `chapter_NN/figures/`. Hand-drawn sketches are recreated as
  hand-authored SVG rather than scanned.

## Math macros

Macros are defined once, in `project.math` in `myst.yml` (ported from the
`notation.sty` of TMA4320): `\RR`, `\mcT`, `\bfn`, `\diam`, `\jump{...}`, ...
Add new ones there rather than per file.

MyST also accepts a `math:` block in a single file's YAML frontmatter, which is
merged over the project-level ones. Use it only for notation local to one
chapter — and note that jupytext strips unknown frontmatter keys on save unless
`math` is listed in `notebook_metadata_filter` in `pyproject.toml`.

**Macros do not render outside MyST.** The VS Code markdown preview, GitHub and
any plain KaTeX/MathJax renderer show them unresolved (see below). If a file
should be readable outside the MyST build, expand the macros into plain LaTeX:

```bash
python scripts/expand_math_macros.py                 # dry run over the project
python scripts/expand_math_macros.py chapter_02 --write
```

The definitions stay in `myst.yml`, so expanded and unexpanded files coexist and
new material can keep using macros. The rewrite is one-way; run it per chapter
once the text has settled, and `jupytext --sync` afterwards on paired files.

## Equations, labels and cross-references

Display math is written as `$$ ... $$`, with the label in parentheses after the
closing delimiter:

```
$$
a(u, v) = l(v) \quad \forall v \in V
$$ (eq:lax-milgram-problem)
```

The `{math}` directive is equivalent MyST but is not rendered by the VS Code
preview; the 136 directives that existed were converted with
`scripts/math_directive_to_dollars.py`, so new material should use `$$` too.

Labels are prefixed by kind — `eq:`, `def:`, `thm:`, `lem:`, `cor:`, `exer:` —
and referenced with `` {eq}`eq:...` `` for equations and
`` {prf:ref}`text<thm-...>` `` for theorem-like blocks. Both `prefix:kebab` and
`prefix-kebab` spellings occur historically; `prefix:kebab` is the majority and
is what new labels should use.

> **A labelled `$$` block needs a blank line before the opening `$$`.**
>
> ```
> Some text.
>
> $$
> a = b
> $$ (eq:my-label)
> ```
>
> Without the blank line the equation still renders, but the label is silently
> dropped and ` (eq:my-label)` leaks into the page as literal text, so every
> reference to it breaks. This is an upstream parser limitation
> ([MyST-Parser #897](https://github.com/executablebooks/MyST-Parser/issues/897)),
> present in `mystmd` too. Unlabelled `$$ ... $$` is unaffected, as is the
> `{math}` directive.

## Directives in use

`{prf:definition}`, `{prf:theorem}`, `{prf:lemma}`, `{prf:corollary}`,
`{prf:remark}`, `{prf:proof}`, `{exercise}`, `{figure}`, `{code-cell}`, and the
admonitions `{note}`, `{hint}`, `{important}`. Hints inside exercises are nested
`{admonition} Hint` blocks with `:class: hint dropdown`.

## Notebooks

Executable pages are paired `.md`/`.ipynb` via jupytext, format `ipynb,md:myst`
(set in `pyproject.toml` under `[tool.jupytext]`). Commit both halves and keep
them in sync with `jupytext --sync`.

## Previewing

Use `myst start` from this directory. It is the only renderer that resolves the
macros, the `{math}` and `{prf:*}` directives, and the cross-references
together.

The VS Code markdown preview is a rough draft view at best:

* the MyST-Markdown extension renders `$...$` with its own KaTeX call that
  ignores `markdown.math.macros`, so macros never resolve there — no setting
  fixes this;
* its directive support covers only admonitions, `code`, `raw`, `image` and
  `figure`, so `{math}` and every `{prf:*}` block stay unrendered.
