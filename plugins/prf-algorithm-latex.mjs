/*
 * Typeset {prf:algorithm} blocks as `algorithm` floats in LaTeX/PDF exports.
 *
 * mystmd's LaTeX writer maps proof kinds to amsthm environments through a
 * hard-coded list (theorem, lemma, definition, remark, ...). "algorithm" is not
 * on it, so the export fails with
 *
 *     Unhandled LaTeX proof environment "algorithm"
 *
 * and the block is dropped from the PDF altogether. This transform rewrites
 * each algorithm into raw LaTeX that brackets the original body:
 *
 *     \begin{algorithm}[H]
 *     \caption{<title>}\label{<identifier>}
 *     <body, converted by mystmd as usual: lists, math, bold, ...>
 *     \end{algorithm}
 *
 * The `algorithm` environment comes from the `algorithm` package (loaded by
 * templates/tma4220-book/template.tex). Only the float, caption and numbering
 * are taken from it; the steps stay the Markdown list they are written as.
 * The label is the directive's `:label:`, so `@alg:...` references resolve to
 * the float's number.
 *
 * Plugin transforms are not told which output is being built, and the same
 * rewrite would break the website, so the build target is read from the
 * command line. See `isLatexBuild` below.
 */

/** Output flags of `myst build` / `jupyter book build`. */
const SITE_FLAGS = ["--site", "--html", "--all"];
const LATEX_FLAGS = ["--pdf", "--tex"];
const OTHER_EXPORT_FLAGS = [
  "--docx",
  "--typst",
  "--jats",
  "--meca",
  "--md",
  "--xml",
  "--cff",
  "--ipynb",
];

/**
 * Whether this process builds LaTeX-based exports and not the website.
 *
 * - `--pdf` / `--tex`: yes.
 * - a bare `build` with no output flag builds the exports declared in
 *   myst.yml, which for this book are PDF and TeX: yes.
 * - `start`, `--site`, `--html`: no, the website renders prf:algorithm itself.
 * - site and LaTeX in one run (`--all`, `--site --pdf`): each output reuses the
 *   same transforms, so they cannot both be served. The website wins, and a
 *   warning says to build the PDF on its own.
 */
let warnedMixedBuild = false;

function isLatexBuild(argv) {
  const buildsSite = argv.includes("start") || argv.some((a) => SITE_FLAGS.includes(a));
  const outputFlags = [...SITE_FLAGS, ...LATEX_FLAGS, ...OTHER_EXPORT_FLAGS];
  const buildsLatex =
    argv.some((a) => LATEX_FLAGS.includes(a)) ||
    (argv.includes("build") && !argv.some((a) => outputFlags.includes(a)));

  if (buildsLatex && buildsSite) {
    // The transform runs once per page; say it once per build.
    if (!warnedMixedBuild) {
      warnedMixedBuild = true;
      console.warn(
        "[prf-algorithm-latex] website and LaTeX built in one run: {prf:algorithm} " +
          "is left as is, so the LaTeX export will report it as unhandled. " +
          "Build the website and the PDF in separate runs, e.g. " +
          "`jupyter book build --html` then `jupyter book build --pdf`.",
      );
    }
    return false;
  }
  return buildsLatex;
}

const raw = (tex) => ({ type: "raw", tex });

/** Replace one algorithm proof node, in place, by its LaTeX rendering. */
function toAlgorithmFloat(node) {
  const children = node.children ?? [];
  const titleNode = children.find((c) => c.type === "admonitionTitle");
  const body = children.filter((c) => c !== titleNode);
  const label = node.identifier ? `\\label{${node.identifier}}` : "";

  // The caption is a paragraph so the title's inline nodes (text, inline math)
  // are converted by mystmd instead of being pasted in as a raw string.
  const caption = {
    type: "paragraph",
    children: [raw("\\caption{"), ...(titleNode?.children ?? []), raw(`}${label}`)],
  };

  for (const key of Object.keys(node)) {
    if (key !== "position") delete node[key];
  }
  Object.assign(node, {
    type: "div",
    children: [raw("\\begin{algorithm}[H]\n"), caption, ...body, raw("\\end{algorithm}\n")],
  });
}

const algorithmTransform = {
  name: "prf-algorithm-latex",
  doc: "Typeset {prf:algorithm} as an `algorithm` float in LaTeX/PDF exports.",
  // "project" runs after cross-references are resolved. At the earlier
  // "document" stage the rewrite would remove the reference target, and
  // `@alg:...` would fall back to being parsed as a citation.
  stage: "project",
  plugin: (_opts, utils) => (tree) => {
    if (!isLatexBuild(process.argv.slice(2))) return;
    utils.selectAll("proof[kind=algorithm]", tree).forEach(toAlgorithmFloat);
  },
};

export default {
  name: "prf:algorithm for LaTeX",
  transforms: [algorithmTransform],
};
