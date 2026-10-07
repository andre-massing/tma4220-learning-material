/*
 * Typeset {exercise} blocks as a numbered `exercise` environment in LaTeX/PDF
 * exports.
 *
 * mystmd's LaTeX writer has no handler for exercise nodes and silently drops
 * them, so neither the exercise nor the target of a reference to it reaches
 * the export, and every `{prf:ref}` to an exercise prints as "??". This
 * transform rewrites each exercise into raw LaTeX that brackets the original
 * body:
 *
 *     \begin{exercise}[<title>]\label{<identifier>}
 *     <body, converted by mystmd as usual: lists, math, bold, ...>
 *     \end{exercise}
 *
 * The `exercise` environment is declared by templates/tma4220-book/template.tex
 * and numbered per chapter: "Exercise 3.2". Hints are separate admonitions
 * following the exercise and are not touched.
 *
 * Plugin transforms are not told which output is being built, and the same
 * rewrite would break the website, so the build target is read from the
 * command line, as in prf-algorithm-latex.mjs.
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
 * Same rules as `isLatexBuild` in prf-algorithm-latex.mjs.
 */
let warnedMixedBuild = false;

function isLatexBuild(argv) {
  const buildsSite = argv.includes("start") || argv.some((a) => SITE_FLAGS.includes(a));
  const outputFlags = [...SITE_FLAGS, ...LATEX_FLAGS, ...OTHER_EXPORT_FLAGS];
  const buildsLatex =
    argv.some((a) => LATEX_FLAGS.includes(a)) ||
    (argv.includes("build") && !argv.some((a) => outputFlags.includes(a)));

  if (buildsLatex && buildsSite) {
    if (!warnedMixedBuild) {
      warnedMixedBuild = true;
      console.warn(
        "[exercise-latex] website and LaTeX built in one run: {exercise} is left " +
          "as is, so the LaTeX export will drop it. Build the website and the PDF " +
          "in separate runs, e.g. `jupyter book build --html` then " +
          "`jupyter book build --pdf`.",
      );
    }
    return false;
  }
  return buildsLatex;
}

const raw = (tex) => ({ type: "raw", tex });

/** Replace one exercise node, in place, by its LaTeX rendering. */
function toExerciseEnvironment(node) {
  const children = node.children ?? [];
  const titleNode = children.find((c) => c.type === "admonitionTitle");
  const body = children.filter((c) => c !== titleNode);
  const label = node.identifier ? `\\label{${node.identifier}}` : "";
  const titleChildren = titleNode?.children ?? [];

  // The opening line is a paragraph so the title's inline nodes (text, inline
  // math) are converted by mystmd instead of being pasted in as a raw string.
  // Braces around the title protect a `]` inside it.
  const opening = {
    type: "paragraph",
    children: titleChildren.length
      ? [raw("\\begin{exercise}[{"), ...titleChildren, raw(`}]${label}`)]
      : [raw(`\\begin{exercise}${label}`)],
  };

  for (const key of Object.keys(node)) {
    if (key !== "position") delete node[key];
  }
  Object.assign(node, {
    type: "div",
    children: [opening, ...body, raw("\\end{exercise}\n")],
  });
}

const exerciseTransform = {
  name: "exercise-latex",
  doc: "Typeset {exercise} as a numbered `exercise` environment in LaTeX/PDF exports.",
  // "project" runs after cross-references are resolved; at the "document"
  // stage the rewrite would remove the reference target.
  stage: "project",
  plugin: (_opts, utils) => (tree) => {
    if (!isLatexBuild(process.argv.slice(2))) return;
    utils.selectAll("exercise", tree).forEach(toExerciseEnvironment);
  },
};

export default {
  name: "Exercises for LaTeX",
  transforms: [exerciseTransform],
};
