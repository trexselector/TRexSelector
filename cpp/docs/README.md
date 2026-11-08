# The C++ user guide

The user guide of the C++ library, thirteen pages with the content of the R vignettes and the Python user
guide, next to the API reference that Doxygen generates from the headers under `../src`.

The code on a page comes from an example program under `examples/`. The build compiles every program
against the installed library, runs it, and captures the output of each documented section into
`examples/output/`. The pages include the code with `\snippet` and the output with the `\output` alias,
so every output shown is the output of the code above it.

## Building

Prerequisites: the library installed under `cpp/install/` (see the README of the library), Doxygen 1.18
or newer, Graphviz for the class diagrams, and `bibtex` with `perl` for the bibliography. The figures are
drawn with matplotlib from the data the examples write, by `figures/make_figures.py` with the style of
the Python guide. `build_docs.sh` runs it with `TREX_DOCS_PYTHON` (default `python3`) when that Python
has matplotlib and numpy, and otherwise keeps the committed PNG files under `pages/figures/`.

```bash
bash cpp/docs/build_docs.sh
open cpp/docs/_build/html/index.html
```

After Doxygen, `tools/add_doi_links.py` appends the DOI of every entry of the citation list from
`references.bib`, since Doxygen's own bibliography style prints none.

Doxygen 1.18.0 as released dies at random with a bus error in its comment scanner, about one run in
three on this guide, because a command without braces makes it read two bytes before the scanner
buffer. Upstream fixed it for 1.19.0 (commit 26e9087). Until that release, build Doxygen 1.18.0 with
`tools/doxygen-1.18.0-commentscan-npos.patch`, with Homebrew through a local tap whose formula
carries the patch.
The script ends with `=== done, exit 0` when the examples compiled, ran, and Doxygen reported no error.
`TREX_INSTALL_PREFIX` points it at another install prefix. The warnings of a run are in
`_build/doxygen_warnings.log`, and a clean run has none.

## Writing a page

- `pages/NN-<name>.md` starts with `# <Title> {#<id>}`, then `[TOC]`, then the sections. Math is
  `\f$ ... \f$` inline and `\f[ ... \f]` displayed, citations are `\cite <key>` with the keys of the
  repository's `references.bib`, figures live under `pages/figures/`.
- `examples/<id>.cpp` holds the code. `DocSections doc("<id>", DOC_OUTPUT_DIR);` at the top of `main()`,
  `doc.begin("<section>");` before every `//! [<section>]` block, and the page shows the block with
  `\snippet{trimleft} <id>.cpp <section>` followed by `\output{<id>.<section>.txt}`.
- The example's id goes into `GUIDE_EXAMPLES` of `examples/CMakeLists.txt`, the page into the lists of
  `pages/index.md` and `pages_nav.js`.
- A figure: the example writes its data into `DOC_OUTPUT_DIR` (`<id>.<name>.csv` and friends), a
  function in `figures/make_figures.py` draws the PNG into `pages/figures/`, and the page shows it with
  `\image html <name>.png "<caption>" width=90%`.
