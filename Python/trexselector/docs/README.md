# The Python user guide

The user guide of the Python package, thirteen pages with the content of the R vignettes, next to the API
reference that Sphinx generates from the docstrings of the installed package. The guide pages are MyST
notebooks in Markdown form. Sphinx executes their code cells when the documentation is built, so every
output shown is the output of the code above it.

## Building

Prerequisites: the package installed in the active environment (`pip install trexselector`, or the
sources of this repository) and the toolchain of `requirements.txt`.

```bash
python -m pip install -r Python/trexselector/docs/requirements.txt
python -m sphinx -b html Python/trexselector/docs Python/trexselector/docs/_build/html
open Python/trexselector/docs/_build/html/index.html
```

A page whose cell raises an error fails the build (`nb_execution_raise_on_error`). The full build with
the execution of all pages takes well under a minute.

## Layout

- `guide/NN-<name>.md`, the thirteen pages, each with a hidden first cell that limits the run to two
  threads, and `guide/figures/` with the static figures. `guide/matplotlibrc` is the style of every
  figure a page draws, read by matplotlib from the working directory, and the C++ guide's figure script
  loads the same file.
- `api/`, the reference in four pages (selectors, T-solvers, ML methods, utilities) from the docstrings.
  `conf.py` shortens the pybind11 signatures and joins the reference lists of the docstrings.
- `references.md`, the site-wide bibliography. The pages cite with ``{cite}`key` ``, the keys come from
  the repository's `references.bib`.
- `_static/custom.css`, the few style rules on top of the furo theme, and `_static/pages_nav.js`, the
  position "n / 13" between the previous and next links of a guide page.
