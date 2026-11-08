# The trexselector Python package

**Fast variable selection in high-dimensional settings with control of the false discovery rate (FDR).**

The Python package of the Terminating-Random Experiments (T-Rex) selector, with bindings to the
[TRexSelector C++ library](https://github.com/trexselector/TRexSelector). The
[R package](https://CRAN.R-project.org/package=TRexSelector) binds the same library. This guide walks
through the package with code that runs when it is built, so every output shown is the output of the
code above it, and the API reference describes every class and function from its docstring. The
[vignettes](https://cran.r-project.org/web/packages/TRexSelector/vignettes/TRexSelector.html) of the
R package and the user guide of the C++ library present the same content page for page.

The method family grows within a series of papers, among which the following three form the theoretical
core of the framework so far:

> J. Machkour, M. Muma, and D. P. Palomar, "The terminating-random experiments selector: Fast
> high-dimensional variable selection with false discovery rate control," *Signal Processing*,
> vol. 231, article 109894, 2025. (<https://doi.org/10.1016/j.sigpro.2025.109894>)
>
> J. Machkour, M. Muma, and D. P. Palomar, "High-dimensional false discovery rate control for
> dependent variables," *Signal Processing*, vol. 234, article 109990, 2025.
> (<https://doi.org/10.1016/j.sigpro.2025.109990>)
>
> T. Koka, J. Machkour, D. P. Palomar, and M. Muma, "Virtual dummies: Enabling scalable FDR-controlled
> variable selection via sequential sampling of null features," arXiv:2604.07464, 2026.
> (<https://arxiv.org/abs/2604.07464>)

```{toctree}
:maxdepth: 1
:numbered: 1
:caption: User guide

Getting started <guide/01-getting-started>
Dimensionality and FDR <guide/02-dimensionality>
The T-Rex selector <guide/03-trex-selector>
T-Rex+DA <guide/04-trex-da>
T-Rex GVS <guide/05-trex-gvs>
Screen-T-Rex <guide/06-screening>
T-Rex SPCA <guide/07-trex-spca>
T-Rex Tikhonov <guide/08-trex-tikhonov>
The T-solvers <guide/09-tsolvers>
The ML methods <guide/10-ml-methods>
Memory mapping <guide/11-memory-mapping>
Without data copies <guide/12-no-copies>
Parallelism <guide/13-parallelism>
```

The pages cite the literature, and the [bibliography](references.md) collects it. The guide is built with
Sphinx {cite}`brandl2026sphinx` and MyST-NB {cite}`executablebooks2026mystnb`, and its figures are drawn
with matplotlib {cite}`hunter2007matplotlib`.

```{toctree}
:maxdepth: 1
:caption: Reference

api/index
references
T-Rex publications <https://github.com/trexselector/TRexSelector/blob/main/PUBLICATIONS.md>
```

## Links

- Project website: <https://trexselector.com>
- Repository of the C++ core and the R and Python packages: <https://github.com/trexselector/TRexSelector>
- Examples in R, Python and C++: <https://github.com/trexselector/TRexSelector_Examples>
- Bug reports: <https://github.com/trexselector/TRexSelector/issues>
- Team T-Rex, the authors and their contact: <https://github.com/trexselector/TRexSelector/blob/main/AUTHORS.md>
