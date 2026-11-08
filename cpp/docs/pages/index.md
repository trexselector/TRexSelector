# The TRexSelector C++ library

**Fast variable selection in high-dimensional settings with control of the false discovery rate (FDR).**

The C++ core library of the Terminating-Random Experiments (T-Rex) selector. The
[R package](https://CRAN.R-project.org/package=TRexSelector) and the
[Python package](https://pypi.org/project/trexselector/) are bindings to it. This guide walks through
the library with example programs that are compiled against the installed library and run when the
documentation is built, so every output shown is the output of the code above it, and the API reference
describes every class and function from its header comments. The
[vignettes](https://CRAN.R-project.org/package=TRexSelector/vignettes/TRexSelector.html) of the
R package and the user guide of the Python package present the same content page for page.

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

<!-- HTML headings: Markdown headings would become sections and show up in the navigation tree
     as if they were pages. -->
<h2>User guide</h2>

1. \subpage getting_started
2. \subpage dimensionality
3. \subpage trex_selector
4. \subpage trex_da
5. \subpage trex_gvs
6. \subpage screening
7. \subpage trex_spca
8. \subpage trex_tikhonov
9. \subpage tsolvers
10. \subpage ml_methods
11. \subpage memory_mapping
12. \subpage no_copies
13. \subpage parallelism

The pages cite the literature, and the \ref citelist "bibliography" collects it. The guide is built with
Doxygen \cite vanheesch2026doxygen, and its figures are drawn with matplotlib \cite hunter2007matplotlib.

<h2>Reference</h2>

- <a href="annotated.html">Classes</a>, the selectors, solvers, ML methods and utilities with their members.
- <a href="namespaces.html">Namespaces</a>, the library's namespace tree.
- <a href="files.html">Files</a>, the headers by directory.
- <a href="https://github.com/trexselector/TRexSelector/blob/main/PUBLICATIONS.md">T-Rex publications</a>,
  every paper of the method family.

<h2>Links</h2>

- Project website: <https://trexselector.com>
- Repository of the C++ core and the R and Python packages: <https://github.com/trexselector/TRexSelector>
- Examples in R, Python and C++: <https://github.com/trexselector/TRexSelector_Examples>
- Bug reports: <https://github.com/trexselector/TRexSelector/issues>
- Team T-Rex, the authors and their contact: <https://github.com/trexselector/TRexSelector/blob/main/AUTHORS.md>
