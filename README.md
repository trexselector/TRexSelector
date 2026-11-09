# T-REX SELECTOR

**FDR-controlled variable selection in high-dimensional regression and classification.
A C++20 library with R and Python bindings.**

---

## A Short Story about the T-Rex Selector

False discoveries plague science and engineering.
In genomics, a single false discovery can wrongly tie a gene to a disease, contributing to the
reproducibility crisis and derailing years of research and clinical effort
 [(Huffman, 2018)](https://doi.org/10.1038/s41467-018-07348-x).
In finance, hundreds of studies claim to have found factors that explain stock returns, yet most of
these findings are likely false [(Harvey, Liu, and Zhu, 2016)](https://doi.org/10.1093/rfs/hhv059).

The T-Rex Selector selects variables in low- and high-dimensional settings alike, and it is tailored
to the high-dimensional case, where the number of variables exceeds the number of observations. It
controls the rate of false discoveries at a level you choose, with rigorous finite-sample guarantees,
and it is fast and memory-efficient enough to run on millions of variables
 [(Machkour, Muma, and Palomar, 2025)](https://www.sciencedirect.com/science/article/pii/S016516842500009X).

The algorithm works by generating random dummy variables that act as a calibration baseline:
during each of *K* independent random experiments, a terminating solver traces solution paths
through both the real predictors and the dummies.
By comparing how often each real predictor is selected relative to the dummies across all
experiments, a voting threshold is derived that provably bounds the FDR.

**Project website:** [trexselector.com](https://trexselector.com)

**Team T-Rex:** [AUTHORS.md](AUTHORS.md) names the members with their contact and the discussions page.

**Publications:** [PUBLICATIONS.md](PUBLICATIONS.md) lists every paper of the method family.

---

## Implemented Variants

| Class                          | Description                                               |
|--------------------------------|-----------------------------------------------------------|
| `TRexSelector`                 | Core FDR-controlled selector for linear regression        |
| `TRexDASelector`               | Dependency-aware variant for correlated predictors        |
| `TRexGVSSelector`              | Grouped variable selection                                |
| `TRexScreeningSelector`        | Fast FDR-controlled variable pre-screening for large *p*  |
| `TRexBiobankScreeningSelector` | High-throughput genomic / biobank screening               |
| `TRexSPCASelector`             | FDR-controlled sparse principal component analysis        |
| `TRexTikhonovSelector`         | Selection with a general Tikhonov penalty matrix          |

---

## Repository Layout

This repository hosts the C++ core and R and Python language packages.
Results are identical across languages for the same data and parameters.
Each component has its own README with installation instructions, examples, and developer
documentation:

* **C++ core library** — [`cpp/`](cpp), documented in [cpp/README.md](cpp/README.md):
  building the library and tests, and reaching it from your own C++ project via
  `find_package(TRexSelector)`.
* **Python package** — [`Python/trexselector/`](Python/trexselector), documented in
  [Python/trexselector/README.md](Python/trexselector/README.md): installation, quick
  start, and development setup. Distributed on PyPI and imported as
  `trexselector`.
* **R package** — [`R/TRexSelector/`](R/TRexSelector), documented in
  [R/TRexSelector/README.md](R/TRexSelector/README.md): installation, quick start, and
  the thirteen vignettes. Distributed on CRAN as `TRexSelector`.

---

## Citing

The file [CITATION.cff](CITATION.cff) carries the citation metadata of the repository, and the R
 package answers `citation("TRexSelector")`.
Every paper the class documentation credits is listed in [references.bib](references.bib), the
 master file from which the reference blocks of the C++, R and Python documentation are rendered
 with the scripts in [tools/references/](tools/references).
The bibliography below is rendered from it as well, and [PUBLICATIONS.md](PUBLICATIONS.md) lists
 every paper of the method family.

## Bibliography

The papers behind the opening story:

<!-- Rendered from references.bib:
     python3 tools/references/render_references.py --format md machkour2025trex huffman2018biobanks
         harvey2016crosssection -->
* J. Machkour, M. Muma, and D. P. Palomar, "The Terminating-Random Experiments Selector: Fast High-Dimensional Variable Selection with False Discovery Rate Control," Signal Process., vol. 231, Art. no. 109894, 2025, [doi: 10.1016/j.sigpro.2025.109894](https://doi.org/10.1016/j.sigpro.2025.109894).
* J. E. Huffman, "Examining the current standards for genetic discovery and replication in the era of mega-biobanks," Nat. Commun., vol. 9, no. 1, Art. no. 5054, 2018, [doi: 10.1038/s41467-018-07348-x](https://doi.org/10.1038/s41467-018-07348-x).
* C. R. Harvey, Y. Liu, and H. Zhu, "… and the Cross-Section of Expected Returns," Rev. Financ. Stud., vol. 29, no. 1, pp. 5-68, 2016, [doi: 10.1093/rfs/hhv059](https://doi.org/10.1093/rfs/hhv059).

The C++ core builds on the following software:

<!-- Rendered from references.bib:
     python3 tools/references/render_references.py --format md guennebaud2010eigen boost2026libraries
         grant2017cereal openmp2021spec -->
* G. Guennebaud, B. Jacob et al., Eigen v3, <https://libeigen.gitlab.io/>, 2010.
* The Boost Community, Boost C++ Libraries, <https://www.boost.org/>, 2026.
* W. S. Grant and R. Voorhies, cereal: A C++11 Library for Serialization, <https://uscilab.github.io/cereal/>, 2017.
* OpenMP Architecture Review Board, OpenMP Application Programming Interface, Version 5.2, <https://www.openmp.org/specifications/>, 2021.

---

## License

GPL ≥ 3 — see [LICENSE](LICENSE) or the [GNU General Public License v3](https://www.gnu.org/licenses/gpl-3.0.html).

```text
This program is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. It is distributed WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
```
