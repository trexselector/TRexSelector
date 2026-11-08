# Screening: Screen-T-Rex and the biobank screening {#screening}

[TOC]

## Screening before selecting

A genomics biobank holds hundreds of thousands of variables and hundreds of phenotypes. Running the
full T-Rex calibration on every phenotype is costly, and for most phenotypes the answer is that
nothing is selected. Screen-T-Rex \cite machkour2023screen is a single-pass version of the random
experiments without the L-loop and the T-loop: \f$K\f$ experiments with \f$p\f$ dummies each, stopped
at the first dummy, and a false discovery rate (FDR) estimate read off the vote. It reports an
estimated FDR instead of controlling a target, so it is a fast first pass. The biobank screening wraps
it in a decision rule per phenotype: when the estimated FDR of the screening falls into an acceptance
window, the screening result is kept, otherwise the full T-Rex selector runs as a fallback.

## Screen-T-Rex

The screening selector takes no target FDR. The headers and names of this page:

\snippet{trimleft} screening.cpp includes

A small helper prints index vectors:

\snippet{trimleft} screening.cpp helpers

The data are generated with a fixed seed, five active variables of coefficient 2 among 300, and the
indices are 0-based.

\snippet{trimleft} screening.cpp data

The screening control nests the base control, so the number of random experiments `K` is set in its
`trex_ctrl`. `select()` returns the base result, and `getScreenResult()` the screening result with
the estimated FDR.

\snippet{trimleft} screening.cpp screen
\output{screening.screen.txt}

`beta_mat` holds the coefficients of the original variables at the stop, averaged over the
experiments, and `dummy_betas` the non-zero coefficients of the dummies collected during the
screening, the material the estimate is read off. The L-loop strategy is forced to `STANDARD` here,
fresh dummies per experiment.

\snippet{trimleft} screening.cpp mats
\output{screening.mats.txt}

### A bootstrap confidence interval

`use_bootstrap_CI = true` in the ScreenTRexControlParameter resamples the experiments `R_boot`
times and reports the upper confidence bound of the estimated FDR instead of the point estimate, at
`confidence_level`. `ci_grid_step` is the resolution of the bound.

\snippet{trimleft} screening.cpp boot
\output{screening.boot.txt}

### Dependent variables

`trex_method` switches the screening to dependency-aware dummies: `ScreenTRexMethod::TREX_DA_AR1`,
`TREX_DA_EQUI` and `TREX_DA_BLOCK_EQUI` draw autoregressive, equi-correlated and block
equi-correlated dummies. `estimated_correlation` reports the correlation the dummies were drawn with.
The design of this page has independent columns, so the estimate is close to zero.

\snippet{trimleft} screening.cpp da
\output{screening.da.txt}

## The biobank screening

BiobankScreenTRex screens one or many phenotypes against one genotype matrix. The constructor takes
`Y` as a vector map for one phenotype or as a matrix map with one column per phenotype. The control
is a BiobankScreenTRexControl: the acceptance window of the estimated FDR is `lower_bound_FDR` to
`upper_bound_FDR`, `target_FDR_trex` is the target of the fallback selector, and the screening
control with its base control is nested in `trex_screen_ctrl`. `screenPhenotypes()` returns one
record per phenotype.

\snippet{trimleft} screening.cpp biobank
\output{screening.biobank.txt}

Each record names the phenotype, the method that produced the result and whether the fallback ran.

\snippet{trimleft} screening.cpp records
\output{screening.records.txt}

The record also keeps the ordinary and the bootstrap screening results
(`selected_indices_screen_ordinary`, `estimated_FDR_screen_bootstrap`, ...), so the decision can be
audited. For a single phenotype `screenPhenotype()` returns one record instead of a vector.

## Large designs

Both screening selectors accept a memory-mapped design (the `getMap()` view of a
MemoryMappedMatrix), which \ref memory_mapping "the memory-mapping page" describes. The combination
of the screening and a mapped design is what solved sparse regression problems with five million
variables on a laptop \cite scheidt2023bigtrex.
