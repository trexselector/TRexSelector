---
file_format: mystnb
kernelspec:
  name: python3
  display_name: Python 3
---

# Screening: Screen-T-Rex and the biobank screening

```{code-cell} python
:tags: [remove-cell]
# The pages build on at most two threads.
import trexselector.utils
trexselector.utils.set_num_threads(min(2, trexselector.utils.get_max_threads()))
```

## Screening before selecting

A genomics biobank holds hundreds of thousands of variables and hundreds of phenotypes. Running the
full T-Rex calibration on every phenotype is costly, and for most phenotypes the answer is that
nothing is selected. Screen-T-Rex {cite}`machkour2023screen` is a single-pass version of the random
experiments without the L-loop and the T-loop: $K$ experiments with $p$ dummies each, stopped at the
first dummy, and a false discovery rate (FDR) estimate read off the vote. It reports an estimated FDR
instead of controlling a target, so it is a fast first pass. The biobank screening wraps it in a
decision rule per phenotype: when the estimated FDR of the screening falls into an acceptance window,
the screening result is kept, otherwise the full T-Rex selector runs as a fallback.

## Screen-T-Rex

The screening selector takes no target FDR. The data are generated with a fixed seed, five active
variables of coefficient 2 among 300, and the indices are 0-based.

```{code-cell} python
import numpy as np
import trexselector as tsn

rng = np.random.default_rng(2024)
n, p = 100, 300
true_support = [2, 49, 119, 199, 279]
beta = np.zeros(p)
beta[true_support] = 2.0
X = rng.standard_normal((n, p))
y = X @ beta + rng.standard_normal(n)
```

```{code-cell} python
ctrl = tsn.TRexControlParameter()
ctrl.K = 20

scr = tsn.TRexScreeningSelector(X, y, trex_control=ctrl, seed=3, verbose=False)
res = scr.select()
print("selected:     ", scr.selected_indices)
print("estimated FDR:", res.estimated_FDR)
```

`beta_mat` holds the coefficients of the original variables at the stop, averaged over the
experiments, and `dummy_betas` the non-zero coefficients of the dummies collected during the
screening, the material the estimate is read off. The L-loop strategy is forced to `STANDARD` here,
fresh dummies per experiment.

```{code-cell} python
print("beta_mat:   ", res.beta_mat.shape)
print("dummy_betas:", res.dummy_betas.shape)
```

### A bootstrap confidence interval

`use_bootstrap_CI = True` in the `ScreenTRexControlParameter` resamples the experiments `R_boot`
times and reports the upper confidence bound of the estimated FDR instead of the point estimate, at
`confidence_level`. `ci_grid_step` is the resolution of the bound.

```{code-cell} python
boot_ctrl = tsn.ScreenTRexControlParameter()
boot_ctrl.use_bootstrap_CI = True
boot_ctrl.R_boot = 200

boot = tsn.TRexScreeningSelector(X, y, screen_control=boot_ctrl, trex_control=ctrl,
                                 seed=3, verbose=False)
boot_res = boot.select()
print("used_bootstrap:  ", boot_res.used_bootstrap)
print("confidence_level:", boot_res.confidence_level)
print("estimated FDR:   ", boot_res.estimated_FDR)
print("selected:        ", boot.selected_indices)
```

### Dependent variables

`trex_method` switches the screening to dependency-aware dummies: `ScreenTRexMethod.TREX_DA_AR1`,
`TREX_DA_EQUI` and `TREX_DA_BLOCK_EQUI` draw autoregressive, equi-correlated and block
equi-correlated dummies. `estimated_correlation` reports the correlation the dummies were drawn with.

```{code-cell} python
ar_ctrl = tsn.ScreenTRexControlParameter()
ar_ctrl.trex_method = tsn.ScreenTRexMethod.TREX_DA_AR1

ar = tsn.TRexScreeningSelector(X, y, screen_control=ar_ctrl, trex_control=ctrl,
                               seed=3, verbose=False)
ar_res = ar.select()
print("estimated correlation:", round(ar_res.estimated_correlation, 3))
```

## The biobank screening

`TRexBiobankScreeningSelector` screens one or many phenotypes against one genotype matrix. `Y` is a
vector for one phenotype or a matrix with one column per phenotype. The control is a
`BiobankScreenTRexControl`: the acceptance window of the estimated FDR is `lower_bound_FDR` to
`upper_bound_FDR`, `target_FDR_trex` is the target of the fallback selector, and the screening
control with its base control is nested in `trex_screen_ctrl`.

```{code-cell} python
Y = np.column_stack([y, X @ (0.5 * beta) + rng.standard_normal(n), rng.standard_normal(n)])
phenotypes = ["strong", "weak", "null"]

screen_ctrl = tsn.ScreenTRexControlParameter()
screen_ctrl.trex_ctrl = ctrl
bio_ctrl = tsn.BiobankScreenTRexControl()
bio_ctrl.lower_bound_FDR = 0.05
bio_ctrl.upper_bound_FDR = 0.15
bio_ctrl.target_FDR_trex = 0.1
bio_ctrl.trex_screen_ctrl = screen_ctrl

bio = tsn.TRexBiobankScreeningSelector(X, Y, bio_ctrl=bio_ctrl, seed=3, verbose=False)
records = bio.select()
len(records)
```

Each record names the phenotype, the method that produced the result and whether the fallback ran.

```{code-cell} python
print(f"{'phenotype':<10}{'selected':>9}{'estimated_FDR':>15}  {'method':<24}{'fallback':>9}")
for r in records:
    print(f"{phenotypes[r.phenotype_index]:<10}{len(r.selected_indices):>9}"
          f"{round(r.estimated_FDR, 3):>15}  {r.method_used:<24}{str(r.used_fallback_trex):>9}")
print("selected for 'strong':", records[0].selected_indices)
```

The record also keeps the ordinary and the bootstrap screening results
(`selected_indices_screen_ordinary`, `estimated_FDR_screen_bootstrap`, ...), so the decision can be
audited. For a single phenotype `select()` returns one record instead of a list.

## Large designs

Both screening selectors accept a memory-mapped design (the `to_numpy()` view of a
`MemoryMappedMatrix`), which [the memory-mapping page](11-memory-mapping.md) describes. The
combination of the screening and a mapped design is what solved sparse regression problems with five
million variables on a laptop {cite}`scheidt2023bigtrex`.

## Examples

The [TRexSelector_Examples](https://github.com/trexselector/TRexSelector_Examples/tree/main)
repository holds runnable demos and Monte Carlo studies of Screen-T-Rex and the biobank screening,
in memory and memory-mapped, with their results and a README per demo, in
[R](https://github.com/trexselector/TRexSelector_Examples/tree/main/R/trex_selector_methods/ad/trex_screening),
[Python](https://github.com/trexselector/TRexSelector_Examples/tree/main/Python/trex_selector_methods/ad/trex_screening)
and
[C++](https://github.com/trexselector/TRexSelector_Examples/tree/main/cpp/trex_selector_methods/ad/trex_screening).