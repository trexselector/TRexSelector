# TRexSelector changelog

## TRexSelector C++ 1.0.1 / R 2.0.2 / Python 1.0.1 (in preparation)

The fix release after R 2.0.1, not published yet. It collects the work since the sources of R 2.0.1
went to CRAN, and it changes results in several places. The dependency-aware selector takes its
penalty from the last step for all methods, and its AR1 window follows the magnitude of the
coefficient. The Tikhonov selector runs on the homotopy solver TIKNET by default. The biobank
screening gives every phenotype its own seed, and the screening selector reports an FDR estimate of 1
for an empty selection. A coordinate-descent solve that has not converged is an error. Beside the
corrections the screening selector gains a nearest-neighbor variant, a read-only memory-mapped design
serves where nothing is written into it, the coordinate-descent solver family is deprecated ahead of
its removal in the next minor releases, and the R package compiles in less time. Versions: C++ SDK
1.0.1, R TRexSelector 2.0.2, Python trexselector 1.0.1. The DESCRIPTION of the R package reads 2.0.1
until the release is prepared.

### 2026-11-15

#### The memory tests measure the live heap and no longer fail under load

MLMemory.RidgeCVLowMemory failed once (0.5007 against a budget of 0.5 times the design) in a parallel test
run beside two builds. Fabian: "failing tests are not accepted. We try it again, otherwise we must refine
the test. No commit before we have not found a solution." Evidence with all series:
notes/handover_2026-11-15/memory_test/.

- The cause. The tests read the peak resident set of the process over one call. On macOS that is no
  constant of the call: a cross-validation runs its fits one after another, and how the allocator
  reuses the freed blocks of one fit for the next differs between runs and with the load. The
  low-memory ridge CV gave 0.20 to 0.54 for identical input. Reproduced in 2 of 20 repetitions of
  `ctest -j 12` and in 3 and 6 of 60 runs with twelve copies at once, never idle. One BLAS thread gives
  the same spread. Under the same load RidgeCV reached 3.95 of 4.0 and ElasticNetCVLowMemory 0.287 of
  0.3. On Ubuntu 24.04 with Eigen 3.4 the value is stable (RidgeCV 4.57 to 4.61).
- The first attempt, wider budgets for the three cross-validation tests, held in every stressed run and
  was dropped: the macOS value has no ceiling, and one fold block that is not reused looks like a copy
  of the design.
- The cure. The tests measure the bytes that are live on the heap, from `malloc_zone_statistics` on
  macOS and `mallinfo2` on glibc 2.33 or later, sampled by a second thread without a pause while the
  call runs. The peak depends on what the call allocates. The sampler can miss a short peak and never
  reports more than the peak, so load lowers the value at most. On every other platform the file holds
  one skipped test, as before on Windows. The tests no longer need a process of their own.
- Measured peaks, macOS with Eigen 5.0 and Ubuntu with Eigen 3.4: ridge CV 2.84 and 4.41, elastic net
  1.09 and 1.09, elastic net CV 1.34 and 1.09, SVD 0.33 and 0.22, ridge 0.27 and 0.08, PCA 0.34 and
  0.22, scaler 0.006 and 0.004, low-memory ridge CV 0.18 and 0.19, low-memory elastic net 0.093 and
  0.091, low-memory elastic net CV 0.094 and 0.091. Idle, 15 runs per test: one value per test, except
  the elastic net CV on macOS, whose peak is very short (1.29 to 1.34).
- The budgets stay as they were. Only ridge CV changes, 4.0 to 3.3 under Eigen 5 and 5.0 to 4.9 under
  Eigen 3.4, half a design above the peaks, because the old ones would let a copy of the design pass
  under the new measure. A copy held over the call raises every test by exactly 1.0 and fails all ten.
- After: twelve copies at once, 456 runs on macOS and 456 on Ubuntu, no failed run and no value above
  the idle value. 20 repetitions of `ctest -j 12` without a failure. C++ 751 of 751.

#### The screening: what two simulations showed, in the documentation

Fabian, 2026-11-15, after the two Monte Carlo runs of the day (notes/handover_2026-11-15/da_nn_sim/ and
global_null_sim/, 1,000 replicates each, with an ANALYSIS.md per study): `TREX_DA_NN` stays in the fix
release (decision D5), the notes go into the documentation, and the bootstrap finding is documented ("the
screening was never published in combination with DA awareness, having some results can be insightful
for more and future work").

- `TREX_DA_NN` against `TREX_DA_AR1` and the plain screening, 34 settings, 375,000 runs. At the default
  threshold 0.02 the paired difference is -0.002 to +0.004 in FDR and 0 to +0.022 in TPR, and 95 % of
  the runs end with the same result. Both keep the expected number of false discoveries below the
  bound of the screening paper, p0 / (p + 1), in 22 of 27 dependent settings and fail in the same five
  (blocks at a correlation of 0.8, 30 actives of 300). The plain screening keeps it in 1 of the 27.
  With larger thresholds the variant buys power with false discoveries. The rule "no penalty for an
  empty group" is right on independent columns and removes the protection only when the threshold
  reaches the correlation of the data. The help of `rho_thr_DA` and of the variant says so.
- The bootstrap rule with a dependency-aware variant ended with an empty selection in 90 % to 100 % of
  the runs, for AR1 as for NN, where the plain screening is empty in 48 % to 73 %. The filter lowers the
  ordinary count that the confidence level has to reach. The help of `use_bootstrap_CI` says so. The
  rule is not changed.
- The equi-correlated regime, recorded: an FDR of 0.80 for the plain screening, AR1 and NN, and
  `TREX_DA_EQUI` selects nothing.
- The global null, 22 settings, 154,000 runs. The T-Rex selector is empty in 98.9 % of the runs, and
  the share of non-empty selections is below the target in all 66 cells (means 0.0001, 0.006 and 0.027
  at the targets 0.05, 0.1 and 0.2). The screening selects something in 24 % to 50 % of the runs with
  0.27 to 0.71 variables on average, which keeps the bound of the proof of Theorem 1. The statement
  "FDR <= 1 / R" holds in the mean and not run by run: with two or more selected variables (13 % of the
  null runs of the plain screening) the estimate is 0.5 or less and the FDP is 1. The window of the
  biobank screening covers the case (3 of 22,000 null runs fall into it). The class text gains the
  section "A response without signal" on the three surfaces, and the help of `estimated_FDR` points
  to it.
- Guide page 06 on the three surfaces carries the three notes in the sections on the bootstrap, on the
  dependent variables and in a new section. Its sentence "The L-loop strategy is forced to STANDARD
  here" was wrong since eea219a and is replaced.
- NEWS of the R package announces `trex_method = "TREX_DA_NN"`, which had no entry.
- Texts only. C++ 751 of 751, R 2,218, Python 798.

#### The guides and the READMEs use the short forms

T-Rex+GVS, T-Rex+Screening, T-Rex+SPCA and T-Rex+Tikhonov replace "T-Rex GVS", "Screen-T-Rex", "T-Rex
SPCA" and "T-Rex Tikhonov" in 95 lines of 28 files: the titles, tables and running text of the guide
pages 01, 02 and 05 to 08 on the three surfaces, the landing pages, the READMEs of R, Python and C++,
the two navigation files and three C++ guide examples. The search index of the R vignettes is rebuilt
(the heading "T-Rex+Screening" has the anchor `t-rexscreening`), and the README of the R package is
knit again. The heading "The elastic net policy" of page 05 stays with the parked page review.

#### The old method names leave the comments, the help texts and the progress messages

Fabian, 2026-11-14: the short forms are T-Rex+DA, T-Rex+GVS and T-Rex+Screening, and the names
"DA-TRex / GVS-TRex / Screen-TRex" change with the step of the core selector.

- 79 lines of the screening sources, the experiment runner, the R glue of the screening and the Python
  bindings. "Screen-TRex" is T-Rex+Screening, "Biobank Screen-TRex" is the biobank screening, and the
  progress messages of the two screening selectors follow. The strings of `method_used` ("Screen-TRex
  (ordinary)", "Screen-TRex (bootstrap-CI)") stay: they are results and belong to the naming pass.
- "the R reference" is "the 1.x versions" in trex_screening.cpp, trex_sd.hpp and the two augmented
  elastic-net solvers. The tests and HISTORY keep the old words.
- The help of `TENETAug_Solver` (R) and the header of the C++ solver named a function `trex_GVS()` of
  the old package. The 1.x versions have no such function. Both name the T-Rex+GVS selector of the 1.x
  versions and give the call, `trex(method = "trex+GVS")` (Fabian: the name is T-Rex+GVS, the lower-case
  string is the value of the argument).
- Comments, help texts and progress messages only. RcppExports.R, the help page, the vendored copy and
  the stub follow.

#### The core selector: the help on the three surfaces, the dummy distributions, the positions of R

Step 3 of the documentation comparison for `TRexSelector`, its control and its result, the last of the
selectors. Findings, probes and the four investigations: notes/handover_2026-11-15/doc_step3/
(CORE_findings.md, 29 findings). Fabian, 2026-11-15: "Let us fix the problems with core according to
CORE_findings.md".

- The result fields. `R_mat` is the number of selected variables per stopping time and voting level and
  has one row less than `FDP_hat_mat` and `Phi_mat`. C++ called it a "selection frequency matrix (T_stop
  × v_len)" and Python "empirical tail probabilities (p x T_stop)". The Python texts of `T_stop`
  ("maximum number of LARS steps run"), `fdp_hat_mat` ("p x |voting_grid|"), `phi_mat`, `phi_prime`
  ("smoothed selection probability") and `voting_grid` ("vote totals") are replaced by the content of
  the C++ texts.
- `T_stop` against the selection. Fabian: "T_stop is the stopping time, alias the number of included
  dummies. The T-solvers have run until T_stop dummies each", and "We must explicitly write down the
  divergence of T_stop vs the final choice, otherwise that confuses a user". The paper's Algorithm 2
  raises T while the estimate at the top voting level is at or below the target and then takes the
  largest selection over T' in {1, ..., T - 1}. The 1.x versions and the library do the same, with the
  tie rule "highest voting level, then latest stopping time". In a probe of 1,080 runs the selection
  sat one step before `T_stop` in most runs, for the path solvers and for the greedy solvers under the
  stagnation stop alike, and up to nine steps before it. The help of `T_stop` and of `select()` on the
  three surfaces and guide page 03 say so. A field for the stopping time of the selection is left to
  the naming pass of 2.1.0. Not changed: after an exit at T_max or by the stagnation stop the last step
  is feasible and dropped all the same, as in the 1.x versions (the paper keeps T_max). It held a
  larger selection in 1 of the 1,080 runs. Fabian: "we keep it, but we should set a comment about 1.x
  version mirror". The comment stands at the dropped row in trex.cpp, and the text of `drop_last_row`
  no longer calls the row "the step that violated the stopping criterion" in every case.
- `x_l2_norms` (Fabian: proposal A). The vector holds the divisor the normalizer applied and needs to
  restore X: the L2 norms under L2, the sample standard deviations under ZSCORE, ones under NONE. The
  values stay, the three help texts say what they are, and the name goes to the naming pass.
- The control. `opt_threshold` is the reference voting level of the L-loop ("Optimization threshold"
  before). `tloop_max_stagnant_steps` counts T steps in a row without a growing selection. A value of
  0 acted as 1 and is refused now (Fabian: "yes a"): the C++ constructor throws, `trex_control()` of R
  stops, and Python raises the ValueError of the core. `use_memory_mapping` concerns the dummy
  matrices. The solver lists of the stagnation default
  are complete in C++ and Python (TIENET, TIENET_AUG, TIKNET and the coordinate-descent family were
  missing). HCONCAT: "same FDR in simulations, not the same selection" where the enum said "same
  result" (3 of 20 data sets differed in a probe). The text of `lloop_strategy` no longer calls fresh
  dummies "essential for unbiased FDR calibration".
- Seeds. A seed of 0 or more makes a run deterministic and reproducible, a negative seed draws from
  the hardware entropy source, and a Monte Carlo study gives every trial a seed of its own or uses -1.
  The note of `trex_control()` in R said that a fixed seed "can inflate the empirical FDR". A probe of
  8 batches of 250 trials gave 0.0865, 0.0848 and 0.0841 for one fixed seed, a seed per trial and -1,
  with the same spread between the batches, so the claim is gone.
- Guide page 03, three surfaces. The stagnation stop is described through the selection (it said "the
  FDP estimate has not improved"). Under NONE only an all-zero column is refused (it said "constant
  columns are rejected in every mode"). The paragraph on the result fields explains `T_stop`, and the
  figure text no longer calls the null variables at `T_stop` the false positives of the selection.
- The dummy distributions. `ConstrainedSparseRademacher` has the default s = 0.1 in C++ and Python,
  which only the R glue had (Fabian: "Default of 0.1 is ok for all three"). R gains `p1` and `p2` of
  `"Mammen"`. The R glue takes both defaults from the core. The R help names every default, says that
  `a` is also the lower bound of the triangle and `b` its mode. The SD family is not reachable from R
  or Python: both hold the `ad` level only, and the sparse Rademacher distribution appears there as a
  dummy distribution of the augmented dummies.
- The positions of R (Fabian: "Yes, R counts from 1, we must fix that"). The errors of the data
  normalizer said "column 4 (0-based)" for the fifth column. `.r_positions()` in R/utils.R raises the
  row and column numbers by one, drops "0-based" and spells the scaling mode as the argument. It wraps
  the constructors of the seven selectors, the run of the SPCA and the biobank selector, the two
  power rules and `trex_normalize_design()`. Every other error passes as it is.
- The Python root selector gains a text on each of its 16 attributes, and `selected_indices` is
  annotated as the list it returns. The text of the core class no longer says that y is used without
  a copy. The comments of the wrapper class say what the getters return.
- "the R reference" is "the 1.x versions" in trex.hpp and trex.cpp. The class text lists the derived
  selectors with their short forms and no longer counts the SPCA selector among them.
- Tests. C++: the default of the sparse Rademacher constructor. R: the Mammen probabilities reach the
  dummy draw, the default s, the positions of a constant column, of an all-zero column under NONE and
  of an entry that is not finite, and an error without a position. The scaling test expects the
  positions of R. Python: the default s. The three surfaces test that a stagnation window of 0 steps
  is refused and that 1 runs.
- C++ 751 of 751, R 2,218, Python 798. `R CMD check --as-cran` Status OK, Sphinx 0 warnings, Doxygen 0
  warning lines.

### 2026-11-14

#### The Tikhonov selector accepts TIKNET alone, on the three surfaces

Fabian, 2026-11-14: the deprecated TCDTIKNET is refused now, in C++, R and Python ("I would prefer for
stability that we accept only Tiknet"). Up to here the selector took it with a deprecation warning, and
its removal was announced for the next minor release.

- C++. The constructor of `TRexTikhonovSelector` refuses every `solver_type` other than TIKNET. For
  TCDTIKNET the message adds that the deprecated solver is no longer accepted by the selector.
  `fold_dummy_coupling = false`, the decoupled dummy block that only TCDTIKNET had, is refused as well
  ("is no longer available"). The field stays until the removal of the family. The solver class
  `TCDTIKNET_Solver`, the enumerator and the dispatch stay as they are, deprecated: no selector reaches
  them any more, and the paths for it in the selector (the group-centering hooks, the flag of the
  decoupled block) are dead until the family is removed.
- R. `.check_tikhonov_solver()` accepts "TIKNET" alone and words the refusal of "TCDTIKNET". The plain
  selector said that "TCDTIKNET" needs a penalty matrix and pointed to `TRexTikhonovSelector`. It says
  now that no selector accepts it. Python: the binding no longer raises the deprecation warning for the
  name (the error of C++ comes instead), and the root check of the plain selector words it as R does.
- The documentation of the three surfaces says that TIKNET is the one solver and that `false` for the
  dummy coupling is an error. The sentences on the hooks of TCDIENET under GROUP_CENTERING are gone,
  because that route needed TCDTIKNET. NEWS of R has the change under "Changed" and the exception in
  "Deprecated".
- Tests, C++. The folded-coupling anchors against the GVS selector (group-mean K, singletons, the demo
  control with L capped at p over 80 cases, the six seeds of the L-growth) run TIKNET against
  GVS-TCDIENET and hold. GROUP_CENTERING is compared with SPARSE_GENERAL on the built matrix, and the
  all-singleton case checks that lambda_2 has no effect. The two tests of the lambda_2 guard run the
  default solver.
- The anchor "K = I equals the TCDENET base run" needed the decoupled dummy block and has no
  counterpart: with the FOLDED coupling a dummy shares the penalty of its variable, so K = I is the
  elastic net in the penalty and not in the run. A probe over 40 seeds (notes/handover_2026-11-14/
  tcdtiknet_refused/) shows it: at lambda_2 = 5 the run differs from the base selector with TENET,
  TENET_AUG and TCDENET in every seed, at lambda_2 = 0.5 in none, and at lambda_2 = 0 it equals them and
  TLASSO in all 40. The new anchor is the last one: at lambda_2 = 0 with i.i.d. dummies the selector
  selects what the plain selector selects with TLASSO.
- The test that a `select()` ending with an error gives X back and retires the selector used an
  unconverged TCDTIKNET solve. It runs the same problem through the GVS selector with TCDIENET and one
  group per variable now (C++ in the restore tests, R and Python in the GVS tests).
- Removed without a replacement in R and Python: the two equivalences of the coordinate-descent solver
  (K = I against TCDENET, group-mean K against GVS-TCDIENET). C++ holds the anchors.
- C++ 751 of 751, R 2,198, Python 797.

#### The GVS selectors of R and Python refuse a solver that does not fit, as C++ does

Finding 1 of TIKHONOV_findings.md, second half (Fabian, 2026-11-14: "we fix now"). His three decisions
of the same day: the base control follows the Tikhonov pattern, TLARS is refused, and on the IEN track
the solver has to be named, with no change in C++.

- C++, unchanged. The GVS control nests the base control, whose solver is TENET. On the EN track the
  constructor derives the solver from `en_solver`: the nested value has to be TENET or that solver, a
  conflicting elastic-net solver and every solver outside the family (TLARS, TOMP) throw. On the IEN
  track the nested solver has to be TIENET, TIENET_AUG or TCDIENET, so `gvs_type = IEN` with an
  untouched control throws.
- Before, in R and Python. The base control came as a second argument with the solver TLARS. The
  bindings set the solver themselves: on EN from `en_solver`, on IEN an IEN solver was kept and
  everything else became TIENET. TLARS passed without a word, another solver gave a warning. So the
  bindings changed the control they were given, and `gvs_type = "IEN"` alone ran TIENET.
- Now. R takes `control = NULL` and Python `trex_control = None` by default, and the base control
  nested in the C++ GVS control applies. A control that is passed replaces it as it is. The glue of
  both surfaces (`merge_gvs_control()`) no longer touches the solver, and the C++ constructor derives
  or refuses it. `TRexGVSSelector` and the GVS power rule follow this rule on both surfaces.
- R words the refusal itself, because the C++ messages name `solver_type` and `trex_ctrl`:
  `control$solver must be "TENET" for gvs_type = "EN"` with the solver it got and the way out, and
  `gvs_type = "IEN" needs a base control that names its solver`. Python raises the ValueError of C++.
- For code written for R 2.0.1 and Python 1.0.0: a base control without a solver is an error for the
  GVS selector, and the IEN track needs a base control with TIENET or TIENET_AUG. NEWS says so. The
  selections do not change, because every run that worked before ran the same solver.
- Tests. R: five solvers outside the track on EN and four on IEN are refused for the selector and the
  rule and leave the design untouched, what fits runs without a warning, `control = NULL` gives the
  run of `trex_control(solver = "TENET")`. Python: the same refusals on the two entry points, the
  nested control applies and is not changed, and IEN runs with the solver named on the nested control.
  R tests 2,199, Python 798.
- Vignette "T-Rex GVS" and guide page 05 of Python name the solver in their base controls and no longer
  call TIENET the default of the IEN track. The C++ page needed no change.

#### The Python record of the power rule flags the upper end of the grid, as C++ does

Found while the R help of the power rule was written from the C++ record. A text of the Python
binding, no behaviour.

- `Lambda2PowerRuleResult.at_lower_end` read "The pick is the first grid value (widen the grid before
  trusting it)" and `at_upper_end` "The pick is the last grid value". The hint stood at the wrong
  flag. C++ had one comment for both flags, and since 771a6a7 it has one per member: the pick on the
  last grid value is the boundary solution that asks for a wider grid, and the first grid value is
  usually 0, the answer "no penalty". The two Python texts are the C++ texts now. 6f9cdce had missed
  them.
- The record of the pilot rule (`Lambda2PilotResult`) had the two texts right.
- The stub is regenerated. Python tests 786, Sphinx 0 warnings.

#### The C++ guide page of the Tikhonov selector follows the class

Guide page 08 of C++ ("T-Rex Tikhonov: a general penalty matrix"), texts only. 771a6a7 had corrected the
documentation of the class and had left the page.

- The cross-validation rules no longer "return lambda_2 on the solver scale". The page gains the
  paragraph on the scale: a supplied `lambda_2` is used as it is, a grid is stated for columns of unit
  norm and runs at the grid times n - 1 under `ScalingMode::ZSCORE`, and every rule reports its pick on
  the scale of a supplied value.
- The page names what the class has: `runPilotOnly()`, `getDummyClusters()`,
  `buildGroupCenteringK()`, `gvs_lambda2_power_rule()` and `default_lambda2_grid()` as the default
  grid of `Lambda2PowerRuleConfig`. A selector runs once.
- The section "The solver" said "It is the default of the control". It says now that TIKNET is the
  solver of the nested base control and that the constructor refuses every other solver, so a base
  control that is assigned to `trex_ctrl` has to name it.
- The compiled example and its outputs are unchanged. Doxygen gives the one warning it gave before, a
  missing example output of page 11, which this change does not touch. The three guide pages 08 now
  agree.

#### The texts of the R Tikhonov class mirror the C++ class

Step 3 of the rework of the R class of the Tikhonov selector, the last one of the class. Texts only.
It is the R side of 6f9cdce.

- The class. The help opens with the sentence of C++, "FDR-controlled variable selection with the
  penalty (lambda_2 / 2) * beta^T K beta, whose positive semi-definite matrix K carries the structure
  among the variables", under the title with the short form T-Rex+Tikhonov. It gains the sections of
  the C++ class: the penalty and its special cases, the solver, the dummies, the L-loop strategies, the
  rules for lambda_2 with the scale, and that a selector runs once. Before, it had one paragraph.
- Wrong texts. The class "drove the Tikhonov-net solvers (TIKNET by default, or the deprecated
  TCDTIKNET)". `lambda_2` and `lambda2_used` were "solver scale". `lambda2_used` was "user-supplied or
  CV-resolved", the default is the pick of the pilot rule. `prior_groups` was "for the dummy
  clustering" and "ignored when cluster_dummies = FALSE", the groups also define the penalty under
  GROUP_CENTERING and the count GROUP_MAJORITY. `corr_max` was the "maximum cross-cluster correlation",
  which holds for single linkage only. `seed` was "Random seed (default: -1)".
- Gaps. `hc_linkage` names why Ward is rejected. `pilot_statistic` says what the three counts are. The
  record `lambda2_pilot` and the result of the power rule are lists of their fields with what each
  holds, among them `curve_soft`, `calibrated`, `group_majority`, `groups_qualifying` and
  `probe_seconds`, which the help had named without a word or not at all.
- The vignette "T-Rex Tikhonov" gains the paragraph on the scale of lambda_2 and of the grids and no
  longer says that the cross-validation rules return lambda_2 "on the solver scale". Its title and the
  titles of the other pages wait for the page review.

#### The R Tikhonov class has what the C++ class has

Step 2 of the rework of the R class of the Tikhonov selector (notes: TIKHONOV_findings.md, findings 3,
5, 6 and 7). It is the R side of 5d06c2c.

- `run_pilot_only()` runs the lambda_2 pilot without the selection and returns its record. It follows
  the rule of C++: it needs the pilot rule, and a selector runs once.
- The field `dummy_clusters` reports the clusters the dummies were drawn for, as a list with the names
  of the fields of the R GVS selector (`groups` numbered from 1, `max_clusters`, `hc_method_used`).
- The pilot record gains `dummy_tag`, which C++ and Python had. The record has 21 fields now, one more
  than `Rcpp::List::create()` takes, so the glue fills the list element by element.
- `trex_gvs_lambda2_power_rule()` is the power rule on the GVS selector. The R class of the GVS
  selector held its solver handling and its length check of `groups` inline. Both moved into internal
  functions that the class and the rule share (`.gvs_base_control()`, `.check_gvs_groups_length()`),
  and the join of the two controls in the glue moved into `merge_gvs_control()`. The class behaves as
  before.
- `trex_tikhonov_lambda2_power_rule()` takes a memory-mapped design, which it had refused. The probes
  work on the file through its Eigen::Map, as a selector does: `"readwrite"` mode, or `"readonly"`
  under `scaling_mode = "NONE"`. The GVS rule takes one as well.
- `trex_default_lambda2_grid()` returns the grid that C++ holds once. `pilot_grid` of
  `trex_tikhonov_control()` and `grid` of the two rules default to the call of this function. R had the
  nine numbers at two places.
- `trex_build_group_centering_k()` builds the centering matrix of a group assignment, which C++ and
  Python had.
- `prior_groups` is checked in the words of R, with the check of `groups` in `trex_gvs_control()`
  (`.check_group_numbers()`). Labels that start at 0 had given "prior_groups must be 0-based
  non-negative integers" and a gap "Cluster ID 5 is missing". A wrong length is reported with the name
  of the argument by the selector and by the power rule.
- No result changes. The vignettes "T-Rex Tikhonov" and "T-Rex GVS" and the README name the new
  functions.
- On the way: the R glue compiles as one unit (`rcpp_bindings_unit.cpp`), so an incremental install
  into a scratch library does not see a changed glue file. The test run of step 1 in the scratch
  library had used the old glue for that reason. The clean build of `R CMD check` had the new one and
  passed with the same 2,102 tests.

#### The R Tikhonov class takes TIKNET only, as C++ does

Step 1 of the rework of the R class of the Tikhonov selector (notes: TIKHONOV_findings.md, finding 1).
It is the R side of what e6705c3 did for Python.

- Before. `TRexTikhonovSelector$new()` and `trex_tikhonov_lambda2_power_rule()` took
  `control = trex_control()`, whose solver is `"TLARS"`. The R class let `"TLARS"` pass, warned for every
  other solver outside the family and set it to `"TLARS"`, and the glue (`merge_tikhonov_control()`)
  turned everything outside TIKNET and TCDTIKNET into TIKNET. So R changed the control it was given, and
  C++ refuses such a solver at construction.
- Now. `control` defaults to `NULL`. Without a base control the one nested in the C++ Tikhonov control
  applies, whose solver is TIKNET. A control that is passed replaces the nested one as it is, and a
  solver other than `"TIKNET"` (or the deprecated `"TCDTIKNET"`) is an error. `"TLARS"` is among them,
  as in C++ and Python. The glue no longer touches the solver.
- The error is raised in R, in the words of R: "control$solver must be "TIKNET" for the Tikhonov
  selector. Got: "TLARS". trex_control() names "TLARS" by default, so pass trex_control(solver =
  "TIKNET", ...), or leave control out to run the default control." The C++ message names `solver_type`,
  a field R does not have.
- A change of the published interface of 2.0.1: a call that passes `trex_control(K = 10)` is an error
  until the solver is named. NEWS says so under "Changed". The selections do not change, because every
  run that worked before ran TIKNET.
- Tests. Five solvers outside the family and the default of `trex_control()` are refused and leave the
  design untouched, `control = NULL` gives the run of `trex_control(solver = "TIKNET")` for the selector
  and for the power rule. The other tests of the selector name the solver. R tests 2,102.
- Vignette "T-Rex Tikhonov". Its calls name the solver, and the section "The solver" no longer says that
  `"TLARS"` maps to `"TIKNET"`.

#### The texts of the Python Tikhonov class mirror the C++ class

Step 5 of the rework of the Python class of the Tikhonov selector, the last one. Texts only.

- The class. The core class and the root class open with the sentence of C++, "T-Rex selector with a
  general Tikhonov penalty (T-Rex+Tikhonov): FDR-controlled variable selection with the penalty
  (lambda_2 / 2) * beta^T K beta, whose positive semi-definite matrix K carries the structure among the
  variables". The core class gains the notes of the C++ class (the penalty and its special cases, the
  solver, the dummies, the L-loop strategies, the rules for lambda_2, the scale, a selector runs once),
  and the root class shows them through the decorator. Before, the core class had a summary and the
  references only.
- Wrong texts. `PilotStatistic.CALIBRATED` was the "soft curve" and `curve_soft` a "calibrated (soft)
  count": the soft curve is the mean number of variables ahead of the t-th dummy, and the calibrated
  count is the T-Rex selection rule on the pilot data. `GROUP_MAJORITY` was a count of groups, it counts
  the members of the groups that reach the fraction. `corr_max` bounded the correlation between
  clusters, which holds for single linkage only. `lambda_2` and `lambda2_used` were "on the solver
  scale", they refer to the scale of the standardized columns. The root class named "the pilot power
  rules".
- Gaps. `hc_linkage` names its values and why Ward is refused. `prior_groups` names its three uses (the
  clusters of the dummies, the penalty under GROUP_CENTERING, the count GROUP_MAJORITY). The pilot
  fields say how they are capped. The fields of the pilot record say what they hold.
- Guide page 08 of Python. Its section "The solver" still said that TLARS "maps to" TIKNET, which
  e6705c3 had ended, an omission of that commit. It now says that the nested base control applies and
  that a passed control has to name TIKNET. The page gains the paragraph on the scale of lambda_2 and of
  the grids, and names `run_pilot_only()` and `dummy_clusters`. Its title and the titles of the other
  pages wait for the page review.
- The reference page shows the notes for the root and the core class. Python tests 786, Sphinx 0
  warnings.

#### The Python Tikhonov class has what the C++ class has

Step 4 of the rework of the Python class of the Tikhonov selector (notes: TIKHONOV_findings.md, findings
3 and 5).

- `run_pilot_only()` of the root class (`runPilotOnly()` of the core class) runs the lambda_2 pilot
  without the selection and returns its record. It follows the rule of C++: it needs the pilot rule, and
  a selector runs once.
- `dummy_clusters` of the root class (`getDummyClusters()` of the core class) reports the clusters the
  dummies were drawn for, as the record `DummyClusters` with the field names of the GVS result.
- `gvs_lambda2_power_rule()` is the power rule on the GVS selector, which C++ has beside the one for
  the Tikhonov selector. The join of the two controls of the Python GVS class moved into one function,
  `merge_gvs_control()`, which the class and the rule share. The class behaves as before.
- `default_lambda2_grid()` returns the grid that C++ holds once. The two root functions of the power
  rule take `grid = None` and read it from there, where the Tikhonov one had the numbers a second time.
- The three new names are exported at the package root and in `trex_selector_methods.ad`.
- Tests. Python pins that the pilot alone gives the pilot of a selection, that it needs the pilot rule
  and spends the selector, the record of the clusters before a run, for supplied groups, for found
  clusters and without cluster dummies, and the GVS rule (the grid, the pick, the same seed, a negative
  seed, the control of the caller unchanged). Python tests 786.

#### The Python Tikhonov class takes TIKNET only, as C++ does

Step 3 of the rework of the Python class of the Tikhonov selector. It corrects 74fef7c, which had kept
one mapping.

- Fabian on 74fef7c: "TLARS can only be a solver if the K matrix is zero. I would prefer for stability
  that we accept only Tiknet. How is it in Cpp?" In C++ the selector accepts TIKNET, the default of the
  base control nested in its control, and the deprecated TCDTIKNET, and refuses every other solver,
  TLARS included.
- 74fef7c had let TLARS stand for "no solver chosen" and give TIKNET, because Python takes the base
  control as a second argument and a fresh `TRexControlParameter()` names TLARS. That mapping is gone:
  the binding changes nothing in a control, and the C++ constructor decides. TLARS raises `ValueError`
  as every solver outside the family does.
- The base control of a run. Without `trex_control` the base control nested in the Tikhonov control
  applies (`tik_control.trex_ctrl`), whose solver is TIKNET, as in C++. Up to now a fresh default
  control replaced it, and settings on the nested control had no effect. A `trex_control` that is passed
  replaces the nested one as it is, and so has to name the solver TIKNET. The root class, the core
  class (`trex_control = None` by default, where it was a default control) and
  `tikhonov_lambda2_power_rule()` follow this rule.
- What changes for code written for 1.0.0. A call that passes a base control built from
  `TRexControlParameter()` without a solver raises `ValueError` and needs
  `solver_type = SolverTypeForTRex.TIKNET`. Guide page 08 did so and sets the solver now.
- Tests. Python pins that TLARS and four other solvers are refused on the three entry points, that the
  nested control applies without `trex_control` and is not changed, and that a passed control replaces
  it. Python tests 782.
- The R class follows with its rework.

#### The penalty matrix of the Python Tikhonov class is sparse

Step 2 of the rework of the Python class of the Tikhonov selector (notes: TIKHONOV_findings.md, finding
4, TIKHONOV_python_rework.md).

- The difference. C++ holds the penalty matrix K sparse. R takes a base matrix or a sparse matrix of the
  Matrix package and returns a sparse one from `trex_gamma_to_k()`. Python took and returned dense
  arrays only: `tikhonov_K`, `gamma_to_k()` and `build_group_centering_k()` converted on the way in and
  on the way out, a scipy sparse matrix raised TypeError, and a read of `tikhonov_K` built a dense copy.
  A dense K for p = 20,000 has 3.2 GB.
- The decision (Fabian: "From an efficiency perspective the scipy inclusion would make sense ... Since
  we use in R also Matrix, that should be fair."). scipy becomes a dependency of the package
  (`scipy>=1.11`), and the matrix is sparse on the Python surface. `tikhonov_K` takes a scipy.sparse
  matrix of any format or a dense array, which is converted, and is read back as a
  `scipy.sparse.csc_matrix`. `gamma_to_k()` takes a sparse or a dense operator and returns a
  `csc_matrix`, and `build_group_centering_k()` returns one.
- What changes for code written for 1.0.0. An assignment of a dense array works as before. A read of
  `tikhonov_K` and the results of the two builders are sparse matrices now, where they were arrays:
  `.toarray()` gives the array.
- The wheels. scipy has binary wheels for every target the package is built for: Python 3.10 to 3.15 on
  Linux x86_64, macOS arm64 and Windows x86_64 (checked with pip download, scipy 1.15.3 for 3.10 and
  1.18.1 for 3.12 to 3.15). The test environments of the wheel build install it with the package.
- Texts. The README names SciPy beside NumPy, guide page 08 says that `gamma_to_k()` returns a sparse
  matrix and prints its type, and the docstrings state the types.
- Tests. Python pins the four input formats, the type that is read back and returned, that a sparse
  and a dense matrix give the same run, and that the centering matrix of 20,000 variables stays sparse
  (a million entries). Python tests 780.

#### The Python Tikhonov class refuses a solver outside the family

Step 1 of the rework of the Python class of the Tikhonov selector (Fabian: "next on Python. Afterwards
R.", notes: TIKHONOV_findings.md, finding 1).

- The difference. C++ refuses a solver other than TIKNET and TCDTIKNET at construction. The Python root
  class warned and ran TIKNET, and the Python core class ran TIKNET without a word:
  `merge_tikhonov_control()` in the binding replaced every solver outside the family. A binding does
  not remap a control (Fabian's rule of 2026-11-14, and on the SPCA selector: "we refuse non-fitting
  solvers").
- Why one value still maps. C++ has one control with the base control nested in it, whose solver
  defaults to TIKNET. Python takes the base control as a second argument, and its solver defaults to
  TLARS. A user who builds a base control to set K therefore passes TLARS without having chosen it
  (guide page 08 does so). TLARS stands for "no solver chosen" and gives TIKNET, without a warning.
- The change. Every other solver is passed on as it is, and the C++ constructor refuses it
  ("solver_type must be TCDTIKNET or TIKNET"): `ValueError` from the root class, from the core class
  and from `tikhonov_lambda2_power_rule()`. The root class and the root function no longer carry their
  own warning. X is left alone by a refused control.
- Tests. Python pins the refusal for four solvers on the three entry points, and that TLARS gives the
  run of TIKNET without a warning. Python tests 777.
- The R class has the same pattern and follows with its rework. So has the GVS selector in R and
  Python (a non-matching solver is replaced with a warning), which is not part of this step.

#### The documentation of the Tikhonov selector sits at the class and says what the code does

Step 6 of the C++ class of the Tikhonov selector (notes: TIKHONOV_cpp_findings.md, findings 9 to 20, 22
and 23). Texts and tests only, no change of behaviour.

- The structure. The method was described in the file header, and the class said "see the file
  header", so the generated class page showed a summary and the references. The description moved to
  the class: the penalty and its special cases, the solver, the dummies, the L-loop strategies, the
  rules for lambda_2 and its scale, and the architecture. The file header keeps two lines.
- The name. The class opens with "T-Rex selector with a general Tikhonov penalty (T-Rex+Tikhonov):
  FDR-controlled variable selection with the penalty (lambda_2 / 2) * beta^T K beta, whose positive
  semi-definite matrix K carries the structure among the variables" (Fabian: yes).
- Wrong texts. The class summary named the "CCD Tikhonov-Net solver", where the default is TIKNET, the
  homotopy solver. `lambda_2 == 0` was the "pure TCDLASSO lasso path", it is the lasso path of the
  solver in use. The constructor gave `tFDR` "in (0, 1]" (it is [0, 1]), asked for a non-empty K (not
  under GROUP_CENTERING) and said that X is "restored on destruction" (it is restored when `select()`
  returns). `prior_groups` was "ignored when cluster_dummies == false", while the groups also define
  the penalty under GROUP_CENTERING and the count GROUP_MAJORITY. `corr_max` bounded the correlation
  between clusters, which holds for single linkage only (the sentence of the GVS selector). The scale
  of `lambda_2` had three names ("solver scale", "working column scale", "unit-L2 solver scale"), it is
  the scale of the standardized columns. The cross-validation rules were "(recommended)", where the
  default is the pilot rule. `counts` left out GROUP_MAJORITY, and `getLambda2Used()` the pilot rule.
- References to a workspace that is not public ("in the Tikhonov T-Rex workspace", "demo_14 in the
  Tikhonov T-Rex workspace", "the rule the Tikhonov study adopted") read "in simulations" (Fabian: yes).
- One comment per member. In `Lambda2PilotResult` and `Lambda2PowerRuleResult` several members shared
  one comment, so the generated documentation described the first of them only. Every member has its
  own text, and the two values of `TikPenaltyStructure` as well.
- Tests that were missing for the C++ class: the count CALIBRATED, the four strategies with cluster
  dummies including SKIPL and the permutations with i.i.d. dummies, and the conversion of a
  cross-validation pick under ZSCORE (n - 1 times the pick under L2).
- Doxygen over the folder gives no warning for these files, and the class page shows the notes.
- The R and the Python class follow with their rework, where these texts are mirrored.

#### The Tikhonov selector reports the clusters of its dummies

Step 5 of the C++ class of the Tikhonov selector (notes: TIKHONOV_cpp_findings.md, finding 21).

- The gap. The selector draws its dummies per cluster of correlated variables. The clusters were kept
  in a private member "for diagnostics" and could not be read. Their number stood in the progress
  lines only. The GVS result reports its groups, their number and their origin.
- `getDummyClusters()` returns a `DummyClusters` record with the field names of the GVS result:
  `groups_vec` (the cluster of each variable, 0-based), `max_clusters` and `hc_method_used` ("prior"
  for supplied groups, else the linkage of the clustering). The record is empty before a run and with
  `cluster_dummies = false`.
- `linkageName()` in trex_cluster_hac.hpp holds the names of the linkages once. The GVS selector takes
  them from there, its results do not change.
- R and Python get the record with the rework of their two classes.
- Tests. C++ pins the record for supplied groups, for found clusters (a partition into non-empty
  clusters, and with single linkage no two variables of different clusters correlate above `corr_max`),
  before a run and without cluster dummies.

#### A lambda_2 grid means the same shrinkage under every scaling mode

Step 4 of the C++ class of the Tikhonov selector (notes: TIKHONOV_cpp_findings.md, finding 8).

- The fact. `scaling_mode` of the base control has L2 (default, columns of unit norm), ZSCORE (columns
  of unit variance, squared norm n - 1) and NONE. The fit term is n - 1 times larger under ZSCORE, so
  the same lambda_2 shrinks n - 1 times less there. The cross-validation rules multiply their pick by
  n - 1 for this reason. The grid of the pilot rule and the grid of the power rule were used as they
  were. Probe (n = 100): the pilot counts 13 13 11 10 at the grid values 0, 1, 3, 10 under L2 came two
  decades later under ZSCORE (13 up to 100, then 11 and 10 at 300 and 1000), and with the grid times
  n - 1 the two profiles were equal. Under ZSCORE the default grid therefore reached a shrinkage of
  3000 / (n - 1) only.
- The decision (Fabian asked for a proposal and agreed, also to the scope). A grid is stated for columns
  of unit norm. The pilot rule and the two power rules (`tikhonov_lambda2_power_rule()`,
  `gvs_lambda2_power_rule()`) run at the grid times `lambda2_scale_factor()`, which is n - 1 under
  ZSCORE and 1 otherwise. Both rules are converted, because the pilot rule is the cheap version of the
  power rule and the two would disagree under ZSCORE otherwise. The pick is reported on the working
  scale, as the cross-validation rules report theirs, so it can be supplied as `lambda_2` and reproduces
  the run. A single supplied `lambda_2` is used as it is, as before.
- The records. `Lambda2PilotResult::grid` and the probes of the power rule hold the values the runs
  were made at, and both results gain `grid_scale`, the factor. R and Python carry the field.
- What changes for a user. Nothing under the default L2 and under NONE. Under ZSCORE a selection that
  resolves lambda_2 by the pilot rule, and a power rule, now see the grid they would see under L2.
- Tests. C++ pins that the pilot counts, the pick and the selection under ZSCORE equal those under L2,
  that the pick supplied as `lambda_2` gives the run, the same for the power rule of the Tikhonov and
  of the GVS selector, and the factor. R and Python pin the same on their surface. C++ tests 747 of
  747, R tests 2,084, Python tests 773.

#### The power rule refuses a negative seed and has the default grid of the pilot rule

Step 3 of the C++ class of the Tikhonov selector (notes: TIKHONOV_cpp_findings.md, findings 6 and 7).

- The seed. `tikhonov_lambda2_power_rule()` and `gvs_lambda2_power_rule()` run one selection per grid
  value with the same seed, so that the probes share their dummies and differ by the penalty alone. A
  negative seed draws fresh dummies for every selector, and C++ accepted it: with seed 5 two calls gave
  the counts 10 12 11 11 twice, with seed -1 they gave 11 10 11 12 and 10 10 14 11. Both functions
  refuse a negative seed with std::invalid_argument (`require_shared_seed()`). R and the Python root
  function had refused it already, the Python core function now does through C++.
- The default grid. `Lambda2PowerRuleConfig::grid` was empty, so a default config failed with "the
  lambda_2 grid is empty", while the pilot rule had 0, 1, 3, 10, 30, 100, 300, 1000, 3000 and R and
  Python had written the same numbers into their bindings. `default_lambda2_grid()` in
  tik_lambda2_pick.hpp holds the grid once, and `pilot_grid`, the config of the power rule and the
  default of the Python core binding take it from there. The defaults in the signatures of the R
  function and of the Python root function hold the same numbers until the two classes are reworked,
  and a test on each surface pins that they agree.
- Tests. C++ pins the refusal for both rules, that the default config runs and that the three defaults
  are one grid. Python pins the refusal of the core function and its default grid, R and Python that
  the power rule and the pilot rule agree.

#### The Tikhonov selector checks its control at construction, and its penalty matrix once

Step 2 of the C++ class of the Tikhonov selector (notes: TIKHONOV_cpp_findings.md, findings 2 to 5).

- `lambda_2`. NaN ended in `select()` with std::logic_error "the pilot rule is resolved by
  runLambda2Pilot(), not here", an internal message, and +infinity was accepted and selected nothing.
  The constructor refuses a value that is not finite.
- Late checks. A K that is not finite, not symmetric or not positive semi-definite was refused in
  `select()` by the solver ("TIKNET_Solver: ..."), an unsupported linkage by the clustering helper
  ("cluster_dummies::runCorrelationHAC: ..."), and wrong sizes of the cross-validation by the tuners
  ("ienet_cv_ccd::fit: ..."). The constructor checks them under the name of the selector: the content
  of K, the linkage when the selector clusters the variables itself, and `cv_n_folds` in [2, n] and the
  grid sizes when a cross-validation rule will run.
- The check of K ran in every solver. Each TIKNET solver validated the matrix (finiteness, symmetry, a
  sparse LDLT factorization), and a default selection builds at least 110 solvers (9 grid values times
  10 pilot experiments, and 20 experiments per step of the L-loop). Measured per solver: 1.2 ms for a
  centering matrix at p = 1,000, 97 ms for a dense K at p = 1,000, 661 ms at p = 2,000, so about 73
  seconds per selection there. The copy of the matrix that each solver also makes costs 0.9 ms at
  p = 2,000 and stays. `SolverOptions` gains `tikhonov_validated` (default false): the solver skips its
  check when the caller has checked the matrix. The selector checks K once in its constructor and sets
  the option for its solvers through `ExperimentRunnerConfig` and `SolverConfig`. The centering matrix
  that the selector builds itself is valid by construction and is not checked. R (`tsolver_options()`)
  and Python (`SolverOptions`) carry the option as well, so the solver options stay the same on the
  three surfaces.
- A supplied K under GROUP_CENTERING. The header said that it "must equal" the centering matrix, and
  nothing compared it. With K = I supplied the solver ran on the groups while a cross-validation rule
  tuned on K = I (lambda_2 0.39 against 1.37 with the matrix left out). The constructor compares the
  supplied matrix with the centering matrix of the groups and refuses another one (Fabian asked for a
  proposal, this is it). The right matrix and no matrix give the same run.
- Texts. The C++ header states the checks at `tikhonov_K`, `lambda_2`, `cv_n_folds`, `hc_linkage` and in
  the list of exceptions of the constructor. The R help and the Python docstrings of the control follow
  for the sentences that changed. Their rework comes with the two classes.
- Tests. C++ pins every refusal at construction with the name of the selector, that unused settings are
  left alone (a linkage with supplied groups, the sizes with a fixed lambda_2), the comparison under
  GROUP_CENTERING and the solver option. R and Python pin the refusals and the option on their
  surface. C++ tests 742 of 742, R tests 2,069, Python tests 770. The pilot probe at lambda_2 = 0 with
  10 experiments and a dense K at p = 600 takes 32 ms where it took 265 ms (23 ms with the identity).

#### A selector runs once: a second select() is an error

- The defect. `select()` gives X back when it ends (`denormalizeAndReleaseX()`). The guard at the start of
  `select()` read a flag that only `restore()` and a failed `select()` set, so a second `select()` on the
  same object ran on the design at its original scale, against dummies of unit norm, and returned a wrong
  selection without an error. Probe: `TRexSelector` selects the 10 true variables of 60 in the first run
  and all 60 in the second, in C++, R and Python. The installed R 2.0.0 does the same, so the published
  R 2.0.1, Python 1.0.0 and C++ 1.0.0 have it. Found while reading the Tikhonov selector, whose
  `runPilotOnly()` has the same hole: `runPilotOnly()` followed by `select()` selected 16 variables at
  lambda_2 = 3000 where a fresh `select()` selects 12 at lambda_2 = 1.
- The fix. `beginRun()` replaces `ensureNotReleased()`. It refuses a selector that gave X back or that
  has run, and marks the selector as run. The three places where a run begins call it: `select()` of the
  base class, `select()` of the screening selector, and `runPilotOnly()` of the Tikhonov selector. A
  second call throws std::runtime_error ("this selector has run already"), under every scaling mode. The
  results of the run stay readable, X is left alone, and a new selector on the same matrix repeats the
  run. The error after `restore()` keeps its text.
- `runPilotOnly()` and `select()` exclude each other in both orders, as the header had said without
  enforcing it. A control without the pilot rule is refused before the run begins, so that selector can
  still select.
- A target FDR that is NaN was accepted by the base class (the check `tFDR < 0 || tFDR > 1` lets NaN
  pass) and selected nothing. The check has the negated form now, as in the SPCA selector.
- Texts. The header, the R help and the Python docstrings of `select()` state the rule, for the base
  class and with one sentence for the derived selectors. Guide page 12 says it on the three surfaces.
- Fabian asked for a proposal (2026-11-14, "What do you propose?"). The alternative, a second `select()`
  that standardizes again and reruns, was not taken: the class holds state of the first run (warm
  starts, the dummy store, memory-mapped files, the results), it was not built for reruns, and the error
  matches the rule that already holds after `restore()`.

#### The SPCA selector checks its bounds at construction, and its help agrees on the three surfaces

- The bounds in C++. The constructor of `TRexSPCA` checked `M > min(n, p)` only. `M = 0` failed in
  `select()` with the message of the SVD and `M = -1` with `std::bad_alloc`. The target FDR was checked by
  the inner selector at the first component, after the PCA. The constructor now refuses `M < 1`, a target
  outside [0, 1] (NaN included) and the four rules of the informed elastic net for `lambda2_method`,
  which the inner selector had refused late and under its own name.
- The range of `tFDR`. R and the Python root class took (0, 1), while C++ and every other selector take
  [0, 1]. They take [0, 1] now. The Python root class no longer checks `M` and `tFDR` itself, the core
  does. R checks the bounds of the core in `$new()`, because its C++ object is built in `select()`, so
  an `M` above `min(n, p)` is refused at construction as well.
- `verbose` of the Python root class. It had no such argument and always ran quiet. It gains
  `verbose = False`, the default of C++ and R. The default stays false on the three surfaces (Fabian).
- Wrong texts. `cumulative_ev` was a "percentage" in R and Python, and it is a fraction of the total
  variance of the centered X. Guide page 07 of R and Python called `adjusted_ev` a share, and it is a
  variance. C++ numbered the two modes "methods 5, 7" and "methods 6, 8". The paper has no numbers, it
  names the T-Rex PCA and the T-Rex Thresholded PCA, and 7 and 8 were variants with the informed
  elastic net, which the constructor refuses. C++ said that `gvs_type` is "forced to EN", where it is
  refused. A supplied `lambda_2` was "in LARS units" on the three surfaces, and it refers to the scale
  of the standardized columns, as corrected for the GVS selector. The Python root class said that X is
  "stored by reference", which holds for a Fortran-ordered float64 array only. The C++ class text
  computed the PCA per component, and it is computed once. The C++ guide page said that
  `lambda_2 = 0` gives a warning, which only the R class and the Python root class give.
- Gaps closed on the three surfaces. A component without a selected variable has a zero column in `V`
  and `Z` and the explained variance 0. In Thresholded mode the non-zero loadings can be other variables
  than `active_sets`, for which the FDR control holds (C++ said so, R and Python did not).
  `adjusted_ev` carries its definition, and `cumulative_ev` the note that it is not the PEV of the
  paper. The class documentation of C++ and Python names the L-loop strategies. Python says that the
  nested `gvs_ctrl.en_solver` and `gvs_ctrl.trex_ctrl.solver_type` are replaced by `en_solver` of the
  control. TENET_AUG is named as the form of the paper and TCDENET as deprecated.
- `cumulative_ev` against the PEV of the paper (Fabian's question). The numerator is the same, the
  adjusted explained variance. The paper divides by the signal EV plus the mixed EV, which needs the
  true support and exceeds 100 % when a method puts variance on null loadings. The library divides by
  the total variance of the centered X, which exists for any data and keeps the value between 0 and 1.
- The name. The class opens with "T-Rex sparse PCA selector (T-Rex+SPCA): sparse principal components
  whose non-zero loadings are selected with FDR control. Each ordinary principal component is the
  response of a T-Rex+GVS run with the elastic net." The short form is Fabian's. The titles of the guide
  pages keep "T-Rex SPCA" until the page review.
- The control of R. `trex_spca_control()` reaches nine parameters, while C++ and Python reach the whole
  GVS control and the base control. Fabian: the missing reach must be repaired, and for 2.0.2 the help
  says that every other parameter of the inner runs stays at its default. The repair belongs to the
  naming pass.
- Tests. C++ pins the refusals at construction (`M`, `tFDR`, the rules, the permutation strategies),
  that the four strategies run, that the nested solver fields have no effect, that `cumulative_ev` is
  the running sum of `adjusted_ev` over the total variance, that Thresholded mode keeps as many
  loadings as variables were selected, and the rule for a component without a selection. R and Python
  pin the same statements on their surface, Python also `verbose` and the copy of an array that is not
  Fortran-ordered. C++ tests 733 of 733, R tests 2,041, Python tests 760.

#### The biobank screening has one parametrization on the three surfaces, and its help is corrected

- The default of `verbose`. R printed by default (`verbose = TRUE`), C++ and Python do not. R takes the
  default of C++, FALSE.
- The target FDR. C++ and Python hold it in the control as `target_FDR_trex`. R took it only as `tFDR`
  in the constructor, and `trex_biobank_control()` had the two bounds. The control of R gains
  `target_FDR_trex` (default 0.1), placed after the bounds so that calls by position keep their
  meaning. `tFDR` of the R constructor now defaults to NULL and, when it is supplied, replaces the
  target of the control, so code written for 2.0.1 keeps running. Whether `tFDR` leaves is a matter of
  the naming pass.
- The screening runs once per phenotype. The C++ class note said that it runs twice, once per variant.
  Both FDR estimates come from one set of K experiments, and the bootstrap variant is derived from the
  coefficients of that run.
- The decision rule is written down on the three surfaces. A variant is accepted when its estimate lies
  in [lower_bound_FDR, upper_bound_FDR]. Of two accepted variants the one with the larger estimate is
  taken, and the bootstrap variant on a tie. With no accepted variant the T-Rex selector runs at
  `target_FDR_trex`. C++ had one line on it, R and Python none.
- `use_bootstrap_CI` of the nested screening control has no effect in the biobank screening, because
  the selector computes both variants for every phenotype. The three surfaces say so.
- Texts. Python called the two bounds the bounds "on the FDR search range", they bound the accepted
  range of the estimate. The three values of `method_used` are named in C++ and R as well. The class
  opens with "Biobank screening with the T-Rex selector: per phenotype the fast screening
  (T-Rex+Screening), and the T-Rex selector as fallback when no FDR estimate of the screening lies in
  the accepted range". The R vignette on the screening sets the target through the control.
- The R help example of `TRexBiobankScreeningSelector` printed `NULL`. `Gauss_data$y` is a matrix with
  one column, which the selector takes as a list of phenotypes, so `res$selected_indices` did not
  exist. The example passes `drop(y)` and prints the selection, and it sets the target through the
  control.
- Tests. C++, R and Python pin the decision rule over six acceptance windows, with the expected method
  computed from the two reported estimates, and that `use_bootstrap_CI` changes nothing. R pins the
  defaults, the target from the control, the override by `tFDR` and a control list of 2.0.1 without
  the target. C++ tests 725 of 725, R tests 2,007, Python tests 738.

#### The screening takes the same strategies on the three surfaces, and its help is corrected

- Rule of Fabian, 2026-11-14: "we aim for standards, so parametrizations between Cpp, Python, and R are
  not allowed to vary, except it is not possible in the respective language. Cpp serves as the reference
  for the texts." A binding passes the control on and lets the core decide.
- The strategies of the screening. The core takes STANDARD (default), SEEDED, PERMUTATION and
  PERMUTATION_SEEDED and refuses HCONCAT and SKIPL, because the screening is one pass at L = p and has
  no L-loop to grow (confirmed by Fabian the same day: "HCONCAT and SKIPL make no sense there.
  PERMUTATION for Screening is valid. default is Standard"). R overwrote the strategy of the control
  with STANDARD in `TRexScreeningSelector` and in `TRexBiobankScreeningSelector`, so SEEDED and the two
  permutations could not be chosen in R, and a chosen one was replaced without a word. The Python root
  class ran HCONCAT as STANDARD without a word and SKIPL with a warning, and wrote STANDARD into the
  control object of the caller. Both now pass the control on. R runs the four strategies, and R and
  Python refuse HCONCAT and SKIPL with the message of the core.
- The variants of the screening are filters. The C++ texts of `ScreenTRexMethod` called them
  "Dependency-Aware T-Rex (...)". Fabian: "Screen TRex does not call Dependency-Aware T-Rex it uses
  dependency aware politics or filters." The texts now say "Screening with the dependency-aware filter
  for ...", and R and Python carry the C++ texts.
- The fields of the result have one text. In R `estimated_correlation` was the "Covariance max value",
  `beta_mat` the "Complete variables parameter space", `dummy_betas` the "Dummy variables parameter
  space", `confidence_level` "Bounds the bootstrap estimates" and `used_bootstrap` "Employs
  bootstrapping". Python gave `beta_mat` the shape p x T_stop, it has one value per variable. Python
  said that the bootstrap serves the FDR estimation. It makes the selection, and the estimate stays one
  over the number of selected variables.
- The values of the base class after a screening run are documented: T_stop 1, L 1, p dummies, the
  voting threshold 0.5, no target FDR (0), a 1 x 1 matrix with the number of selected variables, and a
  1 x 1 matrix with 0 for the FDP estimates. The screening estimates no FDP per voting level, and
  `screen_trex()` of the 1.x versions of the R package returns none either, so the 0 stays.
- The AR(1) window carries the sentence that was corrected for the dependency-aware selector. The class
  opens with "T-Rex screening selector (T-Rex+Screening): fast FDR-controlled variable screening in one
  pass of the random experiments, until one dummy is included (T = 1), with as many dummies as
  variables (L = p)". The class documentation of R and Python names the strategies, as the C++ one did.
- Tests. R and Python pin that the strategy of the control reaches the selector (the permutation runs
  differ from the default run, each strategy equals its seeded twin), that HCONCAT and SKIPL are
  refused, and that Python leaves the control of the caller as it is. R pins the same for the biobank
  screening. C++, R and Python pin the fixed values of the base class and the fields of the result.
  C++ tests 724 of 724, R tests 1,986, Python tests 731.

#### The help of the GVS selector agrees on the three surfaces and with the code

- The help texts of the selector, of its control and of its result were read side by side in C++, R
  and Python and checked against the code, as for the dependency-aware selector before.
- R: the error texts for `groups` spoke the language of the core. Labels that start at 0 gave
  "prior_groups must be 0-based non-negative integers", a gap gave "Cluster ID 1 is missing" for the
  group number 2 of R, and a wrong length named `prior_groups`. `trex_gvs_control()` now checks the
  numbering itself ("groups must number the groups 1 to M without gaps"), and the selector checks the
  length with the name of the R argument.
- `corr_max` was called the largest correlation between variables of different clusters. That holds
  for Single linkage only. On an autoregressive design with `corr_max = 0.5`, Complete, Average and
  WPGMA left correlations up to 0.82 between clusters. The texts now say what the cut at
  1 - corr_max means for Single and for Complete linkage.
- A supplied `lambda_2` was described as "in LARS units" in R and Python. It is used as it is, on the
  scale of the standardized columns, so the same shrinkage needs the value times n - 1 under ZSCORE.
  The three surfaces say so. The C++ test `ScalingModeParity_IEN_FixedLambda2` has pinned this fact
  since it was written.
- `gvs_type` was an "augmentation policy" in R, a "selection policy" in C++ and a "method variant" in
  Python, while all augmentation lives in the solvers. Fabian, 2026-11-14: it is the "penalty of the
  grouped selection: EN (elastic net) or IEN (informed elastic net)".
- The result fields. `lambda2_used` named a "GCV" that the selector does not have, the rules are k-fold
  cross-validations. `max_clusters` counts supplied groups as well. `hc_method_used` is "prior" for
  supplied groups and otherwise the linkage in lower case, which the texts now state. Python called
  the group names "unique group label integers".
- R: `hc_linkage` takes the linkage in lower case, "single", "complete", "average" or "wpgma", in
  `trex_da_control()`, `trex_gvs_control()` and `trex_tikhonov_control()`, as `hc_method_used` returns
  it (Fabian, 2026-11-14: "make them also lower case. For ver 2.0.2"). The controls store the lower-case
  name, and the capitalized names of 2.0.1 are still accepted, so code written for 2.0.1 keeps running.
  One function of the R bindings reads the name for the three selectors. Python and C++ take the
  enumerators of `LinkageMethod` and do not change.
- The control. `tenet_aug_use_lars` also acts on the IEN track with TIENET_AUG, which R and Python did
  not say. `cv_n_lambda` serves the rules CV_*_SVD, CV_*_CCD and CV_*_TIK_SVD, and Python had limited
  it to the "EN-shaped rules". The numbering of the groups without gaps, the users of `hc_linkage`
  and the meaning of a negative `seed` are stated on every surface. Python's `GVSType.EN` no longer
  says that it "uses the TENET solver", because `en_solver` chooses it.
- The class documentation of the three surfaces says which L-loop strategies the selector takes:
  STANDARD, HCONCAT, SKIPL and SEEDED, with SEEDED equal to STANDARD, while PERMUTATION and
  PERMUTATION_SEEDED are refused. Only the C++ file header had said so.
- Names and repairs. The short form is T-Rex+GVS, and the class opens with "T-Rex selector for grouped
  variable selection (T-Rex+GVS): FDR-controlled variable selection for groups of highly correlated
  variables". "The R reference" and its relatives are "the 1.x versions of the R package". The C++
  header lost two sentences on rules that no longer exist, a "Both" after three solver variants and
  three typing errors. Guide page 05 states the numbering of the groups.
- Tests. C++, R and Python pin the values of `hc_method_used`, the number and the names of supplied
  groups and the supplied `lambda_2`. R and Python pin the refusal of the permutation strategies. R
  pins the new error texts and that both spellings of a linkage give the same run for the three
  selectors, also for a control list built by hand. C++ tests 723 of 723, R tests 1,956, Python tests
  727.

#### The Python reference builds without the warning on the penalty formula

- The Sphinx build reported "Undefined substitution referenced" for the penalty formula of the
  dependency-aware selector, `|Phi_T(j) - Phi_T(j')|`, once for the core class and once for the
  selector of the package root. reStructuredText reads a text between two bars as a substitution. The
  hook of `docs/conf.py` protected a single name between bars (`|cor|`, `|A|`) and no expression.
- The hook now protects any text between two bars that has no blank right inside them. The bars of a
  type union such as `A | None` have blanks and stay as they are. The build ends with 0 warnings, and
  the formula shows its bars on both classes.

#### The disk-space tests no longer depend on other processes in the temp folder

- `DiskSpace.TempDirManagerChecksBeforeCreating` failed once in a parallel run of the C++ tests. It
  counted the entries of the temp folder that start with `trex_disk_space_test_` before and after a
  refused create, and two other cases of the file wrote files with that prefix into the same folder.
  ctest runs every case as its own process, so the count could change in between.
  `DiskSpace.MemmapManagerChecksWorstCaseBeforeCreating` counted the `trex_memmap_` folders the same
  way, which every test with a memory-mapped dummy store creates.
- Fabian, 2026-11-14: "I don't accept failing tests, either we can fix it, or we rethink the test."
- Every test of the file that writes now works in a folder of its own with a random name, created and
  removed by a small guard class. For the two manager tests the folder is the base directory of
  `TempDirManager` for the time of the test. The checks are stricter than before: the folder must be
  empty after a refused create. One new test pins that the guard hands the base directory back.
- Evidence. With the old executable the first test failed 5 times in 300 raced runs and the second 137
  times in 300. With the new one both fail 0 times in 300, four executables at once fail 0 times in 600,
  and nothing is left in the temp folder. The full suite passes five times in a row, C++ tests 722 of
  722. No other test of the C++, R or Python trees lists the shared temp folder. No library code
  changes.

#### The reference element of the grid rounds a half up, and the documentation says so

- The reference threshold of BT, NN and the prior groups is element round(0.75 * L) of a grid of length
  L. The library rounds a half up (`std::round`). The 1.x versions of the R package compute
  `round(0.75 * rho_grid_len)` in R, which rounds a half to the even number. The two differ for
  L = 6, 14, 22, 30 and so on, where the 1.x versions take one element lower. The default length 20 is
  not affected.
- Decision of Fabian, 2026-11-14: the library keeps its rule, because the differences are minor, and
  states it. The class documentation of the three surfaces carries one sentence, and a comment stands at
  the three places in the code. No code changes.

#### The selectors of the Python package root carry their documentation

- `trexselector.TRexDASelector`, the class of the guides, showed one line under `help()`. The notes on
  the penalty and the reference list sat on the core class in `trexselector.trex_selector_methods.ad`,
  the binding of the C++ class with its getter methods. All seven selectors were built this way, and
  most of their texts opened with "Python wrapper for the C++ ...".
- Each selector of the package root now opens with what it is, lists the parameters of its constructor
  and ends with the notes and the reference list of its core class. A decorator takes these sections
  from the binding docstring at import, so the rendered reference lists keep one source. `select()` of
  the dependency-aware, the grouped and the screening selector names the result class it returns. The
  parameters of the SPCA selector moved from its constructor to the class text.
- The reference page of the selectors documents the seven classes of the package root first, then the
  controls, results and enums with the core classes. It had said that the classes of the core module
  are "re-exported at the package root", which does not hold for the selectors.
- Left for the comparison of each method: the attributes of the selectors (`T_stop`, `L`, `v_thresh`
  and the others) have no texts yet, and the opening lines are those of the core classes.
- Tests. A new Python test file pins for the seven selectors that the text lists every constructor
  parameter and ends with the sections of the core class. Python tests 720.

#### The FDR limitation of BT, NN and the prior groups stands in the help of the three surfaces

- The file header of the dependency-aware selector in C++ carried a paragraph that R and Python lacked.
  At a very high correlation within the groups (about 0.9 and above) and with no independent predictors,
  the equiangular solver TLARS can exceed the target FDR with BT, NN and the prior groups, because the
  true signal and its collinear shadows become statistically exchangeable. The default rule
  FeasibleOnly and a greedy solver mitigate it.
- The paragraph is now part of the class documentation in C++, a section of the R help of
  `TRexDASelector` and a note of the Python class, and the C++ file header points to the class. It names
  the 1.x versions of the R package where it said "the R reference", so the files of the
  dependency-aware selector no longer use that name. No code changes.

#### The help of the dependency-aware selector agrees on the three surfaces and with the code

- The help texts of the selector, of its control and of its result were read side by side in C++, R and
  Python and checked against the code. One defect in the R binding, one late error and the texts below
  came out of it.
- R: the field `da_method_used` returned "UNKNOWN" for a run with prior groups, because the getter had
  no case for the method. It returns "PRIOR_GROUPS", as the result of C++ and Python does.
- C++ and Python: the method PRIOR_GROUPS without `prior_groups` passed the constructor and failed in
  `select()`. The constructor rejects it. The selector sets the method itself when groups are supplied,
  and the texts of the method say so.
- The reference threshold. BT, NN and the prior groups evaluate every step on the whole grid of
  correlation thresholds. The L-loop and the T-loop read the FDP estimate at one of them, the element at
  three quarters of the grid (element 15 of the default 20), and the selection searches all thresholds.
  The matrix of FDP estimates, the matrix of relative occurrences and Phi' of the selector therefore
  belong to the reference threshold, while the selection counts, the voting threshold and the selection
  belong to `rho_thresh`. No surface said so. The class documentation of the three surfaces and guide
  page 04 state it.
- The result fields have one text on the three surfaces. `rho_thresh` is -2 for AR1 and EQUI and when no
  cell meets the target. `rho_grid` says where its values come from. `cor_coef` says that the methods
  without a coefficient return the value of the control. The R text of `cor_coef` still spoke of a
  "dependency-aware dummy generation" and the R text of `da_method_used` of an "integer code". The
  Python texts called the three arrays "BT method only", and they are filled for BT, NN and the prior
  groups. The C++ header gave `R_array_BT` one row more than it has.
- The control. `hc_grid_length` is the number of thresholds on the grid of BT, NN and the prior groups
  (Python had called it the grid of the hierarchical clustering). `hc_linkage` names its users, BT and
  the prior groups. `bt_selection_mode` acts for BT, NN and the prior groups (Python had named BT).
  `cor_coef` says in R and Python how the coefficient is estimated, and `seed` what a negative value
  does.
- The AR1 window. The texts said that the window holds the lags whose correlation reaches
  `rho_thr_DA`. The half-width is ceil(log(rho_thr_DA) / log(|cor_coef|)), so the last lag of the
  window lies at or below the threshold (0.7 and 0.02 give 11 lags, and 0.7^11 = 0.0198). The help
  texts and guide page 04 give the formula.
- Names. The short form is T-Rex+DA in the help of the three surfaces, where "TRexDA", "DA-TRex" and
  "Dependency Awareness" stood. "The R reference" is now "the 1.x versions of the R package". The C++
  header no longer points to a note outside the repository, and its class details no longer say that
  `select()` is overridden. It is inherited, and the selector acts through the hooks of the base class.
- Tests. C++, R and Python each pin the reference threshold (the matrices of the selector against the
  arrays of the result), the method name of a run with prior groups, and -2 with an empty grid for AR1
  and EQUI. C++ and Python pin the rejection of PRIOR_GROUPS without groups. C++ tests 721 of 721, R
  tests 1,916, Python tests 705.

### 2026-11-13

#### Guide page 05 describes the informed elastic net as its paper does

- The page said that the informed elastic net encodes the groups in the penalty "so that the members of
  a group share a coefficient". That describes a centering penalty. The informed elastic net penalizes
  the squared sum of the coefficients of a group, divided by the group size. The page now says what the
  paper shows: it keeps the grouping effect of the elastic net at a lower computation time, and it is
  the elastic net when every variable is its own group.

#### The guides no longer show the coordinate-descent family

- The coordinate-descent solver family is deprecated and leaves with the next minor releases. The guide
  pages of the three surfaces no longer present it. Pages 05, 08 and 09 lose the examples that ran its
  solvers and the passages that described them: the coordinate-descent choices of the grouped selector,
  the second solver of the Tikhonov selector with its sweep budget, the decoupled dummy block, and the
  row and the options of the family on the T-solver page. Pages 01 and 03 lose the mention of the family,
  and page 01 no longer counts the solvers.
- Page 08 runs the default solver TIKNET throughout, and the T-solver page shows the Tikhonov net with
  TIKNET_Solver.
- The help of the classes, of the solver names and of the options keeps the deprecation notice, and so
  do the READMEs. No code changes.

#### The stand-alone solvers take `lambda2 = 0.1` and `rho_afs = 0.3` by default

- The solver options had `lambda2 = 0` and `rho_afs = 1` as defaults, while a selector hands 0.1 and 0.3
  to the solvers of its random experiments. The two parameters now have one default, 0.1 and 0.3, in
  `SolverOptions` of C++ and Python and in `tsolver_options()` of R.
- What changes for a stand-alone solver that is built without these options. An elastic-net solver
  (TENET, TENETAug, TIENET, TIENETAug, TIKNET and the deprecated coordinate-descent ones) runs with a
  penalty of 0.1 and no longer equals its lasso, and TAFS runs with a shrinkage step of 0.3 and no
  longer as forward stepwise. `lambda2 = 0` and `rho_afs = 1` give the earlier behaviour. The selectors
  set both values themselves and do not change.
- The help texts of the two options on the three surfaces state the default, and the table of guide
  page 09 does.
- Tests. Four C++ tests of "lambda2 = 0 reduces to the lasso" relied on the old default and now set it.
  One Python test pinned the old default. A new test on each surface pins the shared defaults of the two
  parameters and the two tolerances. C++ tests 720 of 720, R tests 1,895, Python tests 699.

#### The two tolerances of the solvers differ on purpose

- A solver has the tolerance `tol` of its internal comparisons, machine epsilon by default. A selector has
  its own `tol`, 1e-6 by default. It is the threshold of the selector's own comparisons, and the selector
  sets it on the solvers of its random experiments. The two serve different purposes, were looked at
  before and stay as they are.
- The documentation of both fields in C++, the help of `trex_control()` and `tsolver_options()` in R and
  the Python texts say which default applies where. No code changes.

#### The R vignettes say why they copy the design with `X + 0`

- The vignettes 03, 09 and 12 keep the design of a run in `X_before <- X + 0`. A plain assignment
  shares the buffer with `X` and would show the working design. A comment beside each of the three
  lines says so.

#### The documentation of the Tikhonov selector says which penalty is which

- The penalty of the informed elastic net is the group-mean matrix, K = sum_m 1_m 1_m^T / p_m. The header
  and the R help of the control said so. Guide page 08 on the three surfaces called the block-diagonal
  centering matrix I - J / p_m "the penalty of the informed elastic net". The page now says that a
  group-mean matrix reproduces that penalty and that a centering matrix asks the members of a group to
  share a coefficient.
- `GROUP_CENTERING` needs `prior_groups`, and the constructor refuses it without them. Page 08, the R
  help and the Python text of the structure said that the groups come "from `prior_groups` or from the
  clustering".
- Page 08 wrote the penalty as lambda_2 ||Gamma beta||^2. The solvers use (lambda_2 / 2) beta^T K beta,
  and the elastic net (lambda_2 / 2) ||beta||^2. Both formulas carry the factor now, so a `lambda_2` of
  the page means what the code does.
- "K = I is the elastic net" and "K = I collapses to the plain elastic net" stood on page 08, in the R
  help and in the header. K = I gives the ridge penalty of the elastic net. The run of the selector
  equals the elastic-net selector only with the decoupled dummy block, which the deprecated TCDTIKNET
  alone has. The three texts name the penalty, and the example of page 08 that ran TCDTIKNET with the
  decoupled dummy block leaves the page on the three surfaces. The page shows the default solver TIKNET.
- The R help called `tikhonov_K` "ignored under GROUP_CENTERING". It is optional there. Left out, the
  selector builds the centering matrix from `prior_groups`, and a supplied matrix is used as it is by
  TIKNET and by the cross-validation tuners.
- Two Python texts of the control and of the pilot result. `L_0` is the number of dummies of the pilot,
  `pilot_L_multiplier` times p, and was called a multiplier. `corr_max` is the maximum correlation
  between variables of different clusters and was called a within-cluster correlation.
- No code changes.

#### The documentation says what each dependency-aware method is made for

- Each method of the dependency-aware selector stands for the dependency model of its paper. The binary
  tree follows the hierarchical graphical model of the T-Rex+DA paper. It is made for disjoint groups of
  highly correlated variables and was studied with one active variable per group. The autoregressive
  method is made for variables whose correlation decays with their distance in the column order, as in a
  first-order autoregressive process, and does not assume one active variable per group. The
  nearest-neighbour method is made for overlapping groups of highly correlated variables, as they are
  characteristic for stock returns. The equi-correlated method is made for variables that all share one
  correlation.
- Guide page 04 says so in a paragraph of "The dependency methods" on the three surfaces, and so do the
  help of `da_method` in R and the texts of `DAMethod` in C++ and Python.

#### The texts of the dummy-multiplier cap and of the stopping-time bound say what they hold

- The Python text of `use_max_T_stop` said "always run solvers to T_stop regardless of early stopping",
  which is not what the option does. With the option on, the default, the T-loop ends at ceiling(n / 2)
  at the latest, and without it the loop can run up to the number of dummies. The Python text of
  `getMaxNumDummies` promised "the maximum number of dummy variables per L level", and the method
  returns the cap of the dummy multiplier.
- The R fields `max_num_dummies` and `max_T_stop` carry names that mislead. The first holds that cap, a
  multiplier, and the second holds whether the bound is on, a logical. Their help says so now, and so
  does the help of `stagnation_check`. The names change with the naming pass of the next minor releases,
  and the Python selector classes gain the three fields then.
- The R help of `max_dummy_multiplier` and `use_max_T_stop`, the C++ documentation of the two getters
  and of the field, and guide page 03 on the three surfaces name the bound. No code changes.

#### The guides say that the design comes back up to rounding

- A selector and a writing T-solver center and scale the design in place and undo it afterwards. The
  undo is not bit-exact, because centering, scaling and their inverse are floating-point operations. The
  difference is of the order of 1e-16 (2.22e-16 in the T-solver example of the C++ guide), so a
  comparison for identity fails while a comparison with a tolerance holds.
- The guide pages print that difference and called it "the original values are back" (pages 03 and 12),
  "the file holds the original values again" (page 11) and "the original values" (page 09). The comments
  and labels now read "restored up to rounding", in the R vignettes, the Python guide pages and the two
  C++ examples. Page 12, the page about the work in place, says it in its text as well. Exact bits would
  need a copy of the design, which no binding makes.

#### The PERMUTATION strategy permutes the dummies

- The strategy tables of the C++, R and Python READMEs and of guide page 03 on the three surfaces
  described `PERMUTATION` as "permuted columns of the design" and its storage as K matrices. The
  strategy draws one base dummy matrix and permutes its rows per experiment, and it stores that one
  matrix. `PERMUTATION_SEEDED` re-derives the base matrix and its permutations from the seed and stores
  nothing. The six tables say so. The code and its header were right.

#### The R controls gain `bt_selection_mode` and `group_labels`

- A comparison of the eight controls on the three surfaces, by parameter and default. Every field of the
  C++ control structs is a Python attribute and the other way round, and R and Python agree on every
  default. R lacked two parameters. `trex_da_control()` gains `bt_selection_mode` ("FeasibleOnly", the
  default, or "RFaithful"), the policy for the cells of the calibration grid that may win the selection
  of BT, NN and prior groups. `trex_gvs_control()` gains `group_labels`, the names of the groups of
  `groups`, which the selector returns in its field `group_labels`.
- The sparse principal component selector of R keeps its one control with the parameters of the wrapper.
  It does not take the controls of the grouped selector and of the T-Rex selector.
- Tests in R: the default, both values, an unknown value at the control and at the parser, a run with
  each policy, the names of the groups returned by the selector and a wrong number of names.

#### The documentation of the dependency-aware selector drops the "dependency-aware dummies"

- There are no dependency-aware dummies. Every method of the dependency-aware selector draws ordinary
  independent dummies, and the dependency enters through the penalty on the relative occurrences. The
  phrase stood in the title of guide page 04 on the three surfaces, in its first heading ("Why dependent
  variables need their own dummies"), in the selector tables of the C++, R and Python READMEs and in the
  description of `trex_da_control()` ("dependency-aware dummy generation"). The page is now "T-Rex+DA:
  dependent variables", the heading "Why dependent variables need a correction", the table row
  "dependency-aware (DA) selection of dependent variables", and the control builds "the dependency-aware
  correction".
- The R help of `hc_linkage` offered "Ward", the R parser accepted it and the constructor rejected it.
  Ward's update is valid for squared Euclidean distances only, and the correlation distance 1 - |cor| is
  not Euclidean. The help no longer offers it and says why, the R parser rejects it with the list of the
  supported linkages, as the parser of the grouped selector does, and the C++ header and the Python
  docstring of `hc_linkage` say the same.
- The R help of the fields `FDP_hat_array_BT`, `Phi_array_BT` and `R_array_BT` was wrong ("Variables x
  Thresholds x K", "variables selected", "mapping original to clustered features"). The arrays hold the
  FDP estimates (stopping times x voting levels x correlation thresholds), the dependency-aware relative
  occurrences (stopping times x variables x correlation thresholds) and the numbers of selected variables
  (the stopping times that can be selected x voting levels x correlation thresholds).
- The Python bindings expanded DA as "Data Augmentation" in six places and called BT a "Block-Toeplitz
  structure". DA is dependency-aware and BT is the binary tree. The texts of the five values of
  `DAMethod`, of four values of `ScreenTRexMethod` and of `method` say what the methods do.
- The search index of the R vignettes is rebuilt. It had missed two headings of earlier commits, "The
  penalty of the last step" and "A read-only file".
- Tests: the R test of the invalid parameters expects the new message for an unknown linkage and for
  Ward, and a Python test expects the rejection of Ward. No result changes.

#### The R help example of `TRexDASelector` runs on dependent variables

- The example ran the selector with the method BT on `Gauss_data`, a design with independent variables,
  and selected nothing, for every one of ten seeds. That is the expected behaviour and no defect. The
  selector penalizes the relative occurrences for a dependency that is not there, and on independent
  variables that costs power without a reason. An example that selects nothing shows the class at its
  worst, so it moves to the design the class is made for.
- The example now draws the design of the vignette "T-Rex+DA", an AR(1) process across the columns with
  coefficient 0.7 (n = 100, p = 100, five active variables). With the method BT it selects four of the
  five active variables for each of ten seeds, with no false variable in eight of them and one in two.
- Only the help changes. `R CMD check --as-cran` runs the example and ends with Status OK.

#### The dependency-aware selector takes its penalty from the last step, for all five methods

- The dependency-aware T-Rex selector works with deflated relative occurrences,
  Phi^DA_t(j) = Psi(j) Phi_t(j) with Psi(j) = 1 / (2 - min |Phi(j) - Phi(j')|) over the variables j'
  that j is compared with. There are two ways to build the history over the steps t = 1, ..., T of
  the T-loop. The penalty per step deflates step t with Psi_t, computed from Phi_t. The penalty of
  the last step deflates every step with Psi_T, computed from Phi_T. BT, NN and PRIOR_GROUPS have
  always used the last step, and AR1 and EQUI used the penalty per step. Since this entry all five
  methods use the penalty of the last step.
- This is a deliberate deviation from the papers (Signal Processing 234, 109990, and the AR1 and
  nearest-neighbour papers). Read literally, they recompute the penalty at every step, and the R
  reference implementation does so for AR1 and EQUI. A comment in the code called the last step of
  the BT-style methods an "R-parity deviation from the paper" and left the per-step penalty open as
  the correct form. It is not, for the reason and with the evidence below.
- The reason. The FDP estimator sums over the increments Phi^DA_t(j) - Phi^DA_(t-1)(j) and is a
  count of expected false discoveries only while no increment is negative. Then every share is in
  [0, 1] and Phi' <= Phi^DA_T <= 1 holds by construction. One factor per variable keeps the
  non-negative increments of Phi_t non-negative. A factor per step does not, because
  Phi^DA_t - Phi^DA_(t-1) = Psi_t (Phi_t - Phi_(t-1)) + Phi_(t-1) (Psi_t - Psi_(t-1)), and the second
  term is negative whenever the neighbours of a variable catch up with it. Three things then happen.
  A step whose summed increment is not positive is dropped. A negative increment meets a negative
  scale factor, Phi' rises above 1, the cap sets it to 1 and the variable counts as a certain true
  discovery. A share exceeds 1. In all three cases the FDP estimate is too small.
- This is the third time the same root shows. The first answer was a clamp of Phi' to [0, 1]
  against negative FDP estimates (documented on 2026-07-05), which came from Phi' > 1 in trex_da.
  On 2026-10-24 the lower clamp was removed, because it made the estimate anti-conservative in long
  T-loops, and the upper cap stayed as a guard for trex_da. The cap never repaired the estimator. It
  hid a per-step history as an estimate that is too small. With the penalty of the last step the cap
  binds in at most 0.009 % of the evaluations, and those come from ordinary occurrences that are not
  monotone near saturation. With the penalty per step it binds in 0.9 % to 26.5 % of them.
- The evidence. A lab build applied either penalty to every method. An exploration ran 150 cells of
  300 runs on seven designs, and a confirmation ran the 69 decisive cells with 1,000 runs on new
  seeds (T-LARS, K = 20, block designs with correlation 0.7 and 0.9, AR(1) chains with 0.3 to 0.7,
  independent designs, targets 0.1 to 0.4, also SKIPL with 10 p dummies).
  - FDR. The penalty of the last step is at or below the target in all 28 cells of the confirmation.
    The penalty per step exceeds a target of 0.1 with NN, 0.152 on the independent design and 0.127
    on the block design, where the last step gives 0.000 and 0.026. No cell has the last step above
    the target and the per-step penalty below it.
  - Power. The penalty per step selects more in eight cells (NN, BT on the independent design,
    PRIOR_GROUPS). Where its history was repaired to a monotone one, the gain disappeared, so it is
    made by capped evaluations. BT on the independent design is the clearest case: 0.03 more TPR at
    an FDR of 0.027 against 0.007, and the numbers of the last step after the repair. Outside these
    cells the penalty per step has at most 0.003 more TPR, and on AR(1) chains the last step has
    0.01 to 0.02 more.
  - AR1, the method that changes its results: the FDR stays within 0.006 on designs with dependence
    and far below the target, the TPR rises by 0.013 to 0.020 on AR(1) chains, and the FDR on the
    independent design falls from 0.044 to 0.008.
- Two statements of equality that were written down before the confirmation run did not hold, and
  the record keeps them. The methods do not all have the same power under both penalties (BT on the
  independent design), and the monotone repair of the per-step penalty does not give the numbers of
  the last step on AR(1) chains, where it stays 0.01 to 0.02 below them in TPR. Neither speaks for
  the penalty per step.
- The change. `TRexDASelector::daCorrect` computes one penalty vector from the occurrences of the
  last step for AR1 and EQUI and divides every column of the history by it. The loop over the steps
  is gone, which also makes it cheaper. The library gives bit-identical results to the last-step
  mode of the lab build on 104 selections (13 design and method pairs, 8 seeds, compared on the
  selection, Phi', the Phi matrix, the FDP matrix and T_stop). Against the commit before, 4 of 40
  AR1 selections differ, and EQUI, BT, NN and PRIOR_GROUPS are bit-identical.
- Documentation, on the three surfaces. The class documentation of `TRexDASelector` in C++, R and
  Python states the construction, the deviation from the papers, the reason and the evidence. The
  guide page 04 has a section "The penalty of the last step" in the three guides. The R news file
  lists the change under 2.0.2.
- Tests. C++: `TRexDATest.Deflation_PenaltyOfTheLastStepForEveryMethod` checks the deflated history
  of AR1, EQUI and the BT-style methods value by value against the formula, with a history in which
  a penalty per step gives a negative increment. R and Python: on a tall AR(1) design the calibrated
  occurrences of AR1 never exceed the deflated occurrences of the last step, over six seeds. The
  penalty per step fails this in 9 of 20 seeds.
- For whoever touches this again. Do not move a method to a penalty per step, and do not add or
  remove a clamp on Phi', without repeating the comparison. The analysis, the lab build, its patch,
  the cell scripts and all results are in the developers' notes under
  notes/handover_2026-11-13/da_per_step (ANALYSIS.md, lab/).
- Limits. T-LARS only, p up to 500, and no proof beyond what the papers claim. A solver that drops
  variables has negative increments of its own, which was not examined.

#### The selector tests of the coordinate-descent family run with a sweep budget of 200,000

- Since a coordinate-descent solve that has not converged is an error (an entry below), a selector run
  of the family can end with that error, depending on the draw of the dummies. Over 200 seeds at the
  default budget of 2,000 sweeps, TCDLASSO did so for 6 seeds on the Python test design (n = 50,
  p = 20), its legacy localization for 3 and TCDIENET in the grouped selector for 1. On a design like
  the C++ core test (n = 150, p = 200) TCDLASSO did so for 8 seeds. TCDENET and the R test designs
  never did. At a budget of 200,000 no seed ends with the error on any of these designs.
- Four Python tests ran such a selector without a seed and failed in about one suite run of twelve.
  They now run with a fixed seed.
- The seeded tests passed on the development machine. The dummies come from the distributions of the
  C++ standard library and differ per platform, so a seed that passes on one platform proves nothing
  for another. Every test that runs a selector of the family to completion now sets `cd_max_sweeps`
  to 200,000, in seven C++ test files, two R test files and three Python test files. The test of the
  error itself keeps the default budget, and the tests that run a solver class directly are unchanged.
- No test result changes and no test runs longer. C++ tests 718 of 718, R tests 1,874, Python tests
  697 in five runs in a row.

#### The coordinate-descent solver family is deprecated

- The coordinate-descent solvers TCDLASSO, TCDENET, TCDIENET and TCDTIKNET will be removed in the next
  round of releases, R 2.1.0 and the next minor releases of Python and C++. They brought no gain over
  the homotopy solvers of the same problems. In a benchmark TCDLASSO and TCDENET were 1.1 to 1.5 times
  slower with equal selections, TCDIENET was 2.8 to 30 times slower and ended three of six runs with
  the convergence error, and TCDTIKNET was 1.8 to 2.3 times slower with one error in six runs. The
  replacements are TLASSO, TENET, TIENET and TIKNET.
- R raises a deprecation warning of the class `deprecatedWarning`, as `.Deprecated()` raises it, where
  one of the four solver classes is constructed and where a solver name of the family reaches a
  selector, through `trex_control(solver = )` or `trex_gvs_control(en_solver = "TCDENET")`. Python
  raises a `FutureWarning` on the same events, in the root selectors, in the selector classes of
  `trex_selector_methods.ad`, in the Tikhonov power rule and in the classes of `tsolvers.ad.cd_based`,
  and attributes it to the line of the caller. The warning comes once per construction, before
  anything is standardized or claimed, and never from inside the random experiments. It names the
  replacement. C++ compiles as before and marks the family with `@deprecated`.
- The options that only the family reads get a note in the help and no warning. These are the six
  `cd_*` knobs, `kappa_dummy`, `fold_dummies`, `group_penalty` and `fold_dummy_coupling` of the
  Tikhonov control. They go with the family.
- The help of the classes, of the solver names and of these options says so on the three surfaces, and
  so do the READMEs and the guide pages 05, 08 and 09. The guide examples that run the family stay
  until the removal and show the warning in R and Python. Two wrong statements leave on the way. The
  R help of `tsolver_options()` called the family's options "not bound to R yet", and the Python help
  of `SolverTypeForTRex.TCDTIKNET` said that the Tikhonov selector was not bound in Python.
- The notice goes out with R 2.0.2, Python 1.0.1 and C++ 1.0.1. The R news file starts its section
  2.0.2 with it.
- Tests: new tests of the warnings in R and Python. They expect one warning per construction with the
  replacement in its text, none from `select()` and none for the homotopy solvers and the defaults.
  In Python a warning turned into an error leaves the design untouched. The existing tests of the
  family run with the warning muffled. C++ tests 718 of 718, R tests 1,874, Python tests 697.

#### A coordinate-descent solve that has not converged is an error

- The coordinate-descent solvers stop a solve at their sweep budget (`cd_max_sweeps`, 2,000) and
  accepted the iterate they had, with one warning. Where that iterate was far from a solution, the
  path and the selection built on it were wrong. That happens with the folded dummy coupling of
  TCDTIKNET and TCDIENET once lambda_2 times K is large, because a variable and its dummy copies share
  one penalty term and the sweeps needed grow in proportion to lambda_2.
- A solve that stops at the budget with a KKT violation above 0.001 lambda_1 now throws. The message
  names the budget, the violation and the way out: a larger `cd_max_sweeps`, or the homotopy solver of
  the same problem. The experiment loops carry the exception out of their parallel region, so R gets
  an error and Python a `RuntimeError`. The warning of the first hit stays.
- A `select()` that ends with an error gives the design back at once and retires the selector. Before,
  the design stayed standardized until the selector was destroyed or garbage collected. The same holds
  for the centering of the sparse principal component selector.
- What the check does in practice: it catches the large lambda_2 case, and it never fired in 120
  ordinary runs of the plain selector with TCDLASSO and TCDENET, whose selections all equalled the
  homotopy's. It also fires on designs driven into saturation, where the accepted iterate was usually
  right: 4 of 80 demo-style runs of TCDTIKNET with up to 10 p dummies end with the error at the default
  budget and none at a budget of 200,000.
- The help of `cd_max_sweeps` and of the two solvers states the rule in C++, R and Python, the grouped
  selector's help calls TCDIENET not finally developed, and the guide pages 05, 08 and 12 follow.
- Tests: the solver test of `precentered` runs TCDLASSO with a budget of 200,000 and stops the two
  ridge-type coordinate-descent solvers after a third of the dummies, because deep in saturation they
  do not converge within a practical budget. The Tikhonov equivalence test with up to 10 p dummies
  runs with a budget of 200,000. New tests: a solver test for the error and for the larger budget that
  reproduces the homotopy, and a selector test on each surface for the error, the design given back
  and the retired selector. C++ tests 718 of 718, R tests 1,842, Python tests 679.

### 2026-11-12

#### The Tikhonov selector runs on the homotopy solver TIKNET by default

- `TRexTikhonovSelector` ran on the coordinate-descent solver TCDTIKNET by default. Once lambda_2 times
  K is large, that solver stops at its sweep budget of 2,000 and accepts the unconverged iterate. On the
  chain Laplacian of the guide (n = 120, p = 60, target level 0.1) a fixed lambda_2 of 1,000 gave an FDP
  of 0.17 and a lambda_2 of 3,000 an FDP of 0.70. The default pilot rule scans a grid up to 3,000, took
  the count of the unconverged probes for the knee and picked 3,000, so the default run selected 33
  variables with an FDP of 0.70 after 142 s. The warnings of the solver were the only sign.
- The default is now TIKNET, the pathwise homotopy solver of the method paper, in C++, R and Python. It
  selects the planted block at every lambda_2 in a fraction of a second and prints no warning, and its
  pilot counts equal those of TCDTIKNET wherever that solver converges. TCDTIKNET stays selectable. The
  help says that it is not finally developed and that it can stop at its sweep budget before it has
  converged, which it reports with a warning.
- Results move for every run of the Tikhonov selector that names no solver. `fold_dummy_coupling = FALSE`,
  a validation switch, exists with TCDTIKNET only, so it now needs that solver named.
- The help of the three surfaces and the guide pages 08 name the default and the state of TCDTIKNET. The
  elastic-net anchor of the guide names TCDTIKNET, and the solver example shows it as the alternative.
  Two sentences of the guide are corrected on the way: the folded dummy coupling exists with both
  solvers, only the decoupled dummy block is limited to TCDTIKNET, and in C++ the solver of a fresh base
  control is not mapped to the default. The R help example of the selector set
  `fold_dummy_coupling = FALSE` without a solver, and it now runs on the defaults.
- A new test on each surface runs the default control on a chain Laplacian and expects the planted block
  and a lambda_2 below the top of the pilot grid. With the old default the same data give 13, 10 and 23
  null variables and a lambda_2 of 3,000. The tests that compare with the coordinate-descent selectors
  name TCDTIKNET, and the folded anchor also checks that TIKNET selects the same set as the grouped
  selector. C++ tests 716 of 716, R tests 1,838, Python tests 678.

#### A read-only memory-mapped design serves where nothing is written into it

- `MemoryMappedMatrix::getMap()` refuses a file opened read-only, and every binding of a selector and
  of a T-solver took its `Eigen::Map` through it. A read-only file therefore served no selector and no
  T-solver, although a selector under `ScalingMode::NONE` and a T-solver with `normalize` and
  `intercept` off never write into the design. In Python `to_numpy()` raised for such a file, so its
  data could not be viewed at all.
- The library now states what writes, in one place per family. `data_normalizer::writesDesign()` is true
  for the scaling modes L2 and ZSCORE, and `SolverOptions::writesDesign()` is true with `normalize` or
  `intercept` on. `MemoryMappedMatrix` gains `getReaderMap()` for code that only reads, `getMapFor(writes)`
  that picks the map from that statement, and `mode()`. `getMap()` keeps its refusal, and its message
  now says why and names the way out.
- R takes the map of every selector and of all 20 T-solvers through `getMapFor()`, so an `mmap_matrix`
  in `"readonly"` mode is accepted exactly where nothing is written. `TRexSPCASelector` centers the
  design in every run and refuses it with a message of its own.
- `X` and `D` reach a T-solver as an `Eigen::Map`, a view without a copy, in the three languages. The
  view refers to a matrix in memory or to a memory-mapped file. In R the dummy matrix `D` could only
  be a matrix in memory, while `X` could be either. `D` may now be a memory-mapped file as well, as
  in Python and C++. One creator per solver takes either kind for both arguments, and the 20 creators
  for a memory-mapped `X` are gone.
- The same gap is closed for the response matrix `Y` of the biobank screening, which R took as a
  matrix in memory only. `Y` may now be a memory-mapped file in R as well. The selector only reads
  `Y`, one column at a time, so a file opened read-only serves under every scaling mode, and Python
  accepts a `Y` that is not writable.
- Python keeps its idiom, the `to_numpy()` view. The view of a read-only file is an array that is not
  writable, and every selector and T-solver has a second constructor that takes a Fortran-contiguous
  float64 array as it is, without a copy. It accepts an array that is not writable where nothing is
  written and raises `ValueError` with the rule everywhere else, where such an array met a `TypeError`
  without an explanation before. The rule reads the write flag, so it holds for every such array, for
  `X` and for `D` of a T-solver.
- An audit ran before the bindings changed, because a write into a memory-mapped file opened
  read-only ends the process instead of raising an error. It runs every case in a child process on a
  design in write-protected memory. 203 cases complete without a write: all solver types of the
  selector on a small and a wide design, the dependency-aware methods, the grouped selector with
  every solver and every rule for lambda_2, the Tikhonov selector with every rule, both solvers and
  both dummy mechanisms, the screening variants with the bootstrap, twelve phenotypes of the biobank
  screening with its fallback, the response matrix of the biobank screening, and each of the 20
  T-solvers as a reader. The three controls that must write are killed.
- The help of the design argument, of the scaling mode, of `mem_matrix()` and of the solver contract
  states the rule in C++, R and Python. The guide pages 11 of the three surfaces gain a section with a
  run on a read-only file, and the pages 07 and 12 state the rule.
- The help of the memory-mapped matrix and the guide pages 11 warn that one file serves one writer at
  a time. Several selectors or T-solvers that standardize one file in place at the same time, through
  a second handle or in several processes, change each other's data and return wrong selections. The
  library itself never does this, and the claim on a buffer refuses a second writer on one matrix or
  handle in one process. To use one design in several processes at once, the design is pre-normalized,
  opened read-only and selected from under `NONE`.
- The Python T-solvers refuse an array they cannot share, for example a C-ordered array, the default
  layout of NumPy, with a `TypeError`. Their docstring and guide page 09 said such an array is copied
  first. Both now state the refusal and name `numpy.asfortranarray()`, and a test pins it.
- The tests use real files opened read-only on the three surfaces. A run on the file equals the run in
  memory, every selector and every solver class completes, the writers are refused, and the file is
  unchanged afterwards. A writer on a memory-mapped `D` standardizes the file in place, claims it and
  restores it. C++ tests 715 of 715, R tests 1,833, Python tests 677.

#### The result of the biobank screening says what its FDR level is

- `estimated_FDR` of a biobank result holds the estimate of the accepted screening variant, or, for a
  fallback result, the target level at which the T-Rex selector controls the FDR. The documentation
  called it an estimate in both cases. An estimate is one over the number of selected variables, so it
  can equal the target, for example 0.1 with ten selected variables, and the number alone does not tell
  the cases apart. The help now says that `used_fallback_trex` does, and that the two screening
  estimates and the two screening index lists are stored for every phenotype, a fallback result
  included. The header said these lists are empty when another method gave the final selection, which
  the code never did.
- The Python text of `method_used` named values that the code does not return. It now names the three
  it returns, "Screen-TRex (ordinary)", "Screen-TRex (bootstrap-CI)" and "T-Rex (fallback)".
- The change is documentation in C++, R and Python and one passage in the guide pages 06 of the three
  surfaces. The code does not change. The tests pin the contract: a fallback result reports the target
  level, carries the flag and keeps both screening estimates. C++ tests 711 of 711, R tests 1,733,
  Python tests 665.

#### The biobank screening gives every phenotype its own seed

- `BiobankScreenTRex` passed its seed to the screening of every phenotype, so all phenotypes competed
  against the same dummies. Phenotype i (0-based) now runs its screening and its fallback selector with
  seed + i, so the phenotypes draw independent dummies, as the R reference does by continuing its random
  stream. The dummy generator mixes the seed, so neighboring seeds give unrelated streams. A negative
  seed is passed on as before, and every selector then draws its own entropy.
- Phenotype 0 and the single-phenotype mode keep the seed, so their results do not change. The results
  of the later phenotypes of a seeded run change. A single-phenotype run with seed + i reproduces
  phenotype i of a run.
- The help of `seed` says this in C++, R and Python, and the guide pages 06 of the three surfaces carry
  one sentence on it.
- The tests screen twelve copies of one noise response, because the screening of such a response selects
  something in about half of the dummy draws. Every copy equals the single run with its seed, and the
  copies do not all agree. Both tests fail when all phenotypes share one seed. C++ tests 711 of 711,
  R tests 1,729, Python tests 664.

#### The screening selector checks its bootstrap parameters and reports 1 for an empty selection

- `ScreenTRexSelector` checks `R_boot` and `ci_grid_step` at construction. `R_boot` is at least 2,
  because the standard error of the bootstrap divides by the number of resamples minus one, and
  `ci_grid_step` lies in (0, 1], because a step of 0 never advanced the grid over the confidence level
  and the search did not end. The check runs whether or not `use_bootstrap_CI` is set, because the
  biobank screening applies the bootstrap to every phenotype. The R glue rejects a negative `R_boot`
  before the conversion to an unsigned count.
- The estimated FDR of an empty selection is 1, as the estimator 1 / max(1, R) of the screening paper
  (Machkour, Muma and Palomar, SSP 2023) and the documentation of the field state it. The code returned
  0 there, like the R reference. An estimate of 0 lay inside every acceptance window that starts at 0,
  so the biobank screening accepted an empty result where the paper falls back to the T-Rex selector.
  With the default window of 0.05 to 0.15 the decision is the same before and after, and only the two
  reported screening estimates of an empty result move from 0 to 1.
- The help of `R_boot`, `ci_grid_step` and `estimated_FDR` says this in C++, R and Python. The guide
  pages 06 of the three surfaces give the estimator and describe the bootstrap variant as what it is, a
  confidence interval of the mean of the non-zero dummy coefficients with a selection by the averaged
  coefficients, where they spoke of resampled experiments and an upper confidence bound.
- The tests use a response that is orthogonal to the design and to the intercept, which leaves every
  variable behind the first dummy, so the selection is empty for certain. They fail on the code before
  the change. C++ tests 709 of 709, R tests 1,655, Python tests 663.

#### The screening selector gains a nearest-neighbor variant

- `TREX_DA_NN` is a new variant of the screening selector. The group of a variable holds the neighbors
  inside the AR(1) window whose sample correlation with the variable reaches `rho_thr_DA`, and the
  relative occurrence is divided by 2 minus the smallest distance to the relative occurrences of that
  group, as in the other dependency-aware variants. A variable without such a neighbor keeps its relative
  occurrence. This rule for an empty group is a choice of the library, because the screening has no
  calibration over the threshold and the penalty of one half of the dependency-aware papers would leave
  nothing above the voting level on uncorrelated data. The window comes from the shared `ar1HalfWidth()`,
  the coefficient from the AR(1) estimate or from `cor_coef`, and the threshold lies in [0, 1). The cost
  is of order n p kappa.
- The variant is bound in R (`trex_method = "TREX_DA_NN"`) and in Python
  (`ScreenTRexMethod.TREX_DA_NN`), with the help of `trex_method` and `rho_thr_DA` on the three surfaces.
  The guide pages 06 of the three surfaces describe the four dependency-aware variants of the screening
  as what they are, a penalty on the relative occurrences with independent dummies, where they spoke of
  dependent dummies, and they name the new variant and the parameters `cor_coef`, `rho_thr_DA` and
  `n_blocks`.
- The C++ code agrees with an independent Python prototype in all 1,200 runs of a probe (independent,
  AR(1), equi-correlated and block designs, three thresholds, 100 runs each). In that probe the variant
  behaves like the AR(1) window on the AR(1) design and gains a little power on blocks of five correlated
  variables, it equals plain screening where no neighbor is correlated, and it does not reach the
  equi-correlated design, whose dependence is not local. The tests pin a zero threshold (equal to
  `TREX_DA_AR1`) and centered orthonormal columns (equal to plain screening, with a selection that is not
  empty), and the rule for `cor_coef` covers the variant. C++ tests 706 of 706, R tests 1,643, Python
  tests 649.

#### The dependency-aware selector and the screening selector share one AR(1) window

- The half-width of the AR(1) window, ceil(log(`rho_thr_DA`) / log(|`cor_coef`|)), is one function now.
  `ar1HalfWidth()` lives in the new header trex_utils/trex_da_window.hpp, and `TRexDASelector` (AR1)
  and `ScreenTRexSelector` (TREX_DA_AR1) both call it. It uses the magnitude of the coefficient, returns
  1 for a zero coefficient and the full window for a zero threshold or a coefficient of magnitude 1,
  never exceeds p, and rejects a threshold outside [0, 1) and inputs that are not numbers. The screening
  selector had a copy of its own, which converted the half-width to an integer without a check.
  `rho_thr_DA = 0` therefore converted infinity, which is undefined behaviour, and a negative threshold
  gave NaN.
- A supplied `cor_coef` outside [-1, 1] gives a warning and is replaced by the estimate, in both
  selectors and for every dependency-aware method. Only the value -2 requested the estimate before,
  and a value such as -1.5 was taken as a coefficient and led to the full window. The value -2 still
  requests the estimate without a warning. `isCorrelationCoefficient()` in the same header holds the rule.
- The screening selector checks its dependency-aware parameters at construction. `rho_thr_DA` lies in
  [0, 1] and in [0, 1) for TREX_DA_AR1, and `n_blocks` lies between 1 and p for TREX_DA_BLOCK_EQUI,
  where 0 divided by zero. `TRexDASelector` rejects a threshold that is not a number as well.
- `trex_screen_control()` of the R package gains `cor_coef`, `rho_thr_DA` and `n_blocks`. C++ and Python
  could set them before and R could not. The defaults are the C++ defaults, so existing calls keep their
  results.
- The help of `cor_coef`, `rho_thr_DA` and `n_blocks` describes the methods in C++, R and Python. The
  Python text of `cor_coef` claimed an estimate for values of -1 and below and now says that setting
  the parameter to -2 estimates it from the data. The guide pages 04 of the three surfaces describe the
  AR1 and EQUI methods as a penalty on the relative occurrences with independent dummies, where they
  spoke of dependent dummies, and name `rho_thr_DA` for both methods and the rule for `cor_coef`.
- The equi-correlated variants are unchanged. A probe of the screening variants on independent, AR(1)
  and equi-correlated designs gives the same output before and after the shared window (nine cells of
  100 runs). The new tests of the checks fail on the code before the change. C++ tests 705 of 705,
  R tests 1,624, Python tests 645, and the install smoke test finds the new header.

#### The help page of the R package lists the authors of DESCRIPTION again

- man/TRexSelector-package.Rd was last generated on 2026-11-05 and kept an author block that
  DESCRIPTION had left since. It still showed an e-mail address of the former contributor, which
  DESCRIPTION no longer carries. roxygen could not repair the page, because it reads `Authors@R` as one
  line. The first `#` comment then hides the rest of the field, the evaluation fails and the page comes
  out without an author block. The four comment lines leave `Authors@R`, and AUTHORS.md keeps the
  grouping of the members. The regenerated page matches DESCRIPTION, and `Config/roxygen2/version`
  records 8.1.1, the version that wrote it. The other Rd pages and NAMESPACE come out byte-identical,
  and the tracked tree holds the address nowhere else.

#### The AR1 window of the dependency-aware selector follows the magnitude of the coefficient

- `TRexDASelector` builds the window of its AR1 method with the half-width
  ceil(log(`rho_thr_DA`) / log(|`cor_coef`|)), as Definition 1 of the paper on autoregressive
  dependency models (Machkour, Muma and Palomar, EUSIPCO 2024) states it. Two inputs left that
  definition. A supplied negative `cor_coef` gave the full window of all p variables, because only a
  positive coefficient entered the formula. `rho_thr_DA = 1` gave a half-width of 0, a window without a
  neighbor, which switched the correction off without a message. The window now follows the magnitude
  of the coefficient, and the AR1 method rejects `rho_thr_DA = 1` with an error, because its range is
  [0, 1) while EQUI keeps [0, 1]. A half-width beyond the design counts as the full window. The
  estimated coefficient was a magnitude before, so runs with the defaults do not change, and `cor_coef`
  in the result stays the coefficient as supplied.
- The help of `cor_coef` and `rho_thr_DA` describes both methods in C++, R and Python (the header,
  `trex_da_control()` with its Rd page, the binding docstrings with the regenerated stub).
- The C++ test pins the half-width for the coefficients 0.5 and -0.5, for a zero and a unit
  coefficient, for a zero threshold and for a window wider than the design. The R and Python tests
  compare the runs with -0.5 and 0.5 and check the rejected threshold. The new tests fail on the code
  before the fix. C++ tests 699 of 699, R tests 1,568, Python tests 633.

### 2026-11-11

#### The solver types register with cereal in their sources

- CRAN's check of 2.0.1 on r-release-macos-x86_64 ended with FAIL and without an error. The log stops
  inside the re-building of the vignettes at 30:03 of wall clock, after an install of 26 minutes, so the
  time limit of the builder ended the check and the compile time was the lever. Eight solver headers
  ended with `CEREAL_REGISTER_TYPE`, two of them with `CEREAL_REGISTER_POLYMORPHIC_RELATION` as well.
  Every source that included a solver header therefore compiled the save and load code of these solvers
  again, in up to 36 of the 48 objects of the R package. Under clang 14 with libc++ 14, the compiler
  generation of that builder, 77 % of the function code in the objects was such repeated code.
- The ten macro lines moved to the end of the eight solver sources (TLASSO, TENETAug, TIENET, TIENETAug,
  TIKNET, TCDENET, TCDIENET and TCDTIKNET). The bindings cover the binary and the portable binary
  archive as before. These are the archives that work with a solver, because the matrices are stored as
  binary blocks and a text archive never compiled. Each header keeps `CEREAL_FORCE_DYNAMIC_INIT` and
  each source carries `CEREAL_REGISTER_DYNAMIC_INIT`. This pair of cereal links the registering object
  out of the static library into every program that includes the header, so a program that reloads a
  solver through a `TSolver_Base` pointer and uses nothing else of its class works as before.
- Measured under clang 14, four files at a time, the sources of the R package take 896 CPU seconds
  before and 351 after, and the largest compiler process needs 1.3 GiB where it needed 1.9 GiB. The C++
  library with its tests builds in 290 CPU seconds where it took 405 (Clang 23, Release), the Python
  extension in about 160 where it took 208, and the manylinux wheel in about 14 CPU minutes where it
  took 17:30. No interface or result changes. The C++ tests pass on macOS and on Linux with GCC 13.3
  (698 of 698), the Python tests on macOS and in the manylinux_2_28 wheel (631 of 631) and the R tests
  (1,562). A consumer of the installed static libraries saves and reloads the augmented solver and a
  solver behind a base pointer with both archives.

#### The Rcpp bindings of the R package compile as one unit

- `bundle_for_cran.sh` writes src/rcpp_bindings_unit.cpp. The generated file includes the twelve binding
  sources at the top of src/, the eleven `rcpp_*.cpp` files and RcppExports.cpp last, and its object
  takes their place in OBJECTS, so the package compiles 37 objects where it compiled 48. The binding
  sources read the same Rcpp, Eigen and selector headers, which the compiler now reads once.
  `--check` compares the unit with the binding sources that exist, so a new binding file without a rerun
  of the script fails the vendor check. `Rcpp::compileAttributes()` writes the same RcppExports files
  with the unit in place. The comment on the debug symbols in the generated Makevars no longer counts
  the translation units. The IntelliSense target of cpp/CMakeLists.txt leaves the unit out.
- Measured on macOS with Apple clang 17, one file at a time as CRAN compiles, the install takes 5:01
  before, 3:32 with the registration in the sources and 2:55 with the unit on top. Under clang 14 the
  unit adds little to the 351 CPU seconds (344 with both). The unit is the largest single compile there
  with 58 seconds and 1.9 GiB. Compiled for Windows with a MinGW-w64 cross compiler it has 3,158
  sections of the 32,767 a COFF object can hold, so `-Wa,-mbig-obj` is not needed. `R CMD check --as-cran`
  of the tarball ends with Status OK on macOS (R 4.6), with an install of 174 seconds. Under R-devel
  with GCC 16.2 in the container it ends with the two NOTEs of the container itself, and the install
  takes 493 seconds where the 2.0.1 tarball took 13 minutes.
- A change in one binding source recompiles the whole unit. An incremental `R CMD INSTALL` inside the
  tree does not notice such a change, because make knows the unit and not the files it includes.
  tools/install.sh removes the objects first and stays the way to reinstall.

### 2026-11-10

#### The R package carries the setup of its pkgdown site

- `_pkgdown.yml` prepares the R part of the documentation site at
  <https://trexselector.github.io/TRexSelector/r/> with the site URL in place of the repository URL,
  Bootstrap 5 with the light switch, the teal #336C68 of the other guides, the footer "Copyright © 2026
  Team T-Rex" and KaTeX for the formulas. MathJax stays out for now, because pkgdown >= 2.1.2 assigns
  `window.MathJax` after the MathJax script and loses the autoloaded commands such as `\boldsymbol`,
  reported as r-lib/pkgdown#3037, and a comment next to the setting says so. pkgdown/extra.js adds the
  page position and the previous/next links of the vignettes to the guide pages. pkgdown/extra.css
  styles them, keeps the white background of the figures in dark mode, where pkgdown inverts them, and
  sets the teal links of the dark theme. .Rbuildignore excludes pkgdown/, so CRAN, `vignette()` and
  RStudio keep the designed vignettes. DESCRIPTION gains the site URL once the site is live.

#### The Python guide encrypts its kernel connection and keeps the property texts

- Since ipykernel 7.4 every notebook of the Sphinx build logs "Kernel is running over TCP without
  encryption". docs/conf.py lets the nbclient kernel manager provision CurveZMQ keys
  (`transport_encryption` "auto"), because MyST-NB has no option for it. "auto" keeps plain TCP for a
  kernelspec without Curve support, and the hook applies only when pyzmq has Curve. requirements.txt pins
  jupyter-client 8.10.0, the version that provides the setting. The docstring hook of the API pages gains
  two fixes. napoleon read a leading "Label: text" of 57 properties as a type and a description, so these
  properties keep their docstring as written, and the pybind11 overload lists such as
  `1. f(self: T, ...) -> R` show as code without self, which also removes the warning about `**kwargs`.

#### A CRAN check script replaces the build chain, and two debugging leftovers leave

- R/TRexSelector/tools/cran_check.sh runs the checks before a CRAN upload in one go. It runs the vendor
  check, `R CMD build` and `R CMD check --as-cran` with the remote part of the incoming check off and
  CRAN's timing thresholds at 2.5, prints the size and the SHA-256 checksum of the tarball, writes into a
  dated folder under the gitignored notes/cran_upload/ and ends with a done line. It replaces
  build_chain.sh, which built and checked without the vendor check and the thresholds, and which nothing
  referenced. Python/trexselector/run_lldb.sh, an lldb wrapper around a test script that no longer
  exists, and R/TRexSelector/tools/debug.R, a devtools loader that the VS Code launch configurations do
  not use, are archived outside the repository.

## TRexSelector R 2.0.1

Patch release of the R package after the public release (CRAN 2026-11-06, GitHub and PyPI
2026-11-07), uploaded to CRAN on 2026-11-10: the C++20 ranges calls that Apple clang 14 rejects are
gone, the R package's configure finds R's OpenMP flags on R 4.6 and R-devel again, the sources compile
without warnings of their own, and the READMEs describe the packages and install from CRAN and PyPI
first. No interface or result changes. The sources of this release carry the C++ SDK and the Python
package as 1.0.1, and neither was published from them. Their 1.0.1 is the release above.

### 2026-11-09

#### The CRAN comments name the win-builder runs of the submitted sources

- cran-comments.md names the win-builder runs of 2026-11-09 on the submitted sources, R-release 4.6.2
  and R-devel r90850, in place of the runs of the day before. Both were Status 1 NOTE, the
  incoming-feasibility NOTE alone, with the day count, the maintainer line of the test copy that routes
  the result mails, and a connection timeout of win-builder towards the GNU site on the README's license
  link, which is reachable and was not flagged on the two days before. The results paragraph says so in
  two sentences, so the reviewer reads the timeout as what it is, and the link stays.

#### The CRAN comments count four days

- cran-comments.md reads "Days since last update: 4", the count of the upload on 2026-11-10. The number
  stays in the text because it is what the reviewer sees in the incoming check, and the sentence
  explains it.

#### The release files of 2.0.1 catch up with the tree

- NEWS.md names what the vignettes gained after its entry was written: the figures of four pages, the
  ggplot2 plots of five, the links to the demos of the Examples repository on the method pages, and
  the package title without the ampersand. cran-comments.md says the same in its summary and names the
  overview vignette, the dimensionality page, the figures, the search box, the title and Team T-Rex.
  Its test environments say that the Actions jobs build and re-build the vignettes and that
  win-builder (R-release 4.6.2 and R-devel, 2026-11-08) built them as well, and the incoming NOTE it
  explains reads "Days since last update: 3", the count of an upload on 2026-11-09. CITATION.cff names
  2026-11-09 as the release day. The tarball built from this tree (366 entries, the fourteen vignettes
  inside) passes `R CMD check --as-cran` locally with Status OK and 1,562 tests.

#### AUTHORS.md separates the former contributor from the members

- The Members table lists the current team alone. A new section "Former contributors" after it names
  Simon Tien with his contribution, contents of the R releases prior to 1.0.0, without an address, the
  form other projects use for their emeritus lists, and the roles paragraph of the R package says the
  same. DESCRIPTION keeps his contributor role but drops the student address, which no longer works; a
  contributor's address is optional for CRAN.

#### The linters follow the 120 ruler

- `.lintr` at the repository root raises lintr's line length to 120, the ruler of the repository,
  and object names to 40 characters, excludes the generated `RcppExports.R`, and turns off four
  linters that contradict the style of the package: object names (parameters such as `tFDR` and
  `R_boot` follow the notation of the papers and the C++ core), indentation (continuation lines are
  indented by a fixed step, not aligned to the opening parenthesis), explicit `return()` at the end
  of a function, and commented code (a header comment lists the control functions). The editor's R
  diagnostics for the package's R folder drop from 462 to none, after the one real finding, a missing
  space in `verbose_ = NULL` of the SPCA selector, was fixed. The tests keep their compound lines
  with semicolons and their one-line function bodies, so the semicolon and brace linters are
  switched off for the tests folder, and their findings drop from 114 to none. lintr walks up the directories from
  the linted file, so the file covers the R package and the R scripts under `tools/`, and it lies
  outside the package, so it never enters the tarball.
- `.markdownlint.yaml` at the root sets the line length to 120 outside tables and code blocks, keeps
  tables column-aligned (MD060, style "aligned") and `<span>` as the only HTML element, and turns the
  undefined-reference rule off, the configuration the editor settings carried so far. Configuration
  files take precedence over the editor setting, which the extension has deprecated. The file also
  forbids trailing whitespace (MD009 with no break spaces) and allows `<br>` next to `<span>`, because
  the editor now trims trailing whitespace on save in Markdown as well, so a hard line break is
  written as `<br>`, never as two trailing spaces. No tracked Markdown page used the two-space form,
  only the license text of the R package, which stays verbatim. The thirteen lines with trailing
  whitespace in the sources, seven in the screening test, three in the memory-mapped matrix source
  and three in two Rcpp wrappers, are clean, and the vendored copy follows.

#### OPENMP.md reads in shorter sentences

- The opening paragraph lists the parallel regions after a colon and continues with one sentence per
  line. The CRAN rule sits in the R section as a subsection, since it concerns the R package alone,
  and the table of the check has aligned columns, which clears the editor's table warnings.

#### The vendored copy of the DA header is in sync again

- The CI job "vendored sources in sync" failed on the bibliography commit: the reference block of
  `trex_da.hpp` gained the pages of the index-tracking paper under `cpp/src`, while the vendored copy
  under `R/TRexSelector/src` kept the old line. `bundle_for_cran.sh` synced the copy, and its
  `--check` passes again. A change under `cpp/src`, even in a comment, needs the bundle script run
  before the commit.

#### PUBLICATIONS.md opens with one sentence

- The page opens with "Every paper of the Terminating-Random Experiments (T-Rex) framework and its
  extensions, newest first." and nothing else before the first part. The sentences on the
  generating script, the two parts and the user guides are gone, since the headings say what the
  parts hold, and an HTML comment at the top of the page, invisible when rendered, names the script
  and says that edits go into the bibliography or the script. The template of
  `tools/bib_to_publications.py` carries the sentence and the comment, so a run reproduces the page.

#### The bibliography gains two papers and the pages of three more

- `references.bib` gains two papers of the method family, both with fields from the local PDF and
  the Crossref record: Vilella, Machkour, Muma and Palomar, "Learning false discovery rate control
  via model-based neural networks," ICASSP 2026, pp. 841-845, and Gargano, Machkour, Nardelli,
  Scilingo and Muma, "Reproducible physiological features in affective computing: A preliminary
  analysis on arousal modeling," MetroXRAINE 2025, pp. 605-610. PUBLICATIONS.md lists them, nineteen
  entries now. They stay out of `inst/CITATION` and CITATION.cff, which name the papers behind the
  implemented methods.
- The file also names the pages of the three papers that lacked them: the EUSIPCO 2026 papers on
  dummy-controlled Gaussian graph learning (pp. 1242-1246) and on FDR-controlled compressed
  beamforming (pp. 2716-2720), read from the footers of the local proceedings PDFs, and the ICASSP
  2025 paper on FDR-controlled index tracking (pp. 1-5) per the Crossref record. Everything rendered
  from the bibliography follows: PUBLICATIONS.md, `inst/CITATION`, the reference block of the
  dependency-aware selector in the C++ header, the roxygen block with its Rd page and the pybind11
  docstring with its stub, and the index-tracking reference of CITATION.cff. The vignette copy of the
  bibliography is the master file again; it had lacked the Huffman 2018 entry of the root README.
  The EUSIPCO 2026 papers have no DOI yet, since the proceedings are not on IEEE Xplore so far.

#### Bare URLs of the README bibliographies stand in angle brackets

- The software entries of the README bibliographies named their home pages as bare URLs, which GitHub
  links by itself while other Markdown renderers do not. `render_references.py --format md` now puts
  a bare URL into angle brackets, so the Eigen, Boost, cereal, OpenMP, pybind11 and CRAN package
  entries of the four READMEs read `<https://...>` and link everywhere. The blocks of the root, R,
  Python and C++ READMEs were re-rendered, and the R package's README.md was knit from README.Rmd
  again.
- The root README wraps its long citation lines, says "for the same data and parameters" about the
  results across languages and splits the Citing paragraph into one sentence per line. Its
  Bibliography opens with "The papers behind the opening story:" and lists exactly those three, the
  T-Rex paper, Huffman 2018 and Harvey, Liu and Zhu 2016; the dependency-aware paper and the
  virtual-dummies preprint left the list, since the story does not cite them and the package READMEs
  and PUBLICATIONS.md carry them. The folder links point at the folders without a trailing slash, the
  form the Markdown link validator of VS Code accepts for every folder (with the slash it reported
  two of them as missing); GitHub resolves both forms alike.
- The R, Python and C++ READMEs open their Bibliography with "The three papers at the theoretical
  core of the framework so far:", in the words of the landing pages, and their Citing sections gained
  the pointer to PUBLICATIONS.md that the root README keeps there. The software lists of the four
  READMEs open with a sentence of the form "The C++ core builds on the following software:".
- The pipe tables of the C++, Python and R READMEs have aligned columns, and the lists of the three
  READMEs use the asterisk throughout, as the root README does. That clears the editor's markdownlint
  warnings for the two files that showed them (table column style, list marker style), about fifty
  each. This file lost its two warnings as well: two member names with trailing underscores in the
  T-MP block of 2026-07 stand in backticks, so they no longer read as emphasis markers, and the second
  of two identical "R <-> Python binding parity" headings of 2026-07-08 names its subject, the
  elastic-net model selection.
- AUTHORS.md opens with "The software centered around the Terminating-Random Experiments (T-Rex)
  selector, the TRexSelector C++ library and its R and Python packages, is developed by Team T-Rex.",
  followed by where the team works and the copyright sentence. The e-mail routing bullet of the Contact
  section is gone, since the Members table carries the addresses, and the table has aligned columns.

### 2026-11-08

#### The repository is cited as software

- GitHub's "Cite this repository" box offered the Signal Processing article of 2025, because
  `CITATION.cff` carried it as `preferred-citation`, and GitHub shows that block instead of the
  software whenever one exists. The block is gone, so the box now cites the software itself, "Team
  T-Rex. (2026). TRexSelector: T-Rex Selector: High-Dimensional Variable Selection and FDR Control
  (Version 1.0.1) [Computer software]", and the article leads the list of references with the other
  method papers. The file also names the release date and the project website. The version stays the
  one of the C++ library and the Python package, whose release tag is `v<version>`, while the R package
  keeps its own version in DESCRIPTION and `citation("TRexSelector")` answers with the package entry.
- The long title reads "High-Dimensional Variable Selection and FDR Control" everywhere now, without
  the ampersand: in CITATION.cff, in the Title field of DESCRIPTION, and in the package entry of
  `inst/CITATION`, which `tools/references/make_inst_citation.py` writes.
- Regenerating `inst/CITATION` from the bibliography also brought the DOIs that `references.bib` gained
  for six conference papers into the R entries, and the same six references of CITATION.cff carry
  them now. Every method paper with a DOI names it in both files.

#### The CRAN links take the canonical form

- win-builder's incoming check of the 2.0.1 tarball (R-release, Status 1 NOTE) listed the link to the
  reference manual PDF in the vignettes' left pane and overview as a CRAN URL not in canonical form.
  The two links now read `https://CRAN.R-project.org/package=TRexSelector/TRexSelector.pdf`, which CRAN
  redirects to the PDF, and the Python and C++ landing pages link the R overview vignette through
  `package=TRexSelector/vignettes/` the same way. The other parts of that NOTE are the two days since
  2.0.0, the temporary maintainer of the uploaded copy, and 404s for AUTHORS.md, OPENMP.md and
  PUBLICATIONS.md, which exist on GitHub only after the push. Tests and the fourteen vignettes passed
  on win-builder, and the install log shows the OpenMP flag on every compile line under GCC 14.3.

#### The four READMEs end with Citing and Bibliography

- The root README folds its References and Acknowledgements sections into one "Bibliography" below
  "Citing", the word the guides use for the literature list: first the papers behind the opening story
  and the variants, then the software the C++ core builds on (Eigen, Boost, cereal, OpenMP), both
  rendered from `references.bib` by `tools/references/render_references.py`, which now keeps a
  corporate author written in braces as written instead of reducing it to initials. The R, Python and
  C++ READMEs get the same pair in place of their Acknowledgements (the C++ README had its Eigen quote
  under Dependencies): Citing with the route of that package, `citation()` or CITATION.cff, and
  Bibliography with the three core papers and the package's own stack, Rcpp, RcppEigen, BH, Rcereal
  and R6 for R, pybind11 and NumPy for Python, Boost, cereal and OpenMP for C++, Eigen on all. The
  Links of the R and Python READMEs point at PUBLICATIONS.md instead of the former References anchor,
  and the R README is knit.

#### The method pages link the demos of the Examples repository

- Pages 3 to 7 and 9 to 11 of the three guides end with an "Examples" section before the bibliography:
  the TRexSelector_Examples repository holds runnable demos and Monte Carlo studies of the page's
  method, with their results and a README per demo, and the section links the demo folders of the
  three surfaces on the main branch of that repository (the selectors under
  `trex_selector_methods/ad/`, the T-solvers under `tsolvers/ad/`, the ML methods and the memory
  mapping). So the guides stay the short walk-through and the comprehensive studies live where they
  run. Pages 2, 8, 12 and 13 have no demo folder yet and keep no section. The search index follows
  (114 headings).

#### A publications page, the software the packages build on, and the documentation toolchain

- `PUBLICATIONS.md` at the repository root lists every paper of the method family, journal articles,
  conference papers, theses and preprints, newest first. `tools/bib_to_publications.py` writes it from
  the "T-Rex framework" section of `references.bib`, so the bibliography stays the single source, and
  `--check` reports a stale page. The section gains the three papers of the library folder it lacked,
  the neighborhood screening paper for Gaussian graphical models (EUSIPCO 2024), the doctoral thesis
  behind the framework (TUprints DOI) and the dummy-controlled graph learning paper (EUSIPCO 2026, no
  DOI yet), fourteen in all, and a second part lists the three papers that carry the framework to
  complex-valued data (the ICASSP 2025 and EUSIPCO 2025 papers with their DOIs, the EUSIPCO 2026 radar
  paper without one yet), which the bibliography keeps in a section of their own. A vignette would
  carry such a page only with a CRAN release, the repository file updates with every paper. The three
  landing pages link it under Reference, as do the R pane and the root README.
- The getting-started pages cite the software the packages build on at the end of their installation
  section: Rcpp, RcppEigen, BH, Rcereal and R6 on the R page, NumPy and pybind11 on the Python page,
  Eigen, Boost with Boost.Interprocess, cereal and OpenMP on the C++ page. The landing pages cite the
  documentation toolchain, knitr, rmarkdown and ggplot2, Sphinx, MyST-NB and matplotlib, Doxygen and
  matplotlib, and the R overview ends with a bibliography like every other page. The R package entries
  come from `citation()` of the installed packages, matplotlib from Crossref, the rest are software
  entries with their home pages. Both copies of `references.bib` carry the sixteen new entries, and
  LinkingTo, Imports and Suggests of the R package are alphabetical.

#### The three surfaces agree on the landing page and on pages 1 and 2

- A comparison tool (`notes/handover_2026-11-07/tools/compare_surfaces.py`) strips a page to its prose,
  normalizes math, citations, links and headings, and diffs the Python and the C++ page against the R
  vignette sentence by sentence. It showed four drifts, now closed: the Python and C++ getting-started
  pages point at the T-Rex selector page for the seed advice, carry the "Page" column of the selector
  table, and say that the T-solvers page describes the families, and the three landing pages open with
  the same three sentences (the package or library and its relation to the others, the code that runs
  at build time with every shown output coming from it, the reference describing every class and
  function) followed by "The pages cite the literature". What remains differs by surface on purpose:
  installation, the shipped data set of R against the drawn one, R6 fields against attributes and
  getters, the registry links, and "Indices are 0-based" on the Python and C++ pages. The paragraph on
  the three result matrices, with its three-line code, joins the Python and C++ getting-started pages
  (a new `matrices` section of the C++ example), and the three landing pages quote the preprint on
  virtual dummies by Koka, Machkour, Palomar and Muma below the two journal articles.

#### The dimensionality pages after the review

- The "Where to go next" section is gone on all three surfaces, since the left pane and the landing
  page hold the map, and no other page had one. The regime map reaches p = 1e8, so the title of the
  highest regime, now broken as "Highest" and "dimensionality", sits inside its region with room to
  spare. The boundary labels run along their lines with the on-screen slope computed after the
  layout (rotating text in data coordinates misleads on logarithmic axes, since the offset is linear),
  "p = 1000 n" below the dashed line and left of the title, and the block of the high regime moved
  up and left, clear of the dashed line. The sentence on false discoveries cites the essays of the
  library for the laboratory (Ioannidis 2005, MacArthur 2012, Uffelmann et al. 2021, Tandy-Connor et
  al. 2018, all read from the local PDFs and verified at Crossref) and Harvey, Liu and Zhu 2016 for the
  investment team. The four essays join the application references of `references.bib`, and the
  vignette copy gains that section with the finance entry. The FDR paragraph says "multiple hypothesis
  testing", cites Benjamini's 2010 review of the rate's origins, Holm 1979 for the family-wise error
  rate (a local PDF, the 1979 volume carries no DOI, so the stable JSTOR link stands in), and the
  low-dimensional section cites Benjamini and Yekutieli 2001 for the dependence conditions of the
  Benjamini–Hochberg procedure, all three entries new in both copies of the bibliography. The high and
  the highest regime are three sentences each: the solvers select without controlling the FDR because
  their stopping rules are tuned for prediction error, the T-Rex selector wraps them and supplies the
  control they lack, and for strongly dependent variables T-Rex+DA restores the exchangeability its
  guarantee rests on; the mechanics of dummies and calibration stay on page 3. The three regime
  sections are "Low dimensionality", "High dimensionality" and "Highest dimensionality", as in the map.
  The R stylesheet justifies the landing page's opening paragraphs and quote, which stand above its
  first section. The page is "Data
  dimensionality and FDR control" (short title "Dimensionality and FDR" in the panes), the title of
  the low regime moved 15 percent of the plot width to the right, the high regime's 10 percent, the
  biobank marker sits on the boundary of the high and the highest regime at p = 8e7 and n = 8e4, the
  scale of the UK Biobank with millions of variants and about a hundred thousand participants, and
  every axis title of the figures on the three surfaces starts with a
  capital letter.

#### The getting-started pages after the review

- On all three surfaces the first paragraph follows Section II of the EUSIPCO 2026 paper in four
  sentences without equations: the inputs and the random experiments in which the real variables compete
  with computer-generated dummies that act as flagged nulls, the forward selection stopped at the T-th
  dummy and the vote, the conservative estimate of the false discoveries and the calibration over
  dummies, stopping time and voting level with its finite-sample guarantee for a broad class of solvers.
  "Linear" is gone from the opening sentence, and the paragraph states the scope instead: the selector
  applies to low- and high-dimensional designs and is tailored to the high-dimensional problem, the
  scheme is not tied to the linear model and scales to millions of variables, while the library
  currently covers the supervised linear model with a response vector and a design matrix. The name is
  written out first, "Terminating-Random Experiments (T-Rex) selector", hyphenated as in the papers,
  on page 1 and the landing pages. The sentence on the predecessor links the last release before the
  rewrite, version 1.0.0 in the CRAN archive on the R page and version 0.6.18 on PyPI on the Python
  page. The selector table says "Independent or modestly dependent variables", since the plain selector
  holds up to an autocorrelation of about 0.7 in an AR(1) design. The
  ruler exempts Doxygen `\image` lines, whose caption stays on the command line, and the checker runs
  the same ruler over the guides, which are untracked until their commit and invisible to the ruler.
  "The package currently has five groups of tools",
  the memory-mapped bullet says that the algorithms are agnostic to the data representation and
  process designs larger than the random-access memory (RAM) out of core, and the utilities bullet names
  the selection metrics, normalizers, dummy distributions and threading functions in words. The table
  section is "T-Rex variants of the library" with "Screen-T-Rex" as the page of the screening rows, and
  the "Map of the pages" and "Links" sections are gone, since the landing page holds both and the left
  pane the map. The R
  page states "R (>= 4.3)", and the registry sentence of the R and Python pages carries a bracketed
  "[package page]" link after "(CRAN)" and "(PyPI)", so the sentence names the registry and the bracket
  leads to the package. The R stylesheet gives the section headings room above them (h1 44 px, h2 32
  px, h3 26 px). The three stylesheets justify the paragraphs and list items of the page content with
  automatic hyphenation (the R navigation script sets `lang="en"` on the page, which html_vignette
  omits, so the hyphenator has a language). The search index is rebuilt (107 headings, all anchors
  verified by the checker).

#### Every page spells out its abbreviations

- Each page of the three guides stands alone, so an abbreviation is written out at its first use on
  that page. A mechanical pass over the sources (headings excepted, names such as T-LARS, R6 or the
  solver classes accepted) found the gaps and the pages close them: "dependency-aware variant T-Rex+DA"
  on page 2, the signal-to-noise ratio (SNR) and the target FDR (tFDR) named before the table of page 3
  (and "signal-to-noise" hyphenated), the false discovery rate and proportion on pages 4, 6, 7 and 8 of
  the R vignettes and page 8 of the other two guides, the singular value decomposition and "sparse
  principal component selector" on page 12. "L2 norm" stays as notation, and "the ML regressions" on
  page 13 stays, since the first paragraph of the ML methods page introduces the term. The checker
  (`notes/handover_2026-11-07/tools/check_vignettes.py`) also verifies links, anchors, images,
  bibliographies and citations of the rendered R pages.

#### DOIs in every bibliography, an OpenMP page, and the regime map with the rule p = 1000 n

- Twenty entries of `references.bib` gain their DOI, each resolved at Crossref and matched by title
  (the six conference papers of the method family, LARS, the lasso and its dual, the elastic net, the
  bridge, ridge, OMP, gradient pursuits, SLINK, Lance and Williams, SimHash, LSH, sparse projections,
  graph kernels). Fu's title reads "Penalized Regressions" as published. Seven entries have no DOI to
  find (JMLR, PMLR, Soviet Math. Dokl., the two French papers, the AFS article, EUSIPCO 2026). The R
  and Python bibliographies printed DOIs already. Doxygen renders its citation list with its own
  BibTeX style, which prints neither DOI nor URL, so `build_docs.sh` runs `tools/add_doi_links.py`
  after Doxygen, which appends the DOI (or the URL) of every entry from `references.bib`.
- Doxygen 1.18.0 as released dies at random with a bus error in `commentscanYYlex`, one run in three
  on this guide. For a command without braces, `commentscan.l` tests `idx > 1` on the `npos` result of
  `find('{')` and reads `yytext[-2]`, two bytes before the scanner buffer. Every Markdown page starts
  with the generated `@page` command at the first byte of a fresh buffer, so the read happens on every
  run and is fatal when malloc places the buffer at the first byte of a region behind reserved address
  space, hence the fault addresses ending in `fffffe`. Upstream fixed it for 1.19.0 (commit 26e9087,
  two weeks after the release, while issue #12326 had been closed as a Homebrew problem).
  `cpp/docs/tools/doxygen-1.18.0-commentscan-npos.patch` holds that one-line fix, the local Doxygen
  is rebuilt with it through a Homebrew tap, and `build_docs.sh` runs Doxygen once instead of retrying.
- `OPENMP.md` at the repository root says how each package gets OpenMP on Windows, Linux and macOS,
  how to check, and where R's rules are written. The Parallelism vignette cites R Installation and
  Administration, Writing R Extensions, CRAN's macOS OpenMP page and that page, so its bibliography
  carries the sources, and the left pane's Reference block and the overview link the page. The
  no-copies pages cite the bindings their contract rests on: Rcpp and RcppEigen on the R page, NumPy
  and pybind11 on the Python page, Eigen and Boost.Interprocess on all three.
- The regime map draws the dashed boundary at p = 1000 n instead of p = 100 n, so the band of high
  dimensionality holds its title clear of the diagonal and the third title stays inside the frame, and
  the three dimensionality pages state the rule accordingly. The left pane of the R vignettes shows its
  entries in the text colour with the current page in teal, as the C++ tree does.

#### The C++ tree reads like the other panes, the regime map breathes, memory mapping credits its roots

- The C++ navigation tree shows the pages with the short titles of the R pane and the Python
  sidebar, hides the section entries under them, and draws the rule before the bibliography, all set
  by the navigation script on the tree that Doxygen builds after the page loads. "Screening" reads
  "Screen-T-Rex" in the page lists of the three guides.
- The regime map names its regimes "Low, High and Highest dimensionality", keeps more room between
  the third title and its keywords, and writes "financial factors" above its marker instead of across
  the boundary. The script draws it into the three figure folders.
- The memory-mapping pages credit Boost.Interprocess, through which the library maps a file, and the
  R page the bigmemory package, on whose `big.matrix` objects the first Big T-Rex was written. A new
  diagram (`tools/figures/memory_agnostic.py`) shows on all three pages that a design in memory and a
  design in a file reach the same implementation through one mapped view. `references.bib` and its
  copy in the vignettes gain the two entries.

#### Framed figures with boxed legends, a search box for the R vignettes, one navigation standard

- The figures of the three guides are fully framed (the matplotlib style turns the top and right
  spines back on, with a 5% margin at the ends of the axes), the legends sit in white boxes, and a
  legend shared by the panels of a figure stands outside them at the right. The TPR panel of the FDR
  study keeps a margin above 1, so its line is no longer cut by the frame.
- The top bar of the R vignettes carries a search box over the titles and headings of the fourteen
  pages. The index lives in `vignette_nav.html` between two marker lines, written by the new
  `tools/build_vignette_search_index.R` from the vignette sources with pandoc's rule for heading
  anchors, so every hit links to its section. Headings only, which keeps the index at a few kilobytes
  per page. NEWS names the box.
- The current page is marked the way the C++ guide does on all three surfaces, teal, bold and a bar
  at the left, without a background. Page entries of the sidebars carry short titles everywhere (the
  Python toctree gets explicit short titles as the R pane has them), and long section entries wrap
  instead of ending in an ellipsis in the C++ tree. The C++ title block gets more room, so the search
  box no longer touches the first tree entry.
- `AUTHORS.md` lists the members with their contact and the roles by package, which differ: for the
  R package Jasin Machkour is the package and IP owner, Fabian Scheidt the main contributor, Simon
  Tien a past contributor; for the Python package Arnau Vilella Pique and Fabian Scheidt are the
  owners, Arnau Vilella Pique proposed the 0.x versions and Fabian Scheidt is the main contributor
  and owner from 1.0.0 on; for the C++ library Fabian Scheidt is the IP owner, package owner and sole
  contributor so far. Daniel P. Palomar and Michael Muma are IP owners of all three.

#### Boxed, sharp and legended figures, and a rule in the sidebars

- The ggplot2 figures of the R pages use theme_bw, so the panels are boxed, and the vignettes render
  their figures at twice the resolution (fig.retina = 2), so the text is sharp on high-resolution
  screens. Every figure of the three guides carries a legend: the two kinds of variables and the
  voting threshold on page 3, the trial mean and the target line in the FDR study, the true support
  and the other entries of the loadings on page 7, the dummy inclusions on page 9 and the CV curve on
  page 10. The sidebars of the Python and the C++ guide draw a rule between the pages of the guide
  and the look-ups, as the left pane of the R vignettes does.

#### The remaining figures follow the pipeline and the guides say Bibliography

- The R pages 3, 7 and 10 draw their figures with ggplot2 like page 9 (the relative occurrences with
  the voting threshold, the FDR and TPR study, the two sparse loading vectors, the ridge CV curve with
  the two choices), so base graphics are gone from the vignettes. The C++ examples of the same pages
  write their figure data into `examples/output/`, and `make_figures.py` draws the four figures with
  the style of the Python guide, so the C++ pages 3, 7 and 10 show what the R and Python pages show.
  The literature list is called "Bibliography" on all three surfaces, in the closing section of the
  thirteen vignettes, in the title of the Python bibliography page and in Doxygen's citation list,
  while "Reference" stays the word for the classes and functions. The READMEs keep their "References"
  headings.

#### The T-solvers page draws the terminated path, the first figure of a pipeline for the three guides

- Page 9 of the three guides shows the coefficient paths of one warm-started T-LARS run: a snapshot
  at T_stop = 1 and the continuation to T_stop = 3 side by side, true predictors in teal, other
  predictors grey, dummies red and dashed, a dotted line at every dummy inclusion and the entering
  variable over each step, after the Part 5 figure of the Examples' T-LARS demo. The R page draws
  with ggplot2 (new in Suggests), the Python page with matplotlib under `docs/guide/matplotlibrc`,
  which gives every figure of the guide one style (brand colour cycle, no top and right spines, a
  light grid, constrained layout, 120 dpi), and the C++ example writes the two snapshots into
  `examples/output/`, from which `cpp/docs/figures/make_figures.py` draws the PNG with the same style.
  `build_docs.sh` runs the script when a Python with matplotlib is found (`TREX_DOCS_PYTHON`, default
  `python3`) and keeps the committed PNG otherwise. The figures of the Python pages 3 and 10 take the
  cycle colours instead of hard-coded red and grey, and the ridge plot gets a legend. The dark block
  of the R stylesheet addresses the code tokens as descendants of the code element, since pandoc wraps
  every line in a span, so the earlier token colours never applied and the dark pages showed pandoc's
  light palette; the dark tokens now use the palette furo uses. The brand in the C++ top bar is white
  again (doxygen-awesome colours links with !important).

#### The three guides share the running bar, centered page numbers and neutral headings

- The R top bar with the brand, the guide's name and "Page n of 13" runs over the Python and C++
  guides too, built by their navigation scripts, with the sidebars and the content starting below it
  (narrow screens keep the themes' own headers). The position "n / 13" of every pager sits in the
  center between the previous and next links on all three surfaces, and the bars' position reads in
  the brand's teal. The R headings drop the teal and use the text colour like the other two guides.
  The left pane of the R vignettes names its blocks "User guide" and "Reference", the latter with the
  reference manual on CRAN, the help pages of the installed package (left out when the page is served
  from CRAN) and the sources on GitHub, and the overview's Reference section lists the same. The C++
  main page writes its three section headings as HTML, so Doxygen's tree no longer lists them as
  pages next to the guide pages. The Python landing page is titled "The trexselector Python package"
  like the R overview and the C++ main page, and the guides' names in the bars read "R package guide",
  "Python package guide" and "C++ library guide".

#### The landing pages cite the dependent-variables paper with its DOI

- The C++ main page gave the article number and DOI of "High-dimensional false discovery rate control
  for dependent variables" as 109982, which resolves to an unrelated paper, and the R overview and the
  Python landing page had taken the line over. All three say 109990 now, as `references.bib`, the
  DESCRIPTION, the CITATION file, the READMEs and the rendered reference blocks do. Every DOI in the
  tree (nineteen, all from `references.bib`) was resolved at Crossref and matches its title, and the
  arXiv identifier of the virtual-dummies paper matches too. The only other DOI-shaped string is the
  example inside a comment of the vendored IEEE citation style.

#### AUTHORS.md names Team T-Rex

- `AUTHORS.md` at the repository root says who Team T-Rex is: the members with their roles as the
  package metadata states them (the R DESCRIPTION, the Python pyproject) and their e-mail addresses,
  the contributor to the R package, and the contact points, the discussions page of the GitHub
  organization, the issue tracker, the organization and the website. The Links of the three landing
  pages and of the READMEs point to it, and the Python landing page gains the two main articles and
  the Links block that the R overview and the C++ main page carry.

#### The three user guides open with an overview and number their pages

- The R vignettes gain an overview, `TRexSelector.Rmd`, which `vignette("TRexSelector")` opens: the
  package in one paragraph, the two main articles, the numbered list of the thirteen pages, the
  reference pointers and the links, as the index pages of the Python and C++ guides have them. The
  navigation header links the overview from the brand of the top bar, lists it above the pages in the
  left pane and as the previous link of page 1, and the overview's own pager leads to page 1. The C++
  pager gets the same previous link to the main page and a next link on the main page. The Python
  guide numbers its toctree, so the sidebar and the previous/next links carry the page numbers 1 to
  13, a script puts the position "n / 13" between the two links as the other two guides show it, the
  heading stays without a number, and the bibliography moves to the Reference group.
  The vignette stylesheet gains a block for a dark system appearance with the palette of the C++
  guide's dark mode, so the R pages follow the reader's setting as the Python and C++ guides do (the
  plots keep their white boxes). `_pkgdown.yml` lists the overview and the current numbering (it
  still had twelve pages in the old order), and the READMEs of the two guides, the Python README and
  the repository README count thirteen pages. The copyright line of the three guides carries the
  year, "Copyright © 2026 Team T-Rex", and the notice at the head of the three LICENSE files follows
  the form of the GPL's own
  appendix, "Copyright (C) 2026 Team T-Rex". NEWS names the overview.

### 2026-11-07

#### A vignette on FDR control and dimensionality opens the R series

- The new vignette `v02-dimensionality.Rmd` follows the getting-started vignette. It states the variable
  selection problem, defines the FDR and the TPR with Benjamini and Hochberg as the source, walks
  through the three regimes (low-dimensional with n at least p, high-dimensional with p above n,
  highest-dimensional with p far above n) and what each allows, shows them in one map of the (p, n)
  plane with the diagonal p = n and the line p = 100 n, and runs the selector in a ten-trial study per
  regime (n = 100, p = 50, 200 and 1,000): the FDR stays below the target in every regime while the TPR
  falls with p. The former vignettes 02 to 12 are 03 to 13, the navigation list, the README's list and
  the map of the first vignette follow, and NEWS names the vignette. The seeded dummy strategy was
  tried for the third regime and taken out again, since it costs four to nine times the standard
  strategy per run at this size on all three surfaces (3.6 to 3.8 s against 0.4 to 0.85 s), a
  performance item for the library.

#### The bibliography gains Benjamini and Hochberg 1995 and a script draws the regime map

- `references.bib` and its copy in the R vignettes gain `benjamini1995fdr`, the 1995 paper that
  introduced the false discovery rate, with the fields from the Crossref record. The map of the three
  regimes that the new dimensionality page shows in all three guides is drawn by
  `tools/figures/dimensionality_map.py` with matplotlib in the package colours, one hue from light to
  dark for the regimes, the boundaries p = n and p = 100 n, and three typical settings as markers.

#### The vignettes end with the copyright line

- The navigation script of the vignettes appends "Copyright © Team T-Rex" under the previous/next links of
  every page, styled like the pager, so the R vignettes close the way the Python and C++ guides do. NEWS
  mentions the line.

#### The READMEs quote the Eigen acknowledgement and the story names both regimes

- The Eigen citation of the Acknowledgements section stands in a block quote in the four READMEs
  (repository, C++, Python, R, the R Markdown knit), as the papers do, so the sentence about the library
  and its dependencies reads on its own. The "Short Story" of the repository README says that the T-Rex
  selector selects variables in low- and high-dimensional settings alike and is tailored to the
  high-dimensional case, where the number of variables exceeds the number of observations.

#### The C++ user guide: thirteen Doxygen pages with compiled examples

- `cpp/docs/` holds the user guide of the C++ library with the content of the R vignettes and the Python
  guide, page for page: getting started, FDR control and dimensionality, the T-Rex selector, T-Rex+DA,
  T-Rex GVS, screening, T-Rex SPCA, T-Rex Tikhonov, the T-solvers, the ML methods, memory mapping,
  working without copies and parallelism. The code on a page comes from an example program under
  `examples/`, which CMake compiles
  against the installed library and runs at build time (`build_docs.sh`), with the output of every
  documented section captured into a text file (`doc_output.hpp`) and shown under the code through the
  `\output` alias, so every output on a page is the output of the code above it. Doxygen 1.18 renders
  the pages as subpages of the main page, the master bibliography through `\cite`, MathJax formulas, the
  two figures, and the API reference of all headers with class diagrams (640 HTML pages, no warning),
  in the doxygen-awesome theme (vendored under `vendor/`, MIT) with the package colours in `custom.css`.
  The version in the header comes from the library's CMake project, the footer template adds
  "Copyright © Team T-Rex" to every page, and a small script (`pages_nav.js`) appends previous/next
  links with the position in the guide to every page, since Doxygen has none for pages. The tree's
  synchronisation toggle is hidden, it was mistaken for a page navigation. The vendored theme is exempt
  from the line ruler as a third-party file. The guide depends on the headers' comments, see the entry
  below on their parameter documentation.

#### The Python user guide: thirteen executed MyST-NB pages under Sphinx

- `Python/trexselector/docs/` holds the user guide of the Python package: the thirteen R vignettes ported
  one to one as MyST notebooks in Markdown form, whose cells Sphinx executes at build time
  (`nb_execution_raise_on_error`, two threads, data generated with seeds, files in temporary folders),
  the API reference in four pages from the docstrings (selectors, T-solvers, ML methods, utilities,
  with `conf.py` hooks that shorten the pybind11 signatures and join the reference lists of the
  docstrings), a site-wide References page from the master bibliography (`sphinxcontrib-bibtex`,
  per-page bibliographies clash on shared keys), and the furo theme in the package colours. The landing
  page links the C++ library, the R package on CRAN and its vignettes, the footer names Team T-Rex. The
  toolchain is pinned in `requirements.txt`, the full build with the execution of all pages takes about
  25 s.

#### The header comments document the current signatures

- Doxygen over `cpp/src` reported 232 warning lines, almost all `@param` entries of the T-solver
  constructors that still described the positional parameters from before the `SolverOptions` struct
  (normalize, intercept, verbose, scaling_mode), plus incomplete parameter lists in the SD family, the
  clustering dispatcher, the elastic-net cross-validation, the data generators and the normalizer, a
  few `#` and `<` sequences Doxygen read as link requests, two `@file` names that did not match their
  files, and two `@enum` lines for enumerators that are not enums. All 45 headers are corrected in
  their comments only (no code line changed), Doxygen reports no warning now, and the vendored copy of
  the R package is re-synced (vendor check OK).

#### The license headers and the citations name Team T-Rex

- The GPL headers of the repository, the C++ library and the Python package named a research group
  with a year as the copyright holder, and the package citations (`inst/CITATION` of the R package and
  `CITATION.cff`) named a "TRexSelector Team". All of them say Team T-Rex now, the name of the group of
  authors across the institutions, without a year that would need updating, and so does the CPack
  vendor of the C++ build. The citation generator `tools/references/make_inst_citation.py` writes the
  new name, and the R package's `LICENSE.md` has no holder line, so it is unchanged.

#### The Python README tables the dummy strategies and takes its badge from shields.io

- The dummy strategies of the L-loop are a table in the Selectors section, as in the C++ README (the six
  strategies with their dummies and storage), after a sentence that introduces K, L and p. The PyPI
  badge comes from shields.io and links to the project page on PyPI, since the badge.fury.io image did
  not render on GitHub and left its alternative text "PyPI version" in the page.

#### The R README tables the dummy strategies and carries the CRAN badge

- The dummy strategies of the L-loop are a table in the Selectors section, as in the C++ README (the six
  strategies with their dummies and storage), after a sentence that introduces K, L and p, and a CRAN
  version badge from r-pkg.org under the title links to the package page on CRAN, as the Python README's
  badge does for PyPI. README.md is knit from README.Rmd.

#### The repository README names every selector and fences the license

- The "Implemented Variants" table of the repository README named `TRexSPCA` for the class every
  surface calls `TRexSPCASelector` and lacked `TRexTikhonovSelector`. Both rows are right now, and
  the GPL paragraph sits in a text fence, as in the three package READMEs.

#### The C++ README describes the library

- The C++ README kept its developer sections (dependencies, building, installing and consuming) and
  gains the opening of the other READMEs, the one-line summary, a lead paragraph and the two main
  journal articles in a block quote, with the repository's references section for the rest. The
  install-and-consume section follows the build section now, so the build flow is contiguous, and a
  new section "What the Library Provides" lists the augmented-dummy selectors with their control
  structs, names the virtual-dummy family with its paper and the sparse-dummy family as future work,
  the twenty T-solvers in four families, the ML methods and the utilities by source directory. The
  dummy strategy and thread sections stay in full, as the reference the bindings rely on, a semicolon
  leaves the thread section, and the GPL paragraph sits in a text fence.

#### The Python README describes the package

- The README of the Python package, the long description on PyPI, opened with the field labels of
  DESCRIPTION as body text, named the current version in prose, attached the OpenMP facts to the source
  build although every wheel bundles the OpenMP runtime of its platform, copied the dummy strategy and
  thread sections of the C++ README verbatim, pointed to "the R vignette page 2" for documentation and
  had no overview of the package. It follows the R README now: the title "trexselector — Python
  Package", the lead paragraph with the two main journal articles in a block quote, installation bullets
  (wheels for three platforms and Python 3.10 to 3.15 with the bundled OpenMP runtime, the source build
  with the toolchain's OpenMP and Homebrew's LLVM on macOS, `main` as the released branch, the
  single-threaded and the development install), the quick start, an overview of the seven selectors
  with their control classes and the dummy strategy, the twenty T-solvers in the four submodules of
  `trexselector.tsolvers.ad`, the ML methods, the utilities and the thread options, a documentation
  section (docstrings, stub files, the Examples repository, the R vignettes), the links (website,
  repository, Examples, references, CITATION.cff, issues) and the GPL paragraph in a text fence.

#### The vignettes share a navigation header

- The `css` option of `html_vignette` replaces rmarkdown's whole default stylesheet, and the package's
  stylesheet held two heading colours, so code blocks, tables and the table of contents of the twelve
  vignettes were unstyled. The stylesheet takes rmarkdown's default as its base and adds the package
  colours, shaded and bordered code blocks as in rendered Markdown, a table of contents that stays on
  the left on wide screens and follows the scroll, and the styles of a shared header. That header,
  `vignette_nav.html`, is included before the body of every page and carries a dark top bar with the
  page position, the list of the twelve pages under the table of contents with the current one marked,
  and previous/next links at the end of the page, all from one page list, so a new vignette is one list
  entry. The format stays `html_vignette`: the `html_document` format of the 1.x vignette embeds
  Bootstrap and weighs 1.2 MB per page, while the twelve pages weigh 995 KB together, 108 KB more than
  before, and build in 15 s.

#### The R README describes the package

- The README of the R package opened with the field labels of DESCRIPTION (`Title`, `Description`,
  `Papers`, `Note`) as body text, carried the quick start as a level-1 heading with the other sections
  below it, and copied the dummy strategy and thread sections of the C++ README verbatim without any
  overview of what the package contains. It is titled like the C++ README now ("TRexSelector — R
  Package") and opens with a lead paragraph and the two main journal articles in a block quote, with
  the repository's references section for the rest. The installation facts are bullets that say where
  OpenMP is available (the Windows binary, the Linux build, a macOS build from source with Homebrew's
  libomp, for which the commands follow, while CRAN's macOS binary runs on one thread), and the
  GitHub passage states that `main` carries the released version and `dev_beyond` the development.
  The quick start is a level-2 section, and a new section lists the seven selectors with their
  control functions and the dummy strategy, the twenty T-solvers in their four families, the ML
  methods, the utilities and the thread options, with pointers to vignettes 2 and 12. The vignette
  list names all twelve, the links carry the project website, the repository, the Examples repository,
  the references section and the issues, the vignette link is absolute (CRAN renders the README outside
  the repository), and the GPL paragraph sits in a text fence. The getting-started vignette carries the
  same branch and OpenMP statements, names the utilities as the fifth group of tools and ends with the
  same links, and a table cell of the T-solvers vignette loses its semicolon. NEWS.md and
  cran-comments.md describe the reorganised README.

#### The submission comments carry the win-builder results of 2.0.1

- `cran-comments.md` lists win-builder on R-release 4.6.2 and on R-devel r90845 among the test
  environments of 2.0.1, both with the incoming-feasibility NOTE alone, which names the one day since
  the last update, and nothing else; the R-devel instance, which returned no result on 2.0.0, answered
  this time.

#### The Python package is 1.0.1

- The 1.0.0 source distribution on PyPI carries the sources before the ranges fix, so a build from it
  fails with Apple clang 14 as the CRAN builders did, while the wheels, built with newer compilers,
  are unaffected. The Python package follows into the patch release 1.0.1: pyproject.toml carries
  the version, and the GitHub release v1.0.1 publishes the wheels and the fixed sdist through the
  existing workflow.

#### The C++ SDK is 1.0.1

- The ranges fix changes the sources every package ships, so the C++ SDK follows the R package into
  a patch release: the CMake project version is 1.0.1, which `find_package(TRexSelector 1.0)` still
  accepts under the SameMajorVersion rule, and CITATION.cff names 1.0.1 as the version of the
  library and the Python package. HISTORY opens the release block above these entries, and the
  1.0.0 block records its release dates.

#### The sources compile without warnings of their own

- CRAN's install logs of 2.0.0 and a local build without OpenMP under `-Wall -Wextra` listed the
  warnings the project's own sources produce: the flag `do_parallel` of the three experiment loops
  is unused when the OpenMP pragma that reads it is compiled out; the Tikhonov selector built its
  group-centering triplets from `Eigen::Index` values, which MSVC narrows to the `int` StorageIndex
  with C4244; three loops of the R bindings compared signed and unsigned sizes; the VD and SD
  selectors take a thread count their sequential build does not use; and VD_Base listed its
  generator before the data views it is declared after. Each site carries the attribute, the cast
  or the order the compiler asked for, with no change in behavior, and the vendored copy of the R
  package follows. The remaining warnings in those logs come from RcppEigen's headers.

#### The configure script reads R's OpenMP flags from Makeconf

- On R 4.6 and R-devel, `R CMD config SHLIB_OPENMP_CXXFLAGS` answers "no information for variable"
  although Makeconf defines it, and the first step of the R package's configure relied on that query,
  so every Linux build on those versions compiled sequentially: CRAN's Fedora builds of 2.0.0 and the
  rocker R-devel container both report "checking for R's OpenMP flags... none", while R 4.3.3 still
  served the value. The configure script reads SHLIB_OPENMP_CXXFLAGS from `${R_HOME}/etc${R_ARCH}/Makeconf`
  and, when it is set, writes the make variable `$(SHLIB_OPENMP_CXXFLAGS)` itself into Makevars, the
  form Writing R Extensions asks for, as Makevars.win already does. Where R has no OpenMP flags, CRAN's
  R on macOS, the steps for Homebrew's omp.h with R's bundled libomp follow unchanged. The configure
  script carries the package version 2.0.1.

#### The R package is 2.0.1

- The ranges fix reaches CRAN as a patch release, since CRAN takes no second upload of a version
  it holds. DESCRIPTION carries 2.0.1, NEWS.md describes the release, inst/CITATION names the
  version, and cran-comments.md states the reason for the update, the test environments, with the
  clang 14 and libc++ 14 container as the closest match to the failing builders, and the check
  results.

#### The solvers compile without the ranges library

- CRAN's daily checks of TRexSelector 2.0.0 reported "Installation failed" on three macOS flavors
  (r-release x86_64, r-oldrel arm64 and x86_64), whose builders run Apple clang 14.0.3 with the
  macOS 11.3 SDK. That compiler accepts `-std=gnu++20` but its standard library has no ranges
  algorithms, so `std::ranges::find` in the column preprocessing and in the inactive-set
  initialization and `std::ranges::sort` of the dropped indices failed with "no member named
  'ranges' in namespace 'std'". The arm64 release builder, on Apple clang 17, passed. The three
  calls are `std::find` and `std::sort` over begin and end now, with identical results, and the
  vendored copy of the R package follows. The tree uses no other C++20 library feature that
  Apple clang 14 lacks.

#### The getting-started vignette installs from CRAN first

- The installation passage of the getting-started vignette still said that CRAN carried the pure-R
  1.x release and sent the reader to GitHub until 2.0.0 was accepted. Like the README since the
  previous commit, it opens with `install.packages("TRexSelector")`, states that CRAN provides the
  Windows and macOS binaries while Linux and every source install compile the C++ backend, and
  presents the GitHub route as the development version, without a current version number.

#### The development versions follow the release

- With TRexSelector 2.0.0 on CRAN and trexselector 1.0.0 on PyPI, the sources carried the released
  version numbers, so an install from the repository called itself the released package. The R
  package is 2.0.0.9000 in DESCRIPTION with a "development version" heading at the top of
  NEWS.md, the R convention after a CRAN release, and the Python package is 1.0.1.dev0 in
  pyproject.toml. Both sort after the released versions and before the next release, whose
  commit sets the final numbers. The C++ SDK version, inst/CITATION and CITATION.cff keep the
  released numbers, since they name what is published.

#### The R README installs from CRAN first

- The installation section of the R package's README opened with the compiler requirement and
  the GitHub route, and closed with a note that CRAN still carried the pure-R 1.x release while
  the 2.0.0 submission was in progress. With 2.0.0 on CRAN the section opens with
  `install.packages("TRexSelector")`, states that CRAN provides the Windows and macOS binaries
  while Linux and every source install compile the C++ backend, and presents the GitHub route as
  the development version. The current version number stays out of the text, so the section
  needs no edit per release. README.md is knit from README.Rmd again, and the root README notes
  that the R package is distributed on CRAN.

## TRexSelector C++ 1.0.0 / R 2.0.0 / Python 1.0.0

Release candidate on `dev_beyond` (2026-09-24); the three packages carry
their own version numbers from here on: the C++ SDK 1.0.0 (CMake project
version, `find_package(TRexSelector 1.0)`), the R package TRexSelector
2.0.0, the Python package trexselector 1.0.0 (the names since 2026-11-04, see
that entry; TRexSelectorNeo before). Released: TRexSelector 2.0.0 on CRAN on 2026-11-06, the GitHub
release v1.0.0 and trexselector 1.0.0 on PyPI on 2026-11-07.

### Breaking changes since 2.0.0 (June 2026)

- L-loop strategies: `ONDEMAND` and `PERMUTATION_ONDEMAND` are `SEEDED` and
  `PERMUTATION_SEEDED` (the seed-derived, nothing-stored strategies); the old strings
  are rejected by `trex_control()` and the Python enum has no alias.
- `SolverTypeForTRex` was renumbered when the CD family and the Tikhonov solvers
  were added; serialized integer codes from before July are not portable.
- CD solver family renamed: `TCCD`, `TCENET`, `TCIENET` are `TCDLASSO`, `TCDENET`,
  `TCDIENET`; the general-Tikhonov mode of `TCIENET` was extracted into
  `TCDTIKNET`, which is the Tikhonov selector's default solver. Checkpoints written
  by the old classes are unreadable by the new ones (formats incompatible).
- The CD family localizes every dummy crossing with certificates by default
  (`cd_certified_localization = TRUE`); the pre-August probe bisection is the
  LEGACY fallback and opt-in.
- `TRexGVSSelector`, `gvs_type = "IEN"`: the default solver is the native pathwise
  `TIENET`; `TIENET_AUG` and `TCDIENET` are the alternatives.
- `TRexTikhonovSelector` defaults: cluster-aware dummies (`cluster_dummies = TRUE`)
  and FOLDED dummy coupling (`fold_dummy_coupling = TRUE`); the `PILOT_MAX` /
  `PILOT_KNEE` lambda_2 rules exist next to the CV tuners (default `CV_1SE_IEN_CCD`).
- The SDK headers install under `<prefix>/include/TRexSelector/` (manual include
  paths change; exported targets do not), see the 2026-09-24 entry.
- The C++ tree has one folder per dummy mechanism: `ad/` (augmented dummies, the
  classical family), `vd/` (virtual dummies), `sd/` (sparse dummies), in both
  `tsolvers/` and `trex_selector_methods/`, and the namespaces follow:
  `trex::tsolvers::ad::linear_model::lars_based::TLARS_Solver`,
  `trex::trex_selector_methods::ad::trex_core::TRexSelector`. Include paths and
  namespaces of every C++ user change; the R and Python surfaces do not. See the
  2026-09-24 entry.
- The selectors reject non-finite entries and constant columns of X (constant
  columns were silently unselectable before), see `ScalingMode::NONE` below.
- Every T-solver is constructed as `(X, D, y, [groups | K,] options)` with one
  options structure; the positional constructors and the per-solver R / Python
  keyword arguments are gone, see below.
- Floating-point classification is bit-based and the Release build no longer uses
  `-ffast-math`: results are IEEE 754 conformant; runs from fast-math builds are
  not bit-reproducible against these.
- The T-Rex selectors mark their centered design `precentered`: on a saturating
  design the solvers stop at n - 1 active variables (the rank of the centered
  design) instead of taking an n-th, numerically collinear step. T-solver
  checkpoints written before 2026-09-27 are unreadable (the base state and T-OLS
  gained fields). See the 2026-09-27 entry.
- T-NCGMP in its line-search variant (the default) stops by the T-MP rules (the
  row bound on the active atoms, 8 x min(max_actives, row bound) steps) instead of running without
  an active-set cap for 20 x (p + L) steps; past saturation its paths are shorter
  and equal T-MP's. The augmented EN wrappers count the intercept in their
  degrees of freedom for lambda2 > 0 (one more than before). See the 2026-09-27
  entry.
- `SolverOptions::tie_seed` defaults to 0 (C++, R `tsolver_options()`, Python
  `SolverOptions`): a solver without an explicit seed breaks exact ties
  reproducibly; a negative seed keeps the former std::random_device behavior. The
  T-Rex selectors are unchanged (seeded selectors derive per-experiment seeds,
  unseeded ones pass -1).
- `TIKOMP`, the Tikhonov-K orthogonal matching pursuit, is removed from the
  C++, R and Python surfaces before the release: the class and its
  `SolverTypeOMPBased` value, `SolverTypeForTRex::TIKOMP`, the R6 class
  `TIKOMP_Solver` and the `"TIKOMP"` solver name, the Python class and enum
  member. The `SolverTypeForTRex` enumerators after it (TGP to TAFS) shift
  down by one, and TIKOMP checkpoints are unreadable. See the 2026-10-26
  entry.
- A T-solver used directly runs its OpenMP regions serially by default
  (`SolverOptions::n_threads = 1`, a whole number of at least 1); before, the
  regions took every core OpenMP offered. `EigenSingleThreadGuard` in
  `tsolver_base.hpp` is `EigenThreadBudgetGuard(n_threads)`. See the
  2026-10-26 entry.
- The selectors run serially by default: the random experiments in sequence
  and every solver on one thread. `TRexControlParameter::max_outer_threads`
  and `max_inner_threads` are honoured as set (defaults 1 and 1, values below
  1 throw) where `autoConfigResources()` used to overwrite them with the core
  count; that method is gone, and so is the `EigenThreadScope` class of the
  experiment runner. In R, `trex_control(use_openmp = )` is
  `parallel_rnd_experiments` (the C++ and Python name), and the two thread
  options now reach the backend. `TRexSDOptions::n_threads` and
  `TRexVDOptions::n_threads` default to 1 instead of "auto" (0 is rejected).
  There is no "every core" value: whoever wants every core reads the count
  and passes it. See the 2026-10-27 entries.
- `LLoopStrategy::SEEDED` is the storage-free twin of `STANDARD` (fresh dummies
  at every L step, bit-identical to STANDARD for the same seed) instead of the
  storage-free twin of HCONCAT in the i.i.d. case. A seeded SEEDED run, and a
  seeded screening run with STANDARD or SEEDED, selects differently from before
  under the same law; unseeded runs are unaffected. See the 2026-10-27 entry.
- T-Stagewise takes fewer steps: a variable its NNLS reduction removes within a
  step is no candidate for that step's move, so the rounding-noise re-entries
  are gone and the path no longer depends on the build. Recorded step counts
  change; coefficients and supports below saturation do not. See the 2026-10-27
  entry.
- R: `MLRidgeCV$fit(num_folds = )` is `n_folds`. C++ `ridge_cv_svd::fit()` and
  Python `RidgeCV.fit()` default to `n_lambda = 100` instead of 1,000 (R's
  default, the glmnet convention). See the 2026-10-27 entry.
- R defaults aligned with C++ and Python: `trex_control()` defaults to
  `lloop_strategy = "STANDARD"` (was HCONCAT), leaves `tloop_stagnation_stop`
  to the C++ solver rule (TIENET and TIENET_AUG no longer stop by default) and
  sets `cd_lambda_rel_tol = 1e-6` (was 1e-5). C++: `TRexGVSSelector` resolves
  the default lambda_2 rule to CV_1SE_IEN_CCD on the IEN track, as the
  bindings did. R selections that relied on these defaults change. See the
  2026-10-27 entry.
- T-Stagewise stops when the support holds T dummies, frozen coefficients
  included; before, only the active dummies counted and the path ran on after a
  removal. Selections and step counts of T-Stagewise runs change. See the
  2026-10-27 entry.
- The FDP estimate keeps negative calibrated occurrences `Phi_prime` (only the
  upper clamp at 1 remains): FDP_hat can only rise, so the L-loop may calibrate
  a larger L, the T-loop may stop earlier and selections can only shrink.
  Results move where a selected variable had Phi' < 0, in long T-loops at
  tFDR 0.3 and above. See the 2026-10-24 entry.
- `TRexSPCASelector` rejects `gvs_type` IEN (`std::invalid_argument`; R
  `trex_spca_control()` stops, Python raises `ValueError`) where the core
  silently ran the elastic net. See the 2026-10-22 entry.
- Two guards throw early with a remedy where a run used to fail late: the
  clustering refuses a RAM buffer beyond the physical memory (`use_mmap`), and
  a memory-mapped file is refused when the disk lacks the space. See the
  2026-10-27 and 2026-10-20 entries.
- The Tikhonov selector's default `lambda2_method` is the pilot knee
  `PILOT_KNEE` (was the CV tuner `CV_1SE_IEN_CCD`), so a Tikhonov selection
  without a user lambda_2 changes. `PILOT_MAX` is removed from
  `TikLambda2Method` (R `trex_tikhonov_control()`, Python
  `TikLambda2Method`), and the power rule's `Lambda2PowerRule` enum with
  `MAX_SELECTION` is gone; the knee with tolerance 0 is the largest
  selection. See the 2026-10-28 entry.

### 2026-11-06

#### The Python README links work on PyPI

- The README of the Python package is the project page on PyPI, where its three relative links to
  the C++ core README, to CITATION.cff and to the package's LICENSE pointed into nothing. They now
  carry the GitHub URLs of those files on the main branch, which work on PyPI and in the
  repository alike.

#### The workflow uploads to TestPyPI on request

- The Python Wheels workflow had no way to rehearse a release: the publish job runs on a GitHub
  release alone, and that upload is final. The manual trigger of the workflow gains a `target`
  input with the choices `none` (the default, a build without upload) and `testpypi`. With
  `testpypi` chosen, a job after the build jobs uploads the same wheels and the sdist to TestPyPI
  through the `testpypi` environment with OIDC trusted publishing, the same mechanism the PyPI job
  uses. No release tag exists in a manual run, so that job skips the tag check and uploads the
  version pyproject.toml carries, which TestPyPI accepts once per version.

#### The publish job checks the release tag against the artifacts

- The publish job of the Python Wheels workflow uploaded whatever the build jobs produced as soon
  as a GitHub release was published, so a release under a mistyped tag or on a commit whose
  pyproject.toml still carried an older version would have reached PyPI. A step before the upload
  derives the version from the release tag, which has the form `v<version>`, and requires every
  wheel and exactly one sdist in `dist/` to carry that version. Any other tag or file fails the
  job and nothing is uploaded.

#### The Python package declares itself stable

- The Development Status classifier of pyproject.toml still said "4 - Beta" from the first
  TRexSelectorNeo wheels, while the package ships as 1.0.0 with the full test suite, the wheels of
  three operating systems and the sdist. The classifier says "5 - Production/Stable", the state
  the version number claims. Nothing else changes.

#### The sdist carries the C++ library

- A source distribution of the Python package could not build, because its CMake reached the
  library through `../../cpp`, which an archive does not contain, so PyPI could have held wheels
  only. The sdist takes `cpp/CMakeLists.txt`, `cpp/cmake`, `cpp/src`, `cpp/LICENSE` and
  `cpp/README.md` from the repository through scikit-build-core's `sdist.force-include` (1.0 is
  the lower bound of the build requirement now) and places them under `cpp/` in the archive; the
  Python CMakeLists.txt uses that copy when present and `../../cpp` otherwise. The Python Wheels
  workflow gains a job that builds the sdist, installs it from source on Ubuntu with the
  header-only dependencies from apt, runs the tests and uploads the archive, and the publish job
  uploads the sdist with the wheels.

#### The Description field cites the journal articles

- The Description field of the R package named no method paper, while 1.0.0 cited the T-Rex
  arXiv preprint, and CRAN's daily check of 1.0.0 now flags that preprint form and asks for a DOI.
  The field names the two Signal Processing articles of 2025, the T-Rex selector and its
  dependency-aware extension, in CRAN's form "Authors (year) <doi:...>", as `inst/CITATION` and
  CITATION.cff already do, and its first sentence no longer says "FDR-controlled" and "controlling
  the FDR" at once.

#### The wheels cover Python 3.14 and 3.15

- The Python Wheels workflow built CPython 3.10 to 3.13 with cibuildwheel 2.21.3 of 2024, while 3.14
  and 3.15 are the current releases and numpy ships wheels for both. The workflow runs cibuildwheel
  4.3.5, builds 3.10 to 3.15, and skips the free-threaded interpreters, which cibuildwheel builds by
  default from 3.14 on and the bindings do not support. cibuildwheel 4 also makes delvewheel the
  default Windows repair, and the explicit repair command with the staged OpenMP runtime overrides
  it as before. The classifiers of pyproject.toml list 3.14 and 3.15.

#### The submission comments carry the R-devel and win-builder results

- `cran-comments.md` of the R package lists the test environments of the 2.0.0 submission: the
  R-devel container (Debian, two NOTEs of the container itself, the hardening flags of its
  Makeconf and the missing V8), win-builder on R-release 4.6.2 (the maintainer NOTE alone), and
  the Windows R-devel job of the workflow. win-builder's R-devel instance returned no result on
  two uploads, which the comments say.

#### The R Package workflow checks the daily R-devel snapshot on Windows

- win-builder's R-devel field stalled twice at "checking CRAN incoming feasibility", with the
  morning's tarball and with the fixed one, while its R-release field passed that step in 28 s, so
  the Windows R-devel half of the CRAN preparation had no result. The R Package workflow gains the
  matrix entry `R windows-latest (R-devel)`: the existing job (install, testthat suite, R CMD check)
  on the daily R-devel snapshot that r-lib's setup-r installs, with the latest Rtools. That is the
  win-builder run without its incoming step, which is platform-independent and passed on
  R-release. Posit's package manager has no R-devel binaries, so that job compiles its dependencies
  from source. R-hub's checks were the alternative and need the package at the root of a GitHub
  repository, which R/TRexSelector is not.

#### The READMEs link the GNU page of the GPL

- win-builder's incoming check of the R package flagged "(possibly) invalid file URI: ../../LICENSE"
  from README.md: the package README linked the repository's LICENSE two directories up, a file the
  tarball does not contain, and the Python README carried the same link, which an sdist cannot
  resolve either. Every README now names the GNU page of GPL-3 next to the local LICENSE where the
  package ships one (root, cpp, Python); the R README and its README.Rmd link the GNU page alone,
  since CRAN packages do not ship the GPL text.

#### The C++ SDK carries its license

- `cmake --install` wrote headers, libraries and the CMake package files and no license text, and
  the CPack packages named the repository's LICENSE as a resource file only, which the TGZ and ZIP
  generators never copy, so the source package held the cpp tree without a license. The cpp tree
  has its own LICENSE (the GPL-3 text, as the Python package has), the SDK install and the binary
  package place it with README.md under `share/doc/TRexSelector`, the CPack resource files point
  at the cpp copies, and the source package leaves the editor and tooling caches out (`.cache`,
  `.clangd`, `.clang-tidy`, `compile_commands.json`).

#### The example of trex_gamma_to_k() runs only with Matrix installed

- Matrix is a suggested package and `trex_gamma_to_k()` stops without it, but its example called
  the function unconditionally, so `R CMD check` fails wherever its example run cannot see Matrix
  (seen under R-devel in a container that hands its libraries to R through `R_LIBS`, which the
  check replaces: two ERRORs, the example and its `--run-donttest` repeat, while tests, vignettes
  and the manual passed). CRAN's machines carry Matrix as a recommended package, and CRAN's policy
  asks for conditional use of suggested packages all the same. The example is wrapped in
  `requireNamespace("Matrix", quietly = TRUE)`, as the vignette and the tests already do.

#### The macOS wheel targets macOS 15

- The macos-15 job of the Python Wheels workflow, the first run after the move from the retired
  macos-14 image, failed in delocate's repair of the cp310 wheel. Homebrew bottles carry the
  runner's macOS as their minimum version, and delocate takes the wheel's floor from
  `MACOSX_DEPLOYMENT_TARGET`, so the image's OpenMP runtime (minimum 15.0) cannot be bundled
  into a wheel for 14.0; the same rule explains the failed repair on the macOS 26 image before.
  The deployment target of the macOS build follows the runner (15.0), the wheel is tagged
  `macosx_15_0_arm64`, and macOS 14 users build from source. A lower floor needs an OpenMP
  runtime built for that version, which Homebrew does not provide.

### 2026-11-05

#### The TIKNET rank-bound tests draw their design from a platform-independent stream

- `TSolverTikTest.TiknetRankBoundDropContinuesPath` failed on every x86_64 CI job (Ubuntu with
  and without OpenMP, Windows) with "rank bound not reached": the path on its design reached 21
  degrees of freedom, one short of the bound 22. The design came from `std::srand` and Eigen's
  `Random`, so it was a different matrix per C library, and the arm64 builds, which contract
  into fused multiply-add, happened to saturate where the x86_64 builds, built for SSE2 without
  FMA, stop one step short. The fixture now fills the design from a `std::mt19937_64` stream
  with 53-bit uniforms assembled by hand, bit-identical everywhere, and the seeds were chosen
  by a scan over 60 seeds on macOS arm64, Linux arm64 and Linux x86_64: seed 45 for the drop
  test and the TIKNET/TIENET match at the bound, seed 12 for the terminal-step test. Both keep
  their path under 40 perturbations of the design at 1e-10 and 1e-8 on FMA and non-FMA
  arithmetic, while 44 of 119 seeds of this saturated 12 x 30 design change their path between
  platforms. The other `srand` sites of the test tree are unchanged.

#### The wheel workflow builds on macos-15 and ubuntu-24.04

- GitHub retired the `macos-14` runner image on 2026-11-02 (actions/runner-images#13518), so the
  macOS wheel job of the Python Wheels workflow failed before its first step with "The
  configuration 'macos-14-arm64-us-default' is not supported". The job runs on `macos-15`, the
  stable arm64 image below macOS 26 (whose image breaks delocate-wheel's repair step); the
  deployment target 14.0 was kept, which the next day's entry corrects. The Linux jobs pin
  `ubuntu-24.04` while `ubuntu-latest` migrates to Ubuntu 26.04 (actions/runner-images#14748,
  until 2026-11-19); the manylinux build runs in its container either way.

#### The C++ CI builds run one compile job per core

- The Build steps of the C++ Tests workflow called `cmake --build --parallel` without a count,
  which hands GNU make an unlimited `-j`: every test unit (about 55 Eigen-heavy compiles)
  started at once on hosted runners with 3 to 4 cores and 7 to 16 GB, and the macOS job of the
  run of 35e3951 was terminated in that phase (exit 143, every compile `Terminated`). The two
  Makefile builds now pass the online processor count (`getconf _NPROCESSORS_ONLN`); the
  Windows build keeps `--parallel`, which MSBuild bounds by the core count itself.

#### An empty coefficient lookup map counts no bucket array in memoryBytes()

- `memoryBytes()` estimates the coefficient lookup map (an `std::unordered_map`) by its bucket
  array and its nodes. A default-constructed map owns no buckets with libc++ and one inside
  the object with libstdc++, but eight on the heap with MSVC's standard library, so the
  Windows wheel reported 64 bytes more than the hand count of the Python test (3,512 instead
  of 3,448 for the TLARS reader at n = 20, P = 20); the R package on Windows builds with
  Rtools' GCC and libstdc++ and was never affected. The estimate counts the bucket array only
  while the map holds entries, so the construction counts agree across standard libraries;
  filled maps are estimated as before.

#### The clustering's verbose test flushes the C runtime on Windows

- `test_verbose_default_is_silent` (Python) reads the backend's allocation line through
  `capfd` after an `fflush` of the C stdio buffer, which it reached with `ctypes.CDLL(None)`,
  the running process on POSIX. Windows has no such handle and raised a `TypeError`, the one
  failure of the Windows wheel job. The test now loads the UCRT (`ucrtbase`) on Windows, the C
  runtime CPython and the MSVC-built extension share, and the process handle elsewhere.

#### The workflows can be started by hand

- The four workflows (C++ Tests, R Package, Python Wheels, Line ruler) carry a `workflow_dispatch`
  trigger next to push and pull request, so a run can be started from the Actions tab for any
  branch, for example after the repository's Actions were switched on later than a push. Jobs
  and path filters are unchanged.

#### The repository addresses point at the GitHub organization trexselector (D15)

- The repository moves from the account G4Lactus to the organization trexselector with its
  history, and every address in the packages follows before the move is published by CRAN or
  PyPI (GitHub forwards the old address, but CRAN's URL check reports a forwarded address as
  "Moved Permanently"): `https://github.com/G4Lactus/TRexSelector` is
  `https://github.com/trexselector/TRexSelector` in DESCRIPTION (URL, BugReports), the R README
  and its `install_github()` line, `_pkgdown.yml`, the getting-started vignette,
  `pyproject.toml` (Homepage, Repository, Bug Tracker), the Python README, CITATION.cff and the
  citation generator. The generated files (the R README.md, `TRexSelector-package.Rd`,
  `inst/CITATION`) are rebuilt from them. No change of code or behaviour.

#### The seeded generators count their draws and engine seedings toward the work threshold

- The work threshold (`kParallelWorkThreshold`, 1e6 entries) was calibrated on plain passes over
  the design, which stream an entry in about 0.2 ns. A random draw costs 15 to 18 ns and the
  seeding of a per-column engine through `std::seed_seq` 5.8 us (a word-seeded engine 0.9 us),
  so `generate_dummies()` and the AR(1) design generator, which counted plain entries, stayed
  serial up to about 15 ms of work (n = 150 with 5,000 columns: 38 ms serial at a budget of 8,
  while 7,000 columns took 7 ms) and the smallest points of a runtime curve generated serially
  where the next ones split. `dummygen::seeded_loop_work(engines, draws, engine_work)` now
  counts such a loop in the threshold's unit with `kDrawWork` = 80, `kSeedSeqEngineWork` =
  32,000 and `kWordSeededEngineWork` = 5,000 entries, so the generators split from about
  0.2 ms of serial work like every other region. The default budget stays 1 and selections are
  unchanged. The thread-budget test pins a seeding-dominated design below the entry threshold
  to a team at a budget of four while the scaler's passes over the same matrix stay serial.
- `calculate_noise_params()` (utils_noisegen.hpp) carries a comment on its variance
  denominator (n - 1 up to n = 100, n above). No change of behaviour.

#### The one-shot kernels take an explicit thread budget: single-threaded by default (WP15 gap)

- The Examples check of 2026-11-05 found the library not single-threaded by default: 28
  OpenMP regions ran without a thread clause and took the runtime's default team (every
  logical core), so a default selection used about 3 cores in Python and 2.5 in R. The data
  normalizer (4 regions), the dummy generator's row permutation (1), the dummy distributions
  (11), the synthetic data generator (1), the DA selector's correlation scan (1), the scalers
  (3) and the hierarchical clustering (7) now carry `num_threads(nt) if (nt > 1)`, with
  `nt = threads_for(budget, work)`: 1 for a budget of 1, below the work threshold
  (`kParallelWorkThreshold`, 1e6 entries, now library-wide in utils_openmp.hpp) or inside a
  caller's team, the budget otherwise. The selectors' kernels take `max_inner_threads`; the
  stand-alone kernels take an explicit `n_threads` (the scalers, the clustering, the data
  normalizer, the generators) or the library's default budget. Selections are unchanged:
  they never depended on the thread count.
- `set_num_threads(m)` is the library's default budget (initially 1) instead of the OpenMP
  runtime's count: `get_num_threads()` reads it, `set_num_threads()` sets it and returns the
  previous value, and refuses values below 1; `get_max_threads()` keeps reporting the
  machine's threads, an upper bound for a budget. The OpenMP runtime's own default team is
  never used by the library.
- The ML regressions (ridge CV, the elastic-net path and CV, the IEN path and CV, the
  Tikhonov CV, the SVD solver, PCA, the ridge solver) run their Eigen products under
  `EigenThreadBudgetGuard` on a budget of their own: `set_num_threads()` / `num_threads()` on
  the classes, an `n_threads` constructor argument for SVD and PCA, a parameter of
  `RidgeSolver::solve()`; the CV classes hand their budget to the fold fits. The guard moved
  from the T-solver base to `utils/openmp/utils_eigen_threads.hpp`. The selectors pass
  `max_inner_threads` to their lambda tuners, SPCA to its scaler, PCA and ridge step.
- R and Python: `trex_get_num_threads()` / `utils.get_num_threads()` are new, the set
  functions return the previous budget, and `n_threads` (`NULL` / `None` = the default
  budget) is an argument of the scalers, `ml_agglomerative_cluster()` /
  `agglomerative_cluster()`, the ML classes, `MLRidge$solve()` / `RidgeSolver.solve()` and
  `trex_normalize_design()` / `normalize_design()`. A team probe
  (`openmp::team_probe`, off by default) records the largest team any budgeted region ran on,
  and the new C++ test pins the default to one thread with the OpenMP runtime's count raised,
  and a budget of four to a larger team. The READMEs and vignette 12 describe the model.

#### NEWS.md and the submission comments of the R package

- `R/TRexSelector/NEWS.md` (new): the 2.0.0 entry for CRAN's news reader, the breaking change from
  the 1.0.0 functional interface to the R6 classes, and the additions. `cran-comments.md` describes
  the update as it is: the rewrite as C++ bindings, the test environments (macOS R 4.6.1, Ubuntu
  24.04 R 4.3.3, the CI jobs), the two Ubuntu NOTEs with their reasons, the Eigen header warnings,
  the OpenMP configure logic and no reverse dependencies; win-builder is marked pending.

### 2026-11-04

#### Eigen's homepage URL follows its move

- `http://eigen.tuxfamily.org` redirects to `https://libeigen.gitlab.io/` (301), which the CRAN
  incoming check reported for the R README and vignette 9 (the one NOTE of the check below). The
  master bibliography carries the new URL and every rendering follows: the Doxygen blocks of the
  SVD, PCA, ridge-CV and Tikhonov-CV headers and their vendored copies, the roxygen blocks and
  man pages, the pybind11 docstrings and stubs (regenerated byte-identical, 629 tests),
  `inst/CITATION`, CITATION.cff and the READMEs. `R CMD check --as-cran` of the tarball with
  CRAN's CPU-to-elapsed thresholds: Status OK.

#### The R and Python packages return to the name TRexSelector (WP19, D13)

- The R package is `TRexSelector` 2.0.0 (was TRexSelectorNeo): it replaces the pure-R 1.x
  releases on CRAN in place, under the same name and maintainer, so DESCRIPTION no longer
  says "alongside which it can be installed". The directory is `R/TRexSelector/`,
  `useDynLib(TRexSelector)`, the types header is `TRexSelector_types.h`, `RcppExports.*`
  carry the `_TRexSelector_` prefix and `R_init_TRexSelector`, `configure` is regenerated
  from `configure.ac`, and `inst/CITATION`, the man pages, the vignettes, the README, the
  tools and the CI workflow follow: `library(TRexSelector)`, `citation("TRexSelector")`.
- The Python package is `trexselector` 1.0.0 (was TRexSelectorNeo / `trex_selector_neo`),
  distributed on PyPI and imported under that one name. 0.6.18 on PyPI is the pure-Python
  predecessor it supersedes, so 1.0.0 continues its numbering. The directory is
  `Python/trexselector/`, the package `src/trexselector/`, and the stubs, the tests, the
  wheels workflow and the README follow.
- The repository metadata carries the new names: CITATION.cff, the root and C++ READMEs,
  the validation scripts, the OrbStack scripts, the line ruler's path list, the VS Code
  tasks and `.gitignore`. Past entries of this file keep the names of their time.
- Verified at this state on macOS: R 1,552 expectations from a scratch install, Python
  629 tests on a fresh wheel with the twelve stubs regenerated byte-identical,
  `configure`, `RcppExports.*` and `inst/CITATION` regenerated identical to the
  substitution, vendor check and line ruler clean. `R CMD check --as-cran` of the tarball
  with CRAN's CPU-to-elapsed thresholds ends with one NOTE that is not from the rename:
  the incoming check reports Eigen's homepage `http://eigen.tuxfamily.org` as moved to
  `https://libeigen.gitlab.io/` (the R README and vignette 9), fixed in the entry above. The
  incoming check accepts 2.0.0 against the pure-R 1.0.0 on CRAN (same maintainer, no "new
  submission" NOTE).

### 2026-11-01

#### SKIPL warns under a large fixed L without the stagnation stop (WP21)

- `TRexSelector` warns when `lloop_strategy` is SKIPL, `max_dummy_multiplier` is 10 or more
  and `tloop_stagnation_stop` is off for a non-greedy solver. No selection moves, and the
  estimate is unchanged.
- The study (200 seeded trials per cell, T-LARS, tFDR 0.1, paired replays): the FDR is 0.136
  at 20 p and 0.103 at 10 p (SNR 5), 0.113 at 20 p (SNR 2), 0.075 at 5 p. The clamps Phi' <= 1
  and FDP_hat in [0, 1] have no share: they never activate, and selections are identical
  without them. At 20 p, Phi' < 0 never occurs either, which is why 415484d did not move that row.
- The cause: the total null inflow is about right (7.0 against 7.0 at 20 p), but it collects
  on a few columns, because X and y are the same in every experiment (4.9 nulls above
  Phi = 0.5, where independent experiments give about 1e-16). After step 20 all gain above 0.5
  is false, while the estimate books 1.9 of the 2.9 realized. This corrects the vignette and the
  reading of 2026-10-27, which blamed the steps without new occurrence. Carrying their inflow
  forward restores control (0.079) but over-reads at low voting levels and under-reads at
  v = 1, so it does not calibrate the estimate.
- Vignette page 2 and the READMEs state the corrected reading. C++, R and Python tests cover
  the warning, R vendored copy re-synced.

### 2026-10-28

#### The Tikhonov selector's lambda_2 rules: k-fold CV, the pilot knee, the full-profile knee (D14)

- The Tikhonov T-Rex offers exactly three ways to set lambda_2 when none is
  given: the k-fold CV tuners (`CV_1SE_IEN_CCD`, `CV_MIN_IEN_CCD`,
  `CV_1SE_TIK_SVD`, `CV_MIN_TIK_SVD`), the pilot knee `PILOT_KNEE` inside the
  selector, and the full-profile knee of `tik_lambda2_power_rule.hpp`, one
  FDR-controlled selection per grid value with common dummies. The knee is
  the smallest lambda_2 whose selection count is within a tolerance of the
  largest; a tolerance of 0 reads the largest selection, ties toward the
  smallest lambda_2.
- `PILOT_KNEE` is the default `lambda2_method` (was `CV_1SE_IEN_CCD`) on the
  three surfaces: in the Tikhonov study it reaches the oracle pin's power
  within 0.02 at an FDR below the full-profile knee's, at 13 to 26 % of its
  cost, while the prediction tuners pick the edges of their own grid on
  structured designs. The pilot defaults (grid 0 to 3,000 in nine steps,
  K 10, T 5, tolerance 1, the majority count) are the study's.
- The argmax variants are gone: `TikLambda2Method::PILOT_MAX` leaves the
  enum (`PILOT_KNEE` takes its integer value), and the power rule's
  `Lambda2PowerRule` enum with `MAX_SELECTION` leaves `tik_lambda2_pick.hpp`;
  `pick_lambda2_index(counts, knee_tolerance)` and
  `Lambda2PowerRuleConfig` carry the tolerance alone. The study found the
  argmax's FDR rising with the SNR (the winner's curse of the largest
  selection) where the knee halves that excess.
- R: `trex_tikhonov_control(lambda2_method = )` defaults to `"PILOT_KNEE"`
  and no longer accepts `"PILOT_MAX"`; `trex_tikhonov_lambda2_power_rule()`
  runs the full-profile knee on the caller's matrix (every probe normalizes
  X in place and hands it back, nothing is copied) and returns the pick,
  every probe's selection and stopping time, the counts and the boundary
  flags. Vignette page 7 shows the three rules. Python:
  `TikLambda2Method.PILOT_MAX` is gone, `tikhonov_lambda2_power_rule()` is
  bound at the root and in `trex_selector_methods.ad` with `Lambda2Probe`
  and `Lambda2PowerRuleResult`; the stubs are regenerated. The C++ header
  gains the `Eigen::Map` overloads the bindings call; the owned-matrix
  overloads copy once and delegate.
- Tests: the pick tests read the argmax as tolerance 0, the pilot tests run
  PILOT_KNEE, R and Python test the full-profile knee against a fixed
  lambda_2 run at the pick; the full suites pass.
- Found on the way: the ruler reflow of 0668de0 had put a `NOLINTNEXTLINE`
  comment between `// [[Rcpp::export]]` and the signature of twelve R
  getters, which `Rcpp::compileAttributes()` reads as the signature; the
  committed RcppExports predated the reflow, so the package still built.
  The suppressions are `NOLINTBEGIN` / `NOLINTEND` blocks around those
  functions now, and the generated files are regenerated.

#### The logistic_model placeholder is gone

- `cpp/src/tsolvers/ad/logistic_model/` held one CMakeLists.txt with no
  sources and no target. The folder and its `add_subdirectory` are removed;
  a logistic T-solver, if it comes, belongs to the virtual-dummy family and
  is not part of this release.

#### OrbStack: the R check runs in a UTF-8 locale the image has

- `tools/orbstack/stage_ubuntu.sh` picks `en_US.UTF-8` when the image has it
  (the apt step generates it) and the image's `C.UTF-8` otherwise. A snapshot
  taken before `locale-gen` fell back to an ASCII session, where `R CMD check`
  ends with a bare ERROR at the DESCRIPTION meta-information step instead of
  the BH NOTE. The Linux run of the dev_beyond batch at 09debc6 is green: C++
  688 tests and the install smoke, R 1,531 tests and the check with two NOTEs
  (installed size, Ubuntu's BH), Python 624 tests, the manylinux cp312 wheel
  with 624 tests.

### 2026-10-27

#### The 120-character ruler holds project wide (D8)

- `tools/check_line_length.py` lists every tracked text line over 120
  characters and exits with status 1 when it finds one; the workflow
  `line-ruler.yml` runs it on every push and pull request, without a path
  filter. Exempt: data and binary files, generated files (the pybind11 stubs,
  the roxygen pages, RcppExports, inst/CITATION, the vendored copy and the
  Makevars of the R package), the clang-tidy check list, Markdown table rows
  and lines that carry a URL or a DOI.
- 394 lines in 52 files were continued on a new line, one mechanical
  rule per kind: pybind11 docstrings split into adjacent literals (the stubs
  stay byte-identical), C++ and R call lines wrapped at commas, the one-line
  Rcpp getters expanded to three lines with `NOLINTNEXTLINE`, Markdown prose
  wrapped at spaces, and a few CMake, YAML, shell and JSON lines by hand. The
  CITATION.cff titles use folded scalars, the reference map holds one key per
  line where an array was long.

#### The READMEs describe the dummy strategies

- The C++, R and Python READMEs gain a "Dummy Strategies" section before their
  Threads section: the six L-loop strategies with their storage, STANDARD as the
  default on every surface, the seeded twins, the measured gain of HCONCAT, the
  SKIPL caveat, the dummy distributions, and where the study lives (vignette
  page 2). The last documentation item of WP18.

#### The control defaults agree on the three surfaces (D11, D12)

- An audit of every control constructor of R against the Python and C++
  structs (identical by construction) found four deviations, all in R or in
  binding logic, and the three surfaces now carry one value each:
  - `trex_control(lloop_strategy = )` defaults to `"STANDARD"` as C++ and
    Python do (was `"HCONCAT"`; the roxygen page and vignette page 2 called
    HCONCAT the default). Decided 2026-10-27 (D11): the fresh-dummy law of the
    T-Rex paper is the default everywhere.
  - `trex_control(tloop_stagnation_stop = NULL)` passes NULL to the C++ core,
    whose solver rule (`isGreedySolver()`) decides on every surface. R used to
    resolve NULL with its own list, which counted TIENET and TIENET_AUG as
    greedy, so an R GVS run on the IEN track stopped its T-loops early where
    Python and C++ ran them to the end (T_stop 7 against 11 in the audit cell).
    The path solvers never stop by default; the greedy ones do.
  - `trex_control(cd_lambda_rel_tol = )` defaults to 1e-6 as C++ and Python do
    (was the former 1e-5).
  - The lambda_2 rule on the IEN track resolves in the C++ GVS constructor
    (D12): the struct default CV_1SE_CCD becomes CV_1SE_IEN_CCD when
    gvs_type = IEN, and the IEN-geometry tuners are refused on the EN track at
    construction. R and the Python wrapper applied that rule in their own code;
    the Python wrapper's copy is gone, R's NULL resolution gives the same value.
    A C++ user who kept the default on the IEN track got the EN-geometry tuner
    (lambda_2 = 29,555 and 2 selected variables in the audit cell against 42).
- The R runtime gap of the Examples is explained on the way: at equal
  parameters the core cell takes the same time in R and Python, and the 17 %
  of the GVS cell is the Eigen version (RcppEigen ships 3.4.0, the C++ tree and
  the wheels use 5.0.1; a wheel built against Eigen 3.4.0 reproduces R's time).
  No library work; the study is archived in the research folder.

#### The root README's reference list is rendered from the master bibliography

- `references.bib` gains an "Application references" section with the two
  papers of the README's opening story (Huffman 2018, Nature Communications,
  a single-author comment article; Harvey, Liu and Zhu 2016, Review of
  Financial Studies), fields verified against their Crossref records.
- `render_references.py` gains `--format md`: one unwrapped Markdown bullet per
  entry with the DOI or the arXiv id as a link. The README's reference list is
  rendered with it (the command stands in an HTML comment above the list) and
  the inline citations of the opening story link the DOIs.

#### R MLRidgeCV aligned with C++ and Python, the clustering's backend line on request, acknowledgements

- `MLRidgeCV$fit()` takes `n_folds` (was `num_folds`; MLElasticNetCV, Python and
  C++ already used that name). The ridge CV grid default is 100 points on every
  surface: C++ `ridge_cv_svd::fit()` and Python `RidgeCV.fit()` come down from
  1,000 to R's 100, the glmnet convention (decided 2026-10-27). The grid stays
  geometric over three decades; the CV curve is smooth in log lambda, so 100
  points locate the minimum and the one-standard-error point within the fold
  noise, while every further point costs a shrinkage pass and a held-out
  prediction per fold. The two queries of the C++ class reach R: `interpolates()` (the smallest
  training fold has p >= n_train, so the lambda -> 0 end of the grid
  interpolates) and `selection_at_grid_edge()` (a pick on the first or last
  grid point). Breaking for callers of `num_folds` (see the release block).
- `ml_agglomerative_cluster()` and `agglomerative_cluster()` gain
  `verbose = FALSE` / `verbose=False`: the backend's allocation line (the size
  of the distance matrix in RAM or in the memory-mapped file) prints on request
  only. The R and Python bindings used to call the dispatcher with its C++
  default of true, so every call wrote to the console and vignette 9 wrapped
  its calls in capture.output(). The C++ default stays true.
- The three READMEs acknowledge cereal, Boost.Interprocess, OpenMP and the
  binding layers (Rcpp and RcppEigen, pybind11) next to Eigen.
- Vignette page 2 carries the SKIPL caveat: under a large fixed L the FDP
  estimate is anti-conservative (an FDR of 0.137 at a target of 0.1 in the
  Examples demo 04).

#### T-Stagewise: the dummies of the support count toward T (D10)

- The early stop at T_stop and the tie budget of the LARS family count a dummy
  while it is active or carries a frozen nonzero coefficient. Before, the count
  followed the active set only: a dummy the NNLS reduction of T-Stagewise
  removed stopped counting although its frozen coefficient kept it selected,
  so the path ran on until T dummies were active again while the removed null
  variables stayed counted. The included nulls at T then exceeded the
  hypergeometric count the FDP estimate assumes (Examples demo_ts_perf_05: a
  mean of 21.35 against 19.60 at T = 20, Kolmogorov distance 0.10 against at
  most 0.014 for the other solver configurations). Now the path stops at the
  first step whose support holds T_stop dummies, the step the T-Rex runner maps
  to t = T_stop, and a frozen dummy that re-enters is not counted twice. The
  lasso-family solvers are unchanged in effect, since a lasso drop zeroes the
  coefficient. T-Stagewise results move (see the release block).
- Test `TStagewiseSelection.FrozenDummiesCountTowardT`: the support of the last
  step holds at least T_stop dummies and the step before fewer, in runs with
  frozen dummies.

#### T-Stagewise: a variable removed within a step is no candidate for that step's move

- The 60 / 56 / 50 steps that T-Stagewise took on bit-identical data in the
  C++, Python and R builds (n = 50, p = L = 1,000, T = 10) were rounding, not
  a tie. A variable the NNLS reduction removes still sits at the active level,
  so its root in the step-size formula is rounding noise over a nonzero
  denominator. A noise root above machine epsilon gave a near-zero step after
  which the variable re-entered as the maximum and was removed again, until the
  noise changed sign. The BLAS and the no-BLAS build of the same machine
  disagreed on 75 of 100 seeds at n = 50 and 100, on 15 at n = 200; T-LARS and
  T-LASSO on none.
- Now the variables removed within a step are excluded from the candidates of
  that step and rejoin with exact correlations after the move, the
  drop-breakpoint rule the lasso-family solvers apply in their recovery step.
  The rule is exact: by the Schur complement the removed variable's correlation
  falls faster than the active level, so it cannot catch up within the step.
- Effect: the two builds agree on every seed; the path equals the old path
  without the pure cycle steps (100 of 100); the coefficients and the supports
  at T are unchanged below saturation (100 of 100 at n = 50, T = 10); the mean
  step count at n = 50 falls from 38.7 to 34.3; runs that the old code
  truncated at the step limit after spending it on cycles now reach T. Genuine
  cycling at saturation remains (the decision of 2026-09-27 stands). Breaking
  for the recorded T-Stagewise step counts only (see the release block).
- Test `TStagewiseSelection.NoRoundingNoiseStepsAfterNNLSRemoval`: every step
  lowers the active level by a relative margin above 1e-12 (the old code
  failed it in most runs, with decreases of exactly zero).

#### The clustering refuses a RAM buffer beyond the physical memory and scans in parallel

- `BlockTiledMatrixPolicy`, the distance matrix of every linkage except single
  under `use_mmap = FALSE`, compares its size with the installed physical memory
  before the zero fill. A buffer the machine cannot hold throws
  `std::runtime_error` with both sizes and the remedy, the memory-mapped
  backend. Before, the fill touched page after page until the kernel killed the
  process (the Examples hac_02 test C). The size comes from the new
  `utils/physical_memory.hpp` (`sysctl` on macOS, `GlobalMemoryStatusEx` on
  Windows, `sysconf` elsewhere), and an unknown size passes. Test
  `AgglomerativeGuardsTest.RamBufferBeyondPhysicalMemoryIsRefused`: the check on
  fixed sizes and a 1,000,000-object construction that throws at once.
- The dispatcher's non-finite scan of the clustering input runs over the
  columns in parallel on the process-wide thread count, like the other one-shot
  kernels. The smallest offending column wins, so the message names the entry a
  serial scan names. The serial scan read a mapped design of 400,000 columns
  once more on a single core (6 to 9 minutes in the Examples).
- TLARS, TIENET, TIKNET, TNCGMP, TAFS and the SD solver base move the updated
  Cholesky factor into place instead of copying it, as TOMP does.
- The TCDTIKNET warnings under `penalty_structure = GROUP_CENTERING` (sweep
  budget, certification fallback) carry the tag `[TCDIENET]` by design: the
  dispatch serves that path with TCDIENET's group-centering hooks. The Tikhonov
  control docs of the three surfaces and vignette page 7 say so.

#### SEEDED is the storage-free twin of STANDARD (WP18)

- `LLoopStrategy::SEEDED` draws fresh dummies at every L step, keyed on
  (experiment k, L step) exactly as `STANDARD` keys its stored matrices, so a
  SEEDED run and a STANDARD run with the same seed draw bit-identical dummies
  and select identically; they differ only in storage. Before, the i.i.d. case
  of SEEDED re-derived its dummies from the prefix-stable tag 0 (the nested
  law of HCONCAT, one draw persisting across the L steps) while its cluster
  case (GVS, Tikhonov) already drew fresh. The nested law stays available as
  HCONCAT, which stores its matrices.
- The decision rests on the paired study of 500 trials in nine cells (Examples
  notes/tools/fresh_vs_nested, recorded in the archive and on vignette page 2):
  fresh and nested dummies give the same FDR, the nested form gains 0.6
  percentage points of TPR on average through a larger calibrated L and longer
  T-loops. No second seeded form and no switch.
- Mechanics: `DummyGenerator::generateSeeded()` / `generateSeededInto()` take the
  L step as their tag, the experiment runner passes `cfg.seed_factor` (the core
  sets it to the L step for SEEDED, as for STANDARD, and leaves 0 for the
  prefix-stable PERMUTATION_SEEDED base), and the screening keys its single
  pass on L step 1 for STANDARD and SEEDED alike, so its pin holds.
- Breaking for reproducibility only (see the release block): a seeded run with
  `lloop_strategy = SEEDED` (and a seeded screening run with STANDARD or
  SEEDED) selects differently from before, from the same law. Unseeded runs
  are unaffected.
- Tests: `LLoopSeeded.SeededMatchesStandard` (three seeds, the support, T_stop,
  L, Phi and FDP-hat), `SeededDrawsAreStandardsDraws` (the generator: the
  SEEDED draw at L step LL equals the stored STANDARD matrix and differs from
  the tag-0 base), `SeededStrategiesAgreeInMemoryAndMapped`, and the same
  identity in the DA and GVS selectors; the Tikhonov and screening pins stay.
  Documentation in the enum, the generator's tag conventions, the R and Python
  control docs.

#### READMEs describe the threading model; the every-core value is gone (WP15 step 5)

- The C++, R and Python READMEs gain a Threads section: the serial defaults, the
  two selector options and the solver budget, the process-wide count of the
  one-shot kernels, determinism, builds without OpenMP.
- The value 0 ("every core") of `max_outer_threads`, `max_inner_threads`,
  `SolverOptions::n_threads` and the SD / VD `n_threads` is removed: every thread
  option is a whole number of at least 1, other values throw, and R and Python
  validate the same way. One keyword that hands the whole machine to a single
  run invites oversubscription next to a caller's own parallel loop or other
  work. A user who wants every core reads the count (`trex_get_max_threads()`,
  `trex_selector_neo.utils.get_max_threads()`, `omp_get_num_procs()`) and
  passes it. The entries of steps 1, 2 and 4 below describe the value as it
  was.
- The work threshold of a million entries per pass is the shipped value, and
  the options carry no caps (decided 2026-10-27). Tests, help pages,
  docstrings, stubs and the two vignettes follow.

#### R and Python: the solvers' thread budget is bound (WP15 step 4)

- R: `tsolver_options(n_threads = 1)` carries the budget of a solver's own passes
  (1 serial, 0 every core), the binding parses it, and the T-solver classes gain
  `$set_num_threads()` and `$get_num_threads()`. Python: `SolverOptions.n_threads`
  with the same meaning, and `setNumThreads()` / `getNumThreads()` on every
  solver; the stubs are regenerated against pybind11 3.1.0. The selector options
  (`parallel_rnd_experiments`, `max_outer_threads`, `max_inner_threads`) were
  bound in step 2.
- Tests: R and Python run T-LARS on a design above the work threshold (400 x
  1,500 predictors and as many dummies) with budgets 1 and 4 and pin the actions
  and the coefficient path bit for bit, resolve 0 to at least one thread, change
  the budget on a live solver and reject a negative one.
- The T-solvers vignette lists the option, the parallelism vignette names it for a
  solver used on its own.

#### Selections are bitwise independent of the thread options: the test (WP15 step 3)

- `test_trex_thread_determinism.cpp` runs every selector (T-Rex with the HCONCAT,
  STANDARD, SEEDED and PERMUTATION_SEEDED strategies, T-Rex+DA, GVS, Screen-T-Rex,
  the biobank screening, SPCA, Tikhonov) on one seeded problem under six thread
  settings: experiments in sequence with a solver budget of 1, 2 and 4, experiments
  in parallel on 2 and 4 threads, and both options at 0. The selected support, the
  stopping time, the voting threshold, the Phi, FDP-hat and R matrices and Phi'
  (the screening's estimate and coefficients, the biobank record, the SPCA supports,
  loadings and scores) equal those of the sequential one-thread run bit for bit.
  The test runs in every build; without OpenMP the settings are inert and the
  identities hold trivially.
- The two ingredients were already in place: experiment k derives its own seed
  and every reduction inside the solvers is summed in a fixed order (the TACGP sum
  since step 1), and the dropped-column lists of the standardization are sorted
  after the per-thread merge in `preprocessMatStatistics()` and in
  `TSolver_Base::preprocess()`, so the plan's "sort after the merge" item needed
  no change.
- Verified in the BLAS build and in a no-BLAS build (where Eigen's own threading
  is live), both with OpenMP.

#### Selectors: the thread options are honoured, serial by default (WP15 step 2)

- One owner of the cores at a time. `parallel_rnd_experiments` (default false)
  runs the K random experiments on `max_outer_threads` threads (default 1, 0 =
  every core, capped at K) and keeps every solver serial inside the team. With
  the experiments in sequence (the default, and always for the lean strategies
  SEEDED, PERMUTATION and PERMUTATION_SEEDED) the solvers get
  `max_inner_threads` as their budget (default 1, 0 = every core), which their
  regions split only above the work threshold of step 1. Nested teams are never
  formed: the K-loop carries an explicit `num_threads` clause, and the solvers'
  regions test `omp_in_parallel()`.
- Gone: `TRexControlParameter::autoConfigResources()`, which overwrote both
  options with the core count (every core for the inner solvers whenever the
  experiments ran in sequence, so a serial selector in an R session took all
  cores above the threshold), the `omp_set_max_active_levels(2)` and the
  per-worker `omp_set_num_threads()` of the K-loops, and the runner's
  `EigenThreadScope`: the solvers pin Eigen's count to their own budget per
  call (step 1), outside a team only. `resolveThreadOptions()` validates and
  resolves the two options at construction.
- The R binding passed `parallel_rnd_experiments` and `use_memory_mapping` to
  the backend but not the two thread counts, so `trex_control(max_outer_threads
  = )` never reached C++. Both are parsed now. `trex_control()` accepts 0 for
  both and rejects negatives. `use_openmp` is renamed `parallel_rnd_experiments`
  (breaking, see the release block). The parallelism vignette describes the
  options; the Python docstrings of the three fields follow.
- SD and VD: `TRexSDOptions::n_threads` and `TRexVDOptions::n_threads` are the
  K-loop threads with the same semantics (1 default, 0 every core, capped at K,
  negatives throw) instead of "0 = auto, every core". `TRexSDOptions` gains
  `n_inner_threads`, the budget of the SD solvers' step-length search, passed
  through `SDTSolver_Base::setNumThreads()`; the O(m) correlation updates and
  the pool refreshes of the SD solvers lost their pragmas (they never paid), the
  step-length search is gated like the ad solvers.
- Tests: `TRexCoreTest.ThreadOptionsResolveAtConstruction` (defaults 1 / 1, 0
  resolves to the core count, the outer count is capped at K, negatives throw at
  resolution and at construction) and `ParallelExperimentsSelectAsSerial` (a
  seeded selector selects the same variables in sequence, with parallel
  experiments and with a solver budget of 2). The thread hygiene tests keep
  their meaning: nothing writes the host's OpenMP or Eigen count any more.
  C++ suite 672 tests.
- R vendored copy re-synced.

### 2026-10-26

#### R: twelve vignettes replace the usage-and-simulations vignette (WP16)

- `vignettes/` holds twelve pages, each a CRAN vignette in `rmarkdown::html_vignette`
  format, seeded, small (every page renders in about one second) and capped at two
  threads by its setup chunk: getting started, the T-Rex selector, T-Rex+DA, T-Rex GVS,
  screening (Screen-T-Rex and the biobank screening), T-Rex SPCA, T-Rex Tikhonov, the
  T-solvers, the ML methods, memory mapping, working without data copies, parallelism.
  The file names carry a `v` and two digits for the order (a name under `inst/doc` must
  start with a letter, or `R CMD build` drops the file), and `_pkgdown.yml` groups them as
  articles.
- The Rd pages of `trex_get_max_threads`, `Gauss_data` and the four `mmap_matrix` S3
  methods carry a `\concept`, so the pkgdown reference groups cover every topic.
- The old `TRexSelector_usage_and_simulations.Rmd` (html_document with a CDN MathJax, a
  prettydoc output outside Suggests, unconditional ggplot2 / patchwork, an unseeded
  selector, the 100-trial simulation at build time) is gone, together with
  `vignettes/debug.R` (a devtools script that shipped in the tarball) and the stale
  `refs.bib`. The pages cite the master `references.bib` of the repository, copied into
  `vignettes/`, with the IEEE style. The plots are base graphics, so ggplot2 and
  patchwork leave Suggests.
- The T-Rex selector page carries the fresh-against-grown dummy comparison of the
  L-loop strategies (the 500-trial study behind decision D3) and the cross-platform note
  on seeded runs. The parallelism page describes the current knobs
  (`trex_control(use_openmp, max_outer_threads, max_inner_threads)`,
  `trex_set_num_threads()`); it is updated when the selector-level thread options land.
- Known cosmetic: pandoc 3.12 prints a deprecation notice for `--mathjax` during the
  build, from rmarkdown, harmless.

#### R CI: a sequential build job; the gate variables of the regions are `[[maybe_unused]]`

- The R workflow gains a fourth matrix entry, Ubuntu with `TREX_OPENMP=no`, the
  sequential build configure falls back to without an `omp.h` (CRAN's macOS R has no
  OpenMP flags). The job asserts that configure reports the sequential build and that no
  compile line carries `-fopenmp`, then runs the testthat suite and R CMD check as the
  other three do. Locally the sequential build passes the suite (1,505 expectations,
  `trex_get_max_threads()` 1). The C++ workflow had this job already.
- In the sequential build the `num_threads(nt) if (nt > 1)` clause is the only reader of
  the gate variable and sits in an ignored pragma, so the 17 sites warned about an
  unused variable. The variable is `[[maybe_unused]]`.
- R vendored copy re-synced.

#### T-solvers: a per-solver thread budget replaces the process-wide OpenMP count (WP15 step 1)

- `SolverOptions::n_threads` (default 1, 0 = every core, negative throws) is the budget
  of the solver's own OpenMP regions. `TSolver_Base::setNumThreads()` (virtual, the
  augmented wrappers forward it to the inner solver together with the budget at
  construction) changes it on a live solver, `getNumThreads()` reads the resolved value.
  The budget is not serialized, so a deserialized solver starts serial and the selector
  dispatch re-applies it to warm-started and loaded solvers.
- Every kept region asks `threadsFor(work)` with the entries of the design the pass
  touches and runs `num_threads(nt) if (nt > 1)`. The answer is 1 below the provisional
  threshold `kParallelWorkThreshold` = 1e6 entries, inside a parallel region of the
  caller (`omp_in_parallel()`, so a solver inside a parallel K-loop never nests a team)
  and with a budget of 1. The threshold comes from the measurement of the inner regions
  in the Examples (notes/tools/wp15_inner_threads): the passes over the design pay from
  n x m ~ 1e6 (x1.5 to x1.9 at two threads, x1.9 to x3.6 at four), below that the
  fork-join (5 to 40 us) costs more than it saves.
- Kept and gated: the correlation passes of the base and the LARS family (initialize,
  refresh inactive, refresh all, the ENET update), the step-length search of TLARS,
  TIENET and TIKNET, the T-OLS norms and downdate, the TACGP dot products with the
  active atoms, the two passes of TCDLASSO, and the standardization and its restore
  (`preprocessMatStatistics` / `restoreMatStatistics` take the thread count as a
  parameter and no longer touch Eigen's global setting).
- Dropped as never paying: the parallel inactive-set scan with per-thread buffers
  (a serial O(m) scan now), the incremental O(m) updates of the step-length search,
  the T-OLS projected-norm scoring loop, the T-GP and TACGP k-loops (direction copy
  and update), and the `if (k > 100)` reductions of TACGP, whose numerator is summed
  in a fixed order so the value does not depend on the thread count.
- `EigenSingleThreadGuard` is `EigenThreadBudgetGuard(n_threads)`: it pins Eigen's
  global count to the budget for the solver call and does nothing inside a caller's
  parallel region, where a write to the global setting would race. Under
  `EIGEN_USE_BLAS` the count is inert either way.
- Behaviour: a T-solver used directly (C++, and through R and Python until the
  bindings expose the budget) runs its regions serially, where before they ran on every
  core OpenMP offered. The selectors are unchanged in this step: the experiment runner
  and the GVS K-loop pass `max_inner_threads` as the budget, so serial random
  experiments keep their inner regions on the cores and parallel experiments keep the
  solvers serial. The selector-level `n_threads` that replaces `max_outer_threads` /
  `max_inner_threads` and `autoConfigResources` follows in step 2.
- Tests: `TSolverInitTest.ThreadBudgetDefaultsToSerial` (default 1, 0 resolves to the
  core count, negative throws, the setter follows) and the typed
  `TSolverPolymorphicTest.PathIsIndependentOfTheThreadBudget` (n 400, p = L = 1,500,
  four steps with budgets 1 and 4, actions and sparse path identical; skipped without
  OpenMP or on one core).
- R vendored copy re-synced.

#### TIKOMP removed: a greedy pursuit cannot give the Tikhonov net's grouping

- The Tikhonov-K OMP was the greedy counterpart of TIKNET and TCDTIKNET under a
  general penalty K. A study on a grouped design (30 equicorrelated groups of 10,
  three active, n = 150, p = 300, the group-centering K of the true groups,
  within-group correlation 0.5 and 0.9, lambda_2 1 and 10, 20 trials; the tool and
  its tables sit in the Examples under notes/tools/tikomp_grouping) shows why it
  cannot work as intended. At the first dummy TIKNET and TCDTIKNET have the active
  groups complete with equal coefficients (the elastic-net grouping), TIKOMP has 7
  to 10 of the 30 true members: its penalized correlation pulls the mates of an
  active variable in at large lambda_2, but the pursuit admits one member per step
  against an orthogonalized residual while the dummies compete at every step. In
  `TRexTikhonovSelector` it reached a TPR of 0.04 to 0.17 at an FDR of 0 (the L-loop
  ran toward its cap) against 0.95 to 1.00 for the other two solvers.
- Removed: `tikomp_solver.hpp` / `.cpp` (archived with the study and pre-removal
  snapshots of the touched files in TRex_Research), the `SolverTypeOMPBased` and
  `SolverTypeForTRex` enumerators, the dispatch case, the Tikhonov selector's gate
  (TCDTIKNET or TIKNET now), the R6 class with its Rcpp creators and getter
  branches, the Python binding and the omp_based re-export, the tests of the solver
  (anchors, augmented cross-check, refit closed form, round trip, validation), the
  selector gate and end-to-end tests, the memory and polymorphic type lists, and the
  references map entry. No demo uses TIKOMP.
- Breaking, see the release block: the R package drops `TIKOMP_Solver` and the
  `"TIKOMP"` value of `trex_control(solver = )`, Python drops
  `tsolvers.ad.omp_based.TIKOMP_Solver` and `SolverTypeForTRex.TIKOMP`.

#### T-solvers: the augmented wrappers forward seed and tolerance to the inner solver; small items

- `TSolver_Base::setTieSeed()` and `setTolerance()` are virtual. `TENETAug_Solver` and
  `TIENETAug_Solver` hid both with forwarding members, so a call through a `TSolver_Base`
  pointer (R's `$set_tie_seed()` and `$set_tolerance()`, which hold the solver as a base
  pointer) seeded the wrapper's unused engine and set the wrapper's unused tolerance. The
  overrides forward to the inner solver. Python holds the concrete type and was not affected.
- The wrappers also built the inner solver's options without `tie_seed`, so the inner
  tie-break engine was always seeded 0 and a seeded selector's per-experiment seed (or -1 for
  an unseeded run) stopped at the wrapper. `tie_seed` is forwarded at construction. A
  selection can only move where the inner path meets an exact dummy tie (identical columns),
  which continuous designs do not produce.
- Tests: `TieSeedDispatch.AugmentedWrappersForwardThroughBasePointer` (two bitwise identical
  dummies and a budget of one; the survivor set through a base pointer or through the options
  equals the one set on the concrete type, and the seed decides it) and the same check in R
  (`test-tsolvers.R`, through `$set_tie_seed()` and `tsolver_options(tie_seed)`). The full
  C++ suite (660 tests) and the R tsolvers testthat files (547 expectations, against a
  scratch install) pass.
- `test_ml_methods_memory.cpp` included `<sys/resource.h>` and called `getrusage` unguarded.
  The Windows CI job builds every test, so the next push would have broken the MSVC job. On
  Windows the file holds one skipped test.
- The dead `.gitignore` line (`TODO before main merge.txt`) is gone. The briefs of the
  `tsolvers/ad.hpp` and `tsolvers/tsolvers.hpp` umbrella headers are one sentence each.
- R vendored copy re-synced.

#### R documentation of the TRexSelector result fields

- The field texts of `T_stop`, `L`, `voting_grid`, `fdp_hat_mat`, `phi_mat`, `r_mat` and `phi_prime` in
  `trex_selector.R` and `TRexSelector.Rd` described the wrong quantities ("probability mass function
  evaluations", "empirical CDF computations", "grid of voters over trials and variables"). They now name the
  relative occurrences after each T step, the number of selected variables and the FDP estimate per T step and
  voting level, the voting levels, the calibrated occurrences (at most 1, possibly negative), the stopping time and
  the dummy multiplier. Documentation only.

### 2026-10-24

#### FDP estimate: `Phi_prime` keeps negative values (lower clamp removed)

- `TRexSelector::computePhiPrime` clamped `Phi_prime` to [0, 1]. The upper clamp
  (Phi' > 1 -> 1) prevents a negative false-discovery mass 1 - Phi' and can only
  raise FDP_hat. The lower clamp (Phi' < 0 -> 0) capped a selected variable's
  false-discovery mass at one, where the T-Rex paper and the R reference let it
  exceed one, and made FDP_hat anti-conservative in long T-loops. In demo 09 of
  the Examples (n 300, p 1,000, 10 unit coefficients, Normal predictors and noise,
  K 20), TLARS at SNR 5 reached an FDR of 0.530 at tFDR 0.5 and 0.415 at tFDR 0.4
  (no stagnation stop, T-loop to T ~ 38), and FDP_hat at fixed T and v -> 1 fell
  below the realized FDP from T = 10 on (0.408 against 0.444 at T = 20). Replaying
  the final selection on the same 200 data sets with the unclipped Phi_prime gives
  0.462 and 0.362. Only the upper clamp remains, together with the [0, 1] clamp of
  `computeFDPHat`, which guards against the negative FDP estimates seen in
  trex_da (they come from Phi' > 1).
- The same change in `TRexSD::Phi_prime_fun_` and `TRexVD::Phi_prime_fun_`.
  trex_da, trex_gvs, trex_screening, trex_tikhonov and trex_spca reuse
  `computePhiPrime` and follow. FDP_hat can only rise, so the L-loop may calibrate
  a larger L, the T-loop may stop earlier and selections can only shrink in size.
  Results move where a selected variable had Phi' < 0, which the replay found at
  tFDR >= 0.3 and not at tFDR 0.1 in the tested cells (TLARS and TOMP, SNR 1
  and 5); the L-loop and T-loop effects need full runs.
- Tests: `TRexCoreTest.PhiPrimeKeepsNegativeValues` builds a Phi history with a
  negative Phi_prime and checks it and FDP_hat against the paper's formula. The
  Phi_prime bound checks of the core, Tikhonov, GVS and VD tests and of the R
  testthat files keep only Phi' <= 1 (`TRexVDTest.PhiPrimeIsCappedAtOne` also
  checks FDP_hat in [0, 1]). The full C++ suite (659 tests), the Python
  selector tests (172) and the R testthat files of the selector, DA, GVS,
  screening and scaling paths pass.
- R vendored copy re-synced (`bundle_for_cran.sh --check` passes).

### 2026-10-22

#### TRex-SPCA accepts the elastic net only; flaky low-memory CV budget

- `TRexSPCASelector` is defined on the elastic net, but the R control
  `trex_spca_control(gvs_type = "IEN")` and the Python docstring offered the
  informed elastic net, and the C++ core silently ran EN instead. An IEN request
  is now an error: the C++ constructor throws `std::invalid_argument` when
  `gvs_ctrl.gvs_type` is not EN (Python raises `ValueError` in `select()`), and
  `trex_spca_control()` stops for any `gvs_type` other than `"EN"`. The IEN
  example of `trex_spca_control()` is gone, the Python docstrings (binding,
  stub, wrapper) say EN only, and the Python test that ran "IEN" as EN now
  checks the error. `cv_dfmax_ratio` stays unexposed in the SPCA controls: it
  only acts on the IEN-geometry tuners, which the EN track rejects.
- `MLMemory.ElasticNetCVLowMemory`: budget 0.3 x X instead of 0.2. One
  low-memory fit needs 0.095 x X in every run, but the six sequential fits of a
  5-fold CV measured 0.12-0.23 x X for identical input (whether the allocator
  reuses the freed paths varies), so about one run in four failed. A copy of X
  would still add 1.0 x X.
- Tests: `test_trex_spca` gains `Exceptions_IENRejected` (16 tests pass),
  `test_ml_methods_memory` passes 12 of 12 repeated runs of the CV case. The R
  test that ran "IEN" as EN now checks that a hand-built IEN control list fails
  in the C++ constructor, next to a check of `trex_spca_control()`. The Python
  SPCA tests (23), the R SPCA testthat file (25 blocks) and the full C++ suite
  (658 tests) pass.
- R vendored copy re-synced (`bundle_for_cran.sh --check` passes).

### 2026-10-21

#### Settable support cap of the IEN-geometry lambda_2 tuners (`cv_dfmax_ratio`)

- The IEN-geometry tuner `ienet_cv_ccd` (`CV_1SE_IEN_CCD` / `CV_MIN_IEN_CCD`)
  stops each lambda_1 path of its reference scan once more than
  floor(`dfmax_ratio` x n_train_min) coefficients are nonzero. `TRexGVSSelector`
  and `TRexTikhonovSelector` called it without the ratio, so the cap was fixed at
  0.8 of the smallest fold's training size. The new control field
  `cv_dfmax_ratio` (default 0.8) sets it: `<= 0` disables the cap, larger values
  allow larger supports at a higher tuner cost, and a non-finite value throws
  `std::invalid_argument` at construction. The selectors now pass the tuner's
  own defaults for `max_iter` (100000) and `tol` (1e-4) explicitly, so the
  default 0.8 gives the same results as before.
- C++: `TRexGVSControlParameter::cv_dfmax_ratio` and
  `TRexTikhonovControlParameter::cv_dfmax_ratio`. Python: the same field on both
  control classes (stubs updated by hand). R: argument `cv_dfmax_ratio` of
  `trex_gvs_control()` and `trex_tikhonov_control()`, which reject a non-finite
  value as well; Rd files updated by hand.
- `lambda2_cv_diagnostics` gains `dfmax_ratio` (the ratio the IEN tuner was
  called with, as passed; 0 when no IEN-geometry tuner ran; `support_cap == -1`
  still marks a disabled cap) and `profile_at_support_cap` (the profiled
  lambda_1 optimum of the chosen lambda_2 sits on the truncation edge, and the
  reference path of that lambda_2 stopped on the support cap rather than on a
  deviance rule or the end of the path). Both are exposed in Python
  (`Lambda2CVDiagnostics`) and in the R `lambda2_cv_diagnostics` list.
- Boundary-pick warning: when the edge is the support cap it says so ("here the
  support cap (more than N nonzero coefficients, cv_dfmax_ratio = r)") and adds
  the remedy "raise the support cap (cv_dfmax_ratio)"; the other wording is
  unchanged.
- Tests: `test_trex_lambda2_guard` (GVS default and explicit 0.8 bit-identical,
  0.95 raises the cap from 43 to 51 at n = 60 with 10 folds, the ratio rides in
  the diagnostics and the GVS result, non-finite values throw; the same
  forwarding for the Tikhonov selector), `test_ml_methods_model_selection` (a
  binding cap is flagged and named in the warning, a disabled cap is not; the
  message contract), Python and R control round trips and diagnostics fields.
  All pass: the full C++ suite (657 tests), the Python GVS/Tikhonov tests (40)
  and the R GVS/Tikhonov testthat files (29 blocks).
- R vendored copy re-synced (`bundle_for_cran.sh --check` passes).

### 2026-10-20

#### Free-space check before memory-mapped files are created

- Memory-mapped files are created sparse, so a full disk only showed at a later
  mapped write and ended the process with SIGBUS (no error, no cleanup). The
  new `utils/disk_space.hpp` compares the bytes a mapping may write with the
  free space of its volume and throws `std::runtime_error` with the required
  and free GiB, the directory and a remedy instead; a volume that cannot be
  queried passes.
- `MemmapManager::initialize()` (the dummy matrices of every T-Rex selector with
  `use_memory_mapping`) checks the worst case before it creates anything: all
  files filled up to the dummy cap, K files (one if shared) of
  n x `max_dummy_multiplier` x p doubles (`requiredBytes()`). The message names
  `max_dummy_multiplier`, K and the temporary directory (TMPDIR; R sets
  `tempdir()` at load; `TempDirManager::setCustomTempDir` in C++).
- `MemoryMappedMatrix` in ReadWrite mode checks the growth of a new or enlarged
  file at its own location before creating it, so user-created designs (R
  `mem_matrix`, Python `MemoryMappedMatrix`) are covered as well.
- `TempDirManager::create` takes the description and the remedy for this
  message, `TempDirManager::baseDir()` returns the directory it creates in; the
  warm-start checkpoints keep their 5 GiB minimum with the new message.
- Test `test_utils_disk_space` (7 cases: helpers, matrix, temp dir and manager
  refuse an unfulfillable size before creating a file or directory, small sizes
  still work). R vendored copy re-synced (`bundle_for_cran.sh --check` passes).
  No result changes.

### 2026-10-08

#### Scientific references in the class documentation; master bibliography; citation files

- `references.bib` at the repository root is the master list of every paper the
  library credits (49 entries: the T-Rex method papers, the origin paper of each
  solver algorithm, the clustering, hashing, cross-validation and ridge
  references, Eigen). `tools/references/render_references.py` renders an entry
  as one plain-text line (initials and surnames, title, venue, volume, pages,
  year, DOI) wrapped for a Doxygen block, a roxygen block or a pybind11
  docstring; `references_map.json` fixes which keys each class cites, the origin
  papers first and the T-Rex paper (the greedy-pursuits paper for the greedy
  solvers) for the terminating form. No citation machinery in the sources: no
  Doxygen `\cite`, no Rdpack, no Sphinx roles.
- C++: every selector (`TRexSelector`, DA, GVS, screening, biobank screening,
  SPCA, Tikhonov), `MemoryMappedMatrix`, all 21 T-solvers, the clustering
  classes (SLINK, NN-chain, generic linkage, the dispatcher, the Lance-Williams
  coefficients, the LSH distance and update policies), the elastic-net, IEN,
  ridge and Tikhonov cross-validation tuners, the EBIC criterion, PCA, the SVD
  solver, the ridge solver and the dummy `Distribution` carry a `References:`
  paragraph in their class block. The hand-written lists of the GVS, TENET,
  TIENET, TIKNET, TMP, TGP, TACGP, SLINK, NN-chain and generic-linkage headers
  and the `@note Reference` lines are replaced by it; TCDLASSO no longer points
  to an internal note that is not in the repository, and two demo pointers are
  gone. Vendored copy in the R package re-synced (`bundle_for_cran.sh --check`
  passes).
- R: `@references` on the seven selectors, `mem_matrix` / `mem_matrix_from`, the
  21 T-solver classes, the six ML classes, the three clustering functions and
  `trex_dummy_distribution`; 40 Rd files regenerated, NAMESPACE unchanged.
  `inst/CITATION` names TRexSelectorNeo, keeps "TRexSelector Team" as the package
  author and lists the eleven T-Rex method papers with text versions rendered
  from the master file (`tools/references/make_inst_citation.py`);
  `citation("TRexSelectorNeo")` prints them.
- Python: the selectors, `DummyDistribution`, `MemoryMappedMatrix`, the T-solver
  classes, the ML classes, the clustering enums and `agglomerative_cluster` end
  their docstring with a numpy-style "References" section. The T-solver and ML
  binding helpers now set a class docstring (the classes had none, only
  `__init__` had one); stubs regenerated with pybind11-stubgen 2.5.5 against a
  pybind11 3.1.0 build (an older pybind11 found by CMake drops the typed enum
  `__eq__` / `__ne__` overloads of the committed stubs; point `pybind11_DIR` at
  3.1.0 when regenerating).
- Repository: `CITATION.cff` at the root for GitHub's "Cite this repository"
  (preferred citation the T-Rex paper, the other method papers and Eigen as
  references); the READMEs acknowledge Eigen and point to the citation files,
  the Python and R READMEs link the journal DOI of the T-Rex paper.

### 2026-09-28

#### OpenMP optional in every build; one OpenMP runtime in the R package on macOS

- C++ and Python: OpenMP is optional, as in R. `find_package(OpenMP)` without
  REQUIRED; when found, `trex_utils` and `trex_tsolvers` link `OpenMP::OpenMP_CXX`
  as before, otherwise CMake reports "OpenMP not found: building single-threaded"
  and the stubs in `utils/openmp/utils_openmp.hpp` take over (pragmas ignored,
  `omp_get_max_threads()` 1). `-DTREX_WITH_OPENMP=OFF` (default ON) forces the
  single-threaded build; Python takes it as `-C cmake.define.TREX_WITH_OPENMP=OFF`.
  A single-threaded build adds `-Wno-unknown-pragmas` (MSVC `/wd4068`), so GCC's
  one warning per inactive pragma is gone. The Homebrew LLVM libomp pin applies
  only to Homebrew's clang and when that libomp exists; MSVC keeps `-openmp:llvm`.
- Installed package config: `find_dependency(OpenMP)` only when the library was
  built with OpenMP; `TRexSelector_WITH_OPENMP` tells consumers which. A
  single-threaded install exports no OpenMP link. The same for BLAS/LAPACK:
  `find_dependency(BLAS/LAPACK)` only for a `TREX_USE_BLAS_LAPACK=ON` build
  (`TRexSelector_WITH_BLAS_LAPACK`); before, every install required them.
- R on macOS: CRAN's R has no OpenMP flags, and configure linked Homebrew's
  libomp by absolute path while `-lomp` resolved to R's bundled runtime, so two
  OpenMP runtimes were mapped. configure now probes (1) R's flags, (2) on macOS
  R's `lib/libomp.dylib` with an omp.h from Homebrew or from the mac.r-project.org
  OpenMP package in /usr/local/include (`-Xpreprocessor -fopenmp`,
  `-L${R_HOME}/lib -lomp`, the runtime data.table and RcppArmadillo load),
  (3) Homebrew's runtime alone, (4) sequential. `TREX_OPENMP=no` forces the
  sequential build and adds `-Wno-unknown-pragmas`; the automatic fallback does
  not, since R CMD check lists that flag as non-portable. configure regenerated
  with autoconf 2.73 (AC_INIT now names version 2.0.0, not 0.2.0); Makevars.win
  unchanged.
- Tests: two OpenMP stub tests failed without OpenMP (status string, the wtime
  stub returns 0 on its first call) and are fixed; `OpenMPTest.SingleThreadedStubs`
  checks the one-thread stubs. CI: `test-no-openmp` (ubuntu-24.04) builds and
  tests with `TREX_WITH_OPENMP=OFF` and fails if any target compiles with -fopenmp,
  then runs the install smoke test with OpenMP disabled for the consumer
  (`TREX_SMOKE_CONSUMER_ARGS`, new in run_install_smoke.sh), so the config of a
  single-threaded install is proven not to ask for OpenMP.
- Verified: C++ 646 of 646 with OpenMP (22 skipped) and without (24 skipped, the
  three OpenMP thread-hygiene tests; no build warnings); R 1,407 with R's runtime
  (12 threads, one libomp mapped) and with TREX_OPENMP=no (1 thread); Python 635
  with and without OpenMP; vendor check OK; install smoke tests OK without
  OpenMP, and without OpenMP, BLAS and LAPACK (consumer with all three disabled).

#### T-solvers: the held heap bytes as `memoryBytes()`, R `$get_memory_bytes()`, Python `getMemoryBytes()`

- Found in the Examples memory demo (demo_ts_perf_03: p = L = 1,000, n up to
  2,000, T = 10): the peak resident set of a child process is no measure there.
  libomp resets the thread count in every forked child on macOS (the
  "single-threaded" children ran on 12 threads), and a spawned R or Python
  interpreter serves the sub-MiB solver allocations from pages it already holds
  (R read about 0 MiB per solver). The library now reports the memory itself, so
  the demo can run in one process per language.
- `TSolver_Base::memoryBytes()` (virtual, const, noexcept) sums the heap bytes of
  the state a solver holds: the y copy, residual, correlations, Cholesky factor,
  sparse path and running coefficients, diagnostics, index sets, the writer's
  means and norms and every reservation (vectors by capacity, Eigen objects by
  size). The classes with own containers add theirs (TLARS, TIENET, TIKNET, TOMP,
  TGP, TACGP, TMP, TOOLS, TIKOMP, TCDLASSO, TCDIENET, TCDTIKNET and the augmented
  wrappers); T-LASSO, T-Stepwise, T-Stagewise, T-ENET, T-NCGMP, T-AFS and T-CDENET
  inherit. Not counted: the borrowed X and D (readers and writers alike), the
  temporaries inside a step and the object itself (sizeof, about 3.3 kB with the
  tie-break engine). The coefficient lookup map is an estimate (bucket array plus
  four words per entry). A call is O(members + steps); nothing runs on the hot
  path, no member is added and checkpoints are unchanged.
- The report counts the up-front reservations as they are: `dummies_at_step_` and
  the lambda sequences (LARS, CD) hold 8 min(max_actives, row bound) + 1 entries
  (T-AFS divides by rho), 70 % of T-LARS at n = 2,000 (256 of 367 kB) and 94 % of
  T-AFS with rho 0.1 (1,279 of 1,361 kB). T-CDLASSO keeps its Gram rows and its
  continuation cache after a step; they are counted too. The augmented EN
  wrappers count their owned augmented copies (8 (n + p + L)(p + L + 1) bytes for
  TENETAug), the diagnostic mirror, the inner solver object and
  `inner_->memoryBytes()`; the inner X and D view the owned copies, so nothing is
  counted twice.
- `solver_utils/memory_bytes.hpp`: `bytesOf()` for vectors (capacity, nested
  heap), dense and sparse Eigen objects, strings and the hash-map estimate; the
  products go through `checked_bytes`.
- Bindings: R `TSolver_Base$get_memory_bytes()` returns a double (a byte count
  overflows R integers at 2 GiB), Python `getMemoryBytes()` exists on every
  solver. RcppExports, the Rd pages and the T-solver stubs are regenerated
  (pybind11 3.1.0, pybind11-stubgen 2.5.5); the vendored copy is re-synced.
- Against the allocator (macOS malloc zones, the 15 Examples demo solvers as
  precentered readers, p = L = 1,000, T = 10): `memoryBytes()` is 88-92 % of the
  bytes in use at n = 50 and 93-98 % at n = 2,000, the rest being the solver
  object and the allocator's block rounding. The peak inside a step is 1.04-1.67x
  the held state at n = 50 and 1.02-1.29x at n = 2,000.
- Tests: C++ 645 of 645. `test_tsolvers_memory` (52): exact hand counts at
  construction for T-LARS, T-OMP, T-OLS, T-AFS (rho 0.5), T-CDLASSO and TENETAug,
  growth bounds after a step, doubling n moves only the two n-vectors, a writer
  holds 16 (p + L) bytes more than a reader, every solver reports a positive and
  growing value, the checkpoint is byte-identical. `test_tsolvers_memory_peak`
  (21): a counting allocator (malloc zones on macOS, `__libc_*` on glibc, skipped
  elsewhere and under sanitizers) bounds peak / held per family; measured up to
  1.45 (LARS), 1.20 (OMP, AFS) and 4.3 (CD, whose certified localization builds
  the restricted Gram for an eigen-solve). Verified on Linux (GCC 13, glibc 2.39)
  as well. R 1407 expectations, Python 635.

### 2026-09-27

#### T-solvers: selection read from the coefficient path, the fixed default tie seed, IEEE Examples

- `cdiagnostics::selected_variables()`: the selection of a fitted solver as the
  T-Rex selector reads it, the variables with |beta_j| > machine epsilon at the
  last recorded step, split into predictors and dummies (active ones in
  active-set order, then the frozen ones). `print_selection()` and
  `print_selection_quality()` use it; they listed the active set, which
  under-reported T-Stagewise (one demo run: 32 active, 44 with a coefficient,
  four true variables near 2.7 among the frozen ones). Frozen variables are
  printed on a line of their own; for the other solvers the output is unchanged.
  The TSTAGEWISE_Solver documentation states the frozen coefficients.
- `SolverOptions::tie_seed` defaults to 0: direct solver runs are reproducible
  also where exact ties occur (T-LASSO at T = 50 changed in 1 of 500 demo runs
  between repetitions with the random_device default). A negative seed keeps the
  nondeterministic seeding; the T-Rex selectors pass their seeds explicitly.
- Cycling at saturation (T-Stagewise through removals and re-entries up to the
  step limit, T-LASSO more mildly) is part of these solvers and is left as is.
- The Examples compile their demos without `-ffast-math`, as the library: the
  header-inline library code they instantiate stays IEEE 754 conformant, and demo
  numbers are bit-comparable with the library builds.
- Solver identities over full paths (Examples cpp/tsolvers/README.md, measured on
  demo_ts_perf_01, 100 runs to T = 100): T-OMP = T-AFS (rho 1) = T-NCGMP (fully
  corrective) bit for bit; T-Stepwise the same steps as T-OMP (coefficients within
  1.1e-12); T-MP = T-NCGMP (line search) in every step (coefficients within
  2.1e-14); T-LASSO = T-LARS bit for bit up to the first T-LASSO drop, T-Stagewise
  within 2.7e-13 up to its first NNLS removal.
- Tests: C++ 572 of 572 (new: `SelectedVariablesReadTheCoefficientPath`,
  `TStagewiseSelection.FrozenCoefficientsAreSelected`,
  `DefaultTieSeedIsReproducible`), R and Python check the default seed.

#### T-solvers: T-NCGMP stopping rules, augmented EN degrees of freedom, install timestamps

- T-NCGMP line search is T-MP step for step (Locatello et al. 2017, Alg. 4 with
  L = 1) but ran without an active-set cap for 20 x (p + L) steps, so at T = 100
  in demo_ts_perf_01 it ended with about 106 selected variables against 54 for
  T-MP. It now stops by T-MP's rules (`effective_n_` distinct atoms,
  8 x min(`max_actives_`, `effective_n_`) steps); the fully corrective variant keeps
  T-OMP's. On the demo setup (100 runs, T up to 100) the line-search actions equal
  T-MP's in every run, the coefficients to 2e-14 (T-NCGMP divides the step by
  `||x_j||^2`, 1 up to rounding on normalized columns).
  `NCGMPVariantsMatchMPAndOMP` pins the full paths of both variants; the
  precentered contract test no longer exempts T-NCGMP from the n - 1 cap.
- TENETAug and TIENETAug reported the degrees of freedom of their inner solver,
  which runs intercept-free on the augmented system for lambda2 > 0, so the
  intercept was missing even with intercept = true. The wrappers add it
  (intercept or precentered), as the pathwise TENET and TIENET count it; the
  TENET / TENETAug and TIENET / TIENETAug parity tests compare the DoF.
- Install timestamps: `file(INSTALL)` gives every installed header the
  modification time of its source. A consumer object compiled after a header
  edit but before the install looked newer than the changed header, was not
  rebuilt, and the changed class layout corrupted its heap (seen in a scratch
  benchmark after the T-OMP / T-OLS members were added). `cmake --install` now
  records the SHA-256 of the header tree before pruning it and stamps new or
  changed files with the install time; unchanged headers keep theirs (verified:
  fresh install 106 files, repeat install 0, one altered header 1; the install
  smoke test passes).
- On a saturating design TCDLASSO can hit its sweep budget in `solveAtLambda`,
  after which its working active set runs ahead of the last recorded step
  (support 16, working set 18, n = 20); the recorded path and its DoF are
  consistent. Noted, not changed.
- R glue: the comment sweep had missed `R/TRexSelectorNeo/src/rcpp_*.cpp` and
  `TRexSelectorNeo_types.h`; 23 plan references, work-package tags and date
  stamps reformulated as the technical statement.

#### T-solvers: T-OLS scoring in one pass per step, the T-OMP refit, and `precentered`

- Found in the Examples performance demo (demo_ts_perf_01: n 100, L = p = 300,
  500 runs): T-OLS took 1.25 ms to T = 1 and 82 ms to T = 100, against 0.17 and
  2.50 ms for T-OMP. The scoring re-orthogonalized every candidate against the
  whole active set at every step (`h_j = X_A^T x_j` and a triangular solve per
  candidate), `O((p + L) n k^2)` per path; 97 % of the T-OLS time.
- T-OLS keeps the squared projected norms `||P_perp x_j||^2` of all candidates and
  downdates them by `(q^T x_j)^2` for each new direction `q = X_A R^{-1} e_k` of the
  active span, fused into the correlation pass: one `O(n (p + L))` pass per step,
  the modified Gram-Schmidt form of OLS (Blumensath & Davies, OMP vs. OLS) without
  an orthogonalized copy of the dictionary. 82 -> 2.9 ms to T = 100, 18.5 ->
  1.2 ms to T = 20, 1.25 -> 0.28 ms to T = 1 (1.3-1.6x T-OMP). The path equals the
  from-scratch score bit for bit up to n - 1 active variables (500 of 500 demo
  runs; pinned by `TOOLSScoring.DowndatedNormsMatchFromScratchProjection`); the
  norms are serialized, so a resumed solver continues bit for bit.
- T-OMP: `updateActiveSet()` no longer copies `R_` as a backup before
  `updateR()`, which never mutates it (T-AFS and TIKOMP already skipped the
  copy), and the refit caches `X_A^T y` in active-set order. 2.50 -> 2.13 ms to
  T = 100 at demo conditions (per core 2.30 -> 1.83 ms); paths bitwise unchanged
  (T-OMP, TOOLS and TIKOMP share the refit).
- `SolverOptions::precentered` (R `tsolver_options(precentered = )`, Python
  `SolverOptions(precentered=)`): the caller centered X, D and y. The solver
  counts the intercept as with `intercept = true` — the active set is capped at
  n - 1, the degrees of freedom and Cp include the intercept — without touching
  the data. Before, a reader on centered data allowed n active variables; at
  saturation the n-th entered as a noise pick through the machine-epsilon
  collinearity check and blew up the refit (100 demo runs at T = 100: T-LARS 98
  runs, coefficients up to 4.4x; T-OMP / T-AFS 76 runs, 7.6x; T-OLS 55 runs,
  coefficients 5e5). The T-Rex dispatch (`buildRunnerConfig`) and the GVS K-loop
  set it in every ScalingMode except NONE; the augmented EN wrappers pass it to
  the inner solver only on the lambda2 = 0 collapse (the augmented rows are not
  centered). The elastic nets with lambda2 > 0 and the line-search T-NCGMP keep
  running without the cap. Contract test for all 18 AD solvers:
  `PrecenteredCountsInterceptLikeIntercept` (same path, DoF and Cp as
  `intercept = true` on exactly centered data).
- With `precentered` left false every other solver path is bitwise unchanged
  (12 solver configurations, 100 demo runs, actions, DoF and coefficients).
- Tests: C++ 569 of 569, R 282 expectations, Python 611.

### 2026-09-26

#### CRAN readiness: the as-cran result, and the Linux installed size

- `R CMD check --as-cran` on macOS (R 4.6): Status 1 NOTE, the "new
  submission" one every first submission gets; tarball 1.7 MB, installed
  11.8 MB, examples (also `--run-donttest`), tests, both manuals and the
  vignette OK; install 323 s, vignette rebuild 158 s.
- The 292 MB installed size seen in the Ubuntu run was debug information,
  not code: R compiles with `-g` on every platform, macOS keeps the DWARF
  in the object files (the shared object is 6.8 MB), ELF platforms link
  the DWARF of the 48 Eigen-templated translation units into the shared
  object. The generated `src/Makevars.in` now strips it after linking on
  Linux (`all: $(SHLIB) strippedLib`, the CRAN-accepted post-link idiom of
  rstan and prophet; no compiler flag is overridden; a no-op elsewhere).
  Measured in the provisioned Ubuntu image: shared object 5.6 MB,
  installed package 7.1 MB. `Makevars.win` is unchanged; the templates in
  `bundle_for_cran.sh` are plain ASCII.
- Still to do before a submission: a precomputed vignette (2.6 min of
  simulations at check time), win-builder for Windows and R-devel, and
  the release from the public main repository.

### 2026-09-25

#### Linux verification of f390bfc (OrbStack): one budget, one locale

- The three Ubuntu jobs on the committed f390bfc: C++ 548 of 549, the R
  suite, `R CMD build`, `R CMD check` (the three known environment NOTEs),
  Python 610, the install smoke, and the cp312 manylinux_2_28 x86_64 wheel
  with its in-container tests — all green, except:
- `MLMemory.RidgeCV` measured 4.50 X on Ubuntu (Eigen 3.4) against the
  4 X budget set on macOS (Eigen 5.0): the Jacobi SVD keeps two internal
  copies of the fold under 3.4 (the scaled copy and the QR preconditioner)
  and one under 5.0. The budget is version-dependent now (4 X / 5 X).
- `R CMD check` reported a WARNING in the syntax check because a bare
  `ubuntu:24.04` image has no `en_US.UTF-8` locale; the stage script
  installs `locales`, generates it and runs the check under it.

#### ML methods: zero-copy harness and the read-only binding boundary (Tier A)

- Ground rule (Fabian): no binding copies X, D or a memory-mapped matrix;
  full copies of the design are avoided wherever they are not needed, and
  functionality is preserved at 100 %. Plan and tiers:
  notes/ZERO_COPY_PLAN.md; inventory: notes/R_PACKAGE_ASSESSMENT.md.
- Harness. `test_ml_methods_reference.cpp` runs every ML core on three
  fixed shapes and compares every output byte for byte against a reference
  recorded outside the repository (`TREX_ML_REF_DIR`; skipped when unset,
  recorded on first use); `tools/zero_copy_reference.{R,py}` do the same
  for the R and Python surfaces (`record` / `compare`). `test_ml_methods_memory.cpp`
  measures the peak-memory growth of each core in units of X on a wide
  500 x 8000 design and pins the current budgets: ridge CV 4.55 X, elastic
  net 1.09 X, elastic net CV 2.17 X, SVD 0.35 X, ridge 0.33 X, PCA 0.35 X,
  scalers 0.03 X (one process per test; ctest does that).
- R binding (Tier A, bit-identical): the eight ML entry points that took
  `const Eigen::MatrixXd&` — and therefore copied X in full on every call
  (RcppEigen converts only to a mutable map) — take the R matrix and build a
  read-only `Eigen::Map<const>` that binds to the cores' `Ref<const>`:
  `ridge_cv_fit`, `enet_fit`, `enet_fit_grid`, `enet_predict`, `enet_cv_fit`,
  `svd_compute`, `pca_transform`, `ridge_solve`. New memory-mapped entry
  points through the const map (read-only files accepted): `enet_fit_mmap`,
  `enet_fit_grid_mmap`, `enet_predict_mmap`, `enet_cv_fit_mmap`,
  `svd_compute_mmap`, `ridge_solve_mmap`, `pca_transform_mmap`, and
  `pca_create_mmap` / `pca_restore_mmap` (readwrite: the file is centered in
  place). The R6 classes route `mmap_matrix` inputs to them instead of
  reading the whole file through `as.matrix()`; `pca_create` pins X in the
  XPtr protection slot. R and Python references identical; new R tests
  compare every mmap path with the in-memory path bit for bit.
- Python (Tier A): `agglomerative_cluster` takes `Ref<const MatrixXd>`
  instead of `const MatrixXd&` (no copy of a Fortran-ordered array).

#### ML cores: copy fusion (Tier B, bit-identical)

- Ridge CV: the gathered training fold is standardized in place instead
  of into a second matrix (the two passes round exactly like the former
  fused expression): peak 4.55 X -> 3.64 X on the wide design; the Jacobi
  SVD's own copies remain.
- Elastic net: the gathered CV fold becomes the fit's working copy
  (`fit_owned` / `preprocess_owned`, centred and standardized in place)
  instead of being copied a second time: elastic-net CV 2.17 X -> 1.37 X.
  A fitted `enet_gaussian` releases its centred working copy after the
  path is solved (`predict` and the accessors read the p x K path only),
  so a fitted object holds p x K, not n x p.
- Ridge: the Cholesky factorization runs in place in the Gram / kernel
  matrix (`LLT<Ref<MatrixXd>>`): one p x p (or n x n) instead of two; the
  rare LDLT fallback rebuilds the matrix from X (a second pass, no copy).
- PCA: `PCA::readOnly(X, M)` fits a const matrix or a read-only
  memory-mapped file without writing (a named factory, so `PCA(X, M)` keeps
  centering and normalizing in place); `restore()` may be called before
  `fit()` (the preprocessing happens at construction; it used to throw),
  is idempotent, and checks that it is given the constructor's matrix; the
  destructor writes X back when `restore()` was never called, and the duty
  travels with a move (the T-solver rule; the Python docstring that
  claimed this is true now). X stays pinned in the R XPtr protection slot
  and by Python's keep_alive.
- Scalers: `fit()` takes a read-only reference, so the statistics can be
  fitted on a read-only memory-mapped file in R and Python (the in-place
  transforms still need a writable mapping).
- Clustering (R): the transpose copy on the in-memory path stays, on
  purpose — the engine clusters columns and the R surface clusters rows,
  and a lazy transpose would make every object a strided row; the comment
  that blamed undefined behaviour on the lazy transpose described an
  NN-chain bug the core fixed in June. The memory-mapped path reads the
  file through the const map (read-only mappings accepted).
- Verified: the 21 C++ reference cases, the R and Python references and
  all three suites are identical / green after every item.

#### ML cores: low-memory engines (Tier C, opt-in; the defaults stay bit-identical)

- Elastic net and its CV (`enet_gaussian(bool low_memory)`,
  `enet_cv_ccd(bool low_memory)`; R `MLElasticNet$new(low_memory = TRUE)`,
  `MLElasticNetCV$new(low_memory = TRUE)`; Python
  `ElasticNet(low_memory=True)`, `ElasticNetCV(low_memory=True)`): the
  coordinate descent runs on the raw X through the caller's view — a
  read-only memory-mapped file included — with the centring carried as a
  scalar residual offset (folded into the residual once per sweep, which
  keeps the training residual mean-zero and the dot products free of the
  cancellation a large column mean would cause) and the standardisation as
  a per-column scale; the CV folds are 0/1 masks and the held-out rows are
  predicted row by row. No allocation of the size of X: 1.09 X -> 0.09 X
  (path) and 1.37 X -> 0.15 X (CV) on the wide design. Same estimator,
  same grid, same deviance rules; the two engines agree to rounding (about
  1e-12 in the coefficients; the tests pin 1e-7), so a CV pick can differ
  only on a flat curve. The default engine is unchanged and bit-identical.
- Ridge CV (`ridge_cv_svd(bool low_memory)`; R `MLRidgeCV$new(low_memory
  = TRUE)`; Python `RidgeCV(low_memory=True)`): the dual (kernel) fold
  engine. The n_tr x n_tr kernel of the standardised training fold and the
  n_test x n_tr cross kernel are accumulated from column blocks of X
  standardised on the fly (read in place), eigen-decomposed, and the ridge
  path is evaluated on the held-out rows as K_te U diag(1/(s^2+lambda)) U^T y
  — no gathered fold, no p-dimensional coefficients, the SVD route's rank
  rule. 3.64 X -> 0.19 X on the wide design; memory O(n^2) instead of
  O(n p), the route for p >> n. Agreement with the SVD route to 1e-7 on
  the CV curve and identical picks in the tests; the kernel carries the
  squared condition number, so on ill-conditioned tall designs the default
  SVD route is the more accurate one.
- SVD dispatch thresholds exposed (`SVDSolver(min_dim_thresh, block_size)`
  in Python, `MLSVD$new(min_dim_thresh, block_size)` in R; defaults
  unchanged): lowering `min_dim_thresh` selects the randomized, block-wise
  path — two passes over X, no copy — for smaller designs; it is the route
  for a large memory-mapped X and is accurate for a decaying spectrum (the
  tests pin 1e-4 on the retained singular values), not for a flat one.
- Ridge `solve` on a memory-mapped X needs no new engine: the Gram or
  kernel matrix is accumulated from the map (Tier A's `ridge_solve_mmap`).
- Peak-memory budgets (`test_ml_methods_memory.cpp`), wide 500 x 8000:
  ridge CV 3.64 X (kernel engine 0.19 X), elastic net 1.09 X (raw engine
  0.09 X), elastic-net CV 1.37 X (raw engine 0.15 X), SVD 0.35 X, ridge
  0.31 X, PCA 0.35 X, scalers 0.03 X.

#### R surface catch-up and small utilities (5.4, 5.5)

- T-solvers (R6 base): `get_residuals()`, `get_beta_path()` (step 0
  included, Python's layout), `get_cp_with_sigma()`, `solver_type_to_string()`,
  `get_dropped_indices()`, `get_active_predictor_indices()`,
  `get_active_dummy_indices()`, `get_num_actives()`, `get_dummy_start_index()`,
  `is_connected()`, `set_tolerance()`, `set_tie_seed()`,
  `set_exchangeable_tie()`, `get_warnings()` (the already-preprocessed
  warning is readable now), `clear_warnings()`. Family getters:
  `TLARS_Solver$get/set_kahan_refresh_interval()`, `TENET_Solver$get_alpha()`,
  `TIKNET_Solver` / `TIKOMP_Solver$get_tikhonov_k()` (a `dgCMatrix`),
  `TCDSolver_Base$get_kkt_violation()`. `reconnect` is not bound in R: an
  R6 solver is bound to its matrices at construction.
- Selectors: active bindings `tFDR`, `K`, `max_num_dummies`, `max_T_stop`,
  `stagnation_check`, `num_selected`, `num_dummies` (new C++ getter
  `TRexSelector::getNumDummies()`, also bound in Python as
  `TRexSelector.num_dummies`); DA `cor_coef`; GVS `hc_method_used`.
- ML: `MLElasticNet$n_nonconverged()`; Python `RidgeCV.interpolates()`,
  `RidgeCV.selection_at_grid_edge()`, `ElasticNetCV.selection_at_grid_edge()`,
  `ElasticNetCV.grid_truncated()`.
- Utils: `trex_get_max_threads()` / `trex_set_num_threads()` (the OpenMP
  thread count of the backend); `mmap_matrix_size` returns a double (an
  int overflowed above 2^31 - 1 elements); `mem_matrix_from()` coerces an
  integer matrix explicitly instead of through a hidden binding copy.

#### Selectors: memory-mapped Biobank and SPCA, restore(), results wrapped once

- Biobank screening and sparse PCA on a memory-mapped X: the R6 classes
  used to read the whole file into RAM through `as.matrix()`. New entry
  points `trex_biobank_screening_{1d,2d}_mmap_create` and
  `trex_spca_select_mmap` map the file directly (readwrite: the inner
  selectors normalize / center it in place and restore it); the R6 classes
  route `mmap_matrix` inputs to them. New R tests compare the mmap and
  in-memory selections exactly.
- `TRexSelector::restore()` (C++, R `$restore()`, Python `restore()`):
  gives X back before `select()` — the constructor normalizes X in place,
  and until now an unselected selector handed it back only from its
  destructor (in R and Python: the garbage collector). Idempotent, a no-op
  after `select()`; a restored selector refuses to run (`select()` throws),
  including the screening selector's own `select()`. Inherited by every R6
  selector class.
- Results wrapped once: the selector's `phi_mat`, `fdp_hat_mat` and `r_mat`
  (T_stop x p) were copied twice per access; the DA `*_array_BT` results
  (rho x T_stop x p) three times (a vector of matrices by value, the wrap,
  `simplify2array`). They are now wrapped once from the const reference,
  the DA arrays built directly as R arrays of the same shape.
  `TRexSPCA::select()` returns a const reference to the stored result, so
  the bindings wrap Z and V once (`auto res = spca.select()` still copies).
- `trex_da_to_base_ptr` (unused from R) pins the owning pointer in its
  protection slot instead of returning an unprotected non-owning view.

#### T-solvers: shared-buffer claim and the in-place writer/reader contract

- Finding (Examples repo, `demo_ts_05_tstagewise` in C++): the checkpoint
  round trip reported "Paths differ" from step 21 on. Not a solver bug: the
  demo constructed the reference solver, the checkpoint solver and the
  resumed solver with `normalize = true, intercept = true` on the SAME
  Eigen maps, so X and D were centered and scaled in place three times. A
  second in-place transformation is not idempotent (it moves the data by
  rounding, RSS difference 9e-13), and the stagewise near-ties downstream
  resolve differently. Pre-scaling once in front and running the solvers
  as readers gives bit-identical paths.
- The contract, now written on `TSolver_Base::restore()` (C++), on every
  Python constructor and on `tsolver_options()` / `TSolver_Base` (R): a
  solver constructed with `normalize = true` or `intercept = true` is a
  WRITER that transforms X and D in place, keeps the statistics, and holds
  a claim on both buffers until `restore()` or destruction, both of which
  write X and D back; a solver constructed with both false is a READER that
  runs as constructed and never touches the buffers; a solver resumed from
  a checkpoint is a reader that expects the buffers in the saved state (so
  the writer that saved it stays alive, or the pipeline pre-scales once in
  front and runs readers); y is never written (the solver centers an
  internal copy).
- The claim is executable: `tsolvers/ad/solver_utils/buffer_registry` (a
  mutex-guarded set of buffer pointers, RAII `BufferClaim`) is shared by the
  T-solvers and the T-Rex selector, which moved its private registry there.
  A second writer on a claimed buffer, T-solver or selector, throws
  `std::runtime_error` (Python `RuntimeError`, R error) with the remedy:
  pre-scale once and construct readers, `restore()` the holder, or let it
  go out of scope. A writer that finds its input already centered and
  unit-norm (every column within 1e-10) records a warning (`getWarnings()`)
  and continues.
- `restore()` is idempotent and a no-op for readers and resumed solvers
  (`holdsPreprocessed()`); the former unconditional inverse scaled and
  shifted the data a second time on a second call. It now also checks that
  the passed X and D are the buffers the solver was constructed on. The
  destructor writes back a writer that was never restored, which also
  covers a derived constructor that throws after the base preprocessed
  (`TENET` with a negative `lambda2` used to leave X normalized). The
  write-back duty travels with a moved solver.
- R binding rewritten zero-copy: the R6 T-solvers used to copy X and D into a
  private holder (`TSolverRcpp`, also for a memory-mapped X, which was read
  into RAM). They now map the R matrix (double storage) or the mmap file
  directly, pin the R objects in the XPtr protection slot, and follow the
  contract: a writer preprocesses the caller's matrix or file in place, the
  new `$restore()` writes it back and releases the claim, the GC finalizer
  writes back a writer that was never restored, `$load()` into a writer
  keeps its write-back duty (the archived statistics are the same). The
  stale copying holder in `TRexSelectorNeo_types.h` (a second definition of
  the same struct, included by `RcppExports.cpp`) is the single definition
  now.
- R `get_actions()` fixed: the C++ solver already records the 1-based
  signed lars codes, and the R6 method shifted them by one more, so
  `r_index = TRUE` returned 2-based codes and `r_index = FALSE` the 1-based
  ones it documented as 0-based. `r_index = TRUE` (default) now returns the
  codes as recorded; `r_index = FALSE` returns the 0-based magnitudes
  (`abs(a) - 1`, sign dropped: a 0-based signed code cannot distinguish the
  addition from the removal of the first variable).
- Tests: `test_tsolvers_buffer_guard.cpp` (registry primitives; second
  writer throws across families; intercept-only and normalize-only claim;
  `restore()` idempotent, readers and resumed solvers never write back;
  destruction writes back, also after a throwing derived constructor; moves
  transfer claim and duty; readers and loaded solvers never claim; the
  already-preprocessed warning; the selector/solver shared registry; the
  pre-scaled shared-buffer checkpoint round trip is bit-identical). The CD
  and Tikhonov checkpoint tests keep the writer alive across the resumed
  run. R: `test-tsolvers-inplace.R` (in-place change visible on the R matrix
  and in the mmap file, second writer errors, readers coexist, `$restore()`
  idempotent and path kept, GC write-back, `$load()` into a writer restores,
  pre-scaled readers bit-identical); `test-tsolvers.R` restores each writer
  and runs the 22-solver checkpoint round trip on pre-scaled readers with
  exact equality. Three Python tests that rebound a writer in a loop now
  release it first.

#### Linux verification (OrbStack): two platform-bound tests, broken Rd pages

- The three Ubuntu CI jobs were reproduced in an `ubuntu:24.04` container
  (GCC 13, Eigen 3.4, R 4.3.3, Python 3.12) on the committed tree, plus the
  install smoke, `R CMD check` and the cibuildwheel manylinux_2_28 x86_64
  wheels with their in-container test run. Everything passed except:
- Two C++ tests calibrated on macOS data: `std::normal_distribution` and
  `std::rand` (behind `Eigen::Random`) are implementation-defined, so the
  same seed draws different data under libc++ / macOS libc and libstdc++ /
  glibc. `EnetCvCcdTest.KeepsGlmnetLambdaRuleAndDetectsGridFloorPin` now
  asserts the wide-versus-pinned ratio with a margin both platforms clear
  (1.2; observed 1.9 on macOS, 1.32 on Linux), and
  `TSolverTikTest.SerializeLoadContinueAcrossDropRecovery_TIKNET` probes a
  small seed range until a design yields a checkpoint mid-recovery instead
  of relying on seed 12. Both pass on both platforms.
- Three R help pages were unparseable: an apostrophe in
  `\code{lambda2 * beta' K beta}` opens an R string literal in Rd (the
  TIKNET, TIKOMP and TCDTIKNET pages; now `t(beta) K beta`), a bare
  `[0.25, 4]` became a broken link and `\ldots` sat inside a `\code{}` on
  the `trex_control` page. CI never saw this because the R job installs
  with `--no-docs` and does not run `R CMD check`; the check step is added
  to the workflow.

#### Python: one subpackage per dummy mechanism, and the Tikhonov selector

- `trex_selector_neo.tsolvers.ad.{lars_based, omp_based, afs_based, cd_based}`
  with `SolverOptions`, `ScalingMode` and `GroupPenaltyVariant` under
  `tsolvers.ad`, and `trex_selector_neo.trex_selector_methods.ad` with every
  selector, mirroring `cpp/src/{tsolvers, trex_selector_methods}/ad`; the
  compiled extension carries the same level (`_core.tsolvers.ad`,
  `_core.trex_selector_methods.ad`). The package root keeps re-exporting the
  AD selectors as the everyday spelling; the `vd` / `sd` families will be
  reachable only through their subpackage. Breaking for
  `from trex_selector_neo.tsolvers.lars_based import ...` and
  `from trex_selector_neo.tsolvers import SolverOptions` (insert `.ad`); R
  needs nothing (flat namespace, group words).
- `TRexTikhonovSelector` is bound: `TRexTikhonovControlParameter` (the p x p
  penalty as a dense `tikhonov_K` array converted to sparse internally, or
  `GROUP_CENTERING`; the lambda_2 rule; the pilot knobs; cluster-aware
  dummies), the enums `TikLambda2Method`, `PilotStatistic`,
  `TikPenaltyStructure`, the records `Lambda2PilotResult` / `Lambda2PickInfo`,
  the builders `gamma_to_k` and `build_group_centering_k`, and the wrapper
  with `lambda2_used`, `lambda2_cv_diagnostics` and `lambda2_pilot`. The base
  control merges as in R (a non-Tikhonov `solver_type` warns and becomes
  TCDTIKNET). Every selector R has is now in Python.

#### GVS K-loop: per-experiment tie-break seed

- The in-memory GVS K-loop now derives a deterministic `tie_seed` per
  experiment from the selector seed, as the base ExperimentRunner does.
  Inert unless the exchangeable-tie band is on (`exch_tie_alpha > 0`); with
  it on, a seeded GVS run is reproducible instead of drawing tie-breaks from
  `random_device`.

### 2026-09-24

#### Voting ties: exact vote counts and a single-division voting grid

- The classical experiment runner accumulated Phi as K additions of 1/K. The
  rounding residue (15 x 0.05 = 0.75 + 1 ulp) decided the strict `Phi > v`
  comparison at the ties the voting lattice makes routine, because the grid
  `0.5, 0.5 + 1/K, ...` lives on the same 1/K lattice by design (R's
  `seq(0.5, 1 - eps, by = 1/K)`). The runner therefore selected a variable
  that sat exactly on a grid point; the GVS K-loop, which divides exact
  counts once, and R's `rowMeans()` did not. One such tie flips the L-loop
  stop, which is how `TRexTikhonovSelector` and `TRexGVSSelector` came to
  stop at different L on the same dummies (the divergence pinned since
  2026-08-25; the dummies and the per-experiment paths were bit-identical).
- Fix: `ExperimentRunner::aggregateResults` accumulates exact counts and
  divides by K once (Phi is the correctly rounded m/K), and
  `createVotingGrid` computes every point as one division (K + 2i)/(2K), so
  a mathematical tie yields the identical double and is excluded, as the
  definition and R do. Affects every selector on the runner (classical, DA,
  screening, SPCA, Tikhonov) only where a variable's vote share equals a
  grid point exactly; such variables are no longer selected at that
  threshold. The equivalence test beyond L = p is enabled (decision 9 in the
  plan is closed: there was no convention difference).

#### Beyond architecture: one folder per dummy mechanism (breaking, C++ only)

- `cpp/src/tsolvers/` and `cpp/src/trex_selector_methods/` each hold three peer
  folders, `ad/` (augmented dummies: the classical T-solvers and every T-Rex
  selector, with `trex_utils`), `vd/` (virtual dummies: `vd/` solvers, `trex_vd`)
  and `sd/` (sparse dummies: `sd/` solvers, `trex_sd`). The families share no
  code, and no class serves more than one mechanism; a fourth mechanism would be a
  fourth sibling. Placement rule: a file lives at the level of its lowest common
  user (`TSolver_Base`, `SolverOptions`, `solver_utils/` and all of `trex_utils/`
  are AD-only today; `VD_Base` and `SDTSolver_Base` are standalone roots).
- Namespaces follow the folders: `trex::tsolvers::ad::{linear_model::…,
  solver_utils::…, TSolver_Base, SolverOptions, ScalingMode, …}`,
  `trex::tsolvers::vd`, `trex::tsolvers::sd` (the SD solvers were declared inside
  the AD subnamespaces before; `sd_calibration` is `sd::calibration`),
  `trex::trex_selector_methods::ad::{trex_core, trex_da, trex_gvs, trex_screening,
  trex_spca, trex_tikhonov, utils::…}`, `…::vd::trex_vd`, `…::sd::trex_sd`.
- Umbrellas: `tsolvers/ad.hpp` (was `lm_tsolvers.hpp`), `tsolvers/vd.hpp`,
  `tsolvers/sd.hpp`; `tsolvers/tsolvers.hpp` still means the classical family.
  Installed header paths move with the tree (`TRexSelector/tsolvers/ad/…`,
  `TRexSelector/trex_selector_methods/ad/…`).
- The AD alias for the solver dispatch is `disp` (was `sd`, which now names the
  sparse-dummy family). Tests and validation programs mirror the layout
  (`test_tsolvers/{ad,sd}/`, `test_trex_selector_methods/{ad,vd}/`,
  `validation/{tsolvers,trex_selector_methods}/ad/`).
- R and Python: the glue includes and namespaces follow; no exported name changes.
  `bundle_for_cran.sh` excludes the `vd` and `sd` folders and umbrellas (four
  patterns, was five).

#### R namespace: one prefix per group (breaking)

- The exported symbols now follow one rule: functions and constant lists carry the
  group prefix (`trex_`, `tsolver_`, `ml_`, `mem_`), R6 classes the group word
  (`TRex*Selector`, `T*_Solver`, `ML*`). Renamed: the ML classes `ZScoreScaler`,
  `LpNormScaler`, `RidgeSolver`, `RidgeCV`, `ElasticNet`, `ElasticNetCV`,
  `SVDSolver`, `PCA` to `MLZScoreScaler`, `MLLpNormScaler`, `MLRidge`, `MLRidgeCV`,
  `MLElasticNet`, `MLElasticNetCV`, `MLSVD`, `MLPCA`; `agglomerative_cluster`,
  `cut_tree`, `LinkageMethod`, `DistanceMetric` to `ml_agglomerative_cluster`,
  `ml_cut_tree`, `ml_linkage_method`, `ml_distance_metric`; the metrics
  `compute_fdp`, `compute_tpp`, `compute_precision`, `compute_recall`,
  `compute_fdp_dense`, `compute_tpp_dense` to `trex_fdp`, `trex_tpp`,
  `trex_precision`, `trex_recall`, `trex_fdp_dense`, `trex_tpp_dense`;
  `gamma_to_k` to `trex_gamma_to_k`; `mmap_matrix` and `convert_to_memory_mapped`
  to `mem_matrix` and `mem_matrix_from` (the S3 class `mmap_matrix` and its methods
  are unchanged). `TSolver_Base`, `TCDSolver_Base` and `get_ptr` are internal now.
- Every help page lists its group ("Other T-solvers: ..."), `?TRexSelectorNeo`
  presents the four groups, and `_pkgdown.yml` indexes the reference by them.

#### The CD family, native TIENET, TIKNET and TIKOMP as R6 and Python solvers

- R: `TCDLASSO_Solver`, `TCDENET_Solver`, `TCDIENET_Solver`, `TCDTIKNET_Solver`
  (sharing the CD getters through `TCDSolver_Base`), `TIENET_Solver`,
  `TIKNET_Solver`, `TIKOMP_Solver`, all on `$new(X, D, y, [groups | tikhonov_K,]
  options)`; a Tikhonov matrix is a `Matrix::dgCMatrix`, another `Matrix` sparse
  class, or a dense matrix.
- Python: the `tsolvers.cd_based` submodule with the four CD solvers;
  `TIENET_Solver` and `TIKNET_Solver` in `lars_based`, `TIKOMP_Solver` in
  `omp_based`; a Tikhonov matrix is a dense array (converted internally, no SciPy
  dependency).

#### C++ install layout (breaking for manual include paths)

- The SDK headers now install under a namespaced root,
  `<prefix>/include/TRexSelector/{ml_methods,trex_selector_methods,tsolvers,utils}`,
  instead of directly under `<prefix>/include/`. Consumers that link the exported
  targets (`find_package(TRexSelector)`, `TRexSelector::trex_selector_methods`, ...)
  need no change; a manual `-I<prefix>/include` becomes `-I<prefix>/include/TRexSelector`.
  The pre-install prune that kept the header tree free of stale files now removes only
  that namespaced root, and only when it carries the `TRexSelector.install-stamp` written
  by a previous install: a foreign directory of the same name is refused, never deleted.
  Previously the prune deleted the generic `include/utils`, `include/ml_methods`, ...
  under the prefix, which in a shared prefix (conda, /usr/local) could belong to
  another package. `cpp/tests/install_smoke/run_install_smoke.sh` installs into a
  scratch prefix and builds a `find_package` consumer against it. Its first run
  found seven headers that the code includes but no install file set carried
  (`lambda2_cv_diagnostics.hpp`, `post_selection_criteria.hpp`,
  `trex_cluster_dummies.hpp`, `trex_cluster_hac.hpp`, three `sd_tsolvers` headers):
  an installed SDK could not compile a consumer since August. Fixed.
- New CMake option `TREX_INSTALL_CPP_SDK` (default ON): headers, static archives and
  the CMake package files. The Python build forces it OFF, so wheels no longer ship
  about 16 MB of C++ SDK into `site-packages/include` and `site-packages/lib`.

#### Hardening

- Agglomerative clustering rejects non-finite input at its entry with a message naming
  the entry (a NaN distance made the nearest-neighbour search index out of range);
  the NN chain also refuses its own sentinel. A zero-norm column keeps its documented
  convention (maximally dissimilar, distance 1).
- `cut_tree` / `cut_tree_by_height` validate a merge matrix before indexing it: row
  count, cluster ids that exist at their step, no self-merge, no id consumed twice,
  finite and (for the height cut) non-decreasing distances. Reachable through the
  public R and Python `cut_tree`.
- Overflow-checked byte counts (`utils/checked_arith.hpp`) for the memory-mapped
  matrix, the tiled distance matrix and the binary Eigen serialization: a wrapping
  `rows * cols * sizeof` throws `std::overflow_error` before any file is touched.
- Eigen's process-global thread count is set once per experiment runner from the
  calling thread and restored on exit, never written inside an OpenMP team. Before,
  the sequential runners left it changed in the host process (an R or Python session)
  and the workers raced on it. Only observable in builds without `EIGEN_USE_BLAS`
  (the R package and the Python wheel); with a BLAS, Eigen's count is inert. The
  thread-hygiene tests watch two globals: the host thread's OpenMP thread count in
  every build (the parallel runners set a per-worker count inside their team, and it
  must not reach the caller), and Eigen's count in no-BLAS builds, where it is live.

#### Scaling mode NONE: the user-managed design

- `ScalingMode::NONE` joins `L2` (default, unchanged) and `ZSCORE`. Under NONE the
  selectors take `X` and `y` exactly as supplied: no centering, no scaling, the
  stored normalization parameters are the identity (means 0, divisors 1, `y_mean` 0)
  and nothing is written back into the caller's buffer. The caller pre-centers and
  pre-scales as the model demands, or supplies an intercept-free dictionary
  (steering vectors, Vandermonde, Fourier atoms) whose columns must not be centered
  and may contain a constant atom. The dummies keep the unit-L2 standard, so a NONE
  design should be unit-L2-commensurate; one warning names how many columns have an
  L2 norm outside [0.25, 4] (z-scored or raw-amplitude data). A design pre-normalized
  with the library's own helper and run under NONE selects exactly what the raw design
  selects under L2 (bitwise, test-enforced in C++, R and Python). The solver-side
  convention and the lambda_2 tuners treat NONE as unit-L2; the correlation HAC of the
  cluster-dummy selectors then measures cosine coherence on the columns as supplied.
- The selectors' preprocessing is callable: `trex_normalize_design(X, y, scaling_mode)`
  in R and `trex_selector_neo.normalize_design(X, y, scaling_mode)` in Python return
  the normalized copies and the parameters. R exposes the mode as
  `trex_control(scaling_mode = "L2" | "ZSCORE" | "NONE")` (it had no scaling switch).
- The normalizer rejects bad input in every mode with the column named (0-based):
  a non-finite entry in `X` or `y`, and a dead column — constant under L2 / ZSCORE,
  all-zero under NONE (`std::invalid_argument`; `ValueError` in Python; an R error).
  Before, a constant column was silently pinned to scale 1 and could never be selected.
  Designs that carry constant columns must drop them now (breaking, deliberate:
  PLAN WP4 4.1 (d)). The matrix is unmodified when the constructor throws.

#### T-solver constructors take one options structure (breaking)

- Every T-solver is constructed as `Solver(X, D, y, options)`, with
  `trex::tsolvers::SolverOptions` (`tsolvers/solver_options.hpp`) carrying every other
  input: the preprocessing flags (`normalize`, `intercept`, `scaling_mode`, `verbose`),
  the tolerance, the tie-break seed and the exchangeable-tie configuration, and the
  family hyperparameters (`lambda2`, `rho_afs`, `ncgmp_variant`, `tenet_aug_use_lars`,
  `kappa_dummy`, `fold_dummies`, `group_penalty`, the `cd_*` knobs). The structural
  inputs that are data, not options, stay explicit: the group vector of the IEN solvers
  (`TIENET`, `TIENETAug`, `TCDIENET`) and the Tikhonov matrix of the sparse-K solvers
  (`TIKNET`, `TIKOMP`, `TCDTIKNET`) follow `y`. Each family reads the fields it needs;
  the defaults are the former solver-level defaults (`lambda2` defaults to 0, the lasso
  path, where it used to be a required argument). The positional constructors with
  their boolean triple and trailing enum are gone; the `NCGMPVariant` and
  `GroupPenaltyVariant` enums moved next to the structure (their old names still
  resolve). The T-Rex layer maps its `SolverHyperparameters` onto the structure in one
  place of the dispatch, which lost its per-family constructor table and its
  post-construction setters; selections are unchanged (the full suite is bitwise the
  same on the anchors).
- R: `tsolver_options()` builds the list every solver's `$new(X, D, y, options = ...)`
  takes (`TIENETAug_Solver$new(X, D, y, groups, options = ...)`); the per-solver
  arguments `normalize`, `intercept`, `verbose`, `lambda2`, `rho`, `variant` and
  `use_lars_inner` are gone (`rho` is `rho_afs`, `variant` is `ncgmp_variant`).
- Python: `trex_selector_neo.tsolvers.SolverOptions` with a keyword constructor
  (`SolverOptions(lambda2=0.5, normalize=False)`), passed as the fourth positional or
  the `options=` keyword; `GroupPenaltyVariant` exported next to it. The former
  keyword arguments of the solver constructors are gone.

#### lambda_2 tuner boundary guard (warning-only)

- The GVS and Tikhonov selectors warn once when a CV tuner returns a
  boundary-dominated lambda_2, and return the pick unchanged. Two modes are
  covered: the returned lambda_2 sits on the first or last point of its grid (the
  grid boundary, not a CV optimum: a flat or still-improving curve, typically low
  SNR), and, for the IEN tuner, the profiled lambda_1 optimum of the chosen
  lambda_2 sits on the truncation edge of the lambda_1 path. Before, both were
  visible only in `n_profile_at_edge`, which is blind to the grid-ceiling mode
  (measured in August on the demo designs: a pick at the ceiling with the counter at
  zero). The warning names the selector, the method, the mode and the remedies
  (1SE over MIN, a fixed lambda_2, a longer lambda_1 path, a finer grid). For the
  EN CCD tuner a last-point pick on a path that glmnet's deviance rules truncated
  is benign and is not flagged.
- The lambda_2 diagnostics record (C++ `lambda2_cv_diagnostics`, R
  `lambda2_cv_diagnostics`, Python `Lambda2CVDiagnostics`) carries the pick:
  `grid_index` / `grid_size`, `at_grid_edge`, `profile_at_edge` and
  `boundary_pick`. Every CV method fills the grid fields now, not only the IEN
  tuner. `ienet_cv_ccd` gained `selection_at_grid_edge()` and a per-lambda_2
  `profile_at_edge(k)`; its `diagnostics()` takes the rule (MIN or 1SE) whose pick
  the caller uses.

### 2026-07-21

#### R package refreshed: vendored core synced, exchangeable tie-break exposed

- The vendored C++ core is byte-identical with `cpp/src` again (tie-break
  additions in `tsolver_base`, the solver dispatcher, TOMP/TAFS, plus the
  cosmetic trex_da drift). `trex_control()` gains `exch_tie_alpha` (default 0)
  and `exch_tie_floor` (default 0.5), parsed into `SolverHyperparameters` by
  the Rcpp layer; man pages regenerated via roxygen. New testthat coverage:
  control defaults and tie-broken TOMP/TAFS selection runs on a correlated
  cluster design. Full testthat suite green.

#### Python package refreshed for the exchangeable tie-break API

- `SolverHyperparameters` now exposes `exch_tie_alpha` / `exch_tie_floor`, and the
  solver wrappers (TOMP/TAFS/TLARS families) gain `setExchangeableTie(alpha, floor=0.5)`;
  type stubs updated accordingly (the PRIOR_GROUPS binding changes shipped with the
  core rewrite). New tests: hyperparameter defaults/assignment, reproducible
  tie-broken paths + alpha = 0 no-op on TOMP/TAFS, data-driven rho grid
  (ascending, rho = 1 anchor) and negative-label rejection on the prior-groups
  DA path. Suite 453/453.

### 2026-07-17

#### TLARS "No valid gamma" warning classified and demoted

- Investigated the frequent `No valid gamma found, using fallback` warnings seen in the
  trex_da demos. Instrumentation on the exact demo DGPs/seeds showed **100 % of observed
  events are the benign path-exhaustion case**: on small designs
  (p + num_dummies <= n - 1, e.g. demo 04's p=25 cells) the T-loop drives the LARS path
  until every real + dummy variable is active, and `computeStepSize()` then runs with an
  empty inactive set — the terminal LS step `Cmax / A_A` is the correct step. A Gaussian
  control reproduced the events at the same rate as the heavy-tailed DGP (173 vs. 166 per
  40 trials), and the events vanish exactly at the exhaustion boundary — not heavy tails,
  not collinearity, not the demo RNG (seeds are `base + mc`, fresh mt19937 per trial).
  Demo 05 (p_total=500) emits none; its earlier attribution was console interleaving.
- `computeStepSize()` now logs the empty-inactive-set case verbose-only
  (`logInfo("Inactive set empty; taking terminal LS step.")`) and keeps a `logWarning`
  (with candidate count) only for the never-observed all-candidates-filtered case.
  Fallback step unchanged. Suite 291/291; R mirror synced. Details in
  `TRex_Research/documentation/CTRexSelector_Open_Tasks.md` §7 (CTLARS should gain the
  same classified logging when synced — it is currently silent in both cases).

#### PRIOR_GROUPS redefined: constrained sub-clustering (deflation-mismatch fix)

- `setupDA_PriorGroups()` no longer uses the user-supplied groups as deflation
  neighbourhoods. The nearest-partner penalty `delta = 2 - min|Phi_j - Phi_k|`
  (Definition 1 of the T-Rex+DA paper) assumes tight neighbourhoods in which
  the closest-Phi member is a genuine shadow; inside a large prior group nearly
  every variable has SOME member with an almost identical Phi, so delta -> 2
  uniformly and the shadow/active ranking is left invariant (the FDR-plateau /
  TPR-cap signature documented in
  `TRex_Research/documentation/Prior_Groups_Deflation_Mismatch_DA_TRex.md`).
- New semantics: the prior groups are hard merge CONSTRAINTS on the dependency
  structure. The finest common refinement of all supplied `prior_groups` levels
  partitions the variables; HAC (`hc_linkage`, correlation distance) runs
  within each constraint group (merges across boundaries are structurally
  impossible); the rho grid (length `hc_grid_length`) is subsampled from the
  pooled within-group dendrogram heights and terminated by the conservative
  rho = 1 singleton anchor (the paper's Psi = 1/2 empty-group penalty). The
  BT-style 3D calibration then operates on tight data-driven neighbourhoods
  INSIDE the known groups — the coarsest grid point recovers at most the whole
  prior group, the finest is provably conservative.
- `rho_grid_labels` removed from `TRexDAControlParameter`, the R
  `trex_da_control()`, and the Python bindings — the grid is now data-driven.
  `prior_groups` labels must be non-negative (validated). Linkage validation
  now covers PRIOR_GROUPS as well as BT.
- New guard tests in `cpp/tests/test_trex_selector_methods/test_trex_da.cpp`:
  boundary containment + rho = 1 anchor, finest-common-refinement semantics,
  full select() smoke run, negative-label rejection.

### 2026-07-15

#### Exchangeable tie-breaking for greedy solvers (DA FDR fix)

- New solver hyperparameters `exch_tie_alpha` (default 0 = off) and
  `exch_tie_floor` (default 0.5) in `SolverHyperparameters`, consumed by
  TOMP/TAFS via `TSolver_Base::applyExchangeableTieBreak()`. When alpha > 0
  and the step's top candidate is a real predictor, inactive real predictors
  with |corr| >= floor to the top candidate whose correlation gap lies within
  `alpha * ||x_j* -/+ x_k|| * ||r|| / sqrt(n)` are treated as statistically
  exchangeable, and ONE member is picked uniformly at random per step (via
  the existing tie-break RNG, so `setTieSeed()` keeps runs reproducible).
  Dummy picks are never randomized — dummy budget and early stopping are
  unchanged.
- Motivation: greedy solvers are winner-take-all — within a T-Rex random
  experiment the in-sample winner of a highly correlated cluster is a
  deterministic function of (X, y) (dummies never perturb the real-variable
  ranking), so occurrence mass concentrates on one cluster member and the
  DA deflation `delta = 2 - min|Phi_j - Phi_k|` degenerates to the identity:
  collinear shadows that beat their active in-sample survive as false
  positives (DA-AR1 FDR ~ intrinsic shadow-win rate, e.g. ~34% per cluster at
  rho = 0.9, SNR = 2, n = 300 — above tFDR = 0.2). Randomizing statistically
  indistinguishable picks restores the within-trial occurrence spread across
  exchangeable cluster members that the DA FDR control relies on (LARS-type
  solvers produce this spread through their path geometry and are unaffected).
- Validated on the demo_trex_da_01 AR(1) cells (rho = 0.9, gap = 100 and
  Random support): alpha = 0.25 restores FDR < tFDR for TOMP/TAFS while
  retaining more power than TLARS; the safe regimes (rho <= 0.7) are
  essentially unchanged. Recommended: `exch_tie_alpha = 0.25` for greedy
  solvers under trex+DA.
- The dispatcher re-applies the parameters on all three solver lifecycle
  paths (fresh, serialized warm start, in-memory warm start), mirroring
  `setTolerance`; the fields are intentionally not serialized.

### 2026-07-08

#### R <-> Python binding parity

- BREAKING (R): `TRexBiobankScreeningSelector$select()` now returns the same
  per-phenotype record shape as the Python binding. A matrix `Y` previously
  returned a list with a `$statistics` data.frame plus parallel top-level
  `$selected_indices` / `$selected_indices_screen_ordinary` /
  `$selected_indices_screen_bootstrap` lists; it now returns a plain list of
  records, one per phenotype, where `res[[i]]` has the identical shape to the
  single-phenotype (vector `Y`) return. No information is lost — every
  `$statistics` column and index list is a field of each record — and the
  summary table is a one-line projection:
  `do.call(rbind, lapply(res, function(r) as.data.frame(r[c("phenotype_index", "estimated_FDR",
  "method_used", "used_fallback_trex")])))`.
  This matches Python, which returns one `BiobankScreenTRexResult` (vector) or a
  list of them (matrix).
- R gains `prior_groups` + `rho_grid_labels` in `trex_da_control()` — the
  prior-groups dependency-aware path (`setupDA_PriorGroups`), which Python
  already exposed; the R glue previously could not route to it. `hc_grid_length`
  is now actually read by the DA glue (it was silently inert before).
- R gains `cv_seed` / `cv_n_folds` / `cv_n_lambda` in `trex_gvs_control()` — the
  cross-validation fold controls for the CV `lambda2_method` paths (`cv_seed`
  enables per-trial fold reproducibility), at parity with Python.
- `TRexGVSSelector` gains read-only `$groups` (1-based per-variable cluster
  labels, length `p`) and `$group_labels`, so cluster diagnostics (e.g. purity)
  can be computed R-side; Python exposes the same via `GVSSelectionResult`.

### 2026-07-07

#### Shared-buffer safety guard (new runtime check)

- `TRexSelector` (and every derived selector) normalizes the input matrix `X`
  in place via an `Eigen::Map` and restores it when `select()` finishes or on
  destruction. Two selectors mapped onto the SAME buffer therefore silently
  corrupted each other: the first to restore left the buffer un-normalized, so
  the second ran on raw data (typically selecting everything). The core now
  claims a buffer when it normalizes it and throws a clear error if a second
  live selector is constructed on the same matrix. Fix each caller by giving
  each selector its own copy of `X`, or by running selectors one at a time
  (construct, `select()`, and let each go out of scope before the next). The
  guard is in the C++ core, so C++, R, and Python all benefit; internal
  sequential orchestrators (SPCA per-PC, biobank screen->fallback) are
  unaffected because `select()` releases the claim before the next selector
  is built.

#### R <-> Python binding parity, elastic-net model selection

- Elastic-net model selection is now exposed in BOTH languages: `ElasticNet`
  (glmnet-style coordinate-descent coefficient path, backed by `enet_gaussian`)
  and `ElasticNetCV` (K-fold CV lambda.min/1se selection, backed by
  `enet_cv_ccd`). `RidgeCV` (SVD-path ridge CV, `ridge_cv_svd`) is added to
  Python, which previously exposed no model-selection CV at all. Verified
  R and Python produce bit-identical elastic-net paths on shared data.
- R gains `da_method = "NN"` in `trex_da_control()` — the nearest-neighbour
  correlation-threshold DA sweep (`DAMethod::NN`), which Python already
  exposed; the R glue previously rejected it. `hc_grid_length = 0` lets the
  core auto-select `min(20, p)`.
- R gains `tenet_aug_use_lars` in `trex_gvs_control()` and
  `trex_spca_control()` — the LARS-vs-LASSO inner-path switch for `TENET_AUG`
  (Python already had it). Documentation clarifies that `en_solver`
  (`TENET` vs `TENET_AUG`) is a solver-EQUIVALENCE knob (both select the same
  variables), NOT the LARS-vs-LASSO distinction.
- `TENETAug_Solver` / `TIENETAug_Solver` are exported from the R package (were
  internal-only), with plain and mmap creators, 1-based group accessors, and
  save/load.

#### GVS / SPCA / clustering fixes

- BREAKING (R): `gvs_type = "IEN"` now derives the `TIENET_AUG` solver
  (matching the 2026-07-05 solver rework and the C++/Python wrappers); the R
  wrapper previously forced the old `TLASSO` solver.
- `TRexGVSSelector` now emits a warning when `control$solver` / `solver_type`
  is set to a value that does not match the solver derived from `gvs_type`
  (the field is silently overridden).
- `lambda_2 == 0` (degenerate no-ridge / pure T-LASSO) now emits a warning —
  the `-1` auto / `0` degenerate / `> 0` fixed sentinel convention had tripped
  in-house callers.
- Ward + `Correlation_LSH_Approx` is now permitted in both the R and Python
  clustering glues (the core dispatcher supports it, running Ward in a 64-D
  SimHash-projected Euclidean space; the glues wrongly rejected every
  non-Euclidean metric under Ward).

#### Python API

- `TRexBiobankScreeningSelector`: the control parameter `biosctrex_ctrl` is
  renamed to `bio_ctrl`, and `BiobankScreenTRexControl` now exposes its nested
  `trex_screen_ctrl` (and that control's nested `trex_ctrl`), so Python can set
  K, solver, screening method, etc., at parity with the R knobs.
- Scalers gain returning, layout-agnostic `transform()` / `fit_transform()` /
  `inverse_transform()` methods alongside the strict zero-copy
  `*_inplace` variants (which require a writeable Fortran-ordered float64
  buffer), removing the sharpest edge in the Python scaler API.

#### Core fixes

- `TENETAug_Solver` / `TIENETAug_Solver` deserialization: loading into an
  already-constructed solver (the R/Python bindings' `load()` path) now rebuilds
  the augmented maps and reconnects the inner solver unconditionally on load,
  instead of gating on a null-map check that only held for fresh-object loads —
  fixing checkpoint resume through the bindings.
- `cmake --install` now prunes the package-owned header trees before copying,
  so a renamed or deleted header can no longer poison the local install tree
  (a stale header that `#include`d a since-renamed file previously broke every
  consumer until the install was cleaned by hand).

### 2026-07-05

#### T-solver review fixes (tsolvers module)

- Full code review of `cpp/src/tsolvers/` with correctness fixes: stale
  correlations after T-Stagewise NNLS drops, missing dummy de-normalization in
  T-ENET's `getBeta`/`getBetaPath`, a data race and a scale-aware
  near-constant-column drop threshold in the preprocessing, guarded
  `getBeta(step)`/`getCp()` accessors, and a consistent
  actives/inactives/dropped partition across all solvers. `restore()` no
  longer corrupts the caller's `y` (the response is copied internally and was
  never centered in place).
- API: `executeStep()` with `T_stop = 0` now computes the full solution path as
  documented; `setTieSeed()` and `getWarnings()`/`clearWarnings()` are public
  and exposed in Python; warnings print unconditionally and are recorded
  per solver; TENETAug mirrors its inner solver's diagnostics, validates
  `lambda2`, and supports `save()`/`load()`.
- Performance: RAII guard restores Eigen's global thread count after solver
  runs; beta-path storage grows geometrically instead of per-step; dead
  Cholesky backup copies removed.
- T-ACGP now applies the conjugate correction across active-set growth steps
  via zero-padded previous directions (Blumensath & Davies 2008).
- Internals deduplicated (correlation refresh, tie scans, save/load,
  cycling-ratio tracking); tie-breaking RNG state, tolerance, and warnings are
  serialized. NOTE: the checkpoint format changed — `.bin` files saved with
  earlier builds cannot be loaded.

#### DA lifecycle refactor and NN-setup parallelisation

- `TRexDASelector` no longer overrides `select()`: the dependency-structure
  setup moved into an `onSelectBegin()` override and the DA selection /
  result widening into a `finalizeSelectionResult()` override, so DA now
  runs the inherited base pipeline like GVS does (behavior-neutral; the
  duplicated scaffolding — voting-grid setup, loop drivers, cleanup,
  X denormalization — existed twice before and could drift).
- The DA-NN neighbourhood setup is OpenMP-parallelised: each variable's
  neighbour list is built independently (race-free row-wise scan instead of
  the symmetric pair loop), with an early exit over the ascending rho grid.
- `corr_max` (GVS) documents that the default 0.5 matches the R reference
  and that the EUSIPCO 2022 paper's GWAS experiments used rho_thr = 1/3.

#### T-Rex+GVS rework: solver-encapsulated augmentation and the TIENETAug solver

- New `TIENETAug_Solver` (tsolvers, LARS family): the Informed Elastic Net
  (Machkour et al., CAMSAP 2023, Theorem 1) as a terminating forward selector.
  Builds the group-informed row-augmented system `[X; B]`, `[D; B tiled per
  dummy layer]`, `[y; 0_M]` (B(m, j) = sqrt(lambda2)/sqrt(p_m) on group-m
  columns) internally from a 0-based group assignment, applies the R-faithful
  post-augmentation column scaling, and runs an inner T-LASSO (or pure T-LARS)
  — full save/load support, registered in the polymorphic solver tests.
- `TRexGVSSelector` now hands every solver the plain (n x p) design, the plain
  (n x L) cluster-MVN dummy block, and the centered response; all augmentation
  lives in the solvers. `GVSType::IEN` routes through `TIENETAug_Solver`
  (new `SolverTypeForTRex::TIENET_AUG`, also in Python; IEN previously
  required TLASSO — breaking config change), and `en_solver = TENET_AUG`
  routes through `TENETAug_Solver`, which finally makes the
  `tenet_aug_use_lars` R-parity flag functional (it was silently ignored).
  The inline augmentation machinery (X_aug_/y_aug_ members, five build
  helpers, the variable effective row count, and the EN-aug mmap
  reinterpretation) is deleted; memory-mapped dummy files now always have n
  rows. EN-aug thereby switches from post-augmentation re-normalization to
  the demo-verified as-constructed convention that reproduces the Gram-based
  TENET path exactly (asserted by a new equivalence test).
- Review fixes: the GVS solver cache is cleared BEFORE the dummy buffers it
  views are reassigned (transient-dangling-view hardening, mirroring the
  base-class contract); EN dummy normalization goes through the shared
  data-normalizer (degenerate columns warn and keep scale 1.0 instead of
  throwing).
- New end-to-end GVS tests (grouped-recovery for EN/TENET, EN/TENETAug,
  EN/TENETAug with pure-LARS inner, IEN/TIENETAug, plus TENET==TENETAug
  selection equivalence). The IEN test pins lambda_2: the ridge-CV value
  (identical to the R reference's, which shares one CV and one * p/2
  conversion between EN and IEN) can make the sqrt(lambda2)-scaled group
  rows dominate the augmented geometry on n > p designs, leaving nothing
  selectable — a documented sensitivity of the IEN track (see
  computeLambda2()).

#### T-Rex selector review fixes (trex_selector_methods module)

- Warm starts now actually engage: `WarmStartManager` never stored solvers or
  checkpoint paths, so every T-loop step silently re-solved from scratch. The
  manager now retains solvers in-memory (co-retaining the dummy matrix for
  the PERMUTATION/DIRECT strategies, whose D_k is a per-call temporary) and
  registers serialized checkpoints for the memory-mapped path; `invalidate()`
  keeps the K slots repopulatable instead of shrinking them away. The T-loop
  now continues each solver's path incrementally, as in the R reference;
  results are unchanged (warm/cold equality is asserted by new tests). The
  checkpoint temp directory is only created for serialized mode.
- DIRECT / PERMUTATION_DIRECT dummy generation honours the user seed: it
  previously derived from the PERMUTATION-only base seed (always 0, so
  different seeds produced identical dummies) and, when unseeded, drew fresh
  `random_device` entropy on every call — changing the dummies between
  T-steps. The base seed is now resolved once per `DummyGenerator` and
  reused, making unseeded runs internally consistent and seeded runs
  reproducible across instances.
- Documentation: the intentional `Phi_prime` [0, 1] clamp (deviation from the
  R reference / paper to avoid negative FDP estimates) now documents its
  conservative/anti-conservative asymmetry; the DA-BT deflation cites
  Definition 1 of the T-Rex+DA paper (empty-group penalty Ψ = 1/2) and the
  deliberate R-parity choice of terminal-Φ deltas vs. the paper's per-step
  penalties (Eq. (10)).
- Tests: warm-start engagement and warm/cold equality (in-memory, DIRECT, and
  an end-to-end memory-mapped vs. in-memory run that asserts the T-loop
  actually continues), DIRECT seed handling; fixed a checkpoint filename race
  between the typed serialization tests under parallel ctest.
- `tloop_stagnation_stop` is now tri-state (`std::optional<bool>`, Python:
  `bool | None`). The default (auto) resolves by solver family: DISABLED for
  the equiangular LARS-path solvers (TLARS/TLASSO/TENET/TENET_AUG), which
  terminate properly on their own — matching the R reference and the T-Rex
  paper, neither of which has a stagnation stop — and ENABLED for greedy
  solvers, which otherwise run into a noise trap and select until the
  maximum number of dummies is reached. An explicit true/false overrides the
  family default.
- Cleanups from the review follow-up: the voting-grid sentinel uses machine
  epsilon instead of the solver tolerance (a large user tol with K >= 100
  could make the grid non-monotone; also matches the R reference grid); the
  AR1 window fallback for an exactly-zero correlation estimate is the
  formula's limit (±1 window) instead of the full window; the
  equi-correlation estimator excludes near-constant columns from the pair
  count; warm-start solvers are invalidated before the dummy matrices they
  view are mutated (STANDARD/HCONCAT), with the buffer-reallocation contract
  documented on the DummyGenerator; dead validation branches removed; the
  Yule-Walker vs. `arima(..., method="ML")` auto-estimation difference for
  the AR1 `cor_coef` is documented on the parameter and covered by a
  known-rho recovery test.

## TRexSelector 2.0.0

### 2026-06-11

#### Architecture

- Complete architectural rework: C++ core library extracted into a standalone
  CMake project (`TRexSelector`) with full `find_package()` support; R and
  Python bindings are now thin language layers on top.
- C++20 standard adopted throughout; Eigen3, Boost, Cereal, and OpenMP are
  managed via vcpkg or system package managers.
- Simulation and demo executables separated into a standalone sibling project
  (`TRexSelector_Simulations`); the core library no longer depends on Matplot++.

#### Selector variants

- `TRexSelector` — base FDR-controlled selector for linear regression.
- `TRexDASelector` — dependency-aware variant for correlated predictors.
- `TRexGVSSelector` — grouped variable selection.
- `TRexScreeningSelector` — fast FDR-controlled variable pre-screening for
  large *p*.
- `TRexBiobankScreeningSelector` — high-throughput genomic / biobank screening.

#### T-algorithm solvers

- LARS-based family: T-LARS, T-LASSO, T-ENET, T-Stagewise, T-Stepwise, TIENET.
- OMP-based family: T-OMP, T-GP, T-ACGP, T-MP, T-NCGMP, T-OLS.
- AFS-based family: T-AFS.
- All solvers support Cereal polymorphic serialization for warm-start workflows.

#### CMake / build system

- Four exported CMake targets:
  `TRexSelector::trex_utils`, `TRexSelector::trex_ml_methods`,
  `TRexSelector::trex_tsolvers`, `TRexSelector::trex_selector_methods`.
- `cmake --install` produces a fully relocatable install tree consumed via
  `find_package(TRexSelector REQUIRED)`.
- `TREX_PORTABLE_BUILD` option disables `-march=native` for portable wheel
  builds (PyPI / CI).
- `TREX_USE_BLAS_LAPACK` option controls BLAS/LAPACK backend (default ON;
  disabled automatically for Python wheel builds).
- CMakePresets.json provides `debug`, `release`, `release-libs`, and
  `vcpkg-release` presets.

#### R package (version 2.0.0)

- New R6 object-oriented API replaces the previous functional `trex()` API.
- C++ backend is compiled via Rcpp/RcppEigen; no longer depends on the `tlars`
  R package.
- R6 classes: `TRexSelector`, `TRexDASelector`, `TRexGVSSelector`,
  `TRexScreeningSelector`, `TRexBiobankScreeningSelector`.
- Helper functions: `compute_fdp()`, `compute_tpp()`.
- Full testthat test suite and a comprehensive vignette added.
- `R CMD check --as-cran` passes with 0 errors, 0 warnings, 1 NOTE (private
  GitHub URL, expected for new submission).
- Compiler flags set to `-O3 -funroll-loops` in Makevars.in / Makevars.win for
  improved runtime performance.

#### Python package (version 0.1.2)

- pybind11 bindings built via scikit-build-core.
- API mirrors the R6 interface; indices are 0-based.
- Classes: `TRexSelector`, `TRexDASelector`, `TRexGVSSelector`,
  `TRexScreeningSelector`, `TRexBiobankScreeningSelector`.
- 355 pytest tests covering all selectors, all solvers, ML utilities, and
  memory-boundary correctness.
- GIL released during long-running C++ computations.
- `pyproject.toml` contains full PyPI metadata (authors, classifiers, URLs).

#### Documentation

- Root `README.md` updated with C++ developer installation guide, full
  `find_package` consume instructions, and published Elsevier paper references.
- R `README.Rmd` and vignette updated with published journal references
  (Signal Processing, 2025).
- `inst/CITATION` lists both published papers.

---

## TRexSelector 0.1.0

### 2026-05-11

First installable release of the C++ TRexSelector library.

#### Package infrastructure

- `cmake --install` produces a fully relocatable install tree:
  `include/`, `lib/`, `lib/cmake/TRexSelector/`.
- Downstream projects can use `find_package(TRexSelector REQUIRED)` and link
  against `TRexSelector::trex_selector_methods`.

#### CMake changes

- All four library targets declare `FILE_SET HEADERS` for correct header
  installation and IDE indexing.
- Single authoritative `install(EXPORT TRexSelectorTargets)` in root
  `CMakeLists.txt`.
- `CMakePresets.json`: added `release-libs`, `vcpkg-release` presets.

#### CPack packaging

- `cmake --build <dir> --target package` produces binary archives.
- `cmake --build <dir> --target package_source` produces source archives.

#### vcpkg dependency manifest

- `vcpkg.json` declares eigen3, boost-headers, boost-iostreams, cereal, and
  libomp (macOS only) as managed dependencies.

#### Licensing

- GNU General Public License v3 added.
