# Parallelism: who owns the cores {#parallelism}

[TOC]

## Two layers of work

A T-Rex run has two layers that could use several cores. The outer layer is the \f$K\f$ random
experiments, which are independent of each other and share the design. The inner layer is the work
inside one T-solver: the passes over the design that refresh the correlations, the search for the
next knot of the path, the standardization. The rule of the library is that one layer owns the cores
at a time. Either the experiments run in parallel and every solver runs serially, or the experiments
run in sequence and a solver may split its passes. Nested teams of threads are never formed.

The default is serial. A selector runs its experiments one after the other, every solver runs on one
thread, and a T-solver used on its own runs serially. Parallel work is opt-in through two fields of
TRexControlParameter: `max_outer_threads` for the experiments and `max_inner_threads` for the
solvers' passes.

## The default thread budget of the library

The library is single-threaded by default on every machine. Every parallel region runs on an
explicit budget, and the methods without an option of their own (the scalers, the clustering, the ML
regressions, the normalizer, the dummy and data generators) take the library's default budget, which
starts at 1. In the namespace `trex::utils::openmp`, `get_num_threads()` reads that default,
`set_num_threads()` sets it for the process and returns the previous value, and the thread argument
or the `set_num_threads()` member of those classes overrides it per object. Inside a selector the
same kernels take `max_inner_threads`. `get_max_threads()` reports the threads OpenMP can run on
this machine, an upper bound for a budget. In a build without OpenMP every setting is a no-op and
the machine count is one. The first section of the example applies the two-thread rule of this
guide, so the first value below is 2.

\snippet{trimleft} parallelism.cpp includes

\snippet{trimleft} parallelism.cpp budget
\output{parallelism.budget.txt}

The selectors and the solvers do not read the default budget. Their budget comes from the control.

## The solvers' passes

`TRexControlParameter::max_inner_threads = m` gives every solver inside the random experiments a
budget of `m` threads for its passes over the design (the correlation refreshes, the step-length
search, the standardization). The default is 1. A pass is split only when it touches at least a
million entries of the design, because the cost of forming a team of threads exceeds the gain below
that, so small problems run serially whatever the budget says. The seeded generators (the dummies,
the synthetic designs) count their work in the same unit with the measured cost of a random draw
(about 80 entries) and of seeding an engine (up to 32,000), so they split from a few thousand
entries. The budget is live when the experiments run in sequence. Inside parallel experiments the
solvers stay serial. A T-solver used on its own takes the same budget through
`SolverOptions::n_threads`, and `setNumThreads()` changes it on a live solver.

## Parallel random experiments

`parallel_rnd_experiments = true` with `max_outer_threads = m` runs the random experiments on `m`
threads (at most \f$K\f$). Inside the parallel experiments every solver runs serially, so the number
of threads in use is `max_outer_threads` and nothing more. Both fields are whole numbers of at least
1, and a value below 1 throws. The lean L-loop strategies (`SEEDED`, `PERMUTATION`,
`PERMUTATION_SEEDED`) run their experiments in sequence in every case and give the cores to the
solvers.

\snippet{trimleft} parallelism.cpp data

\snippet{trimleft} parallelism.cpp parallel
\output{parallelism.parallel.txt}

## Determinism

A seeded selector derives one seed per random experiment from `seed` and the experiment index, so
the dummies of experiment \f$k\f$ do not depend on which thread draws them or on how many threads
there are. The selection of a seeded run is therefore the same with experiments in sequence and in
parallel, and the same for any thread count of the inner passes. Every reduction inside the solvers
is summed in a fixed order for the same reason.

## Combining with your own parallel loops

A Monte Carlo study usually parallelizes over trials, with an OpenMP loop or with threads of the
standard library. Then the trials own the cores, and the library should run serially inside every
worker, which is the default. Two details matter.

- The selectors and the solvers take their budget from their control, so the default control inside
  a worker is serial. The kernels on the library's default budget never form a team inside a
  parallel region of the caller, so a default budget above one does not oversubscribe the machine
  either.
- Every trial constructs its own data and its own selector, and draws its own dummies with
  `seed = -1`, so the trials are independent and no object is shared between the workers.

\snippet{trimleft} parallelism.cpp mc
\output{parallelism.mc.txt}

## Builds without OpenMP

The library takes OpenMP from the toolchain when it is configured, and the CMake option
`TREX_WITH_OPENMP=OFF` forces the single-threaded build where OpenMP exists. `find_package()` sets
`TRexSelector_WITH_OPENMP`, so a consuming project can tell which build it links. A single-threaded
build selects the same variables, it only runs the experiments and the passes in sequence.

## The CRAN rule

CRAN checks the R package on at most two cores. The pages of this guide and the tests of the library
follow the same rule with the idiom of the first section of this page, which caps the default budget
for the duration of the program and restores it at the end.

\snippet{trimleft} parallelism.cpp restore-budget
