# Memory mapping: disk-backed designs {#memory_mapping}

[TOC]

## A design larger than the memory

A memory-mapped matrix lives in a file on disk, and the operating system pages the parts in use into
memory. The library binds such a file exactly as it binds an Eigen matrix, so every selector, every
T-solver, the scalers, the low-memory regression engines and the clustering accept one. That is how
a sparse regression problem with five million variables was solved on a laptop
\cite scheidt2023bigtrex. The file holds the matrix as raw doubles in column-major order, the layout
of Eigen, of R and of Fortran-ordered NumPy arrays, without a header.

The file is mapped through Boost.Interprocess \cite gaztanaga2008interprocess, the library of the
Boost C++ Libraries for memory-mapped files and shared memory. The library binds any buffer, mapped
or not, so the selectors never see where the design lives, as Figure 1 shows.

\image html memory_agnostic.png "Figure 1: A design in memory and a design in a file reach the same implementation through one mapped view." width=90%

## Creating and opening a mapped matrix

The class MemoryMappedMatrix, a template on the scalar type and `MemoryMappedMatrix<double>` here,
creates or opens a file with its dimensions and the access mode. Its `getMap()` is an `Eigen::Map`
of the mapping, a zero-copy view that reads and writes the file, so assigning an Eigen matrix to the
map writes the file. The headers and names of this page:

\snippet{trimleft} memory_mapping.cpp includes

\snippet{trimleft} memory_mapping.cpp create
\output{memory_mapping.create.txt}

The handle supports the element access one needs to inspect and to edit the file, and the map gives
the block access of Eigen. Reading through the map touches the pages of the file, the library never
copies the design. A file opened in `AccessMode::ReadOnly` refuses element writes with an exception.

\snippet{trimleft} memory_mapping.cpp access
\output{memory_mapping.access.txt}

## Streaming a design into the file

A design that does not fit in memory is written column by column. The file is raw doubles in
column-major order, so writing one column after another with an output stream produces it.

\snippet{trimleft} memory_mapping.cpp stream
\output{memory_mapping.stream.txt}

## Selectors on a mapped design

The selectors take the map of the handle where they take `X`. The file must be open in
`AccessMode::ReadWrite`, because the selector centers and scales the design in place and restores it
when `select()` returns (the page \ref no_copies "Working without data copies" has the contract).

\snippet{trimleft} memory_mapping.cpp select
\output{memory_mapping.select.txt}

The same holds for TRexDASelector, TRexGVSSelector, ScreenTRexSelector, BiobankScreenTRex, TRexSPCA
and TRexTikhonovSelector. `TRexControlParameter::use_memory_mapping = true` additionally keeps the
dummy matrices of the random experiments in mapped files under the library's temporary directory
(`trex::utils::get_global_temp_dir()`, which `set_global_temp_dir()` moves), for the case where
\f$K \cdot L \cdot p\f$ dummy columns do not fit either.

## T-solvers and ML methods

A T-solver binds the mapped `X` as it binds a matrix, the scalers transform the file in place, the
low-memory engines of the regressions read it in place, and the clustering accepts it as data.

\snippet{trimleft} memory_mapping.cpp solver
\output{memory_mapping.solver.txt}

## When it pays

A mapped design pays when the design does not fit into memory, when several solvers share one design
and the working set should be loaded once, and when a design is produced by a streaming process and
read afterwards. For a design that fits comfortably in memory an Eigen matrix is faster, since no
paging happens. Mapped files are never read into the process heap by the library, so the memory
footprint of a run is the library's working data, see `memoryBytes()` on the T-solvers. The example
program removes its folder of mapped files at the end, after the mapped objects are gone.

## Examples

The [TRexSelector_Examples](https://github.com/trexselector/TRexSelector_Examples/tree/main)
repository holds runnable demos and Monte Carlo studies of memory-mapped designs, with their results
and a README per demo, in
[R](https://github.com/trexselector/TRexSelector_Examples/tree/main/R/memory_mapping/demo_memory_mapping),
[Python](https://github.com/trexselector/TRexSelector_Examples/tree/main/Python/memory_mapping/demo_memory_mapping)
and
[C++](https://github.com/trexselector/TRexSelector_Examples/tree/main/cpp/memory_mapping/demo_memory_mapping).