# Benchmark Report for *HierarchicalMatrices*

## Job Properties
* Time of benchmark: 19 Mar 2020 - 9:53
* Package commit: dirty
* Julia commit: 2d5741
* Julia command flags: `-O3,--check-bounds=no`
* Environment variables: `JULIA_NUM_THREADS => 2` `OPEN_BLAS_NUM_THREADS => 1`

## Results
Below is a table of this job's results, obtained by running the benchmarks.
The values listed in the `ID` column have the structure `[parent_group, child_group, ..., key]`, and can be used to
index into the BaseBenchmarks suite to retrieve the corresponding benchmarks.
The percentages accompanying time and memory values in the below table are noise tolerances. The "true"
time/memory value for a given benchmark is expected to fall within this percentage of the reported value.
An empty cell means that the value was zero.

| ID                                      | time         | GC time   | memory        | allocations |
|-----------------------------------------|-------------:|----------:|--------------:|------------:|
| `["HMatrix", "assembly", "CPU1"]`       | 3.110 s (5%) | 54.785 ms | 1.29 GiB (1%) |      996137 |
| `["HMatrix", "assembly", "CPUThreads"]` | 1.877 s (5%) | 40.769 ms | 1.29 GiB (1%) |     1038052 |

## Benchmark Group List
Here's a list of all the benchmark groups executed by this job:

- `["HMatrix", "assembly"]`

## Julia versioninfo
```
Julia Version 1.3.1
Commit 2d5741174c (2019-12-30 21:36 UTC)
Platform Info:
  OS: Linux (x86_64-linux-gnu)
      Ubuntu 19.10
  uname: Linux 5.3.0-42-generic #34-Ubuntu SMP Fri Feb 28 05:49:40 UTC 2020 x86_64 x86_64
  CPU: AMD Ryzen 3 2200G with Radeon Vega Graphics: 
              speed         user         nice          sys         idle          irq
       #1  3073 MHz      25531 s        246 s       3018 s    3564076 s          0 s
       #2  3166 MHz      19170 s        712 s       2878 s    3572194 s          0 s
       #3  3182 MHz      21972 s         16 s       2859 s    3571620 s          0 s
       #4  3075 MHz      28288 s        284 s       2968 s    3564911 s          0 s
       
  Memory: 5.798274993896484 GB (727.52734375 MB free)
  Uptime: 36028.0 sec
  Load Avg:  1.6435546875  1.34375  1.00146484375
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-6.0.1 (ORCJIT, znver1)
```