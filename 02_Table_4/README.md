# Table 4: fine-mapping power and FDR

Table 4 of the main text compares single-trait (JAM expanded) and multi-trait
(flashfm) fine-mapping when the effect estimates come from OLS, St, DSh, or
GSR. The design is that of Section 5.1 and Appendix C, with 100 replications
in the *IL2RA* region.

The fine-mapping relies on two pieces of code by Jennifer Asimit, the
flashfmZero package
([jennasimit/flashfmZero](https://github.com/jennasimit/flashfmZero)) and the
JAM and flashfm functions with a user-chosen estimator in `flashfm-savvySh.R`
([jennasimit/flashfm-savvySh](https://github.com/jennasimit/flashfm-savvySh)).
The scripts in this folder install both from GitHub and add the estimators of
the paper, the simulated traits, and the power and FDR. They take the
population of 10,000 individuals as an input file, simulated with HapGen2 as
described in Appendix C. HapGen2 sets its seed from the clock, so every
simulated population is different, and a run with a new population gives
values that differ from those of Table 4.

```
Rscript run_simulation.R 1     # replication 1
Rscript run_simulation.R       # all 100 replications
Rscript make_tables.R
```

| File | What it does |
|---|---|
| `run_simulation.R` | draws the data and fine-maps both traits with each estimator, one CSV per replication in `results/` |
| `make_tables.R` | power and FDR from the replications, written to `tables/` as `Table_4.tex` and `Table_4.csv` |
| `table4_functions.R` | the estimators in the form flashfm-savvySh expects, the trait draws, and power and FDR |

`run_simulation.R` takes replication numbers, for example `1`, `1-20`, or
`1-20 45`. Replications run in parallel on at most `max_workers` workers
(`../R/simulation_tools.R`), one replication per worker.


## Requirements

- Java on the PATH. JAM runs its sampler in Java through the R2BGLiMS code
  in flashfmZero.
- flashfmZero at commit `009fe64` (13 February 2025), the version available
  when Table 4 of the current paper was computed. Commit `feb79cd` (17 April
  2025) switched on, in `flashfmU()`, the flashfm adjustment for the
  correlation between traits, which changes the multi-trait rows.
  `run_simulation.R` installs the package from GitHub into `flashfmZero_lib/`
  on the first run. It has C++ code, so it needs a compiler (Rtools on
  Windows), and its imports lqmm, rlist, gtools, and gplots are installed with
  it when they are missing. When the code is loaded, `table4_functions.R`
  adapts one loop in the internal function `groupmultiU()`, which groups
  correlated variants, to the case where no pair of groups is left to merge.
  The installed package is not changed, and the variant-level MPPs of Table 4
  do not depend on the grouping.
- `flashfm-savvySh.R`, which `run_simulation.R` downloads at commit `ad7a101`
  on the first run, and savvySh, which `flashfm-savvySh.R` loads and which is
  installed from CRAN when it is missing. The estimators themselves are those
  of `../R/All_eight_estimators.R`.
- remotes, future, and future.apply.
- The population file `data/IL2RA_population.rds`.


## Population

As in Appendix C, the population has 10,000 individuals simulated with
HapGen2 from the CEU haplotypes of the 1000 Genomes Phase 3 data in the
*IL2RA* region. `run_simulation.R` reads it from `data/IL2RA_population.rds`
as a matrix of genotype scores (0, 1, 2) with one row per individual and the
rsIDs as column names. It has to contain the causal variants rs61839660 (C1),
rs62626317 (C2), and rs11594656 (C3). Variants with MAF of 0.005 or less are
dropped when the file is read.


## One replication

Replication k sets `set.seed(123 + k)` with L'Ecuyer-CMRG, as in Tables B1-B6.

1. 5,000 of the 10,000 individuals are drawn without replacement.
2. For each trait, the effects of its two causal variants are drawn from
   U(0.15, 0.4). Trait 1 has C1 and C2, and Trait 2 has C1 and C3, so C1 has
   a separate effect on each trait.
3. The errors are bivariate normal with variances 0.20 and correlation 0.40,
   and the traits have no intercept.
4. Both traits are fine-mapped with `FLASHFMwithJAMd_rawX_selreg()` once for
   each of OLS, St, DSh, and GSR. The settings are those of the example in
   `flashfm-savvySh.R`, target odds 1, cumulative posterior probability 0.99,
   models of at most 3 variants (raised only if JAM fails, up to 10), and 2
   million iterations.

JAM takes its sampler seed from R's random number generator, so the generator
is reset with `set.seed(100000 + k)` before each estimator. Every run can then
be repeated on its own, and the four estimators use the same sampler seed
within a replication.

`rep_k.csv` has one row per estimator, trait, and variant, with the
single-trait and multi-trait marginal posterior probabilities (MPP), whether
the variant is causal for the trait, its effect, and the computing time of
the estimator's fine-mapping in seconds.


## Power and FDR

A variant is selected when its MPP is above 0.9. In each replication, power
is the share of the trait's two causal variants that are selected, and FDR is
the share of selected variants that are not causal, or 0 when none is
selected. Both are averaged over the 100 replications, for the single-trait
and multi-trait MPPs separately. In the table, the highest power for a trait
among the single-trait rows, and among the multi-trait rows, is in red when
it is above that of OLS.
