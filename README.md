# AdaFilter-Sidak: simulations and GSE30495 real-data analysis

This repository uses **one R script** to reproduce the simulations and the real-data method comparisons from the manuscript *Uniformly Improved Adaptive Filtering Procedures for Replicability Testing with FWER Control.*.

## File

- `AdaFilter_Sidak_all_analyses.R` — independent normal-mixture simulations, AR(1) dependence simulations, calibrated simulation, and the complete GSE30495 preprocessing + PC-testing workflow.

## Getting started

1. Download this repository and open `AdaFilter_Sidak_all_analyses.R` in RStudio.
2. Check the four `run_*` switches at the top of the script.
3. Run the script (or use `source("AdaFilter_Sidak_all_analyses.R")`).

**By default, only `run_real_data <- TRUE`**; the other three switches are FALSE. This allows checking the manuscript's GSE30495 rejection counts without first running the Monte Carlo simulations. Set all switches to TRUE to run all analyses.

To regenerate the real-data table without rerunning downloads later, keep the downloaded GEO files inside `GSE30495_independent/`; the script reuses available files where supported.

## Dependencies

- Simulation plots: install `ggplot2` and `ggpubr`, e.g. `install.packages(c("ggplot2", "ggpubr"))`.
- Real-data analysis: `GEOquery`, `Biobase`, `affy`, `mouse4302cdf`, and `limma` (Bioconductor). The real-data section attempts to install missing packages using `BiocManager`.
- A connection to GEO and sufficient space for raw CEL files are needed for the first real-data run.

## Simulation settings

- Independent and AR(1) mixture: `m = 500`, `n = 2,3,4`, `mu = 2,3,4,5`, 10,000 repetitions; AR(1) lag-one correlation `rho = -0.8, 0.8`.
- Calibrated simulation: `m = 100`, `n = 5`, `u = 3`, `c = 3.1401031557`, `v = 1e-5, 1, 5, 10, 15`, 10,000 repetitions.

The script keeps the original random seeds and simulation/test calculations. It writes the pre-existing simulation CSV, RDS, and figure filenames to the current working directory.

## Real-data workflow

The real-data section automatically downloads GEO GSE30495, isolates WT/KO samples for each of five treatments, performs **RMA normalisation within each comparison**, calculates two-sided p-values, assembles `p_mat`, and applies eight PC-testing procedures. Its outputs are saved under `GSE30495_independent/`, including `p_mat.rds`, `results_all_methods.csv`, and `results_adafilter_table.csv`.

At the end, the script compares its six AdaFilter rejection-count rows against the manuscript table and reports **MATCH** or the count differences. That check is diagnostic only; it does not change any results.

**Important:** the provided original GSE30495.R calculates **ordinary t-statistics**, not empirical Bayes moderated t-statistics. It uses `limma::lmFit()` without `limma::eBayes()`. This is preserved deliberately to avoid altering results. The manuscript wording should be reconciled with the executed analysis.

## Reproducibility note

The consolidation preserves the calculations of the previously supplied scripts, but the script could not be executed in the preparation environment because R was unavailable. Please run it locally and compare the results before uploading the repository or submitting the paper. Different R/Bioconductor versions may also affect some real-data results.
