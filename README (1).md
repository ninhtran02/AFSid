# AdaFilter-Sidak: simulation and real-data method-comparison code

R scripts accompanying the manuscript *A uniformly improved adaptive filtering procedure for family-wise error rate control*.

## Files

| File | Purpose |
| --- | --- |
| `simulation_1_mixture_and_ar1.R` | Independent two-point normal mixture simulation (main text) and AR(1) dependence simulations (appendix). Produces FWER, power, and power-improvement figures. |
| `simulation_2_calibrated.R` | Calibrated normal simulation with five variance settings; produces the FWER/power figure. |
| `real_data_method_comparison.R` | Function for applying the eight PC testing procedures to an existing matrix of base p-values (one probe set per row, one study per column). |

## Requirements

R and the R packages `ggplot2` and `ggpubr` are needed for the simulation figures. The procedures and the real-data comparison function use base R / the `stats` package.

```r
install.packages(c("ggplot2", "ggpubr"))
```

## Running the simulations

From the directory containing the scripts, execute:

```sh
Rscript simulation_1_mixture_and_ar1.R
Rscript simulation_2_calibrated.R
```

Both scripts use the original `10,000` Monte Carlo repetitions and their existing random seeds. They can be computationally intensive. In `simulation_1_mixture_and_ar1.R`, set `run_independent` or `run_ar1` to `FALSE` near the top to run only the desired part. Generated CSV, RDS, PDF, and PNG files are saved in the current working directory with **the original output filenames**, including those referenced by the manuscript.

For a quick test only, temporarily reduce `n_rep` in the mixture/AR(1) sections or `B` in the calibrated script; restore `10,000` before reproducing manuscript results.

### Produced figures

Independent mixture:
- `adafilter_fwer_power_combined.pdf`
- `adafilter_power_improvement.pdf`

AR(1) dependence:
- `adafilter_ar1_rho_neg08_fwer_power.pdf`
- `adafilter_ar1_rho_pos08_fwer_power.pdf`
- `adafilter_ar1_rho_neg08_power_improvement.pdf`
- `adafilter_ar1_rho_pos08_power_improvement.pdf`

Calibrated simulation:
- `adafilter_calibrated_fwer_power.pdf` (also `.png`)

The scripts also write additional CSV/RDS summaries and, for the independent mixture, separate FWER/power PDFs.

## Using the real-data method-comparison function

The third file expects an existing **numeric** matrix `P_mat` of base p-values in `[0, 1]`. For example, after constructing the matrix using independent study-specific analyses:

```r
source("real_data_method_comparison.R")
result <- compare_pc_methods(P_mat, alpha = 0.05, u_values = 2:5)
print(result, row.names = FALSE)
```

Set `return_decisions = TRUE` to retrieve per-feature rejection indicators. The default CoFilter threshold is `tau_cofilter = 0.1`.

**Important reproducibility limitation:** The uploaded files did not include the raw GSE30495 download, RMA preprocessing, or limma analysis that creates `P_mat`. Consequently, `real_data_method_comparison.R` reproduces the comparison **once the p-value matrix has been prepared**, but these three files alone do not reproduce the manuscript's real-data table starting from GEO raw data. The full preprocessing script and/or the derived p-value matrix would be needed for that.

Dataset: [GSE30495](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30495).

## Notes on preservation

The simulation designs, random seeds, counts, test decisions, reported performance measures, and existing figure/output filenames have not been intentionally modified. The independent and AR(1) simulations are combined as two self-contained blocks, rather than rewriting their numerical implementations into a shared algorithm. The real-data comparison script corrects the display label `AdaFilter-Seq` to `AdaFilter-SeqBon`; this **does not change its rejection calculations**.

The code was reviewed statically in this preparation environment; **R was not available to rerun the 10,000-repetition simulations**, so numerical/graphical equality could not be independently verified by execution. For final archival reproducibility, run the scripts locally and compare their outputs with your current results.
