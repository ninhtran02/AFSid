# =============================================================================
# AdaFilter-Sidak: simulations and real-data analysis (single-file version)
# =============================================================================
# R code accompanying the manuscript:
#   A uniformly improved adaptive filtering procedure for family-wise
#   error rate control
#
# Included analyses:
#   1. Independent two-point normal mixture (main-text Simulation setting 1)
#   2. AR(1) mixture with rho = -0.8 and 0.8 (appendix)
#   3. Calibrated normal model (main-text Simulation setting 2)
#   4. GEO GSE30495: raw CEL download, treatment-specific RMA,
#      five separate WT/KO comparisons, and PC testing (real data)
#
# Running the entire script executes only the sections marked TRUE below.
# For a quick check against the manuscript's GSE30495 table, the default
# enables only real data. Set all four values to TRUE to regenerate everything.
# All original statistical calculations, simulation seeds, and output names
# have been retained; the analysis blocks use local scopes to avoid conflicts.
#
# Requirements:
# - Simulations/plots: ggplot2, ggpubr
# - Real data: Bioconductor GEOquery, Biobase, affy, mouse4302cdf, limma
# - Internet connection needed for GEO download and package installation
# - Run from a directory where output files can be written.
#
# IMPORTANT: The uploaded GSE30495 analysis calculates ordinary t-statistics
# from limma::lmFit, NOT empirical-Bayes moderated t-statistics. This was
# retained intentionally to avoid changing the supplied results.
# =============================================================================

# ---- Select which analyses to execute ----------------------------------------
run_independent <- FALSE    # 10,000 Monte Carlo repetitions per scenario
run_ar1         <- FALSE    # rho = -0.8 and 0.8, 10,000 reps each scenario
run_calibrated  <- FALSE    # 10,000 repetitions per variance value
run_real_data   <- TRUE     # Download GSE30495 and reproduce the table

# To run all analyses, set all four switches above to TRUE.
# To run simulations without downloading GEO data, set run_real_data <- FALSE.

# =============================================================================
# 1. INDEPENDENT TWO-POINT NORMAL MIXTURE
# =============================================================================

if (run_independent) local({
  # ============================================================================
  # Simulation setting 1: independent two-point normal mixture
  #
  # Six (u,n) configurations, eight FWER-controlling procedures.
  # No calibration.
  #
  # Outputs:
  #   1. adafilter_point_mixture_results.csv
  #   2. adafilter_point_mixture_results.rds
  #   3. adafilter_point_mixture_fwer.pdf
  #   4. adafilter_point_mixture_power.pdf
  #   5. adafilter_power_improvement.pdf
  # ============================================================================

  # ---- 1. Settings ----------------------------------------------------------

  m <- 500L
  alpha <- 0.05
  pi1 <- 0.50
  mu_grid <- c(2, 3, 4, 5)

  u_by_n <- list(
    `2` = 2L,
    `3` = 2:3,
    `4` = 2:4
  )

  n_rep <- 10000L
  tau <- 0.10
  base_seed <- 20261020L
  n_cores <- 1L

  method_names <- c(
    "AdaFilter-Bon",
    "AdaFilter-SeqBon",
    "AdaFilter-Bon (Sidak PC)",
    "AdaFilter-SeqBon (Sidak PC)",
    "AdaFilter-Sid",
    "AdaFilter-SeqSid",
    "CoFilter-Holm",
    "Hommel-Fisher"
  )

  config_order <- c(
    "(2,2)", "(2,3)", "(3,3)",
    "(2,4)", "(3,4)", "(4,4)"
  )

  # ---- 2. Fast row-wise sorting (n <= 4) ------------------------------------

  row_sort <- function(x) {
    n <- ncol(x)

    for (pass in seq_len(n - 1L)) {
      for (j in seq_len(n - pass)) {
        a <- x[, j]
        b <- x[, j + 1L]

        x[, j] <- pmin(a, b)
        x[, j + 1L] <- pmax(a, b)
      }
    }

    x
  }

  # ---- 3. Exact AdaFilter cutoff --------------------------------------------

  # N_R(t) = sum_{i notin R} 1{F_i < t}.
  #
  # Bonferroni:
  #   t * N_R(t) <= alpha
  #
  # Sidak:
  #   1 - (1-t)^N_R(t) <= alpha

  af_cutoff <- function(F, rejected, alpha,
                        rule = c("bon", "sidak")) {

    rule <- match.arg(rule)

    f <- sort(F[!rejected])
    L <- length(f)

    k <- 0:L
    lower <- c(0, f)
    upper <- c(f, Inf)

    t_k <- rep(alpha, L + 1L)
    idx <- which(k >= 2L)

    if (length(idx)) {
      if (rule == "bon") {
        t_k[idx] <- alpha / k[idx]
      } else {
        t_k[idx] <- -expm1(log1p(-alpha) / k[idx])
      }
    }

    candidate <- pmin(alpha, upper, t_k)

    max(candidate[candidate >= lower])
  }

  adafilter <- function(S, F, alpha,
                        rule = c("bon", "sidak"),
                        sequential = FALSE) {

    rule <- match.arg(rule)

    rejected <- rep(FALSE, length(S))

    repeat {
      cutoff <- af_cutoff(F, rejected, alpha, rule)

      newly <- !rejected & S < cutoff

      if (!any(newly)) break

      rejected[newly] <- TRUE

      if (!sequential) break
    }

    rejected
  }

  # ---- 4. Partial conjunction statistics ------------------------------------

  pc_stats <- function(P_sorted, u) {

    n <- ncol(P_sorted)
    h <- n - u + 1L

    p_u <- P_sorted[, u]
    p_before <- P_sorted[, u - 1L]

    list(
      S_bon = h * p_u,
      F_bon = h * p_before,

      S_sid = -expm1(h * log1p(-p_u)),
      F_sid = -expm1(h * log1p(-p_before)),

      fisher = pchisq(
        -2 * rowSums(log(P_sorted[, u:n, drop = FALSE])),
        df = 2L * h,
        lower.tail = FALSE
      )
    )
  }

  cofilter_holm <- function(fisher, alpha, tau = 0.1) {

    keep <- which(fisher <= tau)

    rejected <- rep(FALSE, length(fisher))

    if (length(keep)) {
      rejected[keep] <- p.adjust(
        fisher[keep] / tau,
        method = "holm"
      ) <= alpha
    }

    rejected
  }

  compare_methods <- function(P_sorted, u, alpha, tau) {

    ps <- pc_stats(P_sorted, u)

    out <- list(

      "AdaFilter-Bon" = adafilter(
        ps$S_bon, ps$F_bon, alpha,
        rule = "bon"
      ),

      "AdaFilter-SeqBon" = adafilter(
        ps$S_bon, ps$F_bon, alpha,
        rule = "bon",
        sequential = TRUE
      ),

      "AdaFilter-Bon (Sidak PC)" = adafilter(
        ps$S_sid, ps$F_sid, alpha,
        rule = "bon"
      ),

      "AdaFilter-SeqBon (Sidak PC)" = adafilter(
        ps$S_sid, ps$F_sid, alpha,
        rule = "bon",
        sequential = TRUE
      ),

      "AdaFilter-Sid" = adafilter(
        ps$S_sid, ps$F_sid, alpha,
        rule = "sidak"
      ),

      "AdaFilter-SeqSid" = adafilter(
        ps$S_sid, ps$F_sid, alpha,
        rule = "sidak",
        sequential = TRUE
      ),

      "CoFilter-Holm" = cofilter_holm(
        ps$fisher, alpha, tau
      ),

      "Hommel-Fisher" = (
        p.adjust(ps$fisher, method = "hommel") <= alpha
      )
    )

    out[method_names]
  }

  # ---- 5. Simulate one (n, mu) condition ------------------------------------

  run_condition <- function(n, mu) {

    n <- as.integer(n)

    set.seed(
      base_seed + 1000L * n + as.integer(10 * mu)
    )

    u_grid <- u_by_n[[as.character(n)]]

    fwer_count <- matrix(
      0,
      length(u_grid),
      length(method_names)
    )

    power_sum <- matrix(
      0,
      length(u_grid),
      length(method_names)
    )

    for (b in seq_len(n_rep)) {

      # Theta_ij = mu with probability pi1,
      # and Theta_ij = 0 otherwise.

      signal <- matrix(
        rbinom(m * n, size = 1L, prob = pi1),
        m, n
      )

      Z <- matrix(rnorm(m * n), m, n) + mu * signal

      P <- 2 * pnorm(abs(Z), lower.tail = FALSE)

      # Avoid log(0) in Fisher combination.
      P <- pmax(P, .Machine$double.xmin)

      P_sorted <- row_sort(P)

      K <- rowSums(signal)

      # Reuse the same data for all u at this n.
      for (a in seq_along(u_grid)) {

        u <- u_grid[a]

        true_pc <- K < u
        n_alt <- sum(!true_pc)

        decisions <- compare_methods(
          P_sorted, u, alpha, tau
        )

        # Check theoretical rejection-set inclusions.
        if (b == 1L) {
          stopifnot(
            !any(decisions[[1]] & !decisions[[3]]),
            !any(decisions[[3]] & !decisions[[5]]),
            !any(decisions[[2]] & !decisions[[4]]),
            !any(decisions[[4]] & !decisions[[6]])
          )
        }

        for (j in seq_along(method_names)) {

          R <- decisions[[j]]

          fwer_count[a, j] <- fwer_count[a, j] +
            any(R & true_pc)

          # Fraction of false PC nulls rejected.
          # Defined as zero if no false PC null exists.
          if (n_alt > 0L) {
            power_sum[a, j] <- power_sum[a, j] +
              sum(R & !true_pc) / n_alt
          }
        }
      }
    }

    do.call(
      rbind,
      lapply(seq_along(u_grid), function(a) {

        data.frame(
          u = u_grid[a],
          n = n,
          mu = mu,

          configuration = paste0(
            "(", u_grid[a], ",", n, ")"
          ),

          Method = method_names,

          FWER = fwer_count[a, ] / n_rep,
          Power = power_sum[a, ] / n_rep,

          row.names = NULL
        )
      })
    )
  }

  # ---- 6. Run the 12 (n, mu) scenarios -------------------------------------

  scenarios <- expand.grid(
    n = 2:4,
    mu = mu_grid
  )

  jobs <- split(
    scenarios,
    seq_len(nrow(scenarios))
  )

  run_job <- function(job) {

    n <- job$n[1]
    mu <- job$mu[1]

    ans <- run_condition(n, mu)

    message(sprintf(
      "Finished n=%d, mu=%g", n, mu
    ))

    ans
  }

  if (n_cores > 1L && .Platform$OS.type != "windows") {

    pieces <- parallel::mclapply(
      jobs, run_job,
      mc.cores = n_cores
    )

  } else {

    pieces <- lapply(jobs, run_job)

  }

  results <- do.call(rbind, pieces)

  results$Method <- factor(
    results$Method,
    levels = method_names
  )

  results$configuration <- factor(
    results$configuration,
    levels = config_order
  )

  results <- results[
    order(
      results$configuration,
      results$mu,
      results$Method
    ),
  ]

  row.names(results) <- NULL

  # ---- 7. Save simulation results ------------------------------------------

  write.csv(
    results,
    "adafilter_point_mixture_results.csv",
    row.names = FALSE
  )

  saveRDS(
    results,
    "adafilter_point_mixture_results.rds"
  )

  print(results, row.names = FALSE)

  # ---- 8. Publication-style figures ----------------------------------------

  if (requireNamespace("ggplot2", quietly = TRUE)) {

    library(ggplot2)

    # ---- 8.1. Consistent colours and linetypes ------------------------------

    # These named palettes are used in ALL three figures.

    method_colors <- c(
      "AdaFilter-Bon"               = "#999999",
      "AdaFilter-SeqBon"            = "#E69F00",
      "AdaFilter-Bon (Sidak PC)"    = "#009E73",
      "AdaFilter-SeqBon (Sidak PC)" = "#56B4E9",
      "AdaFilter-Sid"               = "#D55E00",
      "AdaFilter-SeqSid"            = "#0072B2",
      "CoFilter-Holm"               = "#CC79A7",
      "Hommel-Fisher"               = "#000000"
    )

    method_linetypes <- c(
      "AdaFilter-Bon"               = "solid",
      "AdaFilter-SeqBon"            = "dashed",
      "AdaFilter-Bon (Sidak PC)"    = "dotted",
      "AdaFilter-SeqBon (Sidak PC)" = "dotdash",
      "AdaFilter-Sid"               = "longdash",
      "AdaFilter-SeqSid"            = "twodash",
      "CoFilter-Holm"               = "solid",
      "Hommel-Fisher"               = "dashed"
    )

    # Reusable manual scales.
    method_scales <- list(
      scale_color_manual(values = method_colors),
      scale_linetype_manual(values = method_linetypes)
    )

    theme_sim <- theme_bw(base_size = 12) +
      theme(
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank()
      )

    # ---- 8.2. Empirical FWER ------------------------------------------------

    g_fwer <- ggplot(
      results,
      aes(
        x = mu,
        y = FWER,
        color = Method,
        linetype = Method
      )
    ) +
      geom_hline(
        yintercept = alpha,
        color = "grey40",
        linetype = "dashed"
      ) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      facet_wrap(
        ~ configuration,
        nrow = 2
      ) +
      scale_x_continuous(
        breaks = mu_grid
      ) +
      method_scales +
      labs(
        x = expression(mu),
        y = "Empirical FWER"
      ) +
      theme_sim

    ggsave(
      "adafilter_point_mixture_fwer.pdf",
      g_fwer,
      width = 11,
      height = 7
    )

    # ---- 8.3. Empirical power -----------------------------------------------

    g_power <- ggplot(
      results,
      aes(
        x = mu,
        y = Power,
        color = Method,
        linetype = Method
      )
    ) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      facet_wrap(
        ~ configuration,
        nrow = 2
      ) +
      scale_x_continuous(
        breaks = mu_grid
      ) +
      scale_y_continuous(
        limits = c(0, 1)
      ) +
      method_scales +
      labs(
        x = expression(mu),
        y = "Empirical power"
      ) +
      theme_sim

    ggsave(
      "adafilter_point_mixture_power.pdf",
      g_power,
      width = 11,
      height = 7
    )

    # ---- 8.4. Power improvement over AdaFilter variants --------------------

    # Reference method: AdaFilter-SeqSid.
    # Compare against the other five AdaFilter procedures.

    competitors <- c(
      "AdaFilter-Bon",
      "AdaFilter-SeqBon",
      "AdaFilter-Bon (Sidak PC)",
      "AdaFilter-SeqBon (Sidak PC)",
      "AdaFilter-Sid"
    )

    reference <- subset(
      results,
      Method == "AdaFilter-SeqSid",
      select = c(u, n, mu, Power)
    )

    names(reference)[4] <- "ReferencePower"

    df <- subset(
      results,
      Method %in% competitors
    )

    df <- merge(
      df,
      reference,
      by = c("u", "n", "mu")
    )

    # Absolute improvement expressed in percentage points.
    #
    # PowerGain = 100 * (Power_SeqSid - Power_Competitor)

    df$PowerGain <- 100 * (
      df$ReferencePower - df$Power
    )

    df$configuration <- factor(
      paste0("(", df$u, ",", df$n, ")"),
      levels = config_order
    )

    df$Method <- factor(
      df$Method,
      levels = competitors
    )

    # Numerical tolerance for simulation estimates.
    if (any(df$PowerGain < -1e-10)) {
      warning(
        "Negative power gains detected. Check the results."
      )
    }

    # ---- 8.5. Power improvement plot ----------------------------------------

    g_gain <- ggplot(
      df,
      aes(
        x = mu,
        y = PowerGain,
        color = Method,
        linetype = Method
      )
    ) +
      geom_hline(
        yintercept = 0,
        color = "grey50",
        linetype = "dashed"
      ) +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      facet_wrap(
        ~ configuration,
        nrow = 2,
        scales = "free_y"
      ) +
      scale_x_continuous(
        breaks = mu_grid
      ) +
      method_scales +
      labs(
        x = expression(mu),
        y = "Power improvement (percentage points)"
      ) +
      theme_sim

    ggsave(
      "adafilter_power_improvement.pdf",
      g_gain,
      width = 11,
      height = 7
    )

    # ---- 8.6. Display all figures -------------------------------------------

    print(g_fwer)
    print(g_power)
    print(g_gain)

  } else {

    message(
      "Install ggplot2 to automatically generate the three PDF figures."
    )

  }

  library(ggpubr)

  g_combined <- ggarrange(
    g_fwer,
    g_power,
    ncol = 1,
    nrow = 2,
    labels = c("(a)", "(b)"),
    common.legend = TRUE,
    legend = "bottom",
    align = "v"
  )

  print(g_combined)

  ggsave(
    "adafilter_fwer_power_combined.pdf",
    g_combined,
    width = 11,
    height = 12
  )

  library(ggpubr)
})

if (run_ar1) local({
  # ============================================================================
  # Simulation setting 1: two-point normal mixture with AR(1) dependence
  # Across features within each study: Sigma[i,i'] = rho^abs(i-i')
  # Across studies: independent
  # Two separate runs: rho = -0.8 and rho = +0.8
  # Each run produces a combined FWER/power figure with one shared legend.
  # ============================================================================

  # ---- 1. Settings ----------------------------------------------------------
  m <- 500L
  alpha <- 0.05
  pi1 <- 0.50
  mu_grid <- c(2, 3, 4, 5)
  u_by_n <- list(`2` = 2L, `3` = 2:3, `4` = 2:4)
  n_rep <- 10000L                 # Change to 100L for a smoke test
  tau <- 0.10
  base_seed <- 20261020L
  n_cores <- 1L                   # > 1 on macOS/Linux if desired
  rho_values <- c(-0.8, 0.8)

  method_names <- c(
    "AdaFilter-Bon",
    "AdaFilter-SeqBon",
    "AdaFilter-Bon (Sidak PC)",
    "AdaFilter-SeqBon (Sidak PC)",
    "AdaFilter-Sid",
    "AdaFilter-SeqSid",
    "CoFilter-Holm",
    "Hommel-Fisher"
  )

  config_order <- c("(2,2)", "(2,3)", "(3,3)",
                    "(2,4)", "(3,4)", "(4,4)")

  # Same method colours and line types as independent Simulation setting 1.
  method_colors <- c(
    "AdaFilter-Bon"                = "#999999",
    "AdaFilter-SeqBon"             = "#E69F00",
    "AdaFilter-Bon (Sidak PC)"     = "#009E73",
    "AdaFilter-SeqBon (Sidak PC)"  = "#56B4E9",
    "AdaFilter-Sid"                = "#D55E00",
    "AdaFilter-SeqSid"             = "#0072B2",
    "CoFilter-Holm"                = "#CC79A7",
    "Hommel-Fisher"                = "#000000"
  )

  method_linetypes <- c(
    "AdaFilter-Bon"                = "solid",
    "AdaFilter-SeqBon"             = "dashed",
    "AdaFilter-Bon (Sidak PC)"     = "dotted",
    "AdaFilter-SeqBon (Sidak PC)"  = "dotdash",
    "AdaFilter-Sid"                = "longdash",
    "AdaFilter-SeqSid"             = "twodash",
    "CoFilter-Holm"                = "solid",
    "Hommel-Fisher"                = "dashed"
  )

  # ---- 2. Row-wise sorting --------------------------------------------------
  row_sort <- function(x) {
    nc <- ncol(x)
    for (pass in seq_len(nc - 1L)) {
      for (j in seq_len(nc - pass)) {
        a <- x[, j]
        b <- x[, j + 1L]
        x[, j] <- pmin(a, b)
        x[, j + 1L] <- pmax(a, b)
      }
    }
    x
  }

  # ---- 3. Stationary AR(1) Gaussian errors across FEATURES -----------------
  # For each study j:
  #   e[1,j] ~ N(0,1)
  #   e[i,j] = rho * e[i-1,j] + sqrt(1-rho^2) * innovation[i,j]
  # The study columns are mutually independent.
  ar1_errors <- function(m, n, rho) {
    stopifnot(abs(rho) < 1)
    errors <- matrix(rnorm(m * n), nrow = m, ncol = n)
    if (rho == 0 || m == 1L) return(errors)
    sd_innov <- sqrt(1 - rho^2)
    for (j in seq_len(n)) {
      errors[-1L, j] <- as.numeric(stats::filter(
        sd_innov * errors[-1L, j],
        filter = rho,
        method = "recursive",
        init = errors[1L, j]
      ))
    }
    errors
  }

  # ---- 4. AdaFilter cutoff and rejection rules -----------------------------
  af_cutoff <- function(F, rejected, alpha, rule = c("bon", "sidak")) {
    rule <- match.arg(rule)
    f <- sort(F[!rejected])
    L <- length(f)
    k <- 0:L
    lower <- c(0, f)
    upper <- c(f, Inf)

    t_k <- rep(alpha, L + 1L)
    idx <- which(k >= 2L)
    if (length(idx)) {
      if (rule == "bon") {
        t_k[idx] <- alpha / k[idx]
      } else {
        t_k[idx] <- -expm1(log1p(-alpha) / k[idx])
      }
    }
    candidate <- pmin(alpha, upper, t_k)
    max(candidate[candidate >= lower])
  }

  adafilter <- function(S, F, alpha, rule = c("bon", "sidak"),
                        sequential = FALSE) {
    rule <- match.arg(rule)
    rejected <- rep(FALSE, length(S))
    repeat {
      cutoff <- af_cutoff(F, rejected, alpha, rule)
      newly <- !rejected & S < cutoff
      if (!any(newly)) break
      rejected[newly] <- TRUE
      if (!sequential) break
    }
    rejected
  }

  # ---- 5. Partial conjunction p-values and methods -------------------------
  pc_stats <- function(P_sorted, u) {
    n <- ncol(P_sorted)
    h <- n - u + 1L
    p_u <- P_sorted[, u]
    p_before <- P_sorted[, u - 1L]

    list(
      S_bon = h * p_u,
      F_bon = h * p_before,
      S_sid = -expm1(h * log1p(-p_u)),
      F_sid = -expm1(h * log1p(-p_before)),
      fisher = pchisq(
        -2 * rowSums(log(P_sorted[, u:n, drop = FALSE])),
        df = 2L * h, lower.tail = FALSE
      )
    )
  }

  cofilter_holm <- function(fisher, alpha, tau = 0.1) {
    keep <- which(fisher <= tau)
    rejected <- rep(FALSE, length(fisher))
    if (length(keep)) {
      rejected[keep] <- p.adjust(fisher[keep] / tau, method = "holm") <= alpha
    }
    rejected
  }

  compare_methods <- function(P_sorted, u, alpha, tau) {
    ps <- pc_stats(P_sorted, u)
    out <- list(
      "AdaFilter-Bon" = adafilter(ps$S_bon, ps$F_bon, alpha,
                                  rule = "bon"),
      "AdaFilter-SeqBon" = adafilter(ps$S_bon, ps$F_bon, alpha,
                                     rule = "bon", sequential = TRUE),
      "AdaFilter-Bon (Sidak PC)" = adafilter(ps$S_sid, ps$F_sid, alpha,
                                             rule = "bon"),
      "AdaFilter-SeqBon (Sidak PC)" = adafilter(ps$S_sid, ps$F_sid, alpha,
                                                rule = "bon", sequential = TRUE),
      "AdaFilter-Sid" = adafilter(ps$S_sid, ps$F_sid, alpha,
                                  rule = "sidak"),
      "AdaFilter-SeqSid" = adafilter(ps$S_sid, ps$F_sid, alpha,
                                     rule = "sidak", sequential = TRUE),
      "CoFilter-Holm" = cofilter_holm(ps$fisher, alpha, tau),
      "Hommel-Fisher" = p.adjust(ps$fisher, method = "hommel") <= alpha
    )
    out[method_names]
  }

  # ---- 6. One (rho, n, mu) condition ----------------------------------------
  run_condition <- function(n, mu, rho) {
    n <- as.integer(n)
    # Same underlying RNG seed for +/- rho at fixed (n,mu), to ease comparison.
    set.seed(base_seed + 1000L * n + as.integer(10 * mu))
    u_grid <- u_by_n[[as.character(n)]]
    fwer_count <- matrix(0, length(u_grid), length(method_names))
    power_sum <- matrix(0, length(u_grid), length(method_names))

    for (b in seq_len(n_rep)) {
      # Independently assigned null and nonnull effects.
      signal <- matrix(rbinom(m * n, 1L, pi1), nrow = m, ncol = n)
      # Dependence is across rows (features), separately within each study.
      Z <- mu * signal + ar1_errors(m, n, rho)
      P <- pmax(2 * pnorm(abs(Z), lower.tail = FALSE),
                .Machine$double.xmin)
      P_sorted <- row_sort(P)
      K <- rowSums(signal)

      for (a in seq_along(u_grid)) {
        u <- u_grid[a]
        true_pc <- K < u
        n_alt <- sum(!true_pc)
        decisions <- compare_methods(P_sorted, u, alpha, tau)

        if (b == 1L) {
          stopifnot(
            !any(decisions[[1]] & !decisions[[3]]),
            !any(decisions[[3]] & !decisions[[5]]),
            !any(decisions[[2]] & !decisions[[4]]),
            !any(decisions[[4]] & !decisions[[6]])
          )
        }

        for (j in seq_along(method_names)) {
          R <- decisions[[j]]
          fwer_count[a, j] <- fwer_count[a, j] + any(R & true_pc)
          if (n_alt > 0L) {
            power_sum[a, j] <- power_sum[a, j] + sum(R & !true_pc) / n_alt
          }
        }
      }
    }

    do.call(rbind, lapply(seq_along(u_grid), function(a) {
      data.frame(
        rho = rho, u = u_grid[a], n = n, mu = mu,
        configuration = paste0("(", u_grid[a], ",", n, ")"),
        Method = method_names,
        FWER = fwer_count[a, ] / n_rep,
        Power = power_sum[a, ] / n_rep,
        row.names = NULL
      )
    }))
  }

  # ---- 7. Run separately at rho = -0.8 and rho = +0.8 -----------------------
  scenarios <- expand.grid(n = 2:4, mu = mu_grid)
  jobs <- split(scenarios, seq_len(nrow(scenarios)))

  run_rho <- function(rho) {
    message(sprintf("Starting AR(1) simulation at rho = %+0.1f", rho))
    run_job <- function(job) {
      n <- job$n[1]
      mu <- job$mu[1]
      ans <- run_condition(n, mu, rho)
      message(sprintf("  Finished rho=%+.1f, n=%d, mu=%g", rho, n, mu))
      ans
    }
    if (n_cores > 1L && .Platform$OS.type != "windows") {
      pieces <- parallel::mclapply(jobs, run_job, mc.cores = n_cores)
    } else {
      pieces <- lapply(jobs, run_job)
    }
    do.call(rbind, pieces)
  }

  results <- do.call(rbind, lapply(rho_values, run_rho))
  results$Method <- factor(results$Method, levels = method_names)
  results$configuration <- factor(results$configuration, levels = config_order)
  results <- results[order(results$rho, results$configuration,
                           results$mu, results$Method), ]
  row.names(results) <- NULL

  write.csv(results, "adafilter_ar1_results.csv", row.names = FALSE)
  saveRDS(results, "adafilter_ar1_results.rds")
  print(results, row.names = FALSE)

  # ---- 8. Separate combined FWER + power plots for each rho ----------------
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggpubr", quietly = TRUE)) {
    message("Install ggplot2 and ggpubr to generate the combined PDF plots.")
  } else {
    library(ggplot2)
    library(ggpubr)

    theme_sim <- theme_bw(base_size = 12) +
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            strip.text = element_text(face = "bold"),
            panel.grid.minor = element_blank())

    method_scales <- list(
      scale_color_manual(values = method_colors),
      scale_linetype_manual(values = method_linetypes)
    )

    for (rho in rho_values) {
      sub <- results[abs(results$rho - rho) < 1e-8, ]
      max_fwer <- max(alpha * 1.1, max(sub$FWER) * 1.1)

      g_fwer <- ggplot(
        sub, aes(x = mu, y = FWER, color = Method, linetype = Method)
      ) +
        geom_hline(yintercept = alpha, color = "grey40", linetype = "dashed") +
        geom_line(linewidth = 0.8) + geom_point(size = 2) +
        facet_wrap(~ configuration, nrow = 2) +
        scale_x_continuous(breaks = mu_grid) +
        coord_cartesian(ylim = c(0, max_fwer)) +
        method_scales +
        labs(x = expression(mu), y = "Empirical FWER") + theme_sim

      g_power <- ggplot(
        sub, aes(x = mu, y = Power, color = Method, linetype = Method)
      ) +
        geom_line(linewidth = 0.8) + geom_point(size = 2) +
        facet_wrap(~ configuration, nrow = 2) +
        scale_x_continuous(breaks = mu_grid) +
        coord_cartesian(ylim = c(0, 1)) +
        method_scales +
        labs(x = expression(mu), y = "Empirical power") + theme_sim

      # Retains your original stacked layout and common legend.
      g_combined <- ggarrange(
        g_fwer, g_power,
        ncol = 1, nrow = 2,
        labels = c("(a)", "(b)"),
        common.legend = TRUE,
        legend = "bottom",
        align = "v"
      )

      rho_tag <- if (rho < 0) "neg08" else "pos08"
      filename <- paste0("adafilter_ar1_rho_", rho_tag, "_fwer_power.pdf")
      ggsave(filename, g_combined, width = 11, height = 12)
      print(g_combined)
      message("Saved: ", filename)
    }
  }

  # Power-improvement figures

  # Power-improvement figures for both AR(1) correlations.
  # Reads the previously saved results CSV; does not rerun Monte Carlo.
  # Requires ggplot2, but does not rerun the Monte Carlo simulation.

  library(ggplot2)
  results <- read.csv("adafilter_ar1_results.csv")

  competitors <- c(
    "AdaFilter-Bon", "AdaFilter-SeqBon",
    "AdaFilter-Bon (Sidak PC)", "AdaFilter-SeqBon (Sidak PC)",
    "AdaFilter-Sid"
  )
  config_order <- c("(2,2)", "(2,3)", "(3,3)",
                    "(2,4)", "(3,4)", "(4,4)")
  mu_grid <- c(2, 3, 4, 5)

  # Exact same colours and linetypes as independent Simulation Setting 1.
  method_colors <- c(
    "AdaFilter-Bon" = "#999999",
    "AdaFilter-SeqBon" = "#E69F00",
    "AdaFilter-Bon (Sidak PC)" = "#009E73",
    "AdaFilter-SeqBon (Sidak PC)" = "#56B4E9",
    "AdaFilter-Sid" = "#D55E00",
    "AdaFilter-SeqSid" = "#0072B2",
    "CoFilter-Holm" = "#CC79A7",
    "Hommel-Fisher" = "#000000"
  )
  method_linetypes <- c(
    "AdaFilter-Bon" = "solid",
    "AdaFilter-SeqBon" = "dashed",
    "AdaFilter-Bon (Sidak PC)" = "dotted",
    "AdaFilter-SeqBon (Sidak PC)" = "dotdash",
    "AdaFilter-Sid" = "longdash",
    "AdaFilter-SeqSid" = "twodash",
    "CoFilter-Holm" = "solid",
    "Hommel-Fisher" = "dashed"
  )

  theme_sim <- theme_bw(base_size = 12) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank())

  for (rho in c(-0.8, 0.8)) {
    sub <- subset(results, abs(results$rho - rho) < 1e-8)
    stopifnot(nrow(sub) > 0)
    reference <- subset(sub, Method == "AdaFilter-SeqSid",
                        select = c(u, n, mu, Power))
    names(reference)[4] <- "ReferencePower"
    gains <- merge(sub[sub$Method %in% competitors, ], reference,
                   by = c("u", "n", "mu"))
    gains$PowerGain <- 100 * (gains$ReferencePower - gains$Power)
    gains$Method <- factor(gains$Method, levels = competitors)
    gains$configuration <- factor(
      paste0("(", gains$u, ",", gains$n, ")"), levels = config_order
    )

    p <- ggplot(gains, aes(mu, PowerGain, color = Method, linetype = Method)) +
      geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
      geom_line(linewidth = 0.8) +
      geom_point(size = 2) +
      facet_wrap(~ configuration, nrow = 2, scales = "free_y") +
      scale_x_continuous(breaks = mu_grid) +
      scale_color_manual(values = method_colors) +
      scale_linetype_manual(values = method_linetypes) +
      labs(x = expression(mu),
           y = "Power improvement (percentage points)") +
      theme_sim

    rho_tag <- if (rho < 0) "neg08" else "pos08"
    filename <- paste0("adafilter_ar1_rho_", rho_tag,
                       "_power_improvement.pdf")
    ggsave(filename, p, width = 11, height = 7)
    print(p)
    message("Saved ", filename)
  }
})

# =============================================================================
# 3. CALIBRATED NORMAL SIMULATION
# =============================================================================

if (run_calibrated) local({
  # =============================================================================
  # Simulation setting 2: calibrated normal model
  # =============================================================================
  # Reproduces the calibrated simulation, with the original random seed,
  # rejection rules, and plot/output filenames unchanged.
  # Run with: Rscript simulation_2_calibrated.R
  # =============================================================================

  # ============================================================================
  # Simulation setting 2: calibrated normal model
  # m = 100, n = 5, u = 3, alpha = 0.05
  # Only the calibrated mean c = 3.1401031557 is considered.
  # Eight FWER-controlling procedures; 10,000 Monte Carlo repetitions.
  #
  # Outputs:
  #   adafilter_calibrated_summary.csv
  #   adafilter_calibrated_replicates.rds
  #   adafilter_calibrated_results.rds
  #   adafilter_calibrated_fwer_power.pdf
  #   adafilter_calibrated_fwer_power.png
  # ============================================================================

  # ---- 1. Settings ----------------------------------------------------------

  set.seed(20261011)

  m <- 100L
  n <- 5L
  u <- 3L
  alpha <- 0.05
  B <- 10000L                  # Use 100L for a smoke test
  c <- 3.1401031557
  v_values <- c(1e-5, 1, 5, 10, 15)
  tau <- 0.10                  # Prespecified CoFilter threshold
  k <- n - u + 1L

  methods <- c(
    "AdaFilter-Bon",
    "AdaFilter-SeqBon",
    "AdaFilter-Bon (Sidak PC)",
    "AdaFilter-SeqBon (Sidak PC)",
    "AdaFilter-Sid",
    "AdaFilter-SeqSid",
    "CoFilter-Holm",
    "Hommel-Fisher"
  )

  # Same colours and linetypes as Simulation setting 1.
  method_colors <- c(
    "AdaFilter-Bon"               = "#999999",
    "AdaFilter-SeqBon"            = "#E69F00",
    "AdaFilter-Bon (Sidak PC)"    = "#009E73",
    "AdaFilter-SeqBon (Sidak PC)" = "#56B4E9",
    "AdaFilter-Sid"               = "#D55E00",
    "AdaFilter-SeqSid"            = "#0072B2",
    "CoFilter-Holm"               = "#CC79A7",
    "Hommel-Fisher"               = "#000000"
  )

  method_linetypes <- c(
    "AdaFilter-Bon"               = "solid",
    "AdaFilter-SeqBon"            = "dashed",
    "AdaFilter-Bon (Sidak PC)"    = "dotted",
    "AdaFilter-SeqBon (Sidak PC)" = "dotdash",
    "AdaFilter-Sid"               = "longdash",
    "AdaFilter-SeqSid"            = "twodash",
    "CoFilter-Holm"               = "solid",
    "Hommel-Fisher"               = "dashed"
  )

  stopifnot(identical(names(method_colors), methods),
            identical(names(method_linetypes), methods))

  # ---- 2. AdaFilter rejection rule ------------------------------------------

  # For continuously distributed base p-values, the following is equivalent
  # to the supremum AdaFilter cutoff. For each currently unrejected feature i,
  # M_i = #{j not yet rejected : F_j <= S_i}.
  # Bonferroni rule: S_i < alpha / M_i.
  # Sidak rule:      S_i < 1 - (1 - alpha)^(1 / M_i).
  # Sequential methods recalculate M_i after each rejection step.

  ada_reject <- function(S, F, alpha,
                         multiplicity = c("bon", "sidak"),
                         sequential = FALSE) {
    multiplicity <- match.arg(multiplicity)
    stopifnot(length(S) == length(F), all(F <= S + 1e-14))

    rejected <- rep(FALSE, length(S))

    repeat {
      remaining <- which(!rejected)
      if (!length(remaining)) break

      counts <- findInterval(S[remaining], sort(F[remaining]))
      stopifnot(all(counts >= 1L))

      cutoff <- if (multiplicity == "bon") {
        alpha / counts
      } else {
        -expm1(log1p(-alpha) / counts)
      }

      newly <- remaining[S[remaining] < cutoff]
      if (!length(newly)) break

      rejected[newly] <- TRUE
      if (!sequential) break
    }

    rejected
  }

  # Optional check against the exact supremum definition.
  ada_cutoff_exact <- function(F, alpha,
                               multiplicity = c("bon", "sidak")) {
    multiplicity <- match.arg(multiplicity)
    ks <- seq_along(F)

    critical <- if (multiplicity == "bon") {
      alpha / ks
    } else {
      -expm1(log1p(-alpha) / ks)
    }

    candidates <- unique(c(alpha, F[F <= alpha], critical))
    candidates <- candidates[candidates >= 0 & candidates <= alpha]

    counts <- findInterval(candidates, sort(F), left.open = TRUE)
    limits <- if (multiplicity == "bon") {
      ifelse(counts == 0L, alpha, alpha / pmax(counts, 1L))
    } else {
      ifelse(counts == 0L, alpha,
             -expm1(log1p(-alpha) / pmax(counts, 1L)))
    }

    max(candidates[candidates <= limits * (1 + 1e-12)])
  }

  # ---- 3. All eight procedures ----------------------------------------------

  run_methods <- function(P, u = 3L, alpha = 0.05, tau = 0.1,
                          check = TRUE) {
    m_local <- nrow(P)
    n_local <- ncol(P)
    h <- n_local - u + 1L
    P_ordered <- t(apply(P, 1L, sort))

    S_bon <- pmin(1, h * P_ordered[, u])
    F_bon <- pmin(1, h * P_ordered[, u - 1L])

    S_sid <- -expm1(h * log1p(-P_ordered[, u]))
    F_sid <- -expm1(h * log1p(-P_ordered[, u - 1L]))

    fisher_stat <- -2 * rowSums(log(P_ordered[, u:n_local, drop = FALSE]))
    S_fisher <- pchisq(fisher_stat, df = 2L * h, lower.tail = FALSE)

    ans <- matrix(FALSE, nrow = m_local, ncol = length(methods),
                  dimnames = list(NULL, methods))

    ans[, 1] <- ada_reject(S_bon, F_bon, alpha, "bon")
    ans[, 2] <- ada_reject(S_bon, F_bon, alpha, "bon", sequential = TRUE)
    ans[, 3] <- ada_reject(S_sid, F_sid, alpha, "bon")
    ans[, 4] <- ada_reject(S_sid, F_sid, alpha, "bon", sequential = TRUE)
    ans[, 5] <- ada_reject(S_sid, F_sid, alpha, "sidak")
    ans[, 6] <- ada_reject(S_sid, F_sid, alpha, "sidak", sequential = TRUE)

    chosen <- which(S_fisher <= tau)
    if (length(chosen)) {
      ans[chosen, 7] <- p.adjust(S_fisher[chosen] / tau,
                                 method = "holm") <= alpha
    }
    ans[, 8] <- p.adjust(S_fisher, method = "hommel") <= alpha

    if (check) {
      nested_pairs <- list(
        c(1, 3), c(3, 5), c(2, 4), c(4, 6),
        c(1, 2), c(3, 4), c(5, 6)
      )
      for (pair in nested_pairs) {
        if (any(ans[, pair[1]] & !ans[, pair[2]])) {
          stop("Unexpected non-nested sets: ",
               methods[pair[1]], " vs. ", methods[pair[2]])
        }
      }
    }

    ans
  }

  check_cutoff_implementation <- function() {
    for (b in seq_len(30L)) {
      P <- matrix(runif(m * n), m, n)
      P_ordered <- t(apply(P, 1L, sort))

      S_b <- pmin(1, k * P_ordered[, u])
      F_b <- pmin(1, k * P_ordered[, u - 1L])
      S_s <- -expm1(k * log1p(-P_ordered[, u]))
      F_s <- -expm1(k * log1p(-P_ordered[, u - 1L]))

      for (setting in list(
        list(S_b, F_b, "bon"),
        list(S_s, F_s, "bon"),
        list(S_s, F_s, "sidak")
      )) {
        S <- setting[[1]]
        F <- setting[[2]]
        correction <- setting[[3]]
        cutoff <- ada_cutoff_exact(F, alpha, correction)
        stopifnot(identical(S < cutoff,
                            ada_reject(S, F, alpha, correction)))
      }
    }
    invisible(TRUE)
  }

  check_cutoff_implementation()

  # ---- 4. Simulate only the calibrated setting ------------------------------

  # Features 1:10 have exactly three non-null base hypotheses.
  # Features 11:100 have all null base hypotheses.
  # Only c = 3.1401031557 is used; v varies over five values.

  summary_results <- vector("list", length(v_values))
  replicate_results <- vector("list", length(v_values))

  for (idx in seq_along(v_values)) {
    v <- v_values[idx]

    theta <- matrix(0, m, n)
    SD <- matrix(1, m, n)

    theta[1:10, 1:2] <- 6
    SD[1:10, 1:2] <- sqrt(0.01)
    theta[1:10, 3] <- c
    SD[1:10, 3] <- sqrt(v)

    true_positive <- matrix(0L, B, length(methods),
                            dimnames = list(NULL, methods))
    false_positive <- matrix(0L, B, length(methods),
                             dimnames = list(NULL, methods))

    for (b in seq_len(B)) {
      Z <- theta + SD * matrix(rnorm(m * n), m, n)
      P <- 2 * pnorm(-abs(Z))  # Reference distribution remains N(0,1).
      P <- pmax(P, .Machine$double.xmin)  # Avoid log(0) in Fisher.

      rejected <- run_methods(P, u, alpha, tau, check = TRUE)

      true_positive[b, ] <- colSums(rejected[1:10, , drop = FALSE])
      false_positive[b, ] <- colSums(rejected[11:m, , drop = FALSE])
    }

    power_per_rep <- true_positive / 10
    fwer_per_rep <- false_positive > 0L
    fwer_hat <- colMeans(fwer_per_rep)

    summary_results[[idx]] <- data.frame(
      c_value = c,
      v = v,
      Method = methods,
      FWER = fwer_hat,
      FWER_SE = sqrt(fwer_hat * (1 - fwer_hat) / B),
      Power = colMeans(power_per_rep),
      Power_SE = apply(power_per_rep, 2, sd) / sqrt(B),
      Mean_TP = colMeans(true_positive),
      Mean_FP = colMeans(false_positive),
      row.names = NULL
    )

    replicate_results[[idx]] <- list(
      c_value = c,
      v = v,
      TP = true_positive,
      FP = false_positive
    )

    message(sprintf("Finished c = %.10f, v = %g (%d repetitions)",
                    c, v, B))
    print(summary_results[[idx]][, c("Method", "FWER", "Power")],
          row.names = FALSE, digits = 4)
  }

  results <- do.call(rbind, summary_results)
  row.names(results) <- NULL

  write.csv(results, "adafilter_calibrated_summary.csv", row.names = FALSE)
  saveRDS(replicate_results, "adafilter_calibrated_replicates.rds")
  saveRDS(
    list(settings = list(m = m, n = n, u = u, alpha = alpha,
                         B = B, v_values = v_values, tau = tau, c = c),
         results = results),
    "adafilter_calibrated_results.rds"
  )

  # ---- 5. Side-by-side FWER and power plots ---------------------------------

  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Install ggplot2 and ggpubr to create the plots; results were saved.")
  }

  library(ggplot2)

  results$Method <- factor(results$Method, levels = methods)
  results$v_plot <- factor(
    results$v,
    levels = v_values,
    labels = c("10^-5", "1", "5", "10", "15")
  )

  method_scales <- list(
    scale_color_manual(values = method_colors, limits = methods),
    scale_linetype_manual(values = method_linetypes, limits = methods)
  )

  theme_sim <- theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    )

  common_layers <- c(
    list(
      geom_line(linewidth = 0.8),
      geom_point(size = 2),
      scale_x_discrete(labels = function(x) parse(text = x))
    ),
    method_scales,
    list(
      labs(x = "Variance v"),
      theme_sim,
      guides(
        color = guide_legend(nrow = 2, byrow = TRUE),
        linetype = guide_legend(nrow = 2, byrow = TRUE)
      )
    )
  )

  fwer_plot <- ggplot(
    results,
    aes(x = v_plot, y = FWER, color = Method,
        linetype = Method, group = Method)
  ) +
    common_layers +
    geom_hline(yintercept = alpha, color = "grey40",
               linetype = "dashed") +
    coord_cartesian(ylim = c(0, alpha * 1.2)) +
    labs(y = "Empirical FWER")

  power_plot <- ggplot(
    results,
    aes(x = v_plot, y = Power, color = Method,
        linetype = Method, group = Method)
  ) +
    common_layers +
    coord_cartesian(ylim = c(0, 1)) +
    labs(y = "Empirical power")

  # FWER on the left (a), power on the right (b); one shared legend.
  combined_plot <- ggpubr::ggarrange(
    fwer_plot,
    power_plot,
    ncol = 2,
    nrow = 1,
    labels = c("(a)", "(b)"),
    align = "h",
    common.legend = TRUE,
    legend = "bottom"
  )

  print(combined_plot)

  ggsave(
    "adafilter_calibrated_fwer_power.pdf",
    combined_plot,
    width = 14,
    height = 5.8
  )

  ggsave(
    "adafilter_calibrated_fwer_power.png",
    combined_plot,
    width = 14,
    height = 5.8,
    dpi = 300
  )
})

# =============================================================================
# 4. GSE30495 REAL-DATA ANALYSIS
# =============================================================================

if (run_real_data) local({
  # =============================================================================
  # GSE30495: raw-data preprocessing and partial conjunction method comparison
  # =============================================================================
  # Fully self-contained script for the real-data analysis in the AdaFilter-Sidak
  # manuscript. Runs five independent WT-vs-KO comparisons, constructs p_mat,
  # and compares the eight PC-testing procedures.
  #
  # Run from the repository directory:
  #   Rscript real_data_gse30495.R
  #
  # Raw data and GEO metadata are downloaded as needed. Outputs are written to
  # ./GSE30495_independent. No machine-specific working directory is required.
  #
  # IMPORTANT: This analysis calculates ordinary (unmoderated) t-statistics
  # from limma::lmFit. It intentionally does NOT call limma::eBayes(), so that
  # the numerical calculations from the original GSE30495.R remain unchanged.
  # =============================================================================

  # ---- Functions implementing the eight PC testing procedures -----------------
  # All rejection calculations are retained from real_data_method_comparison.R.

  compare_pc_methods <- function(p_mat, alpha = 0.05, u_values = NULL,
                                 tau_cofilter = 0.10,
                                 na_action = c("error", "omit"),
                                 return_decisions = FALSE) {
    na_action <- match.arg(na_action)

    if (!is.matrix(p_mat) && !is.data.frame(p_mat)) {
      stop("p_mat must be a numeric matrix or data frame.")
    }
    p_mat <- as.matrix(p_mat)
    if (!is.numeric(p_mat)) stop("All columns of p_mat must be numeric.")
    if (ncol(p_mat) < 2L || nrow(p_mat) < 1L) {
      stop("p_mat must have at least one row and at least two columns.")
    }
    if (length(alpha) != 1L || !is.finite(alpha) ||
        alpha <= 0 || alpha >= 1) {
      stop("alpha must be a single number strictly between 0 and 1.")
    }
    if (length(tau_cofilter) != 1L || !is.finite(tau_cofilter) ||
        tau_cofilter <= 0 || tau_cofilter > 1) {
      stop("tau_cofilter must be in (0, 1].")
    }

    n <- ncol(p_mat)
    if (is.null(u_values)) u_values <- seq.int(2L, n)
    if (!is.numeric(u_values) || !length(u_values) ||
        anyNA(u_values) || any(!is.finite(u_values)) ||
        any(u_values != floor(u_values)) ||
        any(u_values < 2L | u_values > n)) {
      stop("u_values must contain integers from 2 through ncol(p_mat).")
    }
    u_values <- sort(unique(as.integer(u_values)))

    if (na_action == "omit") {
      keep <- stats::complete.cases(p_mat)
      if (!all(keep)) {
        message("Dropping ", sum(!keep), " rows containing missing p-values.")
        p_mat <- p_mat[keep, , drop = FALSE]
      }
    }
    if (!nrow(p_mat)) stop("No complete rows of p-values remain.")
    if (any(!is.finite(p_mat)) || any(p_mat < 0 | p_mat > 1)) {
      stop("All p-values must be finite and in [0, 1]. Set na_action='omit' for missing rows.")
    }
    if (is.null(rownames(p_mat))) {
      rownames(p_mat) <- paste0("Feature_", seq_len(nrow(p_mat)))
    } else if (anyDuplicated(rownames(p_mat))) {
      rownames(p_mat) <- make.unique(rownames(p_mat))
    }

    # Stable Sidak PC transform, valid also at p=0 and p=1.
    sidak_pc <- function(p, k) -expm1(k * log1p(-p))

    # AdaFilter threshold characterization, strictly F_i < t.
    # In sequential mode, remove all current discoveries and repeat.
    adafilter <- function(S, F, multiplicity = c("bon", "sidak"),
                          sequential = FALSE) {
      multiplicity <- match.arg(multiplicity)
      rejected <- rep(FALSE, length(S))
      repeat {
        remaining <- which(!rejected)
        if (!length(remaining)) break
        # Number of remaining filter statistics strictly below S_i.
        nf <- findInterval(S[remaining], sort(F[remaining]),
                           left.open = TRUE)
        score <- if (multiplicity == "bon") {
          nf * S[remaining]
        } else {
          ans <- numeric(length(remaining))
          nz <- nf > 0L
          ans[nz] <- -expm1(nf[nz] * log1p(-S[remaining[nz]]))
          ans
        }
        new <- remaining[S[remaining] <= alpha & score <= alpha]
        if (!length(new)) break
        rejected[new] <- TRUE
        if (!sequential) break
      }
      rejected
    }

    method_names <- c(
      "AdaFilter-Bon", "AdaFilter-SeqBon",
      "AdaFilter-Bon (Sidak PC)", "AdaFilter-SeqBon (Sidak PC)",
      "AdaFilter-Sid", "AdaFilter-SeqSid",
      "CoFilter-Holm", "Hommel-Fisher"
    )

    # Sorting once supports every u at the same n.
    P_sort <- t(apply(p_mat, 1L, sort))
    dimnames(P_sort) <- list(rownames(p_mat), NULL)

    counts <- matrix(0L, nrow = length(method_names),
                     ncol = length(u_values),
                     dimnames = list(method_names, paste0(u_values, "/", n)))
    decisions <- if (return_decisions) vector("list", length(u_values)) else NULL
    if (return_decisions) names(decisions) <- colnames(counts)

    for (u in u_values) {
      k <- n - u + 1L
      S_bon <- k * P_sort[, u]
      F_bon <- k * P_sort[, u - 1L]
      S_sidak <- sidak_pc(P_sort[, u], k)
      F_sidak <- sidak_pc(P_sort[, u - 1L], k)
      S_fisher <- stats::pchisq(
        -2 * rowSums(log(P_sort[, u:n, drop = FALSE])),
        df = 2L * k, lower.tail = FALSE
      )

      # Fisher PC p-values are filtered at a *pre-specified* tau.
      # Selected p-values are converted to conditional p-values by p/tau.
      cofilter_reject <- rep(FALSE, nrow(p_mat))
      selected <- which(S_fisher <= tau_cofilter)
      if (length(selected)) {
        cofilter_reject[selected] <-
          stats::p.adjust(S_fisher[selected] / tau_cofilter,
                          method = "holm") <= alpha
      }

      dec <- cbind(
        "AdaFilter-Bon" = adafilter(S_bon, F_bon, "bon"),
        "AdaFilter-SeqBon" = adafilter(S_bon, F_bon, "bon", TRUE),
        "AdaFilter-Bon (Sidak PC)" = adafilter(S_sidak, F_sidak, "bon"),
        "AdaFilter-SeqBon (Sidak PC)" = adafilter(S_sidak, F_sidak,
                                                  "bon", TRUE),
        "AdaFilter-Sid" = adafilter(S_sidak, F_sidak, "sidak"),
        "AdaFilter-SeqSid" = adafilter(S_sidak, F_sidak,
                                         "sidak", TRUE),
        "CoFilter-Holm" = cofilter_reject,
        "Hommel-Fisher" = stats::p.adjust(S_fisher,
                                          method = "hommel") <= alpha
      )
      rownames(dec) <- rownames(p_mat)
      counts[, paste0(u, "/", n)] <- colSums(dec)
      if (return_decisions) decisions[[paste0(u, "/", n)]] <- dec
    }

    summary <- data.frame(Method = rownames(counts), counts,
                          row.names = NULL, check.names = FALSE)
    if (return_decisions) {
      return(list(summary = summary, decisions = decisions,
                  alpha = alpha, tau_cofilter = tau_cofilter,
                  u_values = u_values, m = nrow(p_mat), n = n))
    }
    summary
  }

  # ============================================================
  # 0. Install and load packages
  # ============================================================

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  pkgs <- c(
    "GEOquery",
    "Biobase",
    "affy",
    "mouse4302cdf",
    "limma"
  )

  missing_pkgs <- pkgs[
    !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_pkgs) > 0L) {
    BiocManager::install(
      missing_pkgs,
      ask = FALSE,
      update = FALSE
    )
  }

  library(GEOquery)
  library(Biobase)
  library(affy)
  library(limma)

  # ============================================================
  # 1. Settings
  # ============================================================

  gse_id <- "GSE30495"
  gpl_id <- "GPL1261"

  conditions <- c(
    "C183TG",
    "C182TG",
    "C181TG",
    "C226TG",
    "WY"
  )

  alpha <- 0.05

  out_dir <- file.path(getwd(), "GSE30495_independent")

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  # ============================================================
  # 2. Load GEO sample metadata
  # ============================================================

  gsets <- getGEO(
    gse_id,
    GSEMatrix = TRUE,
    AnnotGPL = FALSE,
    getGPL = FALSE,
    destdir = out_dir
  )

  idx <- which(
    vapply(gsets, annotation, character(1)) == gpl_id
  )

  stopifnot(length(idx) == 1L)

  gset <- gsets[[idx]]

  metadata <- pData(gset)

  stopifnot("title" %in% colnames(metadata))

  metadata$gsm <- sampleNames(gset)
  metadata$title <- trimws(as.character(metadata$title))

  # Extract groups directly from sample titles
  metadata$group <- sub(
    "_rep[0-9]+$",
    "",
    metadata$title
  )

  # Keep only the five WT-KO comparisons
  groups <- c(
    paste0("WT_", conditions),
    paste0("KO_", conditions)
  )

  metadata <- metadata[
    metadata$group %in% groups,
    ,
    drop = FALSE
  ]

  # Verify sample membership
  group_counts <- table(
    factor(metadata$group, levels = groups)
  )

  print(group_counts)

  stopifnot(
    nrow(metadata) == 47L,
    !anyDuplicated(metadata$gsm),
    all(group_counts >= 3L)
  )

  # ============================================================
  # 3. Download raw CEL files
  # ============================================================

  raw_base <- file.path(out_dir, "raw")

  dir.create(
    raw_base,
    recursive = TRUE,
    showWarnings = FALSE
  )

  tar_file <- file.path(
    raw_base,
    gse_id,
    paste0(gse_id, "_RAW.tar")
  )

  if (!file.exists(tar_file)) {
    
    getGEOSuppFiles(
      gse_id,
      makeDirectory = TRUE,
      baseDir = raw_base,
      filter_regex = "_RAW\\.tar$"
    )
    
  }

  if (!file.exists(tar_file)) {
    stop("Raw CEL archive could not be downloaded.")
  }

  # ============================================================
  # 4. Extract CEL files
  # ============================================================

  cel_dir <- file.path(out_dir, "CEL")

  dir.create(
    cel_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  archive_files <- utils::untar(
    tar_file,
    list = TRUE
  )

  archive_cels <- archive_files[
    grepl("\\.CEL(\\.gz)?$", archive_files, ignore.case = TRUE)
  ]

  if (length(archive_cels) == 0L) {
    stop("No CEL files found in the GEO archive.")
  }

  expected_files <- file.path(cel_dir, archive_cels)

  if (!all(file.exists(expected_files))) {
    
    utils::untar(
      tar_file,
      exdir = cel_dir
    )
    
  }

  cel_files <- list.files(
    cel_dir,
    pattern = "\\.CEL(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  if (length(cel_files) == 0L) {
    stop("No extracted CEL files found.")
  }

  # Match CEL filenames to GSM accessions
  cel_names <- basename(cel_files)

  has_gsm <- grepl(
    "GSM[0-9]+",
    cel_names,
    ignore.case = TRUE
  )

  if (!all(has_gsm)) {
    stop("Some CEL filenames contain no GSM accession.")
  }

  cel_gsm <- toupper(
    sub(
      "^.*?(GSM[0-9]+).*$",
      "\\1",
      cel_names,
      perl = TRUE,
      ignore.case = TRUE
    )
  )

  if (anyDuplicated(cel_gsm)) {
    stop("Duplicate GSM accessions among CEL files.")
  }

  match_idx <- match(metadata$gsm, cel_gsm)

  if (anyNA(match_idx)) {
    stop(
      "Missing CEL files for: ",
      paste(metadata$gsm[is.na(match_idx)], collapse = ", ")
    )
  }

  metadata$cel_file <- cel_files[match_idx]

  stopifnot(
    !anyDuplicated(metadata$cel_file),
    all(file.exists(metadata$cel_file))
  )

  write.csv(
    metadata[, c("gsm", "title", "group", "cel_file")],
    file.path(out_dir, "sample_audit.csv"),
    row.names = FALSE
  )

  # ============================================================
  # 5. Function: independent analysis for one treatment
  # ============================================================

  analyse_condition <- function(cond, metadata) {
    
    cat("\n========================================\n")
    cat("Analysing:", cond, "\n")
    cat("========================================\n")
    
    wt <- paste0("WT_", cond)
    ko <- paste0("KO_", cond)
    
    # Select only the samples for this contrast
    dat <- metadata[
      metadata$group %in% c(wt, ko),
      ,
      drop = FALSE
    ]
    
    dat$group <- factor(
      dat$group,
      levels = c(ko, wt)
    )
    
    print(table(dat$group))
    
    stopifnot(
      all(table(dat$group) >= 3L),
      !anyDuplicated(dat$gsm)
    )
    
    # ----------------------------------------------------------
    # A. Normalize ONLY these samples
    # ----------------------------------------------------------
    
    files <- dat$cel_file
    
    compressed <- grepl(
      "\\.gz$",
      files,
      ignore.case = TRUE
    )
    
    if (!all(compressed) && any(compressed)) {
      stop("Mixed compressed and uncompressed CEL files.")
    }
    
    stopifnot(all(file.exists(files)))
    
    eset <- affy::justRMA(
      filenames = files,
      celfile.path = NULL,
      sampleNames = dat$gsm,
      compress = all(compressed),
      widget = FALSE,
      normalize = TRUE,
      background = TRUE,
      verbose = FALSE
    )
    
    ex <- Biobase::exprs(eset)
    
    stopifnot(
      ncol(ex) == nrow(dat),
      identical(colnames(ex), dat$gsm),
      all(is.finite(ex))
    )
    
    # ----------------------------------------------------------
    # B. Fit two-group linear model
    # ----------------------------------------------------------
    
    design <- model.matrix(~ group, data = dat)
    
    fit <- limma::lmFit(ex, design)
    
    # ----------------------------------------------------------
    # C. Ordinary t-test: NO empirical Bayes moderation
    # ----------------------------------------------------------
    
    beta <- fit$coefficients[, 2]
    
    se <- fit$stdev.unscaled[, 2] * fit$sigma
    
    df <- fit$df.residual
    
    if (any(
      !is.finite(beta) |
      !is.finite(se) |
      se <= 0 |
      !is.finite(df) |
      df <= 0
    )) {
      stop("Invalid t-test quantities in condition: ", cond)
    }
    
    t_stat <- beta / se
    
    p_val <- 2 * stats::pt(
      abs(t_stat),
      df = df,
      lower.tail = FALSE
    )
    
    names(p_val) <- rownames(ex)
    
    stopifnot(
      all(is.finite(p_val)),
      all(p_val >= 0 & p_val <= 1)
    )
    
    cat("Probes:", length(p_val), "\n")
    cat("Minimum p-value:", min(p_val), "\n")
    
    return(p_val)
  }

  # ============================================================
  # 6. Analyse all five conditions independently
  # ============================================================

  p_list <- lapply(
    conditions,
    analyse_condition,
    metadata = metadata
  )

  names(p_list) <- conditions

  # Confirm the same probe IDs and ordering
  reference_ids <- names(p_list[[1]])

  stopifnot(
    all(vapply(
      p_list,
      function(p) identical(names(p), reference_ids),
      logical(1)
    ))
  )

  # ============================================================
  # 7. Construct p-value matrix
  # ============================================================

  p_mat <- do.call(cbind, p_list)

  rownames(p_mat) <- reference_ids
  colnames(p_mat) <- conditions

  stopifnot(
    ncol(p_mat) == 5L,
    !anyDuplicated(rownames(p_mat)),
    all(is.finite(p_mat)),
    all(p_mat >= 0 & p_mat <= 1)
  )

  cat("\nP-value matrix dimensions:\n")
  print(dim(p_mat))

  cat("\nFirst six rows:\n")
  print(head(p_mat))

  # ============================================================
  # 8. Save p-value matrix
  # ============================================================

  saveRDS(
    p_mat,
    file.path(out_dir, "p_mat.rds")
  )

  write.csv(
    p_mat,
    file.path(out_dir, "p_mat.csv")
  )

  # ============================================================
  # 9. Compare partial conjunction methods
  # ============================================================

  # Compare all eight methods, retaining the original testing calculations.
  results <- compare_pc_methods(
    p_mat = p_mat,
    alpha = alpha
  )

  cat("\nAll eight procedures (rejection counts):\n")
  print(results, row.names = FALSE)

  # Preserve the original output name for the full eight-method comparison.
  saveRDS(results, file.path(out_dir, "results.rds"))
  write.csv(results, file.path(out_dir, "results_all_methods.csv"),
            row.names = FALSE)

  # The manuscript's GSE30495 table reports only the six AdaFilter variants.
  adafilter_methods <- c(
    "AdaFilter-Bon",
    "AdaFilter-SeqBon",
    "AdaFilter-Bon (Sidak PC)",
    "AdaFilter-SeqBon (Sidak PC)",
    "AdaFilter-Sid",
    "AdaFilter-SeqSid"
  )
  adafilter_results <- results[
    match(adafilter_methods, results$Method),
    , drop = FALSE
  ]

  cat("\nSix AdaFilter variants (manuscript table):\n")
  print(adafilter_results, row.names = FALSE)
  write.csv(adafilter_results,
            file.path(out_dir, "results_adafilter_table.csv"),
            row.names = FALSE)

  # Compare computed counts with the manuscript table (diagnostic only).
  manuscript_counts <- rbind(
    c(164L, 144L, 129L, 104L),
    c(185L, 175L, 157L, 138L),
    c(164L, 144L, 129L, 104L),
    c(185L, 175L, 157L, 138L),
    c(165L, 145L, 129L, 105L),
    c(185L, 178L, 157L, 139L)
  )
  dimnames(manuscript_counts) <- list(
    adafilter_methods, c("2/5", "3/5", "4/5", "5/5")
  )
  observed_counts <- as.matrix(adafilter_results[, -1L, drop = FALSE])
  storage.mode(observed_counts) <- "integer"
  rownames(observed_counts) <- adafilter_results$Method
  if (identical(observed_counts, manuscript_counts)) {
    message("GSE30495 CHECK: All six AdaFilter rows match the manuscript table.")
  } else {
    message("GSE30495 CHECK: Rejection counts differ from the manuscript table.")
    print(observed_counts - manuscript_counts)
  }
})
