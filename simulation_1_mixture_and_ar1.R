# =============================================================================
# Simulation setting 1: independent and AR(1)-dependent two-point mixtures
# =============================================================================
# Reproduces the main-text independent mixture simulation and the appendix
# AR(1) simulations (rho = -0.8, +0.8), including published plot filenames.
#
# The two original simulations are intentionally isolated in local environments
# so that running them together does not change settings or random-number seeds.
# No simulation formulas, rejection rules, or Monte Carlo counts were changed.
#
# Set either flag to FALSE to skip that simulation and its figures.
# Run with: Rscript simulation_1_mixture_and_ar1.R
# =============================================================================

run_independent <- TRUE
run_ar1 <- TRUE

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
