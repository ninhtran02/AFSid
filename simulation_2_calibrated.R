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
