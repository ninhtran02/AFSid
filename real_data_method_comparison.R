# =============================================================================
# Real-data method comparison: eight FWER-controlling PC procedures
# =============================================================================
# Source this file and call compare_pc_methods(P_mat).
# P_mat must already contain the two-sided base p-values: rows are probe sets,
# columns are studies. This file does NOT download GSE30495, perform RMA
# preprocessing, or compute limma moderated t-test p-values.
# =============================================================================
# Input: m x n matrix, one feature per row, one study per column.
# Defaults: every u=2,...,n, alpha=.05, fixed CoFilter tau=.1.
# Returns a data frame of rejection counts (or detailed results if requested).
# Requires only base R; does not download or preprocess any data.

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

# EXAMPLES (after sourcing this file and constructing P_mat):
# summary <- compare_pc_methods(P_mat)
# print(summary, row.names = FALSE)
# View(summary)
#
# summary_3of5 <- compare_pc_methods(P_mat_5, u_values = 3)
# details <- compare_pc_methods(P_mat, return_decisions = TRUE)
# View(details$summary)
# head(details$decisions[["2/4"]])
