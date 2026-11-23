# helper-sim.R -- simulated data and missingness generators for the mvImputer tests.
#
# Purpose: deterministic, seed-scoped test data with known structure (so accuracy
#   thresholds in later tests have signal) and controlled missingness.
# Estimand: none; these functions generate data, they estimate nothing.
# Usage:  X <- sim_numeric(100, 20, seed = 1)
#         d <- sim_mixed(300, seed = 1, extra = 10)
#         Y <- punch_holes(X, 0.1, "mcar", seed = 2)
#
# Dependencies: withr (scoped seeds, so the caller's .Random.seed and RNGkind are
#   untouched); stats (rnorm, quantile).
#
# The data-generating processes below are pinned: thresholds in the Task 6-8 tests
# were set against them. Changing a constant or the ORDER of random draws changes
# every simulated value for a given seed.

SIM_N_LATENT   <- 3     # latent factors behind sim_numeric()
SIM_NOISE_SD   <- 0.5   # idiosyncratic noise on the log scale in sim_numeric()
SIM_SIGNAL_SCL <- 0.5   # scale of the latent signal on the log scale

# --- sim_numeric: positive, log-normal, low-rank + noise (omics-like) -----------
# Returns an n x p numeric matrix, colnames F1..Fp. Exp of a rank-3 Gaussian signal
# plus noise gives correlated, right-skewed, strictly positive columns.
sim_numeric <- function(n, p, seed) {
  stopifnot(is.numeric(n), length(n) == 1, n >= 1,           # BLOCKING
            is.numeric(p), length(p) == 1, p >= 1)
  withr::with_seed(seed, {
    Z <- matrix(stats::rnorm(n * SIM_N_LATENT), n, SIM_N_LATENT)
    L <- matrix(stats::rnorm(SIM_N_LATENT * p), SIM_N_LATENT, p)
    E <- matrix(stats::rnorm(n * p, sd = SIM_NOISE_SD), n, p)
    X <- exp(SIM_SIGNAL_SCL * Z %*% L + E)
  })
  colnames(X) <- paste0("F", seq_len(p))
  X
}

# --- sim_mixed: mixed-type table whose non-numeric columns are predictable ------
# Latents z1, z2 drive the numeric columns; each non-numeric column is a noisy
# function of numeric columns, so a good imputer can beat the marginal guess.
# e(s) draws fresh normal noise of SD s at every call (draw order is part of the
# pinned process). `extra` adds weakly informative noise columns noise1..
sim_mixed <- function(n, seed, extra = 0) {
  stopifnot(is.numeric(n), length(n) == 1, n >= 4,           # BLOCKING: quartiles need n >= 4
            is.numeric(extra), length(extra) == 1, extra >= 0)
  withr::with_seed(seed, {
    e  <- function(s) stats::rnorm(n, sd = s)
    z1 <- stats::rnorm(n)
    z2 <- stats::rnorm(n)

    num1 <- z1 + e(0.3)
    num2 <- z1 + z2 + e(0.3)
    num3 <- z2 + e(0.3)
    num4 <- z1 - z2 + e(0.3)
    int1 <- as.integer(round(5 + 2 * z1 + e(0.5)))
    # ties.method = "first": ties have probability zero with continuous input, but
    # the default ("random") would consume RNG draws if one ever occurred.
    fac3 <- factor(c("a", "b", "c")[max.col(cbind(num1, num2, num3) + e(0.3),
                                            ties.method = "first")],
                   levels = c("a", "b", "c"))
    ord_driver <- num4 + e(0.3)
    ord4 <- cut(ord_driver, breaks = stats::quantile(ord_driver, probs = 0:4 / 4),
                labels = c("L1", "L2", "L3", "L4"), include.lowest = TRUE,
                ordered_result = TRUE)
    chr2  <- ifelse(num2 + e(0.3) > 0, "x", "y")
    lgl   <- num3 + e(0.3) > 0
    const <- rep(5, n)
    d <- data.frame(num1, num2, num3, num4, int1, fac3, ord4, chr2, lgl, const,
                    stringsAsFactors = FALSE)
    for (k in seq_len(extra)) d[[paste0("noise", k)]] <- 0.5 * z1 + e(1)
  })
  # BLOCKING: no NAs and the ordered factor has all four levels in use.
  stopifnot(!anyNA(d), nlevels(droplevels(d$ord4)) == 4)
  d
}

# --- punch_holes: impose missingness on a matrix or data.frame ------------------
# "mcar": exactly round(frac * cells-in-cols) cells, uniform over cols, any type.
# "left": in each NUMERIC column of cols, the lowest round(frac * n) values become
#   NA (detection-limit pattern: deterministic, so `seed` has no effect there).
# cols: NULL (all), integer positions or character names.
punch_holes <- function(X, frac, mechanism = c("mcar", "left"), seed, cols = NULL) {
  mechanism <- match.arg(mechanism)
  stopifnot(is.matrix(X) || is.data.frame(X),                # BLOCKING
            is.numeric(frac), length(frac) == 1, frac >= 0, frac <= 1)
  n <- nrow(X)
  if (is.null(cols)) cols <- seq_len(ncol(X))
  if (is.character(cols)) {
    cols <- match(cols, colnames(X))
  }
  stopifnot(!anyNA(cols), all(cols >= 1 & cols <= ncol(X)))  # BLOCKING: cols must exist
  is_df <- is.data.frame(X)
  set_na <- function(X, rows, j) {
    if (is_df) X[[j]][rows] <- NA else X[rows, j] <- NA
    X
  }
  col_is_numeric <- function(j) is.numeric(if (is_df) X[[j]] else X[, j])

  withr::with_seed(seed, {
    if (mechanism == "mcar") {
      n_cells <- n * length(cols)
      hit  <- sample.int(n_cells, round(frac * n_cells))     # cell index, column-major in cols
      rows <- (hit - 1L) %% n + 1L
      jj   <- cols[(hit - 1L) %/% n + 1L]
      for (j in unique(jj)) X <- set_na(X, rows[jj == j], j)
    } else {
      k <- round(frac * n)
      for (j in cols[vapply(cols, col_is_numeric, TRUE)]) {
        v <- if (is_df) X[[j]] else X[, j]
        X <- set_na(X, order(v)[seq_len(k)], j)
      }
    }
  })
  X
}
