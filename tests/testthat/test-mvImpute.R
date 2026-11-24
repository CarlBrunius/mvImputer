# test-mvImpute.R -- the exported mvImpute() end to end (spec sections 5.1, 5.5, 6).
#
# Never more than 2 workers except in tests behind skip_on_cran(). Multisession workers
# load the INSTALLED package, so run these with load_package = 'installed'.
#
# Warnings: the package's own final-stage not-converged warning is muffled where it is
# not under test (imp() below). Nothing else is ever muffled here: backend warnings must
# stay visible to the tests that look for them (Review Focus 5).

strip <- function(x) { attr(x, "mvImpute") <- NULL; x }

# The pinned text of the package's own warning; any other warning is left alone.
NC_WARNING <- "imputed cells not converged after maxIter in the final stage"
quiet_nc <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl(NC_WARNING, conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
}
imp <- function(...) quiet_nc(mvImpute(...))

# PIDs of both workers of a 2-worker plan; the sleep keeps the first future busy so the
# second runs on the other worker.
pool_pids <- function() {
  fs <- lapply(1:2, function(i) future::future({ Sys.sleep(0.5); Sys.getpid() }))
  sort(unlist(lapply(fs, future::value)))
}

d0 <- sim_mixed(150, seed = 1)
dh <- punch_holes(d0, 0.08, "mcar", seed = 2, cols = c(1:4, 6:9))

test_that("same seed, same data: identical with 1 and 2 workers on mixed data", {
  a <- imp(dh, seed = 1, workers = 1, maxIter = c(3, 5))
  b <- imp(dh, seed = 1, workers = 2, maxIter = c(3, 5))
  expect_identical(strip(a), strip(b))
  expect_identical(attr(a, "mvImpute")$iterations, attr(b, "mvImpute")$iterations)
  expect_identical(attr(a, "mvImpute")$workers, 1L)
  expect_identical(attr(b, "mvImpute")$workers, 2L)
})
test_that("identical with 4 workers too", {
  skip_on_cran()
  expect_identical(strip(imp(dh, seed = 1, workers = 1, maxIter = c(3, 5))),
                   strip(imp(dh, seed = 1, workers = 4, maxIter = c(3, 5))))
})
test_that("different seeds differ; set.seed reproduces seed = NULL", {
  expect_false(identical(strip(imp(dh, seed = 1, workers = 1, maxIter = c(2, 2))),
                         strip(imp(dh, seed = 2, workers = 1, maxIter = c(2, 2)))))
  withr::local_preserve_seed()
  set.seed(5); a <- imp(dh, workers = 1, maxIter = c(2, 2))
  set.seed(5); b <- imp(dh, workers = 1, maxIter = c(2, 2))
  expect_identical(strip(a), strip(b))
  expect_identical(attr(a, "mvImpute")$seed, attr(b, "mvImpute")$seed)
})
test_that("caller RNG untouched with a seed, one draw without", {
  withr::local_preserve_seed()
  set.seed(3); s0 <- .Random.seed; imp(dh, seed = 1, workers = 1, maxIter = c(2, 2))
  expect_identical(.Random.seed, s0)
  set.seed(3); imp(dh, workers = 1, maxIter = c(2, 2)); s1 <- .Random.seed
  set.seed(3); sample.int(.Machine$integer.max, 1L); expect_identical(s1, .Random.seed)
})
test_that("a non-default caller RNG kind and its seed are restored", {
  withr::local_seed(4, .rng_kind = "Wichmann-Hill")
  kind0 <- RNGkind(); s0 <- .Random.seed
  imp(dh, seed = 1, workers = 1, maxIter = c(2, 2))
  expect_identical(RNGkind(), kind0)
  expect_identical(.Random.seed, s0)
})
test_that("an absent .Random.seed stays absent, kinds unchanged", {
  withr::local_preserve_seed()
  if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv())
  kind0 <- RNGkind()
  imp(dh, seed = 1, workers = 1, maxIter = c(2, 2))
  expect_false(exists(".Random.seed", globalenv()))
  expect_identical(RNGkind(), kind0)
})

test_that("observed values and types come back unchanged; no NA left", {
  out <- imp(dh, seed = 1, workers = 1, maxIter = c(3, 5))
  for (nm in names(dh)) { o <- !is.na(dh[[nm]])
    expect_identical(out[[nm]][o], dh[[nm]][o]) }
  expect_false(anyNA(out)); expect_identical(class(out$ord4), c("ordered", "factor"))
  expect_identical(levels(out$fac3), levels(dh$fac3))
  expect_type(out$chr2, "character"); expect_type(out$lgl, "logical")
  expect_identical(rownames(out), rownames(dh))
})
test_that("an integer column with gaps comes back as double, observed values by value", {
  d <- dh; d$int1[c(2, 9, 30)] <- NA
  out <- imp(d, seed = 1, workers = 1, maxIter = c(3, 5))
  o <- !is.na(d$int1)
  expect_type(out$int1, "double"); expect_equal(out$int1[o], d$int1[o])
  expect_false(anyNA(out$int1))
  # An integer MATRIX with gaps is double as a whole.
  M <- matrix(as.integer(round(100 * sim_numeric(40, 4, seed = 1))), 40, 4)
  M[c(3, 8), 2] <- NA
  expect_type(imp(M, seed = 1, workers = 1, maxIter = c(2, 2)), "double")
})

test_that("exclude is invisible to the models (excluded column first, so W positions shift)", {
  # y is built from num1 so it WOULD inform the models if it leaked in; d2 has unrelated values.
  y1 <- d0$num1 + withr::with_seed(9, rnorm(150, sd = 0.1)); y1[1:10] <- NA
  d1 <- data.frame(y = y1, dh)
  d2 <- data.frame(y = withr::with_seed(10, rnorm(150)), dh); d2$y[1:10] <- NA
  a <- imp(d1, exclude = "y", seed = 1, workers = 1, maxIter = c(3, 5))
  b <- imp(d2, exclude = "y", seed = 1, workers = 1, maxIter = c(3, 5))
  none <- imp(dh, seed = 1, workers = 1, maxIter = c(3, 5))
  expect_identical(a$y, d1$y)
  expect_identical(b$y, d2$y)
  expect_identical(strip(a)[names(dh)], strip(b)[names(dh)])
  expect_identical(strip(a)[names(dh)], strip(none))
})
test_that("exclude is invisible under guess = 'minTotal' when the excluded column holds the global minimum", {
  # If y reached first_guess(), numeric_total_min() would be about -1000 and every guess would move.
  d1 <- data.frame(y = c(rep(NA, 5), d0$num1[-(1:5)] - 1000), dh)
  d2 <- data.frame(y = c(rep(NA, 5), d0$num2[-(1:5)] - 500),  dh)
  a <- imp(d1, exclude = "y", guess = "minTotal", seed = 1, workers = 1, maxIter = c(3, 5))
  b <- imp(d2, exclude = "y", guess = "minTotal", seed = 1, workers = 1, maxIter = c(3, 5))
  none <- imp(dh, guess = "minTotal", seed = 1, workers = 1, maxIter = c(3, 5))
  expect_identical(strip(a)[names(dh)], strip(none))
  expect_identical(strip(b)[names(dh)], strip(none))
  expect_identical(a$y, d1$y)
})

test_that("excluded data never reach a worker", {
  # A 4.8 MB excluded column against a 1 MB cap on what future may send to a worker: had the
  # column (or enc, or X) been exported with a job, future would refuse to start it.
  X <- punch_holes(sim_numeric(60, 6, seed = 1), 0.08, "mcar", seed = 2)
  d <- data.frame(X)
  d$big <- I(lapply(seq_len(nrow(d)), function(i) numeric(1e4)))
  stopifnot(utils::object.size(d$big) > 4 * 1024^2)               # BLOCKING: the probe is big
  withr::local_options(future.globals.maxSize = 1024^2)
  out <- imp(d, exclude = "big", seed = 1, workers = 2, maxIter = c(2, 2))
  expect_identical(attr(out, "mvImpute")$workers, 2L)
  expect_identical(out$big, d$big)
  expect_identical(strip(out)[colnames(X)],
                   strip(imp(data.frame(X), seed = 1, workers = 1, maxIter = c(2, 2))))
})

test_that("caller RNG state and kind are restored at workers = 2; seed = NULL costs one draw", {
  withr::local_seed(4, .rng_kind = "Wichmann-Hill")   # a non-default kind, so a leak shows
  kind0 <- RNGkind(); s0 <- .Random.seed
  imp(dh, seed = 1, workers = 2, maxIter = c(2, 2))
  expect_identical(RNGkind(), kind0)
  expect_identical(.Random.seed, s0)
  s_one <- withr::with_preserve_seed({ sample.int(.Machine$integer.max, 1L); .Random.seed })
  imp(dh, workers = 2, maxIter = c(2, 2))
  expect_identical(.Random.seed, s_one)
  expect_identical(RNGkind(), kind0)
})

test_that("randomForest and pls: identical with 1 and 2 workers", {
  skip_on_cran()
  for (m in c("randomForest", "pls")) {
    if (m == "pls") skip_if_not_installed("mixOmics")
    a <- imp(dh, method = m, seed = 1, workers = 1, maxIter = c(2, 3))
    b <- imp(dh, method = m, seed = 1, workers = 2, maxIter = c(2, 3))
    expect_identical(strip(a), strip(b), info = m)
    expect_identical(attr(a, "mvImpute")$oob_error, attr(b, "mvImpute")$oob_error, info = m)
  }
})

test_that("an empty column name is refused before any work", {
  X <- punch_holes(sim_numeric(40, 4, seed = 1), 0.1, "mcar", seed = 2)
  colnames(X)[3] <- ""
  expect_error(mvImpute(X, seed = 1, workers = 1), "Column names must not be empty or NA")
})

test_that("numeric matrix input: class, dimnames and diagnostics", {
  X <- punch_holes(sim_numeric(60, 8, seed = 1), 0.05, "mcar", seed = 2)
  out <- imp(X, seed = 1, workers = 1, maxIter = c(3, 5))
  expect_true(is.matrix(out)); expect_identical(dimnames(out), dimnames(X)); expect_false(anyNA(out))
  expect_named(attr(out, "mvImpute"), c("seed", "rng_kind", "method", "settings", "workers",
    "iterations", "not_converged", "time", "oob_error", "rows_no_observed"), ignore.order = TRUE)
  dg <- attr(out, "mvImpute")
  expect_identical(dg$rng_kind, "L'Ecuyer-CMRG"); expect_identical(dg$method, "ranger")
  expect_identical(dg$seed, 1L); expect_type(dg$iterations, "integer"); expect_length(dg$iterations, 2L)
  expect_length(dg$time, 2L); expect_type(dg$time, "double")
  expect_identical(dg$rows_no_observed, integer())
  expect_identical(dg$settings$maxIter, c(3, 5))
})
test_that("an unnamed matrix returns unnamed; diagnostics use V<k>", {
  X <- punch_holes(sim_numeric(50, 5, seed = 1), 0.1, "mcar", seed = 2); dimnames(X) <- NULL
  out <- imp(X, seed = 1, workers = 1, tol = c(1e-9, 1e-9), maxIter = c(1, 1))
  expect_null(dimnames(out))
  nc <- attr(out, "mvImpute")$not_converged
  expect_true(all(nc[[2]]$column %in% paste0("V", 1:5))); expect_gt(nrow(nc[[2]]), 0L)
  expect_true(all(names(attr(out, "mvImpute")$oob_error) %in% paste0("V", 1:5)))
})

test_that("forceZero floors numeric imputations at zero and leaves observed values alone", {
  X <- punch_holes(sim_numeric(60, 8, seed = 1), 0.05, "mcar", seed = 2); X <- log(X)  # half negative
  out <- imp(X, seed = 1, workers = 1, maxIter = c(2, 2), forceZero = TRUE)
  expect_true(any(X < 0, na.rm = TRUE))
  expect_true(all(out[is.na(X)] >= 0))
  expect_identical(out[!is.na(X)], X[!is.na(X)])
})
test_that("every guess option runs and observed values stay put", {
  X <- punch_holes(sim_numeric(60, 6, seed = 1), 0.05, "mcar", seed = 2)
  for (g in list("minVar", "median", "minTotal", ifelse(is.na(X), 1, X)))
    expect_identical(imp(X, seed = 1, workers = 1, maxIter = c(2, 2), guess = g)[!is.na(X)], X[!is.na(X)])
})
test_that("a row with no observed value warns and is listed", {
  X <- punch_holes(sim_numeric(40, 5, seed = 1), 0.05, "mcar", seed = 2); X[7, ] <- NA
  expect_warning(out <- imp(X, seed = 1, workers = 1, maxIter = c(2, 2)), "no observed value")
  expect_identical(attr(out, "mvImpute")$rows_no_observed, 7L)
})
test_that("NaN counts as missing; an all-zero-but-one column is fine", {
  X <- sim_numeric(40, 5, seed = 1); X[3, 2] <- NaN
  X[, 5] <- 0; X[1, 5] <- 4; X[10, 5] <- NA
  out <- imp(X, seed = 1, workers = 1, maxIter = c(2, 2))
  expect_false(anyNA(out))
})
test_that("each backend runs on mixed data and reports oob", {
  for (m in c("ranger", "randomForest", "pls")) {
    if (m == "pls") skip_if_not_installed("mixOmics")
    out <- imp(dh, method = m, seed = 1, workers = 1, maxIter = c(2, 3))
    expect_false(anyNA(out))
    oob <- attr(out, "mvImpute")$oob_error
    if (m == "pls") expect_true(all(is.na(oob))) else expect_true(all(is.finite(oob)))
  }
})
test_that("imputation beats column means and modes on simulated truth", {
  out <- imp(dh, seed = 1, workers = 1)
  m <- is.na(dh$num1); mean_err <- mean((mean(dh$num1, na.rm = TRUE) - d0$num1[m])^2)
  expect_lt(mean((out$num1[m] - d0$num1[m])^2), mean_err)
  f <- is.na(dh$fac3); mode_acc <- mean(names(which.max(table(dh$fac3))) == d0$fac3[f])
  expect_gt(mean(out$fac3[f] == d0$fac3[f]), mode_acc)
})

test_that("a tiny table (n = 20) with a factor, a sparse and a constant column: no backend warning", {
  d <- sim_mixed(20, seed = 3)[, c("num1", "num2", "num3", "fac3", "const")]
  d$num1[2] <- NA; d$fac3[5] <- NA; d$const[4] <- NA
  d$sparse <- NA_real_; d$sparse[c(1, 8, 14)] <- c(0.3, -1.2, 0.8)
  for (m in c("ranger", "pls")) {
    if (m == "pls") skip_if_not_installed("mixOmics")
    # Only the package's own "not converged" warning is muffled (quiet_nc); any other
    # warning, from a backend or from us, fails expect_no_warning. The constant-column
    # message is expected and passes through suppressMessages.
    expect_no_warning(out <- suppressMessages(
      quiet_nc(mvImpute(d, method = m, seed = 1, workers = 1, maxIter = c(3, 5)))))
    expect_false(anyNA(out))
    expect_true(all(out$const == 5))
  }
})

test_that("only the final stage warns at the cap", {
  expect_warning(mvImpute(dh, seed = 1, workers = 1, tol = c(1e-9, 1e-9), maxIter = c(1, 1)), "not converged")
  nc <- attr(imp(dh, seed = 1, workers = 1, tol = c(1e-9, 1e-9), maxIter = c(1, 1)), "mvImpute")$not_converged
  expect_length(nc, 2L)
})
test_that("the warning names the count, with the pinned text", {
  w <- tryCatch(mvImpute(dh, seed = 1, workers = 1, tol = c(1e-9, 1e-9), maxIter = c(1, 1)),
                warning = function(w) conditionMessage(w))
  expect_match(w, "^[0-9]+ imputed cells not converged after maxIter in the final stage\\.$")
})
test_that("stage 1 at its cap and stage 2 converged gives no warning", {
  out <- NULL
  expect_no_warning(out <- mvImpute(dh, seed = 1, workers = 1, tol = c(1e-9, 5), maxIter = c(1, 40)),
                    message = "not converged")
  nc <- attr(out, "mvImpute")$not_converged
  expect_gt(nrow(nc[[1]]), 0L)
  expect_identical(nrow(nc[[2]]), 0L)
})
test_that("an empty not_converged is a typed empty data.frame, never NULL", {
  out <- imp(dh, seed = 1, workers = 1, tol = c(1e-9, 5), maxIter = c(1, 40))
  e <- attr(out, "mvImpute")$not_converged[[2]]
  expect_s3_class(e, "data.frame")
  expect_identical(e, data.frame(row = integer(), column = character()))
})

test_that("a constant column with gaps is filled with a message and never stalls a stage", {
  d <- dh; d$const[1:3] <- NA
  expect_message(out <- imp(d, seed = 1, workers = 1), "Column 'const' is constant; its missing values were filled with that value.", fixed = TRUE)
  expect_true(all(out$const == 5))
  nc <- attr(out, "mvImpute")$not_converged[[2]]
  expect_equal(nrow(nc[nc$column == "const", ]), 0)
})
test_that("verbose controls progress messages", {
  expect_no_message(imp(dh, seed = 1, workers = 1, maxIter = c(2, 2)))
  msgs <- capture_messages(imp(dh, seed = 1, workers = 1, maxIter = c(2, 2), verbose = TRUE))
  expect_gte(sum(grepl("^Stage [12] of 2", msgs)), 2L)
  expect_true(any(grepl("Iteration", msgs)))
})
test_that("no missing values: input returned with a message and the diagnostics", {
  expect_message(out <- mvImpute(d0, seed = 1, workers = 1), "no missing")
  expect_identical(strip(out), d0)
  expect_identical(attr(out, "mvImpute")$seed, 1L)
  expect_message(mvImpute(d0, seed = 1, workers = 1), "X has no missing values in the columns to impute; returned unchanged.", fixed = TRUE)
})

test_that("no-missing return: diagnostics are empty, one draw with seed = NULL, gaps only in exclude", {
  dg <- attr(suppressMessages(mvImpute(d0, seed = 1, workers = 1)), "mvImpute")
  expect_true(is.na(dg$workers)); expect_identical(dg$iterations, integer())
  expect_identical(dg$not_converged, list()); expect_identical(dg$time, numeric())
  # seed = NULL costs exactly the draw resolve_seed() makes, and records it.
  withr::local_preserve_seed()
  set.seed(3); drawn <- sample.int(.Machine$integer.max, 1L); s_one <- .Random.seed
  set.seed(3)
  out <- suppressMessages(mvImpute(d0, workers = 1))
  expect_identical(.Random.seed, s_one)
  expect_identical(attr(out, "mvImpute")$seed, drawn)
  # Gaps only in an excluded column: nothing to impute, the gaps come back as they were.
  d <- d0; d$y <- c(NA, rnorm(nrow(d0) - 1L))
  expect_message(o2 <- mvImpute(d, exclude = "y", seed = 1, workers = 1),
                 "no missing values in the columns to impute")
  expect_identical(strip(o2), d)
})
test_that("workers is validated even when there is nothing to impute", {
  expect_error(mvImpute(d0, workers = "x"), "workers must be NULL or a whole number")
  expect_error(mvImpute(dh, workers = 3e9), "workers must be NULL or a whole number")
})
test_that("an error after the pool started restores RNG state, kind and plan", {
  skip_on_cran()
  X <- punch_holes(sim_numeric(40, 5, seed = 1), 0.1, "mcar", seed = 2)
  g <- ifelse(is.na(X), 1, X); g[which(is.na(X))[1]] <- NA   # first_guess() fails, after local_workers(2)
  withr::local_preserve_seed()
  old <- future::plan(future::sequential); withr::defer(future::plan(old))
  set.seed(3); sample.int(.Machine$integer.max, 1L); s_one <- .Random.seed   # state after exactly one draw
  set.seed(3); kind0 <- RNGkind()
  expect_error(mvImpute(X, guess = g, workers = 2, maxIter = c(2, 2)), "guess is NA")
  expect_identical(.Random.seed, s_one)
  expect_identical(RNGkind(), kind0)
  expect_true(inherits(future::plan(), "sequential"))
  expect_equal(future::nbrOfWorkers(), 1)
})

test_that("bad arguments are refused before any work", {
  expect_error(mvImpute(dh, mtry = NULL), "mtry must be a single positive whole number", fixed = TRUE)
  expect_error(mvImpute(dh, mtry = 0), "mtry must be a single positive whole number", fixed = TRUE)
  expect_error(mvImpute(dh, forceZero = NA), "forceZero")
  expect_error(mvImpute(dh, method = "nope"), "should be one of")
  expect_error(mvImpute(dh, ntree = 0), "ntree")
  expect_error(mvImpute(dh, seed = 3e9), "seed")
})

test_that("options and environment are untouched by a call", {
  skip_on_cran()
  imp(dh, seed = 1, workers = 2, maxIter = c(2, 2))   # warm-up: first call may load namespaces
  o0 <- options(); e0 <- Sys.getenv()
  imp(dh, seed = 1, workers = 2, maxIter = c(2, 2))
  imp(dh, seed = 1, workers = 1, maxIter = c(2, 2))
  expect_identical(options(), o0)
  expect_identical(Sys.getenv(), e0)
})

test_that("a user's multisession plan is used with workers = NULL and left in place", {
  skip_on_cran()
  old <- future::plan(future::multisession, workers = 2); withr::defer(future::plan(old))
  before <- pool_pids()
  out <- imp(dh, seed = 1, workers = NULL, maxIter = c(2, 2))
  expect_identical(attr(out, "mvImpute")$workers, 2L)
  expect_true(inherits(future::plan(), "multisession"))
  expect_identical(pool_pids(), before)
  expect_identical(strip(out), strip(imp(dh, seed = 1, workers = 1, maxIter = c(2, 2))))
})
test_that("workers = 1 restores the caller's plan afterwards", {
  skip_on_cran()
  old <- future::plan(future::multisession, workers = 2); withr::defer(future::plan(old))
  imp(dh, seed = 1, workers = 1, maxIter = c(2, 2))
  expect_true(inherits(future::plan(), "multisession"))
  expect_equal(future::nbrOfWorkers(), 2)
})
