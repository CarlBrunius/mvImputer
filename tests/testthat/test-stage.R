# test-stage.R -- stage_jobs() and run_stage(): training rows, frozen cells, seeding by column,
# worker invariance (spec sections 5.4, 5.5). The multisession tests load the INSTALLED
# package in their workers, so run this file with load_package = 'installed'.

test_that("training rows are the originally observed rows, even after freezing", {
  miss <- cbind(a = c(TRUE, TRUE, FALSE, FALSE), b = c(FALSE, TRUE, FALSE, FALSE))
  frozen <- cbind(a = c(TRUE, FALSE, FALSE, FALSE), b = FALSE)
  jobs <- stage_jobs(miss, frozen, targets = 1:2)
  expect_identical(jobs[[1]]$train, 3:4); expect_identical(jobs[[1]]$pred, 2L)
  expect_identical(jobs[[2]]$train, c(1L, 3L, 4L))
})
test_that("observed cells are never changed by a stage", {
  X <- punch_holes(sim_numeric(80, 10, seed = 1), 0.1, "mcar", seed = 2)
  enc <- encode_input(X, NULL); W0 <- first_guess(enc$W, enc, "minVar")
  st <- make_streams(1L, 1L, ncol(W0))[[1]]
  out <- run_stage(W0, enc, col_scales(enc$W, enc$miss, enc$type), st, 0.05, 5L, "ranger", 50, 3, 2, FALSE)
  expect_identical(out$W[!enc$miss], W0[!enc$miss])
})
test_that("a stage gives identical results with 1 and 2 workers", {
  X <- punch_holes(sim_numeric(60, 8, seed = 1), 0.1, "mcar", seed = 2)
  enc <- encode_input(X, NULL); W0 <- first_guess(enc$W, enc, "minVar")
  st <- make_streams(11L, 1L, ncol(W0))[[1]]; s <- col_scales(enc$W, enc$miss, enc$type)
  run <- function(k) {
    old <- future::plan(future::multisession, workers = k); withr::defer(future::plan(old))
    run_stage(W0, enc, s, st, 0.05, 4L, "ranger", 50, 3, 2, FALSE)
  }
  expect_identical(run(1), run(2))
})
test_that("factors with informative predictors converge without oscillating", {
  d <- sim_mixed(300, seed = 4); truth <- d$fac3
  d$fac3[withr::with_seed(5, sample(300, 45))] <- NA
  enc <- encode_input(d, NULL); W0 <- first_guess(enc$W, enc, "minVar")
  st <- make_streams(2L, 1L, ncol(W0))[[1]]
  out <- run_stage(W0, enc, col_scales(enc$W, enc$miss, enc$type), st, 0.025, 20L, "ranger", 100, 5, 2, FALSE)
  # Accuracy is the real check: with sibling indicators as predictors it fell to ~0.58 (spec F3).
  j <- which(colnames(enc$W) == "fac3"); m <- enc$miss[, j]
  expect_gt(mean(enc$levels$fac3[out$W[m, j]] == truth[m]), 0.75)
  # A 2-cycle would run to the cap; converging means stopping below it.
  expect_lt(out$iterations, 20L)
  expect_identical(nrow(out$not_converged), 0L)
})
test_that("oob is reported per imputed column", {
  X <- punch_holes(sim_numeric(80, 6, seed = 1), 0.1, "mcar", seed = 2)
  enc <- encode_input(X, NULL); W0 <- first_guess(enc$W, enc, "minVar")
  out <- run_stage(W0, enc, col_scales(enc$W, enc$miss, enc$type), make_streams(1L, 1L, 6L)[[1]], 0.05, 3L, "ranger", 50, 3, 2, FALSE)
  expect_setequal(names(out$oob), colnames(X)[colSums(is.na(X)) > 0])
})

# --- added beyond the plan's five tests (task rulings T7c and the 4-worker check) ---

test_that("stage_jobs skips non-targets and columns with no unfrozen missing cell", {
  miss   <- cbind(a = c(TRUE, FALSE), b = c(TRUE, FALSE), c = c(FALSE, FALSE), d = c(TRUE, FALSE))
  frozen <- cbind(a = c(TRUE, FALSE), b = FALSE, c = FALSE, d = FALSE)
  jobs <- stage_jobs(miss, frozen, targets = c(1L, 2L, 3L))   # d is not a target
  expect_identical(vapply(jobs, `[[`, 1L, "col"), 2L)
  expect_identical(jobs[[1]]$pred, 1L)
})

test_that("seeds are keyed by W column, not by position in the job list", {
  local_rng_guard()
  X <- sim_numeric(30, 4, seed = 1)
  X[c(2, 5), 2] <- NA; X[c(3, 7, 9), 4] <- NA     # columns 1 and 3 complete: job 1 is column 2
  enc <- encode_input(X, NULL); W0 <- first_guess(enc$W, enc, "minVar")
  s <- col_scales(enc$W, enc$miss, enc$type)
  streams <- make_streams(3L, 2L, ncol(W0))
  # Each prediction is the first draws of whatever RNG state the job runs under. Mocked
  # bindings are invisible to multisession workers, so the plan is sequential here.
  local_mocked_bindings(fit_predict = function(method, y, x_train, x_pred, ...)
    list(pred = stats::runif(nrow(x_pred)), oob = NA_real_))
  old_plan <- future::plan(future::sequential); withr::defer(future::plan(old_plan))
  draw <- function(stream, n) {
    local_rng_guard()
    assign(".Random.seed", stream, envir = globalenv())
    stats::runif(n)
  }
  # tol tiny: no cell converges, so iteration 2 overwrites every cell.
  one <- run_stage(W0, enc, s, streams[[2]], 1e-12, 1L, "ranger", 10, 1, 1, FALSE, crn = TRUE)
  expect_identical(unname(one$W[enc$miss[, 2], 2]),
                   draw(stream_for_iteration(streams[[2]][[2]], 1L, TRUE), 2L))
  expect_identical(unname(one$W[enc$miss[, 4], 4]),
                   draw(stream_for_iteration(streams[[2]][[4]], 1L, TRUE), 3L))
  two <- run_stage(W0, enc, s, streams[[2]], 1e-12, 2L, "ranger", 10, 1, 1, FALSE, crn = FALSE)
  expect_identical(two$iterations, 2L)
  expect_identical(unname(two$W[enc$miss[, 4], 4]),
                   draw(stream_for_iteration(streams[[2]][[4]], 2L, FALSE), 3L))
})

test_that("a constant column never keeps a stage iterating; verbose reports each iteration", {
  X <- punch_holes(sim_numeric(40, 4, seed = 3), 0.1, "mcar", seed = 4)
  X[, 1] <- 5; X[c(1, 2), 1] <- NA
  enc <- encode_input(X, NULL); W0 <- first_guess(enc$W, enc, "minVar")
  s <- col_scales(enc$W, enc$miss, enc$type)          # column 1 has scale 0
  msgs <- capture_messages(
    out <- run_stage(W0, enc, s, make_streams(1L, 1L, 4L)[[1]], 1e-12, 2L, "ranger", 20, 2, 1, FALSE,
                     verbose = TRUE))
  expect_identical(out$iterations, 2L)
  expect_match(msgs, "^Iteration [12]: [0-9]+ cells still iterating\n$")
  expect_length(msgs, 2L)
  expect_identical(unname(out$W[1:2, 1]), c(5, 5))
  expect_false("F1" %in% names(out$oob))
  expect_false("F1" %in% out$not_converged$column)
  expect_named(out$not_converged, c("row", "column"))
  expect_type(out$not_converged$row, "integer")
  expect_null(names(out$not_converged$row))
})

test_that("an error in a fit or in the stopping rule names the column", {
  X <- sim_numeric(30, 4, seed = 1); X[c(2, 5), 3] <- NA          # only F3 is a target
  enc <- encode_input(X, NULL); W0 <- first_guess(enc$W, enc, "minVar")
  s <- col_scales(enc$W, enc$miss, enc$type); st <- make_streams(1L, 1L, 4L)[[1]]
  # Mocked bindings are invisible to multisession workers, so the plan is sequential here.
  old_plan <- future::plan(future::sequential); withr::defer(future::plan(old_plan))
  local_mocked_bindings(fit_predict = function(...) stop("non-finite prediction"))
  # future.apply adds its own "Canceling all iterations" warning to a job error; only that
  # one is muffled.
  expect_error(withCallingHandlers(
    run_stage(W0, enc, s, st, 0.05, 1L, "ranger", 10, 1, 1, FALSE),
    warning = function(w) if (grepl("Canceling all iterations", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")),
    "column 'F3': non-finite prediction", fixed = TRUE)
  local_mocked_bindings(fit_predict = function(method, y, x_train, x_pred, ...)
    list(pred = rep(NaN, nrow(x_pred)), oob = NA_real_))
  expect_error(run_stage(W0, enc, s, st, 0.05, 1L, "ranger", 10, 1, 1, FALSE),
               "column 'F3': non-finite prediction: 2 of 2 cells", fixed = TRUE)
})

test_that("a stage gives identical results with 1 and 4 workers", {
  skip_on_cran()
  X <- punch_holes(sim_numeric(60, 8, seed = 1), 0.1, "mcar", seed = 2)
  enc <- encode_input(X, NULL); W0 <- first_guess(enc$W, enc, "minVar")
  st <- make_streams(11L, 1L, ncol(W0))[[1]]; s <- col_scales(enc$W, enc$miss, enc$type)
  run <- function(k) {
    old <- future::plan(future::multisession, workers = k); withr::defer(future::plan(old))
    run_stage(W0, enc, s, st, 0.05, 4L, "ranger", 50, 3, 2, FALSE)
  }
  expect_identical(run(1), run(4))
})
