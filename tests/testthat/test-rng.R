# test-rng.R -- seeds, stream generation and RNG isolation (spec section 5.5).
#
# These tests guard the safety-critical contract: the caller's .Random.seed (including
# its ABSENCE) and RNGkind() come back untouched, and streams are a pure function of
# (seed, n_stage, n_col). Simulated data only; no data files.

# Put the current RNG kinds back at test exit, so a failing test cannot leak a kind
# into later tests. Registered AFTER withr::local_preserve_seed() in the tests that
# remove .Random.seed, so it runs first (withr defers are last-in, first-out).
local_restore_kind <- function(envir = parent.frame()) {
  k <- RNGkind()
  withr::defer(suppressWarnings(RNGkind(k[1], k[2], k[3])), envir)
}

test_that("make_streams is deterministic, distinct and L'Ecuyer", {
  s <- make_streams(42L, 2L, 3L)
  expect_identical(s, make_streams(42L, 2L, 3L))
  expect_length(s, 2); expect_length(s[[1]], 3)
  flat <- unlist(lapply(s, function(x) vapply(x, paste, "", collapse = ",")))
  expect_false(anyDuplicated(flat) > 0)
  expect_identical(s[[1]][[1]][1], 10407L)
  expect_true(all(vapply(unlist(s, recursive = FALSE),
                         function(v) is.integer(v) && length(v) == 7L, TRUE)))
})

test_that("make_streams leaves the caller's RNG exactly as it was", {
  set.seed(1); s0 <- .Random.seed
  make_streams(42L, 2L, 3L)
  expect_identical(.Random.seed, s0); expect_identical(RNGkind()[1], "Mersenne-Twister")
})

test_that("make_streams restores a non-default caller kind with a seed present", {
  withr::local_preserve_seed()
  local_restore_kind()
  suppressWarnings(RNGkind("Wichmann-Hill", "Box-Muller", "Rounding"))
  set.seed(3); k0 <- RNGkind(); s0 <- .Random.seed
  make_streams(42L, 1L, 2L)
  expect_identical(RNGkind(), k0)
  expect_identical(.Random.seed, s0)
})

test_that("an absent .Random.seed stays absent", {
  withr::local_preserve_seed()
  local_restore_kind()
  # Start from the default kinds, whatever an earlier test left behind.
  suppressWarnings(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))
  if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv())
  k0 <- RNGkind()
  make_streams(1L, 1L, 1L)
  expect_false(exists(".Random.seed", globalenv()))
  # The case a withr-only guard gets wrong: the seed stays absent but the kind leaks.
  expect_identical(RNGkind(), k0)
})

test_that("an absent .Random.seed with a non-default kind keeps that kind", {
  withr::local_preserve_seed()
  local_restore_kind()
  suppressWarnings(RNGkind("Wichmann-Hill", "Box-Muller", "Rounding"))
  rm(".Random.seed", envir = globalenv())
  k0 <- RNGkind()
  make_streams(1L, 1L, 1L)
  expect_false(exists(".Random.seed", globalenv()))
  expect_identical(RNGkind(), k0)
})

test_that("local_rng_guard restores seed and kind on an error exit", {
  set.seed(5); s0 <- .Random.seed; k0 <- RNGkind()
  f <- function() {
    local_rng_guard()
    suppressWarnings(RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"))
    stats::runif(3)
    stop("boom")
  }
  expect_error(f(), "boom")
  expect_identical(.Random.seed, s0)
  expect_identical(RNGkind(), k0)
})

test_that("make_streams is a pure function of its arguments", {
  a <- make_streams(9L, 2L, 2L)
  withr::local_preserve_seed()
  local_restore_kind()
  suppressWarnings(RNGkind("Wichmann-Hill", "Box-Muller", "Rounding"))
  set.seed(99)
  expect_identical(make_streams(9L, 2L, 2L), a)
  expect_false(identical(make_streams(10L, 2L, 2L), a))
})

test_that("a stage-1 stream depends only on its position, not on how many columns follow", {
  small <- make_streams(5L, 2L, 2L)
  big   <- make_streams(5L, 2L, 3L)
  # Stage 1 only: stream k is the same whatever n_col is (keyed by position in W). Later
  # stages start n_col streams further along the chain, so they depend on n_col; harmless
  # because n_col = ncol(W) is fixed for a data set.
  expect_identical(small[[1]][[1]], big[[1]][[1]])
  expect_identical(small[[1]][[2]], big[[1]][[2]])
})

test_that("make_streams validates its arguments", {
  expect_error(make_streams(NULL, 1L, 1L), "seed")
  expect_error(make_streams(1.5, 1L, 1L), "seed")
  expect_error(make_streams(1L, 0L, 1L), "n_stage")
  expect_error(make_streams(1L, 1L, 0L), "n_col")
  expect_error(make_streams(1L, 1.5, 1L), "n_stage")
})

test_that("resolve_seed(NULL) follows set.seed and costs exactly one draw", {
  set.seed(7); a <- resolve_seed(NULL); s1 <- .Random.seed
  set.seed(7); b <- resolve_seed(NULL)
  expect_identical(a, b)
  expect_true(is.integer(a) && length(a) == 1L)
  set.seed(7); sample.int(.Machine$integer.max, 1L)
  expect_identical(s1, .Random.seed)
})

test_that("resolve_seed returns a valid seed as integer and rejects the rest", {
  expect_identical(resolve_seed(42), 42L)
  expect_identical(resolve_seed(-3L), -3L)
  for (bad in list("a", c(1, 2), NA, NA_integer_, 1.5, Inf, -Inf, NaN, 3e9, TRUE, integer(0))) {
    expect_error(resolve_seed(bad), "seed must be NULL or a single finite whole number",
                 info = deparse(bad))
  }
})

test_that("crn = FALSE gives a different stream each iteration, crn = TRUE the same", {
  s <- make_streams(1L, 1L, 1L)[[1]][[1]]
  expect_identical(stream_for_iteration(s, 3L, TRUE), s)
  expect_false(identical(stream_for_iteration(s, 3L, FALSE), stream_for_iteration(s, 2L, FALSE)))
})

test_that("stream_for_iteration walks sub-streams from the first iteration", {
  s <- make_streams(1L, 1L, 1L)[[1]][[1]]
  expect_identical(stream_for_iteration(s, 1L, FALSE), s)
  expect_identical(stream_for_iteration(s, 2L, FALSE), parallel::nextRNGSubStream(s))
  expect_identical(stream_for_iteration(s, 3L, FALSE),
                   parallel::nextRNGSubStream(parallel::nextRNGSubStream(s)))
  expect_error(stream_for_iteration(s, 0L, FALSE), "iter")
})
