# test-workers.R -- worker plan handling (spec section 5.1). Never more than 2 workers.
#
# The multisession workers only run base functions here (Sys.getpid), so they do not
# need an installed mvImputer. Each test that sets a plan resets it to sequential at
# exit, so one failing test cannot leave workers running for the next.

# Sorted PIDs of both workers of a 2-worker plan. The sleep keeps the first future busy,
# so the second is scheduled on the other worker: the pair identifies the worker pool.
worker_pids <- function() {
  fs <- lapply(1:2, function(i) future::future({ Sys.sleep(0.5); Sys.getpid() }))
  sort(unlist(lapply(fs, future::value)))
}

test_that("workers = 2 sets multisession for the call and restores the old plan", {
  skip_on_cran()
  future::plan(future::sequential); withr::defer(future::plan(future::sequential))
  f <- function() { n <- local_workers(2L); c(n, future::nbrOfWorkers()) }
  expect_equal(f(), c(2, 2))
  expect_true(inherits(future::plan(), "sequential"))
})

test_that("an active non-sequential plan is used and left in place", {
  skip_on_cran()
  future::plan(future::multisession, workers = 2); withr::defer(future::plan(future::sequential))
  before <- worker_pids()
  f <- function() local_workers(NULL)
  expect_equal(f(), 2)
  expect_true(inherits(future::plan(), "multisession"))
  # Same worker processes: the plan was not replaced by a new one (or by a restored copy).
  expect_identical(worker_pids(), before)
})

test_that("workers = 1 runs sequentially and restores the previous plan", {
  skip_on_cran()
  future::plan(future::multisession, workers = 2); withr::defer(future::plan(future::sequential))
  f <- function() { n <- local_workers(1L); c(n, future::nbrOfWorkers(),
                                              inherits(future::plan(), "sequential")) }
  expect_equal(f(), c(1, 1, 1))
  expect_true(inherits(future::plan(), "multisession"))
  expect_equal(future::nbrOfWorkers(), 2)
})

test_that("workers = 1 from a sequential plan stays sequential", {
  future::plan(future::sequential); withr::defer(future::plan(future::sequential))
  f <- function() { local_workers(1L); future::nbrOfWorkers() }
  expect_equal(f(), 1)
  expect_true(inherits(future::plan(), "sequential"))
})

test_that("workers = NULL under a sequential plan sizes the pool from availableCores()", {
  skip_on_cran()
  future::plan(future::sequential); withr::defer(future::plan(future::sequential))
  f <- function() { n <- local_workers(NULL); c(n, future::nbrOfWorkers()) }
  # 3 cores -> 2 workers (one core left for the main process).
  testthat::local_mocked_bindings(availableCores = function(...) 3L, .package = "future")
  expect_equal(f(), c(2, 2))
  expect_true(inherits(future::plan(), "sequential"))
})

test_that("workers = NULL with one available core stays sequential (nested call)", {
  future::plan(future::sequential); withr::defer(future::plan(future::sequential))
  f <- function() { n <- local_workers(NULL); c(n, future::nbrOfWorkers()) }
  testthat::local_mocked_bindings(availableCores = function(...) 1L, .package = "future")
  expect_equal(f(), c(1, 1))
  expect_true(inherits(future::plan(), "sequential"))
})

test_that("an error exit restores the plan, .Random.seed and RNGkind()", {
  skip_on_cran()
  future::plan(future::sequential); withr::defer(future::plan(future::sequential))
  set.seed(11); s0 <- .Random.seed; k0 <- RNGkind()
  f <- function() {
    local_rng_guard()
    local_workers(2L)
    suppressWarnings(RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection"))
    stats::runif(2)
    stop("boom")
  }
  expect_error(f(), "boom")
  expect_true(inherits(future::plan(), "sequential"))
  expect_identical(.Random.seed, s0)
  expect_identical(RNGkind(), k0)
})

test_that("local_workers rejects anything but NULL or a whole number >= 1", {
  # 3e9 is a whole number but not an integer: as.integer() would give NA with a warning.
  for (bad in list(0, -1, 1.5, NA, NA_integer_, Inf, "2", c(1, 2), TRUE, 3e9)) {
    expect_error(local_workers(bad), "workers must be NULL or a whole number >= 1",
                 info = deparse(bad))
  }
})
