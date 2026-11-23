test_that("punch_holes mcar hits the requested fraction and leaves other cells alone", {
  X <- sim_numeric(100, 20, seed = 1)
  Y <- punch_holes(X, 0.1, "mcar", seed = 2)
  expect_equal(mean(is.na(Y)), 0.1, tolerance = 0.01)
  expect_identical(Y[!is.na(Y)], X[!is.na(Y)])
})
test_that("punch_holes left removes the lowest values of each column", {
  X <- sim_numeric(100, 5, seed = 1)
  Y <- punch_holes(X, 0.2, "left", seed = 2)
  for (j in 1:5) expect_lt(max(X[is.na(Y[, j]), j]), min(X[!is.na(Y[, j]), j]))
})
test_that("sim_mixed has the documented columns and classes", {
  d <- sim_mixed(300, seed = 1)
  expect_identical(vapply(d, function(x) class(x)[1], ""),
    c(num1 = "numeric", num2 = "numeric", num3 = "numeric", num4 = "numeric", int1 = "integer",
      fac3 = "factor", ord4 = "ordered", chr2 = "character", lgl = "logical", const = "numeric"))
  expect_false(anyNA(d))
})
test_that("helpers are reproducible and do not touch the caller's RNG", {
  withr::local_seed(9); s <- .Random.seed
  expect_identical(sim_mixed(50, seed = 3), sim_mixed(50, seed = 3))
  expect_identical(s, .Random.seed)
})
test_that("helpers give the same data whatever RNG kind the caller has active", {
  X0 <- sim_numeric(30, 4, seed = 3)
  D0 <- sim_mixed(50, seed = 3)
  H0 <- punch_holes(D0, 0.2, "mcar", seed = 5)
  kind_default <- RNGkind()
  kind_other   <- c("L'Ecuyer-CMRG", "Box-Muller", "Rounding")
  # Cleanups run last-in-first-out: kind is reset first, then the old seed is put back.
  withr::local_preserve_seed()
  withr::defer(suppressWarnings(RNGkind(kind_default[1], kind_default[2], kind_default[3])))
  suppressWarnings(RNGkind(kind_other[1], kind_other[2], kind_other[3]))
  set.seed(9); s <- .Random.seed
  expect_identical(sim_numeric(30, 4, seed = 3), X0)
  expect_identical(sim_mixed(50, seed = 3), D0)
  expect_identical(punch_holes(D0, 0.2, "mcar", seed = 5), H0)
  expect_identical(RNGkind(), kind_other)      # the caller's kind is restored ...
  expect_identical(s, .Random.seed)            # ... and so is the caller's seed
  # A caller with no .Random.seed at all must keep its kind too.
  rm(".Random.seed", envir = globalenv())
  expect_identical(sim_mixed(50, seed = 3), D0)
  expect_identical(RNGkind(), kind_other)
  expect_false(exists(".Random.seed", envir = globalenv()))
})
