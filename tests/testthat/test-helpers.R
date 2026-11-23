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
  set.seed(9); s <- .Random.seed
  expect_identical(sim_mixed(50, seed = 3), sim_mixed(50, seed = 3))
  expect_identical(s, .Random.seed)
})
