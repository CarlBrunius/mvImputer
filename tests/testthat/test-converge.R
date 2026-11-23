test_that("a negative cell that doubles is not converged (old rule said it was)", {
  expect_false(converge_cells(-1, -2, FALSE, 1, FALSE, 0.05)$converged)
})
test_that("a cell at zero can converge", {
  expect_true(converge_cells(0, 0.001, FALSE, 1, FALSE, 0.05)$converged)
})
test_that("a change exactly at tol is not converged", {
  expect_false(converge_cells(0, 0.5, FALSE, 10, FALSE, 0.05)$converged)
})
test_that("the rule is invariant to shifting the data", {
  old <- c(-3, 0, 2, 7); new <- old + c(0.01, -0.2, 0.03, 0.5); s <- rep(2, 4)
  expect_identical(converge_cells(old, new, rep(FALSE, 4), s, rep(FALSE, 4), 0.05)$converged,
                   converge_cells(old + 1000, new + 1000, rep(FALSE, 4), s, rep(FALSE, 4), 0.05)$converged)
})
test_that("a factor cell needs two consecutive unchanged predictions", {
  expect_false(converge_cells(2, 2, FALSE, NA, TRUE, 0.05)$converged)
  expect_true(converge_cells(2, 2, TRUE, NA, TRUE, 0.05)$converged)
  expect_false(converge_cells(2, 3, TRUE, NA, TRUE, 0.05)$converged)
})
test_that("non-finite predictions are an error, never silently converged", {
  expect_error(converge_cells(1, NaN, FALSE, 1, FALSE, 0.05), "non-finite prediction")
})
test_that("col_scales uses observed values only", {
  W <- cbind(a = c(1, 2, 3, 100), b = c(5, 5, 5, NA), f = c(1, 2, 1, 2))
  miss <- cbind(a = c(FALSE, FALSE, FALSE, TRUE), b = c(FALSE, FALSE, FALSE, TRUE), f = FALSE)
  expect_equal(col_scales(W, miss, c("numeric", "numeric", "nominal")), c(1, 0, NA))
})

# Additions beyond the brief: mixed vectors, and the zero-scale / NA-scale edges.
test_that("mixed numeric and factor cells are judged per type, never NA", {
  res <- converge_cells(old = c(1, 2, 3), new = c(1.001, 2, 4),
                        prev_same = c(FALSE, TRUE, TRUE),
                        scale = c(1, NA, 1), is_factor = c(FALSE, TRUE, FALSE), tol = 0.05)
  expect_identical(res$converged, c(TRUE, TRUE, FALSE))
  expect_identical(res$same, c(FALSE, TRUE, FALSE))
  expect_false(anyNA(res$converged))
})
test_that("a zero scale converges only when nothing moved", {
  expect_true(converge_cells(5, 5, FALSE, 0, FALSE, 0.05)$converged)
  expect_false(converge_cells(5, 5.1, FALSE, 0, FALSE, 0.05)$converged)
})
test_that("a numeric cell with an NA scale is an error, not a silent NA", {
  expect_error(converge_cells(1, 2, FALSE, NA, FALSE, 0.05), "missing scale")
})
