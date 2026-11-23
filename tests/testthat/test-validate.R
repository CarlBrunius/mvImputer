# test-validate.R -- check_input() and rows_without_observed() (spec section 5.1).
# One test per row of the input-check table, plus the ruled additions (nComp, ntree,
# mtry, guess) and the helper. Inputs come from sim_mixed(40, seed = 1) (a data
# frame with no gaps, 10 columns) or small literal matrices.

# ci(): check_input() with valid defaults, so each test changes exactly one thing.
ci <- function(X, exclude = NULL, tol = c(0.05, 0.025), maxIter = c(15L, 60L),
               method = "ranger", nComp = 2L, ...) {
  check_input(X, exclude, tol, maxIter, method, nComp, ...)
}

d40 <- function() sim_mixed(40, seed = 1)

test_that("X that is neither data.frame nor matrix is an error", {
  expect_error(ci(list(a = 1:5, b = 1:5)), "data.frame or matrix")
  expect_error(ci(1:10), "data.frame or matrix")
})

test_that("a non-numeric matrix is an error", {
  M <- matrix(c("a", "b", "c", "d", "e", "f"), 3, 2)
  expect_error(ci(M), "numeric matrix")
  expect_error(ci(matrix(c(TRUE, FALSE, TRUE, TRUE), 2, 2)), "numeric matrix")
})

test_that("exclude names not in colnames(X) are an error, also for an unnamed matrix", {
  expect_error(ci(d40(), exclude = "nope"), "unknown column")
  expect_error(ci(matrix(rnorm(12), 4, 3), exclude = "V1"), "unknown column")
  expect_error(ci(d40(), exclude = 1), "character vector")   # not silently coerced
})

test_that("fewer than 2 non-excluded columns is an error", {
  d <- d40()
  expect_error(ci(d[, 1, drop = FALSE]), "at least 2 columns")
  expect_error(ci(d, exclude = names(d)[-(1:1)]), "at least 2 columns")
  expect_error(ci(d, exclude = names(d)), "at least 2 columns")
})

test_that("an unsupported column type is an error naming the column", {
  d <- d40()
  d$when <- as.Date("2020-01-01") + 0:39
  expect_error(ci(d), "unsupported column type")
  expect_error(ci(d), "when")
  expect_silent(ci(d, exclude = "when"))        # excluded columns are not checked
})

test_that("Inf and -Inf are an error that says how to fix it", {
  d <- d40(); d$num1[3] <- Inf
  expect_error(ci(d), "Inf")
  expect_error(ci(d), "set them to NA")
  d <- d40(); d$num2[5] <- -Inf
  expect_error(ci(d), "set them to NA")
  d <- d40(); d$num1[3] <- NaN                  # NaN is missing, not an error
  expect_silent(ci(d))
})

test_that("a column with gaps and fewer than 2 observed values is an error", {
  d <- d40(); d$num1[-7] <- NA
  expect_error(ci(d), "fewer than 2 observed")
  d <- d40(); d$num1[] <- NA                    # all missing
  expect_error(ci(d), "fewer than 2 observed")
  d <- d40(); d$num1[-(7:8)] <- NA              # exactly 2 observed is allowed
  expect_silent(ci(d))
  M <- matrix(rnorm(12), 4, 3); M[-1, 2] <- NA  # unnamed matrix: label V2
  expect_error(ci(M), "V2")
})

test_that("a character or factor column with all-distinct observed values is an identifier", {
  d <- d40(); d$id <- as.character(1:40)
  expect_error(ci(d), "identifier")
  expect_error(ci(d), "exclude")
  d <- d40(); d$id <- factor(1:40)
  expect_error(ci(d), "identifier")
  d <- d40(); d$id <- as.character(1:40); d$id[c(2, 9)] <- NA   # gaps do not hide it
  expect_error(ci(d), "identifier")
  d <- d40(); d$few <- c("a", "b", NA, rep("c", 37))             # repeats: fine
  expect_silent(ci(d))
  d <- d40(); d$few <- c("a", "b", rep(NA, 38))                  # < 3 observed
  expect_silent(ci(d))
})

test_that("length(tol) != length(maxIter) is an error", {
  expect_error(ci(d40(), tol = c(0.05, 0.025, 0.01)), "one value per stage")
  expect_error(ci(d40(), tol = numeric(0), maxIter = integer(0)), "one value per stage")
})

test_that("tol <= 0 and maxIter < 1 are errors, as are missing or fractional values", {
  expect_error(ci(d40(), tol = c(0.05, 0)), "tol")
  expect_error(ci(d40(), tol = c(-1, 0.05)), "tol")
  expect_error(ci(d40(), tol = c(0.05, NA)), "tol")
  expect_error(ci(d40(), maxIter = c(15L, 0L)), "maxIter")
  expect_error(ci(d40(), maxIter = c(15, NA)), "maxIter")
  expect_error(ci(d40(), maxIter = c(15, 2.5)), "maxIter")
  expect_error(ci(d40(), maxIter = c(15, Inf)), "maxIter")
})

test_that("method pls without mixOmics is an error", {
  testthat::local_mocked_bindings(has_pkg = function(pkg) FALSE)
  expect_error(ci(d40(), method = "pls"), "mixOmics")
  expect_silent(ci(d40(), method = "ranger"))   # the package is only needed for pls
})

test_that("a row with no observed value is a warning", {
  d <- d40(); d[3, ] <- NA
  expect_warning(ci(d), "no observed value")
  M <- matrix(rnorm(20), 5, 4); M[2, ] <- NA
  expect_warning(ci(M), "no observed value")
})

test_that("numeric 0/1 columns forming a one-hot block are a warning", {
  d <- d40(); g <- rep(1:3, length.out = 40)
  d$ha <- as.numeric(g == 1); d$hb <- as.numeric(g == 2); d$hc <- as.numeric(g == 3)
  expect_warning(ci(d), "one-hot.*factor")
  d$hb[4] <- NA                                  # a gap in one column: still detected
  expect_warning(ci(d), "one-hot.*factor")
  d$hc[g == 3] <- 0                              # now rows of group 3 sum to 0
  expect_silent(ci(d))
})

test_that("independent 0/1 columns and a constant 0/1 column do not trigger the one-hot warning", {
  d <- d40()
  d$b1 <- rep(c(0, 1), 20); d$b2 <- rep(c(0, 0, 1, 1), 10)    # co-occur in some rows
  expect_silent(ci(d))
  d$zero <- 0                                    # constant columns never join a block
  expect_silent(ci(d))
})

test_that("a valid mixed table passes silently", {
  expect_silent(check_input(punch_holes(sim_mixed(100, 1), 0.1, "mcar", 2, cols = 1:9),
                            NULL, c(0.05, 0.025), c(15L, 60L), "ranger", 2L))
})

test_that("an excluded identifier column passes", {
  d <- d40(); d$id <- as.character(1:40)
  expect_silent(ci(d, exclude = "id"))
  d$id[5] <- NA                                  # excluded gaps are not checked either
  expect_silent(ci(d, exclude = "id"))
})

test_that("nComp, ntree and mtry must be positive whole numbers", {
  expect_error(ci(d40(), nComp = 0), "nComp")
  expect_error(ci(d40(), nComp = 1.5), "nComp")
  expect_error(ci(d40(), nComp = NA), "nComp")
  expect_error(ci(d40(), nComp = c(1, 2)), "nComp")
  expect_error(ci(d40(), ntree = 0), "ntree")
  expect_error(ci(d40(), ntree = "100"), "ntree")
  expect_error(ci(d40(), mtry = -2), "mtry")
  expect_error(ci(d40(), mtry = 2.5), "mtry")
  expect_silent(ci(d40(), ntree = 100, mtry = 5, nComp = 2))
  expect_silent(ci(d40(), mtry = NULL))          # NULL means not supplied
})

test_that("guess is one of the three strings or a table of X's dimensions", {
  d <- d40()
  expect_error(ci(d, guess = "mean"), "guess")
  expect_error(ci(d, guess = c("median", "minVar")), "guess")
  expect_error(ci(d, guess = NA_character_), "guess")
  expect_error(ci(d, guess = 3), "guess")
  expect_error(ci(d, guess = matrix(0, 3, 3)), "guess")
  expect_error(ci(d, guess = d[-1, ]), "guess")
  for (g in c("minVar", "median", "minTotal")) expect_silent(ci(d, guess = g))
  expect_silent(ci(d, guess = d))                # data frame of X's shape
  M <- matrix(rnorm(40), 10, 4)
  expect_silent(ci(M, guess = M))
})

test_that("rows_without_observed() lists rows empty in the non-excluded columns", {
  d <- d40(); d[c(3, 8), ] <- NA
  expect_identical(rows_without_observed(d), c(3L, 8L))
  expect_identical(rows_without_observed(d, exclude = "num1"), c(3L, 8L))
  d2 <- d40(); d2[4, names(d2) != "num1"] <- NA  # only the excluded column is observed
  expect_identical(rows_without_observed(d2), integer(0))
  expect_identical(rows_without_observed(d2, exclude = "num1"), 4L)
  M <- matrix(c(1, NA, 3, NA, NA, 6), 3, 2)
  expect_identical(rows_without_observed(M), 2L)
  expect_identical(rows_without_observed(sim_mixed(10, 1)), integer(0))
})
