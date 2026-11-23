for (m in c("ranger", "randomForest", "pls")) {
  test_that(paste(m, "predicts a numeric target"), {
    if (m == "pls") skip_if_not_installed("mixOmics")
    x <- withr::with_seed(1, matrix(rnorm(200 * 4), 200)); colnames(x) <- paste0("x", 1:4)
    y <- x[, 1] + withr::with_seed(2, rnorm(200, sd = 0.2))
    r <- withr::with_seed(3, fit_predict(m, y[1:150], x[1:150, ], x[151:200, ], FALSE, NA, sd(y), 100, 2, 2))
    expect_gt(cor(r$pred, y[151:200]), 0.9)
  })
  test_that(paste(m, "predicts a factor target as a valid code"), {
    if (m == "pls") skip_if_not_installed("mixOmics")
    x <- withr::with_seed(1, matrix(rnorm(300 * 3), 300)); colnames(x) <- paste0("x", 1:3)
    y <- max.col(x)                          # class = which column is largest
    r <- withr::with_seed(3, fit_predict(m, y[1:200], x[1:200, ], x[201:300, ], TRUE, 4L, NA, 100, 2, 2))
    expect_true(all(r$pred %in% 1:4)); expect_gt(mean(r$pred == y[201:300]), 0.8)
  })
}
test_that("a constant predictor is dropped without backend warnings", {
  skip_if_not_installed("mixOmics")
  x <- cbind(withr::with_seed(1, matrix(rnorm(60 * 3), 60)), k = 1); colnames(x)[1:3] <- paste0("x", 1:3)
  for (m in c("ranger", "randomForest", "pls"))
    expect_no_warning(withr::with_seed(2, fit_predict(m, x[1:50, 1], x[1:50, -1], x[51:60, -1], FALSE, NA, 1, 50, 2, 2)))
})
test_that("droplevels mapping: a level absent from training is never predicted", {
  x <- withr::with_seed(1, matrix(rnorm(100 * 2), 100)); colnames(x) <- c("x1", "x2")
  y <- ifelse(x[, 1] > 0, 1L, 3L)                    # level 2 of 3 never observed
  r <- withr::with_seed(2, fit_predict("randomForest", y[1:80], x[1:80, ], x[81:100, ], TRUE, 3L, NA, 50, 1, 2))
  expect_true(all(r$pred %in% c(1L, 3L)))
})
test_that("oob: near 0 for an exact function, near 1 for noise, NA for pls", {
  x <- withr::with_seed(1, matrix(rnorm(300 * 3), 300)); colnames(x) <- paste0("x", 1:3)
  exact <- withr::with_seed(2, fit_predict("ranger", 2 * x[, 1], x, x[1:2, ], FALSE, NA, sd(2 * x[, 1]), 200, 3, 2))
  noise <- withr::with_seed(2, fit_predict("ranger", rnorm(300), x, x[1:2, ], FALSE, NA, 1, 200, 3, 2))
  expect_lt(exact$oob, 0.3); expect_gt(noise$oob, 0.8)
  skip_if_not_installed("mixOmics")
  expect_true(is.na(fit_predict("pls", x[, 1], x[, -1], x[1:2, -1], FALSE, NA, 1, 10, 2, 2)$oob))
})
test_that("one prediction row, one predictor, and no surviving predictor all work", {
  x <- withr::with_seed(1, matrix(rnorm(50 * 2), 50)); colnames(x) <- c("x1", "x2"); y <- x[, 1]
  for (m in c("ranger", "randomForest"))
    expect_length(withr::with_seed(2, fit_predict(m, y[1:49], x[1:49, , drop = FALSE], x[50, , drop = FALSE], FALSE, NA, 1, 50, 2, 2))$pred, 1)
  expect_length(withr::with_seed(2, fit_predict("ranger", y, x[, 1, drop = FALSE], x[1:3, 1, drop = FALSE], FALSE, NA, 1, 50, 2, 2))$pred, 3)
  k <- cbind(k = rep(1, 50))
  r <- fit_predict("ranger", y, k, k[1:2, , drop = FALSE], FALSE, NA, 1, 50, 2, 2)
  expect_equal(r$pred, rep(mean(y), 2)); expect_true(is.na(r$oob))
})
test_that("randomForest regression on a 0/1 numeric target gives no warning", {
  x <- withr::with_seed(1, matrix(rnorm(80 * 2), 80)); colnames(x) <- c("x1", "x2")
  expect_no_warning(withr::with_seed(2, fit_predict("randomForest", as.numeric(x[, 1] > 0), x, x[1:3, ], FALSE, NA, 0.5, 50, 1, 2)))
})
test_that("same RNG state gives identical ranger and randomForest output", {
  x <- withr::with_seed(1, matrix(rnorm(100 * 3), 100)); colnames(x) <- paste0("x", 1:3); y <- x[, 1]
  for (m in c("ranger", "randomForest"))
    expect_identical(withr::with_seed(5, fit_predict(m, y, x, x[1:5, ], FALSE, NA, 1, 50, 2, 2)),
                     withr::with_seed(5, fit_predict(m, y, x, x[1:5, ], FALSE, NA, 1, 50, 2, 2)))
})

# --- additions (pre-flight rulings T6a, T6b, T6c and the audit floor) -----------------

test_that("a factor target with exactly one prediction row returns one valid code", {
  x <- withr::with_seed(1, matrix(rnorm(120 * 3), 120)); colnames(x) <- paste0("x", 1:3)
  y <- max.col(x)
  for (m in c("ranger", "randomForest", "pls")) {
    if (m == "pls") skip_if_not_installed("mixOmics")
    r <- withr::with_seed(2, fit_predict(m, y[1:119], x[1:119, ], x[120, , drop = FALSE], TRUE, 3L, NA, 50, 2, 2))
    expect_length(r$pred, 1); expect_true(r$pred %in% 1:3)
    expect_true(is.integer(r$pred))
  }
})

test_that("class argmax sends ties to the first level and draws no random numbers", {
  scores <- rbind(c(0.4, 0.4, 0.2), c(0.1, 0.45, 0.45), c(0.2, 0.3, 0.5), c(1, 1, 1))
  withr::with_seed(7, {
    before <- .Random.seed
    got <- argmax_first(scores)
    expect_identical(.Random.seed, before)
  })
  expect_identical(got, c(1L, 2L, 3L, 1L))
  expect_error(argmax_first(rbind(c(0.5, NaN))), "non-finite prediction")
  expect_error(argmax_first(rbind(c(0.5, NA))), "non-finite prediction")
})

test_that("pls consumes no random numbers; ranger and randomForest do", {
  x <- withr::with_seed(1, matrix(rnorm(60 * 3), 60)); colnames(x) <- paste0("x", 1:3); y <- x[, 1]
  for (m in c("ranger", "randomForest", "pls")) {
    if (m == "pls") skip_if_not_installed("mixOmics")
    withr::with_seed(9, {
      before <- .Random.seed
      fit_predict(m, y, x, x[1:2, ], FALSE, NA, 1, 20, 2, 2)
      if (m == "pls") expect_identical(.Random.seed, before)
      else expect_false(identical(.Random.seed, before))
    })
  }
})

test_that("predictions are unnamed even when the matrices carry row names", {
  x <- withr::with_seed(1, matrix(rnorm(60 * 3), 60, dimnames = list(paste0("r", 1:60), paste0("x", 1:3))))
  y <- x[, 1]; yc <- max.col(x)
  for (m in c("ranger", "randomForest", "pls")) {
    if (m == "pls") skip_if_not_installed("mixOmics")
    expect_null(names(withr::with_seed(2, fit_predict(m, y[1:50], x[1:50, ], x[51:60, ], FALSE, NA, 1, 20, 2, 2))$pred))
    expect_null(names(withr::with_seed(2, fit_predict(m, yc[1:50], x[1:50, ], x[51:60, ], TRUE, 3L, NA, 20, 2, 2))$pred))
  }
})

test_that("mtry and nComp above the number of predictors are capped, duplicated names are harmless", {
  x <- withr::with_seed(1, matrix(rnorm(60 * 2), 60)); colnames(x) <- c("a", "a"); y <- x[, 1] + x[, 2]
  for (m in c("ranger", "randomForest", "pls")) {
    if (m == "pls") skip_if_not_installed("mixOmics")
    r <- withr::with_seed(2, fit_predict(m, y[1:50], x[1:50, ], x[51:60, ], FALSE, NA, 1, 20, 5, 5))
    expect_length(r$pred, 10)
  }
})

test_that("a factor target with one observed training level takes the mode fallback", {
  x <- withr::with_seed(1, matrix(rnorm(30 * 2), 30)); colnames(x) <- c("x1", "x2")
  r <- fit_predict("ranger", rep(2L, 30), x, x[1:4, ], TRUE, 3L, NA, 20, 2, 2)
  expect_identical(r$pred, rep(2L, 4)); expect_true(is.na(r$oob))
})

test_that("factor oob is the misclassification rate: low when separable, high for noise, in [0, 1]", {
  x <- withr::with_seed(1, matrix(rnorm(300 * 3), 300)); colnames(x) <- paste0("x", 1:3)
  sep <- max.col(x)
  noise <- withr::with_seed(4, sample(1:3, 300, replace = TRUE))
  for (m in c("ranger", "randomForest")) {
    lo <- withr::with_seed(2, fit_predict(m, sep, x, x[1:2, ], TRUE, 3L, NA, 200, 3, 2))$oob
    hi <- withr::with_seed(2, fit_predict(m, noise, x, x[1:2, ], TRUE, 3L, NA, 200, 3, 2))$oob
    expect_true(lo >= 0 && lo < 0.25); expect_true(hi > 0.5 && hi <= 1)
  }
})

test_that("ranger prediction.error of a probability forest is a Brier-type score, not the misclassification rate", {
  # Evidence for the header note in R/backends.R and for computing the factor oob from
  # fit$predictions. If a ranger update changes this, the note must be revisited.
  x <- withr::with_seed(1, matrix(rnorm(300 * 3), 300)); colnames(x) <- paste0("x", 1:3)
  y <- factor(max.col(x))
  fit <- withr::with_seed(2, ranger::ranger(x = x, y = y, num.trees = 200, mtry = 2,
                                            probability = TRUE, num.threads = 1, seed = NULL))
  p <- fit$predictions
  p_true <- p[cbind(seq_along(y), as.integer(y))]
  misclass <- mean(max.col(p, ties.method = "first") != as.integer(y))
  expect_equal(fit$prediction.error, mean((1 - p_true)^2))
  expect_gt(abs(fit$prediction.error - misclass), 0.01)
})

test_that("non-finite predictor values are refused (blocking input check)", {
  x <- withr::with_seed(1, matrix(rnorm(40 * 2), 40)); colnames(x) <- c("x1", "x2")
  expect_error(fit_predict("ranger", x[, 1], x, rbind(c(NA_real_, 1)), FALSE, NA, 1, 20, 2, 2))
})
