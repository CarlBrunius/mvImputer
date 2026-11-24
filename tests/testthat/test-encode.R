# test-encode.R -- encode_input(), design_matrix(), first_guess() and decode_output(): the
# internal data model and its lossless round trip (spec sections 5.1, 5.2). No workers.

test_that("encode then decode returns the input unchanged, gaps included", {
  d <- sim_mixed(200, seed = 1)
  d$int1 <- as.numeric(d$int1)              # integer-with-gaps returns double by design
  dh <- punch_holes(d, 0.1, "mcar", seed = 2, cols = 1:9)
  enc <- encode_input(dh, exclude = NULL)
  expect_identical(decode_output(enc$W, enc), dh)
})
test_that("an integer column without gaps stays integer", {
  d <- sim_mixed(50, seed = 1); d$num1[3] <- NA
  expect_type(decode_output(encode_input(d, NULL)$W, encode_input(d, NULL))$int1, "integer")
})
test_that("unused levels survive and are never imputed", {
  d <- sim_mixed(100, seed = 1); d$fac3 <- factor(d$fac3, levels = c("a", "b", "c", "zz"))
  d$fac3[1:5] <- NA
  enc <- encode_input(d, NULL); W <- first_guess(enc$W, enc, "minVar")
  out <- decode_output(W, enc)
  expect_identical(levels(out$fac3), c("a", "b", "c", "zz")); expect_false(any(out$fac3 == "zz"))
})
test_that("design matrix: one-hot for nominal, integer code for ordered, own columns mapped", {
  d <- sim_mixed(30, seed = 1); enc <- encode_input(d, NULL); dm <- design_matrix(enc$W, enc)
  expect_true(all(c("fac3_a", "fac3_b", "fac3_c", "chr2_x", "chr2_y", "lgl_FALSE", "lgl_TRUE") %in% colnames(dm$D)))
  expect_identical(colnames(dm$D)[dm$map[["fac3"]]], c("fac3_a", "fac3_b", "fac3_c"))
  expect_identical(colnames(dm$D)[dm$map[["ord4"]]], "ord4")
  expect_setequal(unique(dm$D[, "ord4"]), 1:4)
})
test_that("excluded columns are absent from W and D but returned untouched", {
  d <- sim_mixed(30, seed = 1); d$y <- c(NA, rnorm(29)); enc <- encode_input(d, exclude = "y")
  expect_false("y" %in% colnames(enc$W))
  expect_identical(decode_output(first_guess(enc$W, enc, "minVar"), enc)$y, d$y)
})
test_that("first guesses: minVar, median, minTotal, factor mode, user matrix at missing cells only", {
  m <- cbind(a = c(1, NA, 3, 10), b = c(-5, 2, NA, 4)); enc <- encode_input(m, NULL)
  expect_equal(first_guess(enc$W, enc, "minVar")[2, "a"], 1)
  expect_equal(first_guess(enc$W, enc, "median")[3, "b"], 2)
  expect_equal(first_guess(enc$W, enc, "minTotal")[2, "a"], -5)
  g <- m; g[] <- 99; expect_equal(first_guess(enc$W, enc, g)[1, "a"], 1)
  expect_error(first_guess(enc$W, enc, m), "guess")
})
test_that("unnamed matrix: dimnames stay absent, internal names are V1..Vp", {
  m <- matrix(c(1, NA, 3, 4, 5, 6), 3); enc <- encode_input(m, NULL)
  expect_identical(enc$colnames, c("V1", "V2"))
  expect_null(dimnames(decode_output(first_guess(enc$W, enc, "minVar"), enc)))
})
test_that("constant columns are detected", {
  expect_identical(encode_input(sim_mixed(30, seed = 1), NULL)$constant, "const")
})

# --- beyond the plan: data-model labels, levels, ties, boundaries ------------------

test_that("type and class labels are the ones the stage engine keys on", {
  enc <- encode_input(sim_mixed(30, seed = 1), NULL)
  expect_identical(enc$type, c(rep("numeric", 5), "nominal", "ordered", "nominal",
                               "nominal", "numeric"))
  expect_identical(enc$class, c(rep("numeric", 4), "integer", "factor", "ordered",
                                "character", "logical", "numeric"))
  expect_identical(enc$levels$ord4, c("L1", "L2", "L3", "L4"))
})

test_that("level sets: logical is always FALSE/TRUE, character is radix-sorted", {
  d <- data.frame(l = c(TRUE, TRUE, NA, TRUE, TRUE), s = c("b", "B", "a", NA, "b"),
                  v = c(1, 2, 3, 4, 5), stringsAsFactors = FALSE)
  enc <- encode_input(d, NULL)
  expect_identical(enc$levels$l, c("FALSE", "TRUE"))
  expect_identical(enc$levels$s, c("B", "a", "b"))       # C order: upper case first
  expect_identical(unname(enc$W[, "l"]), c(2, 2, NA, 2, 2))
  out <- decode_output(enc$W, enc)
  expect_identical(out$l, rep(TRUE, 5))                   # constant (only TRUE observed): filled
  expect_identical(out$s, d$s)
})

test_that("observed values come back bit-identical from first_guess + decode", {
  d <- sim_mixed(100, seed = 3)
  dh <- punch_holes(d, 0.2, "mcar", seed = 4, cols = 1:9)
  enc <- encode_input(dh, NULL)
  out <- decode_output(first_guess(enc$W, enc, "median"), enc)
  for (j in names(dh)) {
    obs <- !is.na(dh[[j]])
    if (j == "int1") expect_equal(as.numeric(out[[j]][obs]), as.numeric(dh[[j]][obs]))
    else expect_identical(out[[j]][obs], dh[[j]][obs])
  }
  expect_false(anyNA(out))
  expect_type(out$int1, "double")                         # integer with gaps -> double
})

test_that("factor first guess is the most frequent observed code, ties to the lowest", {
  d <- data.frame(f = factor(c("b", "a", "b", "a", NA), levels = c("a", "b", "c")),
                  x = c(1, 2, 3, 4, 5))
  enc <- encode_input(d, NULL)
  expect_identical(decode_output(first_guess(enc$W, enc, "minVar"), enc)$f[5],
                   factor("a", levels = c("a", "b", "c")))
  d$f[1:2] <- c("b", "b")                                # b is now the mode
  enc <- encode_input(d, NULL)
  expect_identical(as.character(decode_output(first_guess(enc$W, enc, "minVar"), enc)$f[5]), "b")
})

test_that("constant columns are filled with their constant, whatever guess says", {
  d <- sim_mixed(30, seed = 1); d$const[2] <- NA
  enc <- encode_input(d, NULL)
  g <- matrix(99, nrow(d), ncol(d))
  expect_identical(decode_output(enc$W, enc)$const[2], 5)
  expect_equal(first_guess(enc$W, enc, g)[2, "const"], 5)
})

test_that("a user guess is read only at missing cells of non-constant numeric columns", {
  d <- data.frame(a = c(1, NA, 3), b = c(NA, 2, 2), f = factor(c("u", NA, "u")),
                  y = c(7, 8, 9))
  enc <- encode_input(d, exclude = "y")
  g <- data.frame(a = c(NA, 50, NA), b = c(NA, NA, NA), f = c("junk", "junk", "junk"),
                  y = c(NA, NA, NA))
  W <- first_guess(enc$W, enc, g)
  expect_equal(W[2, "a"], 50)                            # read at its missing cell
  expect_equal(W[1, "b"], 2)                             # constant column ignores the table
  expect_equal(W[2, "f"], 1)                             # factor column ignores the table
  expect_error(first_guess(enc$W, enc, g[1:2, ]), "dimensions")
  g$a[2] <- Inf
  expect_error(first_guess(enc$W, enc, g), "Inf")
  g$a[2] <- NA
  expect_error(first_guess(enc$W, enc, g), "guess is NA")
})

test_that("guess table is located by position in X when a column is excluded", {
  m <- cbind(id = c(100, 200, 300), a = c(1, NA, 3), b = c(4, 5, 6))
  enc <- encode_input(m, exclude = "id")
  g <- m; g[] <- 0; g[2, "a"] <- 42                       # a is column 2 of X, column 1 of W
  expect_equal(first_guess(enc$W, enc, g)[2, "a"], 42)
})

test_that("exclude is invisible: changing the excluded column changes nothing in W", {
  d <- sim_mixed(30, seed = 1); d$y <- c(NA, rnorm(29))
  d2 <- d; d2$y <- c(5, rep(NA, 29))
  e1 <- encode_input(d, "y"); e2 <- encode_input(d2, "y")
  expect_identical(e1$W, e2$W); expect_identical(e1$miss, e2$miss)
  expect_identical(e1$constant, e2$constant)
  expect_false("y" %in% colnames(design_matrix(e1$W, e1)$D))
  expect_identical(decode_output(e2$W, e2)$y, d2$y)
})

test_that("unnamed matrix: messages call columns V<k>", {
  m <- matrix(c(1, NA, 3, 4, 5, 6), 3); enc <- encode_input(m, NULL)
  expect_error(first_guess(enc$W, enc, m), "V1")
  m2 <- matrix(c(1, 2, 3, 4, NA, 6), 3); enc2 <- encode_input(m2, NULL)
  expect_error(first_guess(enc2$W, enc2, m2), "V2")
})

test_that("an integer matrix stays integer without gaps and becomes double with gaps", {
  mi <- matrix(1:6, 3, dimnames = list(NULL, c("a", "b")))
  e <- encode_input(mi, NULL)
  expect_identical(decode_output(e$W, e), mi)
  mi[2, 1] <- NA
  e <- encode_input(mi, NULL)
  out <- decode_output(first_guess(e$W, e, "median"), e)
  expect_type(out, "double"); expect_equal(unname(out[2, "a"]), 2)
})

test_that("design_matrix refuses unfilled factor codes", {
  d <- sim_mixed(30, seed = 1); d$fac3[1] <- NA
  enc <- encode_input(d, NULL)
  expect_error(design_matrix(enc$W, enc), "fac3")
  expect_no_error(design_matrix(first_guess(enc$W, enc, "minVar"), enc))
})

test_that("design matrix keeps numeric NA before the first guess and has one block per column", {
  d <- sim_mixed(30, seed = 1); d$num1[2] <- NA
  enc <- encode_input(d, NULL); dm <- design_matrix(enc$W, enc)
  expect_true(is.na(dm$D[2, "num1"]))
  expect_identical(names(dm$map), enc$colnames)
  expect_identical(sort(unlist(dm$map, use.names = FALSE)), seq_len(ncol(dm$D)))
})

# --- fix round 1 ----------------------------------------------------------------

test_that("a constant is written into missing cells only: 0 and -0 keep their bits", {
  m <- cbind(a = c(0, -0, NA, 0), b = c(1, 2, 3, 4))
  enc <- encode_input(m, NULL)
  expect_identical(enc$constant, "a")
  out <- decode_output(first_guess(enc$W, enc, "minVar"), enc)
  expect_identical(1 / out[1:2, "a"], c(Inf, -Inf))        # sign of zero is observed data
  expect_identical(1 / unname(out[4, "a"]), Inf)
  expect_identical(unname(out[3, "a"]), 0)
})

test_that("an observed constant column is not rewritten when another column is constant too", {
  m <- cbind(a = c(1, NA, 3, 4), b = c(5, 5, NA, 5), c = c(1, 2, 3, 5))
  enc <- encode_input(m, NULL)
  out <- decode_output(first_guess(enc$W, enc, "median"), enc)
  expect_identical(unname(out[, "a"]), c(1, 3, 3, 4))
  expect_identical(unname(out[, "b"]), c(5, 5, 5, 5))
})

test_that("decode ignores a W that was altered at an observed cell", {
  d <- sim_mixed(30, seed = 1); d$num1[2] <- NA
  enc <- encode_input(d, NULL); W <- first_guess(enc$W, enc, "minVar")
  W[1, "num1"] <- 1e6; W[1, "fac3"] <- 3; W[1, "const"] <- -1
  expect_identical(decode_output(W, enc)[1, ], d[1, ])
})

test_that("minTotal ignores an excluded column holding the smallest value", {
  m <- cbind(y = c(-1000, 1, 2, 3), a = c(1, NA, 3, 10), b = c(-5, 2, NA, 4))
  enc <- encode_input(m, exclude = "y")
  expect_equal(unname(first_guess(enc$W, enc, "minTotal")[2, "a"]), -5)
})

test_that("design matrix refuses level codes outside 1..K and has unique column names", {
  d <- sim_mixed(30, seed = 1); enc <- encode_input(d, NULL)
  W <- enc$W; W[1, "fac3"] <- 4
  expect_error(design_matrix(W, enc), "fac3.*valid level codes")
  W <- enc$W; W[1, "chr2"] <- 1.5
  expect_error(design_matrix(W, enc), "chr2.*valid level codes")

  d2 <- data.frame(f = factor(c("a", "b", "a", "b", "a")), f_a = c(0.1, 0.2, 0.3, 0.4, 0.5))
  e2 <- encode_input(d2, NULL); dm <- design_matrix(e2$W, e2)
  expect_false(anyDuplicated(colnames(dm$D)) > 0)
  expect_identical(unname(dm$D[, dm$map[["f_a"]]]), d2$f_a)           # map still points right
  expect_identical(unname(dm$D[, dm$map[["f"]][1]]), c(1, 0, 1, 0, 1))
})
