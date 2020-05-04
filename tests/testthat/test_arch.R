context("Testing motif architectures")

f <- system.file("files/meme.xml", package = "motiftools")
x <- readMEME(f)
s <- getMotifsBySeq(x)

test_that("getMotifsBySeq works as expected", {
  expect_is(s, "list")
  expect_identical(length(s), 4L)
  expect_identical(names(s), c("RASH_MOUSE", "RASK_MOUSE", "RASN_MOUSE", "RASM_HUMAN"))
  expect_identical(s[["RASH_MOUSE"]], c("3", "1", "2", "4", "5", "6"))
})

test_that("getMotifArchString works as expected", {
  z <- getMotifArchString(x)
  expect_is(z, "character")
  expect_true(length(z) == 4)
  expect_identical(names(z), c("RASH_MOUSE", "RASK_MOUSE", "RASN_MOUSE", "RASM_HUMAN"))
  expect_equal(z[1], c("RASH_MOUSE" = "3-1-2-4-5-6"))
  expect_equal(z["RASH_MOUSE"], c("RASH_MOUSE" = "3-1-2-4-5-6"))
  
  z <- getMotifArchString(x, return.unique = TRUE)
  expect_is(z, "character")
  expect_true(length(z) == 3)
  expect_identical(names(z), NULL)
  expect_equal(z[1], "3-1-2-4-5-6")
})

test_that("getMotifSimilarity works as expected", {
  z <- getMotifSimilarity(x)
  expect_is(z, "matrix")
  expect_equal(dim(z), c(4, 4))
  expect_identical(rownames(z), c("RASH_MOUSE", "RASK_MOUSE", "RASN_MOUSE", "RASM_HUMAN"))
  expect_identical(colnames(z), c("RASH_MOUSE", "RASK_MOUSE", "RASN_MOUSE", "RASM_HUMAN"))
})

test_that("getMotifArchSimilarity works as expected", {
  z <- getMotifArchSimilarity(x)
  expect_is(z, "matrix")
  expect_equal(dim(z), c(4, 4))
  expect_identical(rownames(z), c("RASH_MOUSE", "RASK_MOUSE", "RASN_MOUSE", "RASM_HUMAN"))
  expect_identical(colnames(z), c("RASH_MOUSE", "RASK_MOUSE", "RASN_MOUSE", "RASM_HUMAN"))
})

test_that("plotMotifArchSimilarity works as expected", {
  z <- getMotifArchSimilarity(x)
  p <- plotMotifArchSimilarity(z)
  expect_is(p, "gg")
  expect_is(p$data, "data.frame")
  expect_equal(nrow(p$data), 16)
  expect_equal(as.character(p$data[, 1]), rep(c("RASH_MOUSE", "RASN_MOUSE", "RASK_MOUSE", "RASM_HUMAN"), 4))
  expect_equal(as.character(p$data[, 2]), rep(c("RASH_MOUSE", "RASN_MOUSE", "RASK_MOUSE", "RASM_HUMAN"), each = 4))
})