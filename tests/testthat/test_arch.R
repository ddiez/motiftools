context("Testing motif architectures")

f <- system.file("files/meme_ras/meme.xml", package = "motiftools")
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
  
  z <- getMotifArchString(x, convert.to.letter = TRUE, return.unique = TRUE)
  expect_is(z, "character")
  expect_true(length(z) == 3)
  expect_identical(names(z), NULL)
  expect_equal(z[1], "CABDEF")
})

test_that("getMotifArchSimilarity works as expected", {
  z <- getMotifArchSimilarity(x)
  expect_is(z, "matrix")
  expect_equal(dim(z), c(4, 4))
})