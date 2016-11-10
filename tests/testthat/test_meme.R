context("Testing readMEME")

f <- system.file("files/meme_ras/meme.xml", package = "motiftools")
x <- readMEME(f)
m <- getMotifMatrix(x)

test_that("slots", {
  expect_identical(nmotif(x), 9L)
  expect_identical(nseq(x), 4L)
  expect_identical(nmotif(x[, 1]), 1L)
  expect_identical(nseq(x[1, ]), 1L)
  expect_equal(sequenceNames(x), c("RASH_MOUSE", "RASK_MOUSE", "RASN_MOUSE", "RASM_HUMAN"))
  expect_equal(motifNames(x), as.character(1:9))
})

test_that("accessors", {
  expect_is(m, "matrix")
  expect_type(m, "integer")
  expect_identical(m[1, 1], 1L)
  expect_identical(m[1, 9], 0L)
})