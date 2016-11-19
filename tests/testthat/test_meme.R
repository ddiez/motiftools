context("Testing readMEME")

f <- system.file("files/meme_ras/meme.xml", package = "motiftools")
x <- readMEME(f)
m <- getMotifMatrix(x)

test_that("class", {
  expect_true(class(x) == "MotifSearchResult")
})

test_that("slots", {
  expect_identical(nmotif(x), 9L)
  expect_identical(nseq(x), 4L)
  expect_identical(nmotif(x[, 1]), 1L)
  expect_identical(nseq(x[1, ]), 1L)
  expect_equal(sequenceNames(x), c("RASH_MOUSE", "RASK_MOUSE", "RASN_MOUSE", "RASM_HUMAN"))
  expect_equal(motifNames(x), as.character(1:9))
})

test_that("accessors", {
  expect_true(class(m) == "matrix")
  expect_is(m, "matrix")
  expect_type(m, "integer")
  expect_identical(m[1, 1], 1L)
  expect_identical(m[1, 9], 0L)
})

test_that("scores returns the correct values", {
  s <- scores(x)
  expect_is(s, "list")
  expect_is(s[[1]], "matrix")
  expect_identical(s[[1]][1, 1], -260L)
})

test_that("pwm returns the correct values", {
  p <- pwm(x)
  expect_is(p, "list")
  expect_is(p[[1]], "matrix")
  expect_identical(p[[1]][1, 1], 0)
})

test_that("plotMotifMatrix works", {
  p <- plotMotifMatrix(x)
  expect_is(p, "gtable")
})