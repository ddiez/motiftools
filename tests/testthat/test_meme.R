context("Testing readMEME")

f <- system.file("files/meme_ras/meme.xml", package = "motiftools")
x <- readMEME(f)

test_that("readMEME returns object of correct class", {
  expect_true(class(x) == "MotifSearchResult")
  expect_message(show(x), "MotifSearchResult object.")
  expect_message(show(x), "Motif tool: MEME")
  expect_message(show(x), "Number of sequences: 4")
  expect_message(show(x), "Number of motifs: 9")
})

test_that("methods return correct values", {
  expect_identical(nmotif(x), 9L)
  expect_identical(nseq(x), 4L)
  expect_identical(nmotif(x[, 1]), 1L)
  expect_identical(nseq(x[1, ]), 1L)
  expect_equal(sequenceNames(x), c("RASH_MOUSE", "RASK_MOUSE", "RASN_MOUSE", "RASM_HUMAN"))
  expect_equal(motifNames(x), as.character(1:9))
})

test_that("scores returns the correct values", {
  s <- scores(x)
  expect_is(s, "list")
  expect_true(inherits(s[[1]], "matrix"))
  expect_identical(s[[1]][1, 1], -260L)
})

test_that("pwm returns the correct values", {
  p <- pwm(x)
  expect_is(p, "list")
  expect_true(inherits(p[[1]], "matrix"))
  expect_identical(p[[1]][1, 1], 0)
})

test_that("getMotifMatrix works as expected", {
  m <- getMotifMatrix(x)
  expect_true(inherits(m, "matrix"))
  expect_type(m, "integer")
  expect_identical(m[1, 1], 1L)
  expect_identical(m[1, 9], 0L)
})

test_that("plotMotifMatrix works as expected", {
  # default plot.
  p <- plotMotifMatrix(x)
  expect_is(p, "gtable")
  
  # with tree
  library(ape)
  tree <- read.tree(system.file("files/ras.tree", package = "motiftools"))
  p <- plotMotifMatrix(x, tree = tree)
  expect_is(p, "gtable")
  
  # matrix.
  p <- plotMotifMatrix(getMotifMatrix(x))
  expect_is(p, "gtable")
  
  # add some options.
  annot <- matrix(1, nrow = nseq(x), ncol = 1, dimnames = list(sequenceNames(x)))
  annot["RASH_MOUSE", ] <- 2
  annot["RASK_MOUSE", ] <- 3
  p <- plotMotifMatrix(list("MEME motifs" = x), fill = c("white", "grey"), col = "black", annot = list(annot), annot.fill = list(c("white", "red", "blue")), show.tips = TRUE)
  expect_is(p, "gtable")
})