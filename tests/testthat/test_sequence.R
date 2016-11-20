context("Testing sequence functions")

f <- system.file("files/ras.faln", package = "motiftools")

test_that("conservationMatrix works", {
  aln <- readAAMultipleAlignment(f)
  cons <- conservationMatrix(aln)
  
  expect_type(cons, "integer")
  expect_true(cons[1,1] == 1L)
  expect_true(cons[1,14] == 5L)
})

test_that("plotConservationMatrix works", {
  aln <- readAAMultipleAlignment(f)
  BiocGenerics::rownames(aln) <- sub(" .*", "", BiocGenerics::rownames(aln)) # need to fix sequence names.
  
  library(ape)
  tree <- read.tree(system.file("files/ras.tree", package = "motiftools"))
  
  p <- plotConservationMatrix(aln, tree, show.tips = TRUE)
  expect_is(p, "gtable")
  
  cons <- conservationMatrix(aln)
  p <- plotConservationMatrix(cons, tree)
  expect_is(p, "gtable")
})
