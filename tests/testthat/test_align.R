library(motiftools)

test_that("alignment", {
  expect_equal(align("ALA", "ALA", score.matrix = "BLOSUM62", type = "global")$alignment[1, 1], "A")
  expect_equal(align("ALA", "ALA", score.matrix = "BLOSUM62", type = "local")$alignment[1, 1], "A")
})