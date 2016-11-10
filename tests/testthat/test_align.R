library(motiftools)

test_that("alignment", {
  expect_equal(as.character(align("ALA", "ALA", score.matrix = "BLOSUM62", type = "global")$alignment[1, 1]), "A")
})