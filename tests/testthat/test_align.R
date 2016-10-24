library(motifTools)

test_that("alignment", {
  expect_equal(align("FOO", "FOO")$alignment[1, 1], "F")
})