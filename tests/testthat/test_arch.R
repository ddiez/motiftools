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