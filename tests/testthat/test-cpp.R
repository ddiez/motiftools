context("Testing C++ code not exported")
test_that("Catch unit tests pass", {
    expect_cpp_tests_pass("motiftools")
})
