#include <testthat.h>
#include "misc.h"

context("Testing misc.cpp") {
  CharacterVector x(2);
  x(0) = "A"; x(1) = "B";
  CharacterVector y(3);
  y(0) = "A"; y(1) = "C"; y(2) = "D";
  
  test_that("getUniqueLetters returns correct dimension and values") {
    expect_true(getUniqueLetters(x, y).length() == 4);
  }
  
  test_that("getScoreMatrix returns correct dimension and values") {
    expect_true(getScoreMatrix(x, y).nrow() == 4);
    expect_true(getScoreMatrix(x, y).ncol() == 4);
    expect_true(getScoreMatrix(x, y)(1,1) == 1);
    expect_true(getScoreMatrix(x, y)(1,2) == -1);
  }
}
