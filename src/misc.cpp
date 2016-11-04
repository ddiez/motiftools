#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// Get a CharacterVector with the unique characters in x and y.
CharacterVector getUniqueLetters(CharacterVector x, CharacterVector y) {
  CharacterVector z(x.length() + y.length());
  // fill x.
  for (int k = 0; k < x.length(); k++) {
    z(k) = x(k);
  }
  // fill y.
  for (int k = 0; k < y.length(); k++) {
    z(x.length() + k) = y(k);
  }
  return(unique(z));
}

// Get score matrix from the from input sequences using the provided scores.
IntegerMatrix getScoreMatrix(CharacterVector x, CharacterVector y, int match_score = 1, int mismatch_score = -1) {
  // get unique character vector.
  CharacterVector z = getUniqueLetters(x, y);
  
  // create matrix.
  IntegerMatrix sm(z.length(), z.length());
  for (int i = 0; i < sm.nrow(); i++) {
    for (int j = 0; j < sm.ncol(); j++) {
      if (i == j)
        sm(i, j) = match_score;
      else
        sm(i, j) = mismatch_score;
    }
  }
  rownames(sm) = z;
  colnames(sm) = z;
  return(sm);
}
