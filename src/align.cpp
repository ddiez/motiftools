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

// Perform Smith-Waterman local alignment.
// [[Rcpp::export]]
List sw(CharacterVector x, CharacterVector y, Rcpp::Nullable<IntegerMatrix> score_matrix = R_NilValue, int gap_score = -1, int debug = 0) {

  IntegerMatrix sm;
  if (score_matrix.isNull()) {
    sm = getScoreMatrix(x, y);
  } else {
    sm = score_matrix.get();
  }
  // initialize variables.
  CharacterVector schar = as<CharacterVector>(rownames(sm));
  IntegerMatrix m(y.length() + 1, x.length() + 1);
  int nrow = m.nrow();
  int ncol = m.ncol();
  
  CharacterVector rn(nrow);
  for (int k = 0; k < rn.length(); k++) {
    if (k == 0) {
      rn[k] = "";
    } else {
      rn[k] = y[k - 1];
    }
  }
  rownames(m) = rn;
  
  CharacterVector cn(ncol);
  for (int k = 0; k < cn.length(); k++) {
    if (k == 0) {
      cn[k] = "";
    } else {
      cn[k] = x[k - 1];
    }
  }
  colnames(m) = cn;
  
  // initialize the borders.
  CharacterMatrix rtype(nrow, ncol);
  rtype(0, 0) = "r";
  for (int i = 1; i < nrow; i++) {
    m(i, 0) = 0;
    rtype(i, 0) = "v";
  }
  
  for (int j = 1; j < ncol; j++) {
    m(0, j) = 0;
    rtype(0, j) = "h";
  }
  
  // iterate and fill.
  for (int i = 1; i < nrow; i++) { // iterate over y.
    for (int j = 1; j < ncol; j++) { // iterate over x.
      // diagonal.
      int xi = 0;
      int yi = 0;
      for (int k = 0; k < schar.length(); k++) {
        if (schar[k] == x[j - 1]) {
          xi = k;
        }
        if (schar[k] == y[i - 1]) {
          yi = k;
        }
      }
      int d = m(i - 1, j - 1) + sm(xi, yi);
      
      // horizonal
      int h = m(i, j - 1) + gap_score;
      
      // vertical
      int v = m(i - 1, j) + gap_score;
      
      // check max.
      if (d > h && d > v) {
        rtype(i, j) = "d";
        m(i, j) = d;
      } else if (h > d && h > v) {
        rtype(i, j) = "h";
        m(i, j) = h;
      } else if (v > d && v > h) {
        rtype(i, j) = "v";
        m(i, j) = v;
      } else {
        // all equal, choose any.
        rtype(i, j) = "d"; // this?
        m(i, j) = d;
      }
    }
  }
  
  // find indexes for max score.
  int mv = 0;
  int mi = 0;
  int mj = 0;
  for (int ii = 0; ii < nrow; ii++) {
    for (int jj = 0; jj < ncol; jj++) {
      if (m(ii, jj) >= mv) {
        mv = m(ii, jj);
        mi = ii;
        mj = jj;
      }
    }
  }

  // backtrace.
  int i = mi;
  int j = mj;
  int s = 0;
  int rsize = 0;
  while (m(i, j) != 0) { // maybe rtype(i, j) == "r"?
    
    if (rtype(i, j) == "d") {
      s += m(i, j);
      i = i - 1;
      j = j - 1;
      rsize++;
      continue;
    }

    if (rtype(i, j) == "h") {
      j = j - 1;
      rsize++;
      continue;
    }

    if (rtype(i, j) == "v") {
      i = i - 1;
      rsize++;
      continue;
    }
  }
  
  // generate alignment.
  CharacterMatrix aln(2, rsize); // matrix of proper size.
  int k = rsize - 1; // begin to fill from the end. Why?
  i = mi;
  j = mj;
  while (m(i, j) != 0) { // maybe rtype(i, j) == "r"?

    if (rtype(i, j) == "d") {
      aln(0, k) = x[j - 1];
      aln(1, k) = y[i - 1];
      i = i - 1;
      j = j - 1;
      k--;
      continue;
    }
    
    if (rtype(i, j) == "h") {
      aln(0, k) = x[j - 1];
      aln(1, k) = "-";
      k--;
      j = j - 1;
      continue;
    }
    
    if (rtype(i, j) == "v") {
      aln(0, k) = "-";
      aln(1, k) = y[i - 1];
      k--;
      i = i - 1;
      continue;
    }
  }
  
  rownames(aln) = CharacterVector::create("x", "y");
  CharacterVector cn2(aln.ncol());
  for (int k = 0; k < cn2.length(); k++) {
    cn2[k] = k + 1;
  }
  colnames(aln) = cn2;
  DataFrame alndf = aln; // for compatibility with previous way- may change.
  
  List ret; ret["alignment"] = alndf; ret["score"] = s;
  if (debug) ret["scores"] = m;
  return(ret);
}

// Perform Needleman-Wunsch global alignment.
// [[Rcpp::export]]
List nw(CharacterVector x, CharacterVector y, Rcpp::Nullable<IntegerMatrix> score_matrix = R_NilValue, int gap_score = -1, int debug = 0) {
  
  IntegerMatrix sm;
  if (score_matrix.isNull()) {
    sm = getScoreMatrix(x, y);
  } else {
    sm = score_matrix.get();
  }
  // initialize variables.
  CharacterVector schar = as<CharacterVector>(rownames(sm));
  IntegerMatrix m(y.length() + 1, x.length() + 1);
  int nrow = m.nrow();
  int ncol = m.ncol();
  
  CharacterVector rn(nrow);
  for (int k = 0; k < rn.length(); k++) {
    if (k == 0) {
      rn[k] = "";
    } else {
      rn[k] = y[k - 1];
    }
  }
  rownames(m) = rn;
  
  CharacterVector cn(ncol);
  for (int k = 0; k < cn.length(); k++) {
    if (k == 0) {
      cn[k] = "";
    } else {
      cn[k] = x[k - 1];
    }
  }
  colnames(m) = cn;
  
  // initialize the borders.
  CharacterMatrix rtype(nrow, ncol);
  rtype(0, 0) = "r";
  for (int i = 1; i < nrow; i++) {
    m(i, 0) = -1 * i;
    rtype(i, 0) = "v";
  }
  
  for (int j = 1; j < ncol; j++) {
    m(0, j) = -1 * j;
    rtype(0, j) = "h";
  }
  
  // iterate and fill.
  for (int i = 1; i < nrow; i++) { // iterate over y.
    for (int j = 1; j < ncol; j++) { // iterate over x.
      // diagonal.
      int xi = 0;
      int yi = 0;
      for (int k = 0; k < schar.length(); k++) {
        if (schar[k] == x[j - 1]) {
          xi = k;
        }
        if (schar[k] == y[i - 1]) {
          yi = k;
        }
      }
      int d = m(i - 1, j - 1) + sm(xi, yi);
      
      // horizonal
      int h = m(i, j - 1) + gap_score;
      
      // vertical
      int v = m(i - 1, j) + gap_score;
      
      // check max.
      if (d > h && d > v) {
        rtype(i, j) = "d";
        m(i, j) = d;
      } else if (h > d && h > v) {
        rtype(i, j) = "h";
        m(i, j) = h;
      } else if (v > d && v > h) {
        rtype(i, j) = "v";
        m(i, j) = v;
      } else {
        // all equal, choose any.
        rtype(i, j) = "d"; // this?
        m(i, j) = d;
      }
    }
  }
  
  // find indexes for max score.
  int mv = 0;
  int mi = 0;
  int mj = 0;
  for (int ii = 0; ii < nrow; ii++) {
    for (int jj = 0; jj < ncol; jj++) {
      if (m(ii, jj) >= mv) {
        mv = m(ii, jj);
        mi = ii;
        mj = jj;
      }
    }
  }
  
  // backtrace.
  int i = mi;
  int j = mj;
  int s = 0;
  int rsize = 0;
  while (m(i, j) != 0) { // maybe rtype(i, j) == "r"?
    
    if (rtype(i, j) == "d") {
      s += m(i, j);
      i = i - 1;
      j = j - 1;
      rsize++;
      continue;
    }
    
    if (rtype(i, j) == "h") {
      j = j - 1;
      rsize++;
      continue;
    }
    
    if (rtype(i, j) == "v") {
      i = i - 1;
      rsize++;
      continue;
    }
  }
  
  // generate alignment.
  CharacterMatrix aln(2, rsize); // matrix of proper size.
  int k = rsize - 1; // begin to fill from the end. Why?
  i = mi;
  j = mj;
  while (m(i, j) != 0) { // maybe rtype(i, j) == "r"?
    
    if (rtype(i, j) == "d") {
      aln(0, k) = x[j - 1];
      aln(1, k) = y[i - 1];
      i = i - 1;
      j = j - 1;
      k--;
      continue;
    }
    
    if (rtype(i, j) == "h") {
      aln(0, k) = x[j - 1];
      aln(1, k) = "-";
      k--;
      j = j - 1;
      continue;
    }
    
    if (rtype(i, j) == "v") {
      aln(0, k) = "-";
      aln(1, k) = y[i - 1];
      k--;
      i = i - 1;
      continue;
    }
  }
  
  rownames(aln) = CharacterVector::create("x", "y");
  CharacterVector cn2(aln.ncol());
  for (int k = 0; k < cn2.length(); k++) {
    cn2[k] = k + 1;
  }
  colnames(aln) = cn2;
  DataFrame alndf = aln; // for compatibility with previous way- may change.
  
  List ret; ret["alignment"] = alndf; ret["score"] = s;
  if (debug) ret["scores"] = m;
  return(ret);
}

/*** R
#sw(c("A", "L", "D"), c("A", "R", "L", "E"))
library(Biostrings)
data("BLOSUM62")
# local:
sw(c("L", "D"), c("A", "R", "L", "E"), score_matrix = BLOSUM62, debug = FALSE)
motiftools:::.sw(c("L", "D"), c("A", "R", "L", "E"), score.matrix = BLOSUM62, debug = FALSE)
# global:
nw(c("L", "D"), c("A", "R", "L", "E"), score_matrix = BLOSUM62, debug = FALSE)
motiftools:::.nw(c("L", "D"), c("A", "R", "L", "E"), score.matrix = BLOSUM62, debug = FALSE)
*/