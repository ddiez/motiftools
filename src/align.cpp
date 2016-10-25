#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int whichChar(CharacterVector x, char q) {
  int idx = 0;
  for (int k = 0; k < x.length(); k++) {
    if (x[k] == q) {
      idx = k;
    }
  }
  return(idx);
}

//' @export
// [[Rcpp::export]]
List sw(StringVector x, StringVector y, IntegerMatrix score_matrix, int gap_score = -1) {

  // initialize variables.
  CharacterVector srow = as<CharacterVector>(rownames(score_matrix));
  CharacterVector scol = as<CharacterVector>(colnames(score_matrix));
  IntegerMatrix m(x.length() + 1, y.length() + 1);
  int nrow = m.nrow();
  int ncol = m.ncol();
  int rsize = nrow * ncol;
  IntegerVector ri(rsize);
  IntegerVector rj(rsize);
  //IntegerVector rii(rsize);
  //IntegerVector rjj(rsize);
  CharacterVector rtype(rsize);
  
  // initialize the borders.
  for (int i = 0; i < nrow; i++) {
    m(i, 0) = 0;
    int idx = 0 * nrow + i;
    ri[idx] = i;
    rj[idx] = 0;
    //rii[idx] = i - 1;
    //rjj[idx] = 0;
    rtype[idx] = "v";
  }
  for (int j = 0; j < ncol; j++) {
    m(0, j) = 0;
    int idx = j * nrow + 0;
    ri[idx] = 0;
    rj[idx] = j;
    //rii[idx] = 0;
    //rjj[idx] = j - 1;
    rtype[idx] = "h";
  }
  
  // iterate and fill.
  for (int i = 1; i < m.nrow(); i++) {
    for (int j = 1; j < m.ncol(); j++) {
      // diagonal.
      //int xi = whichChar(srow, x[j - 1]);
      //int yi = whichChar(srow, y[i - 1]);
      int xi = 0;
      for (int k = 0; k < x.length(); k++) {
        if (srow[k] == x[j - 1]) {
          xi = k;
        }
      }
      int yi = 0;
      for (int k = 0; k < x.length(); k++) {
        if (scol[k] == y[j - 1]) {
          yi = k;
        }
      }
      int score = score_matrix(xi, yi);
      //cout << s << endl;
      int d = m(i - 1, j - 1) + score; // has to fix assigning real score above.
      // horizonal
      int h = m(i, j - 1) + gap_score;
      // vertical
      int v = m(i - 1, j) + gap_score;
      
      // check max.
      int idx = j * nrow + i;
      if (d > h && d > v) {
        ri[idx] = i;
        rj[idx] = j;
        // rii[idx] = i - 1;
        // rjj[idx] = j - 1;
        rtype[idx] = "d";
        m(i, j) = d;
      } else if (h > d && h > v) {
        ri[idx] = i;
        rj[idx] = j;
        // rii[idx] = i;
        // rjj[idx] = j - 1;
        rtype[idx] = "h";
        m(i, j) = h;
      } else if (v > d && v > h) {
        ri[idx] = i;
        rj[idx] = j;
        // rii[idx] = i - 1;
        // rjj[idx] = j;
        rtype[idx] = "v";
        m(i, j) = v;
      }
    }
  }
  
  // find indexes for max score.
  int mv = 0;
  int i = 0;
  int j = 0;
  for (int ii = 0; ii < nrow; ii++) {
    for (int jj = 0; jj < ncol; jj++) {
      if (m(ii, jj) >= mv) {
        //int k = j * nrow + i;
        mv = m(ii, jj);
        i = ii;
        j = jj;
      }
    }
  }
  cout << "max_i: " << i << endl;
  cout << "max_j: " << j << endl;
  cout << "max_v: " << mv << endl;
  cout << "rtype: " << rtype << endl;
  
  // backtrace.
  int s = 0;
  while (m(i, j) != 0) {
    int idx = j * nrow + i;

    if (rtype[idx] == "d") {
      i = i - 1;
      j = j - 1;
      s += m(i, j);
    }

    if (rtype[idx] == "h") {
      j = j - 1;
    }

    if (rtype[idx] == "v") {
      i = i - 1;
    }
  }

  List ret; ret["total_score"] = s; ret["scores"] = m;
  return(ret);
}


/*** R
library(Biostrings)
data("BLOSUM62")
sw(c("A", "A"), c("A", "R", "A"), score_matrix = BLOSUM62)
#motifTools:::.sw(c("A", "A"), c("A", "R"), score.matrix = BLOSUM62)
*/