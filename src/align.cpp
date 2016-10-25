#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
List sw(CharacterVector x, CharacterVector y, IntegerMatrix score_matrix, int gap_score = -1) {

  // initialize variables.
  CharacterVector schar = as<CharacterVector>(rownames(score_matrix));
  IntegerMatrix m(y.length() + 1, x.length() + 1);
  int nrow = m.nrow();
  int ncol = m.ncol();
  CharacterMatrix rtype(nrow, ncol);
  
  // initialize the borders.
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
       cout << "x: " << x[j - 1] << endl;
       cout << "y: " << y[i - 1] << endl;
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
      cout << xi << yi << endl;
      
      int score = score_matrix(xi, yi);
      cout << "score : " << score << endl;
      int d = m(i - 1, j - 1) + score; // has to fix assigning real score above.
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
  cout << "score matrix: " << endl << m << endl;
  int mv = 0;
  int i = 0;
  int j = 0;
  for (int ii = 0; ii < nrow; ii++) {
    for (int jj = 0; jj < ncol; jj++) {
      if (m(ii, jj) >= mv) {
        mv = m(ii, jj);
        i = ii;
        j = jj;
      }
    }
  }
  cout << "max_i: " << i << endl;
  cout << "max_j: " << j << endl;
  cout << "max_v: " << mv << endl;
  cout << "rtype: " << rtype(i, j) << endl;
  
  // backtrace.
  int s = 0;
  while (m(i, j) != 0) { // maybe rtype(i, j) == "r"?
    cout << "i: " << i << endl;
    cout << "j: " << j << endl;
    cout << "score: " << m(i, j) << endl;
    cout << "type: " << rtype(i, j) << endl;
    
    if (rtype(i, j) == "d") {
      s += m(i, j);
      i = i - 1;
      j = j - 1;
      continue;
    }

    if (rtype(i, j) == "h") {
      j = j - 1;
      continue;
    }

    if (rtype(i, j) == "v") {
      i = i - 1;
      continue;
    }
  }

  List ret; ret["total_score"] = s; ret["scores"] = m;
  return(ret);
}


/*** R
library(Biostrings)
data("BLOSUM62")
sw(c("A", "L"), c("A", "R", "L"), score_matrix = BLOSUM62)
#motifTools:::.sw(c("A", "L"), c("A", "R", "L"), score.matrix = BLOSUM62, debug = TRUE)
*/