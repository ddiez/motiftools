#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int foo() {
  IntegerMatrix x(2, 2);
  int msize = x.nrow() * x.ncol();
  IntegerVector idx(msize);
  for (int k = 0; k < msize; k ++) {
    idx[k] = k;
  }
  x[0] = 0;
  x[1] = 2;
  x[2] = 50;
  x[3] = 50;
  cout << x << endl;
  
  int mv = 0;
  int mi = 0;
  for (int k = 0; k < msize; k ++) {
    if (x[k] > mv) {
      mv = x[k];
      mi = k;
    }
  }
  cout << "max: " << mv << endl;
  cout << "idx: " << mi << endl;
  
  mv = 0;
  mi = 0;
  for (int i = 0; i < x.nrow(); i++) {
    for (int j = 0; j < x.ncol(); j++) {
      if (x(i, j) > mv) {
        int k = j * x.nrow() + i;
        mv = x(i, j);
        mi = k;
      }
    }
  }
  cout << "max: " << mv << endl;
  cout << "idx: " << mi << endl;
  
  return 0;
}

/*** R
foo()
*/