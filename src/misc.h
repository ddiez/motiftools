#include <Rcpp.h>
using namespace Rcpp;

CharacterVector getUniqueLetters(CharacterVector x, CharacterVector y);
IntegerMatrix getScoreMatrix(CharacterVector x, CharacterVector y, int match_score = 1, int mismatch_score = -1);