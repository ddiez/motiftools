// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// sw
List sw(CharacterVector x, CharacterVector y, IntegerMatrix score_matrix, int gap_score);
RcppExport SEXP motifTools_sw(SEXP xSEXP, SEXP ySEXP, SEXP score_matrixSEXP, SEXP gap_scoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type score_matrix(score_matrixSEXP);
    Rcpp::traits::input_parameter< int >::type gap_score(gap_scoreSEXP);
    rcpp_result_gen = Rcpp::wrap(sw(x, y, score_matrix, gap_score));
    return rcpp_result_gen;
END_RCPP
}
// goo
void goo();
RcppExport SEXP motifTools_goo() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    goo();
    return R_NilValue;
END_RCPP
}
// foo
int foo();
RcppExport SEXP motifTools_foo() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(foo());
    return rcpp_result_gen;
END_RCPP
}
