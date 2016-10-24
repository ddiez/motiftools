// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// whichChar
int whichChar(CharacterVector x, char q);
RcppExport SEXP motifTools_whichChar(SEXP xSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< char >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(whichChar(x, q));
    return rcpp_result_gen;
END_RCPP
}
// alignc
List alignc(StringVector x, StringVector y, IntegerMatrix score_matrix, int gap_score);
RcppExport SEXP motifTools_alignc(SEXP xSEXP, SEXP ySEXP, SEXP score_matrixSEXP, SEXP gap_scoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< StringVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type score_matrix(score_matrixSEXP);
    Rcpp::traits::input_parameter< int >::type gap_score(gap_scoreSEXP);
    rcpp_result_gen = Rcpp::wrap(alignc(x, y, score_matrix, gap_score));
    return rcpp_result_gen;
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
