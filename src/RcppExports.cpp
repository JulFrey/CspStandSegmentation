// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fast_unlist
IntegerVector fast_unlist(List list, int l);
RcppExport SEXP _CspStandSegmentation_fast_unlist(SEXP listSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type list(listSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_unlist(list, l));
    return rcpp_result_gen;
END_RCPP
}
// fast_unlist_dist
NumericVector fast_unlist_dist(List list, int l);
RcppExport SEXP _CspStandSegmentation_fast_unlist_dist(SEXP listSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type list(listSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    rcpp_result_gen = Rcpp::wrap(fast_unlist_dist(list, l));
    return rcpp_result_gen;
END_RCPP
}
// eigen_decomposition
NumericMatrix eigen_decomposition(S4 las, int k, int ncpu);
RcppExport SEXP _CspStandSegmentation_eigen_decomposition(SEXP lasSEXP, SEXP kSEXP, SEXP ncpuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type las(lasSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type ncpu(ncpuSEXP);
    rcpp_result_gen = Rcpp::wrap(eigen_decomposition(las, k, ncpu));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CspStandSegmentation_fast_unlist", (DL_FUNC) &_CspStandSegmentation_fast_unlist, 2},
    {"_CspStandSegmentation_fast_unlist_dist", (DL_FUNC) &_CspStandSegmentation_fast_unlist_dist, 2},
    {"_CspStandSegmentation_eigen_decomposition", (DL_FUNC) &_CspStandSegmentation_eigen_decomposition, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_CspStandSegmentation(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
