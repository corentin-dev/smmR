// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// setSeed
void setSeed(unsigned int seed);
RcppExport SEXP _SMM_setSeed(SEXP seedSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type seed(seedSEXP);
    setSeed(seed);
    return R_NilValue;
END_RCPP
}
// getChain
arma::vec getChain(arma::vec J, arma::vec T);
RcppExport SEXP _SMM_getChain(SEXP JSEXP, SEXP TSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type J(JSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type T(TSEXP);
    rcpp_result_gen = Rcpp::wrap(getChain(J, T));
    return rcpp_result_gen;
END_RCPP
}
// getProcesses
List getProcesses(List sequences);
RcppExport SEXP _SMM_getProcesses(SEXP sequencesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type sequences(sequencesSEXP);
    rcpp_result_gen = Rcpp::wrap(getProcesses(sequences));
    return rcpp_result_gen;
END_RCPP
}
// getCountingProcesses
List getCountingProcesses(List& Jm, List& Lm, arma::uword& s, arma::uword& kmax);
RcppExport SEXP _SMM_getCountingProcesses(SEXP JmSEXP, SEXP LmSEXP, SEXP sSEXP, SEXP kmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List& >::type Jm(JmSEXP);
    Rcpp::traits::input_parameter< List& >::type Lm(LmSEXP);
    Rcpp::traits::input_parameter< arma::uword& >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::uword& >::type kmax(kmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(getCountingProcesses(Jm, Lm, s, kmax));
    return rcpp_result_gen;
END_RCPP
}
// getCountingNiuj
arma::cube getCountingNiuj(List& Ym, List& Um, arma::uword& s, arma::uword& kmax);
RcppExport SEXP _SMM_getCountingNiuj(SEXP YmSEXP, SEXP UmSEXP, SEXP sSEXP, SEXP kmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List& >::type Ym(YmSEXP);
    Rcpp::traits::input_parameter< List& >::type Um(UmSEXP);
    Rcpp::traits::input_parameter< arma::uword& >::type s(sSEXP);
    Rcpp::traits::input_parameter< arma::uword& >::type kmax(kmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(getCountingNiuj(Ym, Um, s, kmax));
    return rcpp_result_gen;
END_RCPP
}
// computeKernelNonParamEndcensoring
arma::cube computeKernelNonParamEndcensoring(arma::cube& p);
RcppExport SEXP _SMM_computeKernelNonParamEndcensoring(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(computeKernelNonParamEndcensoring(p));
    return rcpp_result_gen;
END_RCPP
}
// C_rdweibull
double C_rdweibull(double& q, double& beta);
RcppExport SEXP _SMM_C_rdweibull(SEXP qSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double& >::type q(qSEXP);
    Rcpp::traits::input_parameter< double& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(C_rdweibull(q, beta));
    return rcpp_result_gen;
END_RCPP
}
// simulateParam
List simulateParam(unsigned int& seed, arma::Col<arma::uword>& nsim, arma::vec& init, arma::mat& ptrans, arma::Mat<int>& distr, arma::mat& param1, arma::mat& param2, bool censBeg, bool censEnd);
RcppExport SEXP _SMM_simulateParam(SEXP seedSEXP, SEXP nsimSEXP, SEXP initSEXP, SEXP ptransSEXP, SEXP distrSEXP, SEXP param1SEXP, SEXP param2SEXP, SEXP censBegSEXP, SEXP censEndSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< arma::Col<arma::uword>& >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type init(initSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ptrans(ptransSEXP);
    Rcpp::traits::input_parameter< arma::Mat<int>& >::type distr(distrSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type param1(param1SEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type param2(param2SEXP);
    Rcpp::traits::input_parameter< bool >::type censBeg(censBegSEXP);
    Rcpp::traits::input_parameter< bool >::type censEnd(censEndSEXP);
    rcpp_result_gen = Rcpp::wrap(simulateParam(seed, nsim, init, ptrans, distr, param1, param2, censBeg, censEnd));
    return rcpp_result_gen;
END_RCPP
}
// simulateNonParam
List simulateNonParam(unsigned int& seed, arma::Col<arma::uword>& nsim, arma::vec& init, arma::mat& ptrans, arma::cube& distr, bool censBeg, bool censEnd);
RcppExport SEXP _SMM_simulateNonParam(SEXP seedSEXP, SEXP nsimSEXP, SEXP initSEXP, SEXP ptransSEXP, SEXP distrSEXP, SEXP censBegSEXP, SEXP censEndSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int& >::type seed(seedSEXP);
    Rcpp::traits::input_parameter< arma::Col<arma::uword>& >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type init(initSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type ptrans(ptransSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type distr(distrSEXP);
    Rcpp::traits::input_parameter< bool >::type censBeg(censBegSEXP);
    Rcpp::traits::input_parameter< bool >::type censEnd(censEndSEXP);
    rcpp_result_gen = Rcpp::wrap(simulateNonParam(seed, nsim, init, ptrans, distr, censBeg, censEnd));
    return rcpp_result_gen;
END_RCPP
}
// convolution
arma::vec convolution(arma::vec& f, arma::vec& g);
RcppExport SEXP _SMM_convolution(SEXP fSEXP, SEXP gSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type f(fSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type g(gSEXP);
    rcpp_result_gen = Rcpp::wrap(convolution(f, g));
    return rcpp_result_gen;
END_RCPP
}
// matrixConvolution
arma::cube matrixConvolution(arma::cube& A, arma::cube& B);
RcppExport SEXP _SMM_matrixConvolution(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube& >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(matrixConvolution(A, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SMM_setSeed", (DL_FUNC) &_SMM_setSeed, 1},
    {"_SMM_getChain", (DL_FUNC) &_SMM_getChain, 2},
    {"_SMM_getProcesses", (DL_FUNC) &_SMM_getProcesses, 1},
    {"_SMM_getCountingProcesses", (DL_FUNC) &_SMM_getCountingProcesses, 4},
    {"_SMM_getCountingNiuj", (DL_FUNC) &_SMM_getCountingNiuj, 4},
    {"_SMM_computeKernelNonParamEndcensoring", (DL_FUNC) &_SMM_computeKernelNonParamEndcensoring, 1},
    {"_SMM_C_rdweibull", (DL_FUNC) &_SMM_C_rdweibull, 2},
    {"_SMM_simulateParam", (DL_FUNC) &_SMM_simulateParam, 9},
    {"_SMM_simulateNonParam", (DL_FUNC) &_SMM_simulateNonParam, 7},
    {"_SMM_convolution", (DL_FUNC) &_SMM_convolution, 2},
    {"_SMM_matrixConvolution", (DL_FUNC) &_SMM_matrixConvolution, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_SMM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
