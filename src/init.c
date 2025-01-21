#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP baseline_reject(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compute_A(SEXP, SEXP, SEXP, SEXP);
extern SEXP compute_A_deriv(SEXP, SEXP, SEXP, SEXP);
extern SEXP compute_alphaC_deriv(SEXP);
extern SEXP compute_d2C_deriv(SEXP);
extern SEXP compute_deriv_G(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compute_derivatives(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compute_G(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compute_gamma(SEXP, SEXP, SEXP, SEXP);
extern SEXP compute_gamma_deriv(SEXP, SEXP, SEXP, SEXP);
extern SEXP compute_intensity_integrals(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP compute_Vd_int(SEXP, SEXP, SEXP);
extern SEXP compute_vector_v(SEXP, SEXP, SEXP);
extern SEXP dbeta_V(SEXP, SEXP);
extern SEXP dbeta_vector_v(SEXP, SEXP, SEXP);
extern SEXP generate_spin_off_events(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP multiply_covariates(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"baseline_reject",             (DL_FUNC) &baseline_reject,             7},
  {"compute_A",                   (DL_FUNC) &compute_A,                   4},
  {"compute_A_deriv",             (DL_FUNC) &compute_A_deriv,             4},
  {"compute_alphaC_deriv",        (DL_FUNC) &compute_alphaC_deriv,        1},
  {"compute_d2C_deriv",           (DL_FUNC) &compute_d2C_deriv,           1},
  {"compute_deriv_G",             (DL_FUNC) &compute_deriv_G,             5},
  {"compute_derivatives",         (DL_FUNC) &compute_derivatives,         8},
  {"compute_G",                   (DL_FUNC) &compute_G,                   6},
  {"compute_gamma",               (DL_FUNC) &compute_gamma,               4},
  {"compute_gamma_deriv",         (DL_FUNC) &compute_gamma_deriv,         4},
  {"compute_intensity_integrals", (DL_FUNC) &compute_intensity_integrals, 5},
  {"compute_Vd_int",              (DL_FUNC) &compute_Vd_int,              3},
  {"compute_vector_v",            (DL_FUNC) &compute_vector_v,            3},
  {"dbeta_V",                     (DL_FUNC) &dbeta_V,                     2},
  {"dbeta_vector_v",              (DL_FUNC) &dbeta_vector_v,              3},
  {"generate_spin_off_events",    (DL_FUNC) &generate_spin_off_events,    6},
  {"multiply_covariates",         (DL_FUNC) &multiply_covariates,         2},
  {NULL, NULL, 0}
};

void R_init_Hausal(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
