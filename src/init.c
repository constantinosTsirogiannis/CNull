#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP CNull_communities_individual_based_sampling_alpha(SEXP, SEXP);
extern SEXP CNull_communities_individual_based_sampling_beta(SEXP, SEXP);
extern SEXP CNull_communities_individual_based_sampling_beta_interleaved_matrices(SEXP, SEXP);
extern SEXP CNull_communities_permutation_sampling_alpha(SEXP, SEXP);
extern SEXP CNull_communities_permutation_sampling_beta(SEXP, SEXP);
extern SEXP CNull_communities_permutation_sampling_beta_interleaved_matrices(SEXP, SEXP);
extern SEXP CNull_compute_pvalues(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"CNull_communities_individual_based_sampling_alpha",                     (DL_FUNC) &CNull_communities_individual_based_sampling_alpha,                     2},
    {"CNull_communities_individual_based_sampling_beta",                      (DL_FUNC) &CNull_communities_individual_based_sampling_beta,                      2},
    {"CNull_communities_individual_based_sampling_beta_interleaved_matrices", (DL_FUNC) &CNull_communities_individual_based_sampling_beta_interleaved_matrices, 2},
    {"CNull_communities_permutation_sampling_alpha",                          (DL_FUNC) &CNull_communities_permutation_sampling_alpha,                          2},
    {"CNull_communities_permutation_sampling_beta",                           (DL_FUNC) &CNull_communities_permutation_sampling_beta,                           2},
    {"CNull_communities_permutation_sampling_beta_interleaved_matrices",      (DL_FUNC) &CNull_communities_permutation_sampling_beta_interleaved_matrices,      2},
    {"CNull_compute_pvalues",                                                 (DL_FUNC) &CNull_compute_pvalues,                                                 2},
    {NULL, NULL, 0}
};

void R_init_CNull(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
