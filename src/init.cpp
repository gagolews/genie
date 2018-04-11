// #include <R.h>
// #include <Rinternals.h>
// #include <Rversion.h>
// #include <cstdlib>
// #include <R_ext/Rdynload.h>
//
//
// extern "C" SEXP _genie_hclust2_gini(SEXP distanceSEXP, SEXP objectsSEXP, SEXP controlSEXP);
//
// static const R_CallMethodDef cCallMethods[] = {
//     {"genie_hclust2_gini", (DL_FUNC) &_genie_hclust2_gini, 3},
//     {NULL, NULL, 0}
// };
//
// extern "C" void R_init_genie(DllInfo *dll)
// {
//    R_registerRoutines(dll, NULL, cCallMethods, NULL, NULL);
//    R_useDynamicSymbols(dll, (Rboolean)FALSE);
// #if defined(R_VERSION) && R_VERSION >= R_Version(3, 0, 0)
//    R_forceSymbols(dll, (Rboolean)TRUE);
// #endif
// }
//
