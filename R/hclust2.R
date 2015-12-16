#' @title
#' Faster Hierarchical Clustering in Pseudometric Spaces
#'
#' @param method a single string,
#'        one of: \code{single}, \code{gini}
#' @param d an object of class \code{\link{dist}},
#' \code{NULL}, an R function, or a single string
#' @param objects \code{NULL}, numeric matrix, a list, or a character vector
#' @param ... internal tuning parameters,
#' e.g., \code{control=list(thresholdGini=0.3)} for the \code{gini} method.
#'
#' @details
#' For compatibility with \code{\link{hclust}}, \code{d} may be an object
#' of class \code{dist}. In such a case, the \code{objects} argument is ignored.
#'
#' If \code{d} is an R function, then \code{objects} should be an R list.
#' Here, \code{d(objects[[i]], objects[[j]])} gives the value of a pseudometric.
#'
#' If \code{objects} is a character vector or a list, then \code{d}
#' should be a single string, one of: \code{levenshtein} (or \code{NULL}),
#' \code{dinu}, \code{hamming}, \code{euclinf}. Note that the list must consist either
#' of integer or of numeric vectors only (depending on the pseudometric of choice).
#' Each string must be in ASCII, but you can always convert it to UTF-32 with
#' \link[stringi]{stri_enc_toutf32}.
#'
#' Otherwise, if \code{objects} is a numeric matrix, then \code{d} should be
#' a single string, one of: \code{euclidean} (or \code{NULL}),
#' \code{manhattan}, \code{maximum}, or \code{hamming}.
#'
#' @return
#' A named list of class \code{hclust}, see \code{\link{hclust}},
#' with additional components:
#' \itemize{
#'      \item stats   performance statistics
#'      \item control internal tuning parameters used
#' }
#'
#' @export
hclust2 <- function(
      d=NULL,
      method=c("single", "gini", "exemplar", "exemplar2",
               "exemplar_approx", "exemplar_naive"), # , "complete", "exemplar", "exemplar_naive"
      objects=NULL,
      ...)
{
   method <- match.arg(method)
   result <- switch(method,
      single=.hclust2_single(d, objects, ...),
      gini=.hclust2_gini(d, objects, ...),
      single_approx=.hclust2_single_approx(d, objects, ...),
#       complete=.hclust2_complete(d, objects, ...),
      exemplar=.hclust2_exemplar(d, objects, ...),
      exemplar2=.hclust2_exemplar2(d, objects, ...),
      exemplar_approx=.hclust2_exemplar_approx(d, objects, ...),
      exemplar_naive=.hclust2_exemplar_naive(d, objects, ...)
   )
   result[["call"]] <- match.call()
   result[["method"]] <- method

   if (any(result[["height"]]<0)) {
      nonNegative <- which(result[["height"]]>=0)
      lastNonNegative <- nonNegative[length(nonNegative)]
      result[["height"]][1:lastNonNegative] <-
         approx(nonNegative,
            result[["height"]][nonNegative],
            1:lastNonNegative)$y
      result[["height"]][result[["height"]] < 0] <- cummax(-result[["height"]][result[["height"]] < 0])
   }
   result
}
