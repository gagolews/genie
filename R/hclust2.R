#' @title
#' Faster Hierarchical Clustering in Pseudometric Spaces
#'
#' @param method a single string,
#'        one of \code{single}, \code{complete}, or \code{exemplar}
#' @param d an object of class \code{\link{dist}},
#' \code{NULL}, an R function or a single string
#' @param objects \code{NULL}, numeric matrix or a list
#' @param ... internal tuning parameters
#'
#' @details
#' For compatibility with \code{\link{hclust}}, \code{d}
#' may be an object of class \code{dist}. In such a case, the \code{objects}
#' argument is ignored.
#'
#' If \code{d} is an R function, then \code{objects} should be an R list.
#' Here, \code{d(objects[[i]], objects[[j]])} gives the value of a pseudometric.
#'
#' If \code{objects} is a character vector, then \code{d} is a single string,
#' one of: \code{levenshtein} (or \code{NULL}), \code{dinu}.
#'
#' Otherwise, if \code{objects} is a numeric matrix,
#' then \code{d} is a single string, one of: \code{euclidean} (or \code{NULL}),
#' \code{manhattan}, \code{maximum}, or \code{hamming}.
#'
#' @return
#' A named list of class \code{hclust}, see \code{\link{hclust}},
#' with additional components:
#' \itemize{
#' \item stats ...
#' \item control ...
#' }
#'
#' @export
hclust2 <- function(d=NULL, method=c("single", "complete", "exemplar", "exemplar_naive"), objects=NULL, ...) {
   # TO DO: `single` has list control arg........
   method <- match.arg(method)
   result <- switch(method,
      single=.hclust2_single(d, objects, ...),
      complete=.hclust2_complete(d, objects, ...),
      exemplar=.hclust2_exemplar(d, objects, ...),
      exemplar_naive=.hclust2_exemplar_naive(d, objects, ...)
   )
   result[["call"]] <- match.call()
   result
}
