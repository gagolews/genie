#' @title
#' Faster Hierarchical Clustering in Pseudometric Spaces
#'
#' @param method a single string,
#'        one of \code{single}, \code{complete}
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
#' Otherwise, if \code{objects} is a numeric matrix,
#' then \code{d} is a single string, one of: \code{euclidean} (or \code{NULL}),
#' \code{manhattan}, or \code{maximum}.
#'
#' @return
#' A named list of class \code{hclust}, see \code{\link{hclust}}.
#'
#' @export
hclust2 <- function(d=NULL, method=c("single", "complete"), objects=NULL, ...) {
   method <- match.arg(method)
   result <- switch(method,
      single=.hclust2_single(d, objects, ...),
      complete=.hclust2_complete(d, objects, ...)
   )
   result[["call"]] <- match.call()
   result
}
