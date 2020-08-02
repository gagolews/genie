#' @title
#' Fast Hierarchical Clustering in Spaces Equipped With
#' a Dissimilarity Measure
#'
#' @description
#' The reference implementation of the fast, robust and outlier resistant
#' Genie algorithm described in (Gagolewski, Bartoszuk, Cena, 2016).
#' Note that the \code{genie} package has been superseded by \code{genieclust},
#' see \code{\link[genieclust]{gclust}} and \code{\link[genieclust]{genie}}
#' for more details.
#'
#' @param d an object of class \code{\link[stats]{dist}},
#' \code{NULL}, or a single string, see below
#' @param objects \code{NULL}, numeric matrix, a list, or a character vector
#' @param thresholdGini single numeric value in [0,1],
#' threshold for the Gini index, 1 gives the standard single linkage algorithm
#' @param useVpTree single logical value, whether to use a vantage-point tree
#' to speed up nearest neighbour searching in low-dimensional spaces
#' @param ... internal parameters used to tune up the algorithm
#'
#' @details
#' The time needed to apply a hierarchical clustering algorithm
#' is most often dominated by the number of computations of a pairwise
#' dissimilarity measure. Such a constraint, for larger data sets,
#' puts  at a disadvantage the use of all the classical linkage
#' criteria but the single linkage one. However, it is known that the single
#' linkage clustering algorithm is very sensitive to outliers, produces highly
#' skewed dendrograms, and therefore usually does not reflect the true
#' underlying data structure -- unless the clusters are well-separated.
#'
#' To overcome its limitations, in (Gagolewski, Bartoszuk, Cena, 2016)
#' we proposed a new hierarchical clustering linkage
#' criterion. Namely, our algorithm links two clusters in such a way that a chosen
#' economic inequity measure (here, the Gini index) of the cluster
#' sizes does not increase drastically above a given threshold. The
#' benchmarks indicate a high practical usefulness of the introduced method:
#' it most often outperforms the Ward or average linkage in terms of the
#' clustering quality while retaining the single linkage speed.
#' The algorithm can be run in parallel (via OpenMP) on multiple threads
#' to speed up its execution further on.
#' Its memory overhead is small: there is no need to precompute the complete
#' distance matrix to perform the computations in order to obtain a desired
#' clustering.
#'
#' For compatibility with \code{\link[stats]{hclust}}, \code{d} may be an object
#' of class \code{\link[stats]{dist}}. In such a case, the \code{objects}
#' argument is ignored. Note that such an object requires ca. \emph{8n(n-1)/2}
#' bytes of computer's memory, where \emph{n} is the number of objects to cluster,
#' and therefore this setting can be used to analyse data sets of sizes
#' up to about 10,000-50,000.
#'
#' If \code{objects} is a character vector or a list, then \code{d}
#' should be a single string, one of: \code{levenshtein} (or \code{NULL}),
#' \code{hamming}, \code{dinu} (Dinu, Sgarro, 2006),
#' or \code{euclinf} (Cena et al., 2015).
#' Note that the list must consist
#' either of integer or of numeric vectors only (depending on the dissimilarity
#' measure of choice). On the other hand, each string must be in ASCII,
#' but you can always convert it to UTF-32 with
#' \code{\link[stringi]{stri_enc_toutf32}}.
#'
#' Otherwise, if \code{objects} is a numeric matrix (here, each row
#' denotes a distinct observation), then \code{d} should be
#' a single string, one of: \code{euclidean_squared} (or \code{NULL}),
#' \code{euclidean} (which yields the same results as \code{euclidean_squared})
#' \code{manhattan}, \code{maximum}, or \code{hamming}.
#'
#' If \code{useVpTree} is \code{FALSE}, then the dissimilarity measure
#' of choice is guaranteed to be computed for each unique pair of \code{objects}
#' only once.
#'
#' @return
#' A named list of class \code{hclust}, see \code{\link[stats]{hclust}},
#' with additional components:
#' \itemize{
#'      \item \code{stats} -   performance statistics
#'      \item \code{control} - internal parameters used
#' }
#'
#' @examples
#' library("datasets")
#' data("iris")
#' h <- hclust2(objects=as.matrix(iris[,2:3]), thresholdGini=0.2)
#' plot(iris[,2], iris[,3], col=cutree(h, 3), pch=as.integer(iris[,5]), asp=1, las=1)
#'
#' @references
#' Cena A., Gagolewski M., Mesiar R., Problems and challenges of information
#' resources producers' clustering, \emph{Journal of Informetrics} 9(2), 2015,
#' pp. 273-284.
#'
#' Dinu L.P., Sgarro A., A Low-complexity Distance for DNA Strings,
#' \emph{Fundamenta Informaticae} 73(3), 2006, pp. 361-372.
#'
#' Gagolewski M., Bartoszuk M., Cena A.,
#' Genie: A new, fast, and outlier-resistant hierarchical clustering algorithm,
#' \emph{Information Sciences} 363, 2016, pp. 8-23.
#'
#' Gagolewski M., Cena A., Bartoszuk M.
#' \emph{Hierarchical clustering via penalty-based aggregation and the Genie
#' approach}, In: Torra V. et al. (Eds.), \emph{Modeling Decisions for
#' Artificial Intelligence} (\emph{Lecture Notes in Artificial Intelligence}
#' 9880), Springer, 2016.
#'
#' @importFrom stats approx
#' @importFrom genieclust gclust
#' @importFrom genieclust genie
#' @export
hclust2 <- function(d=NULL, objects=NULL, thresholdGini=0.3, useVpTree=FALSE, ...)
{
   opts <- list(thresholdGini=thresholdGini, useVpTree=useVpTree, ...)
   result <- .hclust2_gini(d, objects, opts)
   result[["call"]] <- match.call()
   result[["method"]] <- "gini"

   if (any(result[["height"]]<0)) {
      # corrections for departures from ultrametricity
      # negative heights denote force Genie merges
      # we could just use have used cummax, but then we'd get multiple
      # merges at the same level; instead we'll linearly interpolate
      # between the points
      nonNegative <- which(result[["height"]]>=0)
      lastNonNegative <- nonNegative[length(nonNegative)]
      result[["height"]][1:lastNonNegative] <-
         approx(nonNegative, # linear interpolation
            result[["height"]][nonNegative],
            1:lastNonNegative)$y
      result[["height"]][result[["height"]] < 0] <- cummax(-result[["height"]][result[["height"]] < 0])
   }
   result
}
