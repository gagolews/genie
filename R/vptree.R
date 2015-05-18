#' @title VpTree
#'
#' @description
#' VpTree is an implementation of vantage-point tree, which is a metric tree.
#'
#' @rdname vptree
#' @name VpTree
invisible(NULL)

functionSuffix <- "MetricFunction"
vptreeSuffix <- "Tree"
pointsSuffix <- "Points"
extension <- "serialized"

#' @rdname vptree
#' @details
#' \code{vptree_save} saves the tree into a group of 3 files: Tree itself,
#' points which are R objects and metric function. Only prefix of file names 
#' is needed, which can be absolute or relative path.
#'
#' @return
#' \code{vptree_save} does not return anything interesting.
#'
#' @param filenamePrefix the prefix of names of group of 3 files. Can be 
#' absolute or relative path.
vptree_save <- function(vptree, filenamePrefix)
{
   folderName <- dirname(filenamePrefix)
   filenamePrefix <- basename(filenamePrefix)
   
   items <- vptree_getItems(vptree)
   metricFunction <- vptree_getFunction(vptree)
   saveRDS(items, file.path(folderName,paste(filenamePrefix, "_", pointsSuffix, ".", extension, sep="")))
   saveRDS(metricFunction, file.path(folderName,paste(filenamePrefix, "_", functionSuffix, ".", extension, sep="")))
   vptree_serialize(vptree, file.path(folderName,paste(filenamePrefix, "_", vptreeSuffix, ".", extension, sep="")))
}

#' @rdname vptree
#' @details
#' \code{vptree_load} loads the tree from a group of 3 files: Tree itself,
#' points which are R objects and metric function. Only prefix of file names 
#' is needed, which can be absolute or relative path.
#'
#' @return
#' \code{vptree_load} returns a loaded vptree.
#'
#' @param filenamePrefix the prefix of names of group of 3 files. Can be 
#' absolute or relative path.
vptree_load <- function(filenamePrefix)
{
   folderName <- dirname(filenamePrefix)
   filenamePrefix <- basename(filenamePrefix)
   
   items <- readRDS(file.path(folderName,paste(filenamePrefix, "_", pointsSuffix, ".", extension, sep="")))
   metricFunction <- readRDS(file.path(folderName,paste(filenamePrefix, "_", functionSuffix, ".", extension, sep="")))
   vptree <- vptree_read(file.path(folderName,paste(filenamePrefix, "_", vptreeSuffix, ".", extension, sep="")))
   vptree_setItems(vptree, items)
   vptree_setMetricFunction(vptree, metricFunction)
   vptree
}