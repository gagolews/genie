functionSuffix <- "MetricFunction"
vptreeSuffix <- "Tree"
pointsSuffix <- "Points"
extension <- "serialized"

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