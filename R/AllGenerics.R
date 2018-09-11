
#' The scfind main class object
#' @export
setClass("SCFind", representation(index = "Rcpp_EliasFanoDB", datasets = "character", serialized = "raw"))

#' @examples TODO
#' 
#' @export
setGeneric(name = "buildCellTypeIndex",
           def = function(sce,
                          dataset.name = '',
                          assay.name = 'logcounts',
                          cell.type.label = 'cell_type1')
           {
               standardGeneric("buildCellTypeIndex")
           })


#' 
#' @examples TODO
#' @export 
setGeneric(name = "mergeDataset", def = function(object, new.object) {
    standardGeneric("mergeDataset")
})
 
#' @examples TODO
#' @export 
setGeneric(name = "mergeSCE", def = function(object, sce, dataset.name) {
    standardGeneric("mergeSCE")
})


#' @export
#' 
#' @examples TODO
#'
#' 
setGeneric(name = "queryGene", def = function(object, gene, datasets) {
    standardGeneric("queryGene")
})

#' @export
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
#' # this is needed to calculate dropout rate for feature selection
#' # important: normcounts have the same zeros as raw counts (fpkm)
#' counts(sce) <- normcounts(sce)
#' logcounts(sce) <- log2(normcounts(sce) + 1)
#' # use gene names as feature symbols
#' rowData(sce)$feature_symbol <- rownames(sce)
#' isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' index <- buildCellIndex(sce)
#' res <- findCell(index, genelist = c('SOX6', 'SNAI3'))
#' 
setGeneric(name = "findCellTypes", function(object, gene.list, datasets) {
    standardGeneric("findCellTypes")
})


#' @export
#'
setGeneric(name = "cellTypeMarkers" ,  function(object,
                                               cell.types,
                                               background.cell.types,
                                               top.k = 5,
                                               sort.field = 'f1'){
    standardGeneric("cellTypeMarkers")

})

#' @export
#'
setGeneric(name = "cellTypeNames", function(object){
   standardGeneric("cellTypeNames")
})

#' @export
#'
setGeneric(name = "evaluateMarkers", function(object, gene.list,
                                              cell.types,
                                              background.cell.types,
                                              sort.field = 'f1'){
    standardGeneric("evaluateMarkers")
})


#' @export
setGeneric(name = "scfindShiny", function(object) {
    standardGeneric("scfindShiny")
})


#' Generic to be used instead of readRDS
#' @export
setGeneric(name = "loadObject", function(filename){
    standardGeneric("loadObject")
})

#' Generic to be used instead of saveRDS
#' @export
setGeneric(name = "saveObject", function(object, file){
    standardGeneric("saveObject")
})

#' @export
setGeneric(name = "markerGenes", function(object, gene.list, datasets)
{
    standardGeneric("markerGenes")
})



