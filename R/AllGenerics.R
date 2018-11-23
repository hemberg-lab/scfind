#' The scfind main class object
#' @export
setClass("SCFind",
         representation(
             index = "Rcpp_EliasFanoDB",
             datasets = "character",
             serialized = "raw"))


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


#' @export 
setGeneric(name = "mergeDataset", def = function(object, new.object) {
    standardGeneric("mergeDataset")
})
 

#' @export 
setGeneric(name = "mergeSCE", def = function(object, sce, dataset.name) {
    standardGeneric("mergeSCE")
})


#' queries cells that contain the genes from the list#
#' @export 
setGeneric(name = "queryGene", def = function(object, gene, datasets) {
    standardGeneric("queryGene")
})

#' queries cells that contain all the genes from the list
#' @export
#' 
setGeneric(name = "findCellTypes", function(object, gene.list, datasets) {
    standardGeneric("findCellTypes")
})


#' return all the gene markers for a specified cell.type
#'
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
#' 
#' @export
setGeneric(name = "saveObject", function(object, file){
    standardGeneric("saveObject")
})

#' @export
setGeneric(name = "hyperQueryCellTypes", function(object,
                                                  gene.list,
                                                  datasets){
    standardGeneric("hyperQueryCellTypes")

})

#' Performs query optimization and return the best candidate gene sets
#'
#' @export
setGeneric(name = "markerGenes", function(object, gene.list, datasets)
{
    standardGeneric("markerGenes")
})


#' @export
setGeneric(name = "scfindGenes", function(object){
    standardGeneric("scfindGenes")
})

#' @export
setGeneric(name = "scfindShinyServer", function(object){
    standardGeneric("scfindShinyServer")
})

