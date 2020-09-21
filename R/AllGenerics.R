#' The scfind main class object
#'
#' @export
setClass("SCFind",
         representation(
             index = "ANY", #Index is type of Rcpp_EliasFanoDB but this removes the warning
             datasets = "character",
             serialized = "raw",
             metadata = "list"
             ))


#' 
#' @export
setGeneric(name = "buildCellTypeIndex",
           def = function(sce,
                          dataset.name = '',
                          assay.name = 'logcounts',
                          cell.type.label = 'cell_type1',
                          qb = 2)
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


#' queries cells that contain all the genes from the list
#' @export
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
                                               sort.field = 'f1'
                                               ){
    standardGeneric("cellTypeMarkers")
})


#' @export
#'
setGeneric(name = "cellTypeNames", function(object, datasets){
   standardGeneric("cellTypeNames")
})

#' @export
#'
setGeneric(name = "evaluateMarkers", function(object,
                                              gene.list,
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

setGeneric(name = "markerGenes", function(object, gene.list, datasets, exhaustive = FALSE, support.cutoff = -1)
{
    standardGeneric("markerGenes")
})


#' @export
setGeneric(name = "scfindGenes", function(object){
    standardGeneric("scfindGenes")
})


#' Find out how many cell-types each gene is found
#' 
#' @export
setGeneric(name = "findCellTypeSpecificities", function(object, 
                                                        gene.list,
                                                        datasets,
                                                        min.cells = 10, 
                                                        min.fraction = .25){
    standardGeneric("findCellTypeSpecificities")
})

#' Find out how many tissues each gene is found
#' 
#' @export
setGeneric(name = "findTissueSpecificities", function(object, 
                                                        gene.list,
                                                        min.cells = 10){
    standardGeneric("findTissueSpecificities")
})

#' Find the set of genes that are ubiquitously expressed in a query of cell types
#' 
#' @export
setGeneric(name = "findHouseKeepingGenes", function(object, 
                                                    cell.types,
                                                    min.recall = .5, 
                                                    max.genes = 1000){
    standardGeneric("findHouseKeepingGenes")
})

#' Find the signature genes for a cell-type
#' 
#' @export
setGeneric(name = "findGeneSignatures", function(object, 
                                                 cell.types,
                                                 max.genes = 1000,
                                                 min.cells = 10,
                                                 max.pval = 0){
    standardGeneric("findGeneSignatures")
})

#' Look at all other genes and rank them based on the similarity of their expression pattern to the pattern defined by the gene query
#' 
#' @export
setGeneric(name = "findSimilarGenes", function(object, gene.list, datasets, top.k = 5){
    standardGeneric("findSimilarGenes")
})




#' @export
setGeneric(name = "scfindShinyServer", function(object){
    standardGeneric("scfindShinyServer")
})

