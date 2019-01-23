#' Builds an \code{SCFind} object from a \code{SingleCellExperiment} object
#'
#' This function will index a \code{SingleCellExperiment} as an SCFind index.
#'
#' @param sce object of SingleCellExperiment class
#' @param dataset.name name of the dataset that will be prepended in each cell_type
#' @param assay.name name of the SingleCellExperiment assay that will be considered for the generation of the index
#' @param cell.type.label the cell.type metadata of the colData SingleCellExperiment that will be used for the index
#' 
#' @name buildCellTypeIndex
#'
#' @return an SCFind object
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment rowData rowData<- colData colData<- assayNames assays
#' @importFrom hash hash
#' @importFrom methods new
#' 
#' @importFrom Rcpp cpp_object_initializer
#' @useDynLib scfind 
#' 
buildCellTypeIndex.SCESet <- function(sce, dataset.name, assay.name = 'counts', cell.type.label = 'cell_type1')
{

    if (grepl(dataset.name,'.'))
    {
        stop("The dataset name should not contain any dots")
    }
    
    
    cell.types.all <- as.factor("[["(colData(sce), cell.type.label))
    cell.types <- levels(cell.types.all)
    new.cell.types <- hash(keys = cell.types, values = paste0(dataset.name, '.', cell.types))
    genenames <- unique(rowData(sce)$feature_symbol)
    
    if (length(cell.types) > 0)
    {
        non.zero.cell.types <- c()
        index <- hash()
        ## print(paste("Found", length(cell.types), "clusters on", ncol(sce), "cells"))
        if( ! assay.name %in% assayNames(sce))
        {
            stop(paste('Assay name', assay.name, 'not found in the SingleCellExperiment'))
        }
        else
        {
            message(paste("Generating index for", dataset.name, "from '", assay.name, "' assay"))
        }
        exprs <- "[["(sce@assays$data, assay.name)

        ef <- new(EliasFanoDB)
        for (cell.type in cell.types) {
            inds.cell <- which(cell.type == cell.types.all)
            if(length(inds.cell) < 2)
            {
                ## print(paste('Skipping', cell.type))
                next
            }
            non.zero.cell.types <- c(non.zero.cell.types, cell.type)
            message(paste("\tIndexing", cell.type, "as", new.cell.types[[cell.type]], " with ", length(inds.cell), " cells."))
            cell.type.exp <- exprs[,inds.cell]
            if(is.matrix(exprs))
            {
                ef$indexMatrix(new.cell.types[[cell.type]], cell.type.exp)
            }
            else
            {
                ef$indexMatrix(new.cell.types[[cell.type]], as.matrix(cell.type.exp))
            }
        }
    }
    
    index <- new("SCFind", index = ef, datasets = dataset.name, metadata = list())
    return(index)
}

#' @rdname buildCellTypeIndex
#' @aliases buildCellTypeIndex buildIndex
setMethod("buildCellTypeIndex",
          signature(sce = "SingleCellExperiment"),
          buildCellTypeIndex.SCESet)

#' This function serializes the DB and save the object as an rds file
#'
#' This function can be used to enable the user save the loaded file in a database
#' to avoid re-indexing and re-merging individual assays.
#'
#' After serializing and saving it clears the redundant bytestream from memory
#' because the memory is already loaded in memory
#' @param object an SCFind object
#' @param file the target filename that the object is going to be stored
#'
#' @return the \code{SCFind} object
#' @name saveObject
save.serialized.object <- function(object, file){
    object@serialized <- object@index$getByteStream()
    a <- saveRDS(object, file)
    # Clear the serialized stream
    object@serialized <- raw()
    gc()
    return(object)
}

#' @rdname saveObject
#' @aliases saveObject
setMethod("saveObject",  definition = save.serialized.object)


#' This function loads a saved \code{SCFind} object and deserializes
#' the object and loads it into an in-memory database.
#'
#' After loading the database it clears the loaded bytestream from the memory.
#'
#' @param filename the filepath of a specialized serialized scfind object
#' 
#' @return an \code{SCFind} object
#' @name loadObject
#'
#' @useDynLib scfind
load.serialized.object <- function(filename){
    object <-  readRDS(filename)
    # Deserialize object
    object@index <-  new(EliasFanoDB)
    success <- object@index$loadByteStream(object@serialized)
    object@serialized <- raw()
    gc()
    ## Dirty hack so we do not have to rebuild again every scfind index
    if(is.null(object@metadata))
    {
        object@metadata <- list()
    }
    return(object)
}

#' @rdname loadObject
#' @aliases loadObject
setMethod("loadObject",  definition = load.serialized.object)



#' Merges an external index into the existing object
#'
#' This function is useful to merge \code{SCFind} indices.
#' After this operation object that was merged can be discarded.
#' 
#' The only semantic limitation for merging two databases is to
#' have different dataset names in the two different indices.
#' If that is not case user may run into problems masking datasets
#' from the different datasets while there is a possibility of having
#' different cell types under the same name. This will most likely cause
#' undefined behavior during queries.
#' 
#' @param object the root scfind object
#' @param new.object external scfind object to be merged
#'
#' @name mergeDataset
#' @return the new extended object
#' 
merge.dataset.from.object <- function(object, new.object)
{
    common.datasets <- intersect(new.object@datasets, object@datasets)
    
    message(paste('Merging', new.object@datasets))
    if(length(common.datasets) != 0)
    {
        warning("Common dataset names exist, undefined merging behavior, please fix this...")
    }
    
    object@index$mergeDB(new.object@index)
    object@datasets <- c(object@datasets, new.object@datasets)
    return(object)
}

#' Used to merge multiple eliasfanoDB
#'
#' 
#' @rdname mergeDataset
#' @aliases mergeDataset mergeObjects
setMethod("mergeDataset",
          signature(
              object = "SCFind",
              new.object = "SCFind"
          ),
          merge.dataset.from.object)

#' Merges a SingleCellExperiment object into the SCFind index
#'
#' It creates an \code{SCFind} for the individual assay and then invokes
#' the \code{mergeDataset} method obeying the same semantic rules.
#'
#' @param object the root scfind object
#' @param sce the \code{SingleCellExperiment} object to be merged
#' @param dataset.name a dataset name for the assay
#' @name mergeSCE
#' @return the new object with the sce object merged
merge.dataset.from.sce <- function(object, sce, dataset.name)
{
    object.to.merge <- buildCellTypeIndex(sce, dataset.name)
    return(mergeDataset(object, object.to.merge))
}
#' @rdname mergeSCE
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @aliases mergeSCE
setMethod("mergeSCE",
          signature(
              object = "SCFind",
              sce = "SingleCellExperiment",
              dataset.name = "character"
          ),
          merge.dataset.from.sce)


#' Query Optimization Function for SCFind objects.
#'
#' This function can be used with quite long gene lists
#' that otherwise would have no cell hits in the database
#' 
#' @name markerGenes
#' @param object SCFind object
#' @param gene.list A list of Genes existing in the database
#' @param datasets the datasets of the objects to be considered
#' 
#' @return hierarchical list of queries and their respective scores
find.marker.genes <-  function(object, gene.list, datasets, message = 0)
{
    datasets <- select.datasets(object, datasets)
    results <- object@index$findMarkerGenes(as.character(caseCorrect(object, gene.list)), as.character(datasets), 5, message)
    
    return(results)
}


#' @rdname markerGenes
#' @aliases markerGenes
setMethod("markerGenes",
          signature(
              object = "SCFind",
              gene.list = "character"),
          find.marker.genes)

#' Find marker genes for a specific cell type
#'
#' @name cellTypeMarkers
#' 
#' @param object SCFind object
#' @param cell.types the cell types that we want to extract the marker genes
#' @param background.cell.types the universe of cell.types to consider
#' @param top.k how many genes to retrieve
#' @param sort.field the dataframe will be sorted according to this field
#'
#' @return a data.frame that each row represent a gene score for a specific cell type 
cell.type.marker <- function(object, cell.types, background.cell.types, top.k, sort.field)
{
    if (missing(background.cell.types))
    {
        message("Considering the whole DB..")
        background.cell.types <- cellTypeNames(object)
    }
    all.cell.types <- object@index$cellTypeMarkers(cell.types, background.cell.types)
    if (!(sort.field %in% colnames(all.cell.types)))
    {
        message(paste("Column", sort.field, "not found"))
        sort.field <- 'f1'
    }
    all.cell.types <- all.cell.types[order(all.cell.types[[sort.field]], decreasing = T)[1:top.k],]
    return(all.cell.types)
}


#' @rdname cellTypeMarkers
#' @aliases cellTypeMarkers
setMethod("cellTypeMarkers",
          signature(
              object = "SCFind",
              cell.types = "character"
          ),
          cell.type.marker)


#' Return a vector with all existing cell type names in the database
#' 
#' @name cellTypeNames
#' @param object SCFind object
#'
#' @return a character list
get.cell.types.names <- function(object)
{
    return(object@index$getCellTypes())
}
#' @rdname cellTypeMarkers
#' @aliases cellTypeMarkers
setMethod("cellTypeNames",
          signature(
              object = "SCFind"),
          get.cell.types.names)


#' Evaluate a user specific query by calculating the precision recall metrics
#'
#' @name evaluateMarkers
#' @param object the \code{SCFind} object
#' @param gene.list the list of genes to be evaluated
#' @param cell.types a list of cell types for the list to evaluated
#' @param background.cell.types the universe of cell.types to consider
#' @param sort.field the dataframe will be sorted according to this field
#'
#' @return a DataFrame that each row represent a gene score for a specific cell type
#'
evaluate.cell.type.markers <- function(object, gene.list, cell.types, background.cell.types, sort.field){
    if(missing(background.cell.types))
    {
        message("Considering the whole DB..")
        background.cell.types <- cellTypeNames(object)
    }
    all.cell.types <- object@index$evaluateCellTypeMarkers(cell.types, caseCorrect(object, gene.list), background.cell.types)

    if(!(sort.field %in% colnames(all.cell.types)))
    {
        message(paste("Column", sort.field, "not found"))
        sort.field <- 'f1'
    }
    all.cell.types <- all.cell.types[order(all.cell.types[[sort.field]]),]
    return(all.cell.types)
    
}

#' @rdname evaluateMarkers
#' @aliases evaluateMarkers
setMethod("evaluateMarkers",
          signature(
              object = "SCFind",
              gene.list = "character"
              ),
          evaluate.cell.type.markers)
              


#' Runs a query and performs the hypergeometric test for the retrieved cell types
#'
#' @name hyperQueryCellTypes
#' @param object the \code{SCFind} object
#' @param gene.list the list of genes to be queried
#' @param datasets the datasets vector that will be tested as background for the hypergeometric test
#'
#' @return a DataFrame that contains all cell types with the respective cell cardinality and the hypergeometric test
cell.types.phyper.test <- function(object, gene.list, datasets)
{
    result <- findCellTypes.geneList(object, caseCorrect(object, gene.list), datasets)
    
    return(phyper.test(object, result, datasets))
    
}

#' @rdname hyperQueryCellTypes
#' @aliases hyperQueryCellTypes
#' 
setMethod("hyperQueryCellTypes",
          signature(object = "SCFind",
                    gene.list = "character"),
          cell.types.phyper.test)


#' Find cell types associated with a given gene list. All cells
#' returned express all of the genes in the given gene list
#' 
#' @param object the \code{SCFind} object
#' @param gene.list genes to be searched in the gene.index
#' @param datasets the datasets that will be considered
#' 
#' @name findCellTypes
#' @return a named numeric vector containing p-values
findCellTypes.geneList <- function(object, gene.list, datasets)
{
    
    datasets <- select.datasets(object, datasets)
    return(object@index$findCellTypes(caseCorrect(object, gene.list), datasets))
    
}

#' @rdname findCellTypes
#' @aliases findCellTypes
setMethod("findCellTypes", 
          signature(object = "SCFind",
                    gene.list = "character",
                    dataset = "character"), 
          findCellTypes.geneList)

#' Get all genes in the database
#'
#' @name scfindGenes
#' 
#' @param object the \code{scfind} object
#'
#' @return the list of genes present in the database
scfind.get.genes.in.db <- function(object)
{
    
    return(object@index$genes())

}


#' @rdname scfindGenes
#' @aliases scfindGenes
setMethod("scfindGenes", signature(object = "SCFind"), scfind.get.genes.in.db)

