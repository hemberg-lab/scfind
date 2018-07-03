#' Builds an scfind index from a SingleCellExperiment object
#' 
#'
#' @param sce object of SingleCellExperiment class
#' @param dataset.name name of the dataset that will be prepended in each cell_type
#' @param assay.name name of the SingleCellExperiment assay that will be considered for the generation of the index
#' @param cell.type.label the cell.type metadata of the colData SingleCellExperiment that will be used for the index
#' 
#' @name buildCellTypeIndex
#'
#' @return a `data.frame` containing calculated gene index
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment rowData rowData<- colData colData<- assayNames assays
#' @importFrom hash hash
#' 
#' @importFrom Rcpp sourceCpp
#' @useDynLib scfind 
#' 
buildCellTypeIndex.SCESet <- function(sce, dataset.name, assay.name, cell.type.label)
{

    if(grepl(dataset.name,'.'))
    {
        error("The dataset name should not contain any dots")
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
            error(paste('Assay name', assay.name, 'not found in the SingleCellExperiment'))
        }
        else
        {
            message(paste("Generating index for", dataset.name, "from '", assay.name, "' assay"))
        }
        exprs <- "[["(sce@assays$data, assay.name)

        loadModule('EliasFanoDB')
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
    
    index <- new("SCFind", index = ef, datasets = dataset.name)
    return(index)
}

#' @rdname buildCellTypeIndex
#' @aliases buildCellTypeIndex
setMethod("buildCellTypeIndex",
          signature(sce = "SingleCellExperiment"),
          buildCellTypeIndex.SCESet)

#' Load a binary file and deserializes the database
#'
#' @param filename the compatible file
#'
#' @name loadSerializedObject
#' @return the loaded database
load.from.serialized.object <- function(filename)
{
    loadModule('EliasFanoDB')
    ef <-  new(EliasFanoDB)
    ef$loadFromFile(filename)
    index <-  new("SCFind", index = ef, datasets = filename)
    return(index)

}

#' @rdname loadFromFile
#' @aliases loadFromFile
setMethod("loadFromFile",  definition = load.from.serialized.object)


#' Merges external index to existing object
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

#' @rdname mergeDataset
#' @aliases mergeDataset
setMethod("mergeDataset",
          signature(
              object = "SCFind",
              new.object = "SCFind"
          ),
          merge.dataset.from.object)

#' Merges a SingleCellExperiment object into the SCFind index
#'
#' @param object the root scfind object
#' @param sce
#' @param dataset.name
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


#' Find cell types associated with a given gene list
#' 
#' Calculates p-values of a log-likelihood of a list of genes to be associated
#' with each cell type. Log-likelihood is based on gene expression values.
#'
#' @param gene_index a data.frame with cell types in columns and genes in rows
#' @param gene_list genes that need to be searched in the gene_index
#' 
#' @name findCellTypes
#'
#' @return a named numeric vector containing p-values
findCellTypes.geneList <- function(object, gene.list)
{
    if (is.null(object))
    {
        stop("Please define a scfind object using the `object` parameter!")
    }
    if (is.null(gene.list))
    {
        stop("Please define a list of genes using the `gene.list` parameter!")
    }
    
    return(object@index$findCellTypes(gene.list))
}

#' @rdname findCellType
#' @aliases findCellType
setMethod("findCellTypes", signature(object = "SCFind", gene.list = "character"), findCellTypes.geneList)

