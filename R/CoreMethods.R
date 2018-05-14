#' Build a cell type Index
#' 
#' Calculates a fraction of expressed cells per gene per cell type
#'
#' @param sce object of SingleCellExperiment class
#' @param dataset.name name of the dataset that will be prepended in each cell_type
#' 
#' @name buildCellTypeIndex
#'
#' @return a `data.frame` containing calculated gene index
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom hash hash
#' @importFrom bit as.bit
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
    
    print(paste("Reading", dataset.name))
    d <- sce
    cell.types.all <- as.factor("[["(colData(sce), cell.type.label))
    cell.types <- levels(cell.types.all)
    new.cell.types <- hash(keys = cell.types, values = paste0(dataset.name, '.', cell.types))
    genenames <- unique(rowData(d)$feature_symbol)
    
    if (length(cell.types) > 0)
    {
        non.zero.cell.types <- c()
        object <- hash()
        ## print(paste("Found", length(cell.types), "clusters on", ncol(sce), "cells"))
        if( ! assay.name %in% names(assays(sce)))
        {
            error(paste('Assay name', assay.name, 'not found in the SingleCellExperiment'))
        }
        exprs <- "[["(d@assays$data, assay.name)
        ## Check if we have a sparse represantation
        if(!is.matrix(exprs))
        {
            ## Cast the matrix, expensive operation
            exprs <- as.matrix(exprs)
        }
        genes.nonzero <- which(rowSums(exprs) > 0)
        if(length(genes.nonzero) == 0)
        {
            return(hash())
        }
        
        print(paste0("Non zero genes ", length(genes.nonzero) ))
        for (cell.type in cell.types) {
            inds.cell <- which(cell.type == cell.types.all)
            if(length(inds.cell) < 2)
            {
                ## print(paste('Skipping', cell.type))
                next
            }
            non.zero.cell.types <- c(non.zero.cell.types, cell.type)
            print(paste("Indexing", cell.type," to ", new.cell.types[[cell.type]], " with ", length(inds.cell), " cells."))
            ## Calculate the baseline probability that a gene will be expressed in a cell
            object[new.cell.types[[cell.type]]] <- hash(
                keys = genenames[genes.nonzero],
                values = apply(exprs[genes.nonzero, inds.cell], 1,
                               function (x)
                               {
                                   ef.gene <- eliasFanoCodingCpp(x)
                                   if(length(ef.gene) == 0){
                                       return(ef.gene)
                                   }
                                   return(list(
                                       H = as.bit(ef.gene$H),
                                       L = as.bit(ef.gene$L),
                                       l = ef.gene$l
                                   ))
                               }))
            
              }
    }
    
    message('Finalizing index...')
    index.value.model <- hash()
    for (cell.type in non.zero.cell.types)
    {
        index.value.model[[ as.character(new.cell.types[[cell.type]]) ]] <- list()
    }
    
    new.obj <- hash()
    for( gene in genenames)
    {
        new.obj[[gene]] <- index.value.model
    }


    for (cell.type in non.zero.cell.types)
    {
        new.cell.type <- new.cell.types[[cell.type]]
        for (gene in keys(object[[new.cell.type]]))
        {
            ## This if clause has to be consistent with the Rcpp side
            if (length(object[[new.cell.type]][[gene]]) != 0)
            {
                new.obj[[gene]][[new.cell.type]] <- object[[new.cell.type]][[gene]]
            }
        }
    }


    index <- new("SCFind", index = new.obj, datasets = c(dataset.name))
    return(index)
}

#' @rdname buildCellTypeIndex
#' @aliases buildCellTypeIndex
setMethod("buildCellTypeIndex",
          signature(sce = "SingleCellExperiment"),
          buildCellTypeIndex.SCESet)


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
    if(len(common.datasets) != 0)
    {
        warning("Common dataset names exist, undefined merging behavior, please fix this...")
    }
    
    object@index <- mergeIndices(object@index, new.object@index)
    object@datasets <- c(object@datasets, new.object@datasets)
    return(object)
}
#' @rdname mergeDataset
#' @aliases mergeDataset
setMethod("mergeDataset",
          signature(object = "SCFind",
                    new.object = "SCFind"),
          merge.dataset.from.object)

#' Merges another sce object
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
          signature(object = "SCFind",
                    sce = "SingleCellExperiment",
                    dataset.name = "character"), merge.dataset.from.sce)


#' Retrieves all relative celltypes with their correspodent cell matches
#'
#' @param object an scfind object
#' @param gene an scfind object
#'
#' @name queryGene
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib scfind
#' 
#' @return nada
query.gene <- function(object, gene)
{
    efdb <- object@index
    if(is.null(efdb[[gene]]))
    {
        warning(paste('Requested gene', gene, 'not available in the index'))
        return(hash())
    }
    else
    {
        return(hash(keys = keys(efdb[[gene]]),
                    values = lapply(keys(efdb[[gene]]),
                                    FUN = function(cell.type)
                                    {   
                                        v <- efdb[[gene]][[cell.type]]
                                        return(
                                            eliasFanoDecodingCpp(
                                                as.logical(v$H),
                                                as.logical(v$L),
                                                v$l))
                                    })))
    }
}

#' @rdname queryGene
#' @aliases queryGene
setMethod("queryGene", signature(object = "SCFind", gene = "character"), query.gene)

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
#'
#' @importFrom hash hash keys
#' @importFrom stats pchisq
#' @importFrom methods is
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

    gene.results <- hash()
    for(gene in gene.list)
    {
        gene.results[[gene]] <- query.gene(object, gene)
    }

    genes <-  keys(gene.results)
    genes.queried <-  genes[1]
    
    query.results <- gene.results[[genes[1]]] # cold start operator
    genes <- tail(genes, -1) # pop first element from gene list
    
    for(gene in genes)
    {
        query.results <- and.operator(query.results, gene.results[[gene]])
        existing.cell.types <- keys(query.results)
        genes.queried <- c(genes.queried, gene)
        message(paste("Genes queried (", cat(genes.queried),") with", length(existing.cell.types)))
        if(length(existing.cell.types) == 0)
        {
            warning("Empty set, breaking operation")
        }        
    }
    return(query.results)
}

#' @rdname findCellType
#' @aliases findCellType
setMethod("findCellTypes", signature(object = "SCFind", gene.list = "character"), findCellTypes.geneList)

