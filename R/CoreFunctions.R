


#' Gets two list of gene hits results, intersects the cell.types and applies an operator
#' on the common cell types cell list
#' 
cell.type.intersect.operator <-  function(hash.index1, hash.index2, operator.function)
{    cell.types <- intersect(keys(hash.index1), keys(hash.index2))
    if (length(cell.types) == 0)
    {
        message('No cell types found')
        return(hash())
    }
    else
    {
        return.set <- hash()
        for( cell.type in cell.types)
        {
            return.set[[cell.type]] <- operator.function(hash.index1[[cell.type]], hash.index2[[cell.type]])
        }
        return(return.set)
    }
}
#' Gets two list of gene hits results, merges cell.types and applies an operator
#' on the common cell types cell list
#'  
cell.type.union.operator <-  function(hash.index1, hash.index2, operator.function)
{   cell.types.common <- union(keys(hash.index1), keys(hash.index2))
    cell.types.diff <- setdiff(keys(hash.index2),keys(hash.index2))
    return.set <- hash.index1
    for( cell.type in cell.types.diff)
    {
        return.set[[cell.type]] <- hash.index2[[cell.type]]
    }
    
    for( cell.type in cell.types)
    {
        res <-  operator.function(hash.index1[[cell.type]], hash.index2[[cell.type]])
        if(length(res) != 0)
        {
            return.set[[cell.type]] <- res
        }
    }
    return(return.set)
}



#' Applies the intersect operator on the cell types and on the cells
#' that participate on common cell.types
#'
#' @param hash.index1 result index to be intersected with hash.index2
#' @param hash.index2 result index to be intersected with hash.index1
#'
#' @return the intersected set
and.operator <- function(hash.index1, hash.index2)
{
    return(operator(hash.index1, hash.index2, intersect, intersect))
}


#' Applies the union operator on two result data structures
#'
#' @param hash.index1 result index to be unionized with hash.index2
#' @param hash.index2 result index to be unionized with hash.index1
#'
#' @return the unionized set
or.operator <- function(hash.index1, hash.index2)
{
    return(operator(hash.index1, hash.index2, union, union))
}

#' Remove cells that have some specific gene hits that we do not want in our list
#'
#' @param hash.index the result index that cell types will be filtered
#' @param hash.diff the cell types to be removed in case there is an overlap
#'
#' @return an new filtered index
not.operator <-  function(hash.index, hash.diff)
{
    return(operator(hash.index, hash.diff, union, setdiff))
}

#' Prototype function to read the mouse cell atlas
readMCA <- function()
{
    filenames <- list.files('/nfs/team218/MCA')
    file.paths <- paste0('/nfs/team218/MCA/', filenames)
    i <- 1
    efdb <- NULL
    for ( rds.path in file.paths)
    {
        print(rds.path)
        sce <- readRDS(rds.path)
        efdb.part <- parseDataset(sce , gsub("_.*", "", filenames[i]), assay.name = 'counts')
        if(is.null(efdb))
        {
            efdb <- efdb.part
        }
        else
        {
            efdb <- mergeDatasets(efdb, efdb.part)
            print(paste("DB size:",object_size(efdb)))
            ## return(efdb)
        }
        i <- i + 1
    }

    return(efdb)
}

#' merge two elias fano indices
#'
#' @param efdb.root the root index
#' @param efdb the index to be merged
#' 
mergeIndices <- function(efdb.root, efdb)
{
    genes.to.merge <- intersect(keys(efdb.root), keys(efdb))
    new.genes <-  setdiff(keys(efdb), keys(efdb.root))
    
    ## print(paste(length(new.genes),length(genes.to.merge)))
    for(gene in new.genes)
    {
        efdb.root[[gene]] <- efdb[[gene]]
        
    }
    for ( gene in genes.to.merge)
    {
        ## maybe we want to merge??
        efdb.root[[gene]][keys(efdb[[gene]])] <- values(efdb[[gene]], simplify = F)
    }
    return(efdb.root)
    
}
