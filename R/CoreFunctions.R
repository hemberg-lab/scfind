#' Gets two list of gene hits results, intersects the cell.types and applies an operator
#' on the common cell types cell list
#'
#' @return nothing
cell.type.intersect.operator <- function(hash.index1, hash.index2, operator.function)
{
    cell.types <- intersect(keys(hash.index1), keys(hash.index2))
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
#' @name and
#' @param hash.index1 result index to be intersected with hash.index2
#' @param hash.index2 result index to be intersected with hash.index1
#'
#' @return the intersected set
and.operator <- function(hash.index1, hash.index2)
{
    return(cell.type.intersect.operator(hash.index1, hash.index2, intersect))
}


#' Applies the union operator on two result data structures
#'
#' @name or
#' @param hash.index1 result index to be unionized with hash.index2
#' @param hash.index2 result index to be unionized with hash.index1
#'
#' @return the unionized set
or.operator <- function(hash.index1, hash.index2)
{
    return(cell.type.union.operator(hash.index1, hash.index2, union, union))
}

#' Remove cells that have some specific gene hits that we do not want in our list
#'
#' @name not
#' @param hash.index the result index that cell types will be filtered
#' @param hash.diff the cell types to be removed in case there is an overlap
#'
#' @return an new filtered index
not.operator <-  function(hash.index, hash.diff)
{
    return(cell.type.union.operator(hash.index, hash.diff, setdiff))
}

#' merge two elias fano indices
#'
#' @param efdb.root the root index
#' @param efdb the index to be merged
#'
#' @importFrom hash values keys
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


contigency.table <- function(query.results)
{
    data <- as.data.frame(lapply(values(query.results), length))
    names(data) <-  keys(query.results)
    return(data)
    
}


phyper.test <- function(object, result, datasets)
{
    df <- query.result.as.dataframe(result)
    datasets <- select.datasets(object, datasets)
    
    cell.types.df <- aggregate(cell_id ~ cell_type, df, FUN = length)
    colnames(cell.types.df)[colnames(cell.types.df) == 'cell_id'] <- 'cell_hits'
    cell.types.df$total_cells<- object@index$getCellTypeSupport(cell.types.df$cell_type)
    query.hits <- nrow(df)
    print(df)
    
    cell.types.df$pval <- p.adjust(1 - phyper(cell.types.df$cell_hits, # total observed successes ( query.hits for cell type)
                                 cell.types.df$total_cells, # total successes ( cell type size )
                                 sum(cell.types.df$total_cells) - cell.types.df$total_cells, # total failures( total cells excluding cell type)
                                 query.hits # sample size 
                                 ), n = object@index$numberOfCellTypes(datasets))


    return(cell.types.df)

}

query.result.as.dataframe <- function(query.result)
{
    if (is.data.frame(query.result))
    {
        return(query.result)
    }
    if (length(query.result) == 0)
    {
        return(data.frame(cell_type = c() , cell_id = c()))
    }
    else
    {
        result <- setNames(unlist(query.result, use.names=F), rep(names(query.result), lengths(query.result)))
        return(data.frame(cell_type = names(result), cell_id = result))
    }            
    
}

select.datasets <- function(object, datasets)
{
    
    if (missing(datasets))
    {
        ## Select all available datasets
        datasets <- object@datasets
    }
    else
    {
        ## datasets should not be a superset of the data
        if(length(setdiff(datasets, object@datasets)) != 0)
        {
            stop(paste("Dataset", setdiff(datasets,object@datasets), "does not exist in the database"))
        }
    }
    return(datasets)

}


