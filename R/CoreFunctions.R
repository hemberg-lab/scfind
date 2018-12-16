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

#' @importFrom stats aggregate p.adjust phyper setNames
phyper.test <- function(object, result, datasets)
{
    df <- query.result.as.dataframe(result)
    datasets <- select.datasets(object, datasets)
    
    cell.types.df <- aggregate(cell_id ~ cell_type, df, FUN = length)
    colnames(cell.types.df)[colnames(cell.types.df) == 'cell_id'] <- 'cell_hits'
    cell.types.df$total_cells<- object@index$getCellTypeSupport(cell.types.df$cell_type)
    query.hits <- nrow(df)
    
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



scfind.get.genes.in.db <- function(object){
    
    return(object@index$genes())

}


