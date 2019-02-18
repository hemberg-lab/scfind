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

caseCorrect <- function(object, gene.list)
{
    gene.list <- gene.list[gene.list != ""]

    if(length(gene.list) != 0)
    {
        gene.corr <- object@index$genes()[match(tolower(gene.list), tolower(object@index$genes()), nomatch = 0)]
        
        if(length(setdiff(tolower(gene.list), tolower(gene.corr))) != 0) message(paste0("Ignored ", toString(setdiff(gene.list, gene.list[match(tolower(gene.corr), tolower(gene.list), nomatch=0)])), ". Not found in the index"))

        return(unique(gene.corr))
    }
    else
    {
        return(c())
    }
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

pair.id <- function(cell.list = list()){
    if(length(cell.list) == 0) 
    {
        return(c())
    } 
    else
    {
        pair.vec <- stack(cell.list)
        return (paste0(pair.vec$ind, "#",pair.vec$values))
    }
    
}

find.signature <- function(object, cell.type, max.genes=1000, min.cells=10, max.pval=0) {
    # Use this method to find a gene signature for a cell-type. 
    # We do this by ranking genes by recall and then adding genes to the query until we exceed a target p-value threshold or until a minimum number of cells is returned from the query
    df <- cellTypeMarkers(object, cell.type, top.k=max.genes, sort.field="recall", message=F)
    genes <- as.character(df$genes)
    genes.list <- c()
    thres = max(c(min.cells, object@index$getCellTypeMeta(cell.type)$total_cells))
    for (j in 1:dim(df)[1]) {
        res <- hyperQueryCellTypes(object, c(genes.list, genes[j]))
        if (dim(res)[1]>0) {
            ind <- which(res[,1]==cell.type)
            if (length(ind)==0) {
                break
            }
            else {
                if (res[ind,4]>max.pval | res[ind,2]<thres) {
                    break
                }
            }
        }
        genes.list <- c(genes.list, genes[j])
    }
    return( genes.list )
}
