


#' Gets a list of gene hits results
#' 
operator <-  function(hash.index1, hash.index2, operator.function)
{    cell.types <- operator.function(keys(hash.index1), keys(hash.index2))
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

#' Applies the and operator on two result data structures
#'
#' @param hash.index1 result index to be intersected with hash.index2
#' @param hash.index2 result index to be intersected with hash.index1
#'
#' @return the intersected set
and.operator <- function(hash.index1, hash.index2)
{
    return(operator(hash.index1, hash.index2, intersect))
}


#' Applies the union operator on two result data structures
#'
#' @param hash.index1 result index to be unionized with hash.index2
#' @param hash.index2 result index to be unionized with hash.index1
#'
#' @return the unionized set
or.operator <- function(hash.index1, hash.index2)
{
    return(operator(hash.index1, hash.index2, union))
}

#' Remove cells that have some specific gene hits that we do not want in our list
#'
#' @param hash.index the result index that cell types will be filtered
#' @param hash.diff the cell types to be removed in case there is an overlap
#'
#' @return an new filtered index
not.operator <-  function(hash.index, hash.diff)
{
    return(operator(hash.index, hash.diff, setdiff))
}
