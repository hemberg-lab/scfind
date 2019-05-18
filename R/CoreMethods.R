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
    index <- new("SCFind", index = ef, datasets = dataset.name)
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
#' @param message turn on/off messsage
#'
#' @return a data.frame that each row represent a gene score for a specific cell type 
cell.type.marker <- function(object, cell.types, background.cell.types, top.k, sort.field, message = T)
{
    if (missing(background.cell.types))
    {
        if(message == T) message("Considering the whole DB..")
        background.cell.types <- cellTypeNames(object)
    }
    all.cell.types <- object@index$cellTypeMarkers(cell.types, background.cell.types)
    if (!(sort.field %in% colnames(all.cell.types)))
    {
        if(message == T) message(paste("Column", sort.field, "not found"))
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
get.cell.types.names <- function(object, datasets)
{
    if(missing(datasets)) 
    {
        return(object@index$getCellTypes())
    }
    else
    {
        return(object@index$getCellTypes()[lapply(strsplit(object@index$getCellTypes(), "\\."), `[[`, 1) %in% datasets])
    }
    
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
#' @param gene.list genes to be searched in the gene.index 
#' (Operators: "-gene" to exclude a gene | "*gene" either gene is expressed
#' "*-gene" either gene is expressed to be excluded)
#' @param datasets the datasets vector that will be tested as background for the hypergeometric test
#'
#' @return a DataFrame that contains all cell types with the respective cell cardinality and the hypergeometric test
cell.types.phyper.test <- function(object, gene.list, datasets)
{
    
    result <- findCellTypes.geneList(object, gene.list, datasets)
    
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
#' (Operators: "-gene" to exclude a gene | "*gene" either gene is expressed
#' "*-gene" either gene is expressed to be excluded)
#' @param datasets the datasets that will be considered
#' 
#' @name findCellTypes
#' @return a named numeric vector containing p-values
findCellTypes.geneList <- function(object, gene.list, datasets)
{
    datasets <- if(missing(datasets)) object@datasets else select.datasets(object, datasets)
    
    if(length(grep("^-|^\\*", gene.list)) == 0)
    {
        return(object@index$findCellTypes(caseCorrect(object, gene.list), datasets))
    }
    else
    {
        pos <- caseCorrect(object, grep("^[^-\\*]", gene.list, value = T))
        excl.or <- grep("^-\\*|^\\*-", gene.list, value = T)
        or <- caseCorrect(object, sub("\\*", "", setdiff(grep("^\\*", gene.list, value = T), excl.or)))
        excl <- caseCorrect(object, sub("-", "", setdiff(grep("^-", gene.list, value = T), excl.or)))
        excl.or <- caseCorrect(object, sub("\\*-||-\\*", "", grep("^-\\*|^\\*-", gene.list, value = T)))
        
        if(length(c(intersect(pos, or), intersect(pos, excl), intersect(pos, excl.or), intersect(or, excl), intersect(or, excl.or), intersect(excl, excl.or))) != 0)
        {
            message ("Warning: Same gene labeled with different operators!") 
            message ("There is a priority to handle operators:")
            message (paste("Cells with", paste(pos, collapse=" ∧ "),"expression will be included.", 
                           if(length(or) != 0) "Then cells with", paste(or, collapse=" ∨ "), "expression will be included."))
            message (paste("The result will be excluded by", paste(excl, collapse=" ∧ "), 
                           if(length(excl.or != 0)) paste("and further be excluded by", paste(excl.or, collapse=" ∨ "))))
            cat('\n')
        }

        
        cell.to.id  <- NULL
        
        # Using pair.id to create unique variable for each cell by pairing cell types to cell ID
        
        if(length(pos) == 0 && length(or) == 0 && (length(excl) != 0 || length(excl.or) != 0)) 
        {
            # When no positive selection, include all cells
            cell.to.id <- lapply(as.list(object@index$getCellTypeSupport(cellTypeNames(object, datasets))), seq)
            names(cell.to.id) <- cellTypeNames(object, datasets)
            cell.to.id <- pair.id(cell.to.id)
        } 


        if(length(or) != 0)
        {
            # Include any cell expresses gene in OR condition
            gene.or <- c()
            for(i in 1: length(or))
            {
                tmp.id <- pair.id(object@index$findCellTypes(c(pos, or[i]), datasets))
                if(length(pos) != 0 && !is.null(tmp.id)) message(paste("Found", length(tmp.id), if(length(tmp.id) > 1)"cells" else "cell", "co-expressing", paste(c(pos, or[i]), collapse=" ∧ ") ))
                if(!is.null(tmp.id))
                {
                    cell.to.id <- unique(c(cell.to.id, tmp.id))
                    # Store used query
                    gene.or <- c(gene.or, or[i])
                }
                else
                {
                    cell.to.id <- cell.to.id
                }
            }
                if( length(pos) == 0 && length(gene.or) != 0) message(paste("Found", length(cell.to.id), if(length(cell.to.id) > 1) "cells" else "cell", "expressing", paste(gene.or, collapse=" ∨ ")))
        }
        else
        {
            cell.to.id  <- if(length(pos) != 0) pair.id(object@index$findCellTypes(pos, datasets)) else cell.to.id
            if(length(pos) != 0) message(paste("Found", length(cell.to.id), if(length(pos) > 1) "cells co-expressing" else "cell expressing", paste(pos, collapse = " ∧ ")))
        }
        
            count.cell <- length(cell.to.id)
            gene.excl <- NULL

            if(length(excl.or) != 0)
            {
                # Negative select cell in OR condition
                for(i in 1: length(excl.or))
                {
                    ex.tmp.id <- pair.id(object@index$findCellTypes(c(excl, excl.or[i]), datasets))
                    
                    message(paste("Excluded", sum(cell.to.id %in% ex.tmp.id),
                                  if(sum(cell.to.id %in% ex.tmp.id) > 1)"cells" else "cell",
                                  if(length(excl) != 0) paste("co-expressing", paste( c(excl, excl.or[i]), collapse=" ∧ ")) else paste("expressing", excl.or[i]) ))
                    
                    if(!is.null(ex.tmp.id))
                    {
                        cell.to.id <- setdiff(cell.to.id, ex.tmp.id)
                        gene.excl <- c(gene.excl, excl.or[i])
                    }
                    else
                    {
                        cell.to.id <- cell.to.id
                    }
                }
                count.cell <- count.cell - length(cell.to.id)
                if(count.cell > 0 && length(gene.excl) == 0) message("Excluded", count.cell, if(count.cell > 1) "cells" else "cell", "expressing", paste(excl, collapse=" ∧ "))    
            }
            else
            {
                if(length(excl) != 0)
                {
                    # Negative selection
                    cell.to.id <- setdiff(cell.to.id, pair.id(object@index$findCellTypes(excl, datasets)))
                    count.cell <- count.cell - length(cell.to.id)
                    if(count.cell > 0) message(paste("Excluded", count.cell, if(count.cell > 1) "cells" else "cell", if(length(excl) > 1) "co-expressing" else "expressing", paste(excl, collapse = " ∧ "))) else message("No Cell Is Excluded!")
                }
            }
            
            # Generate a new list
            df <- do.call(rbind, strsplit(as.character(cell.to.id), "#"))
            if(!is.null(df))
            {
                result <- as.list(setNames(as.numeric(split(df[,2], seq(nrow(df)))), df[,1]))
                if(length(unique(df[,1])) == nrow(df)) 
                {
                    return(result)
                } 
                else 
                {

                    if(length(unique(names(result))) == 1) 
                    {
                        tmp <-list(stack(result)$values)
                        names(tmp) <- unique(names(result))
                        return(tmp)
                    }  
                    else 
                    {
                        return(unstack(stack(result)))
                    }

                }
            }
            else
            {
                message("No Cell Is Found!")
                return(list())
            }
        }
}

#' @rdname findCellTypes
#' @aliases findCellTypes
setMethod("findCellTypes", 
          signature(object = "SCFind",
                    gene.list = "character"), 
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


#' Find out how many cell-types each gene is found
#' 
#' @param object the \code{SCFind} object
#' @param gene.list genes to be searched in the gene.index
#' @param datasets the datasets that will be considered
#' @param min.cells threshold of cell hit of a cell type
#' @param min.fraction portion of total cell as threshold
#' 
#' @name findCellTypeSpecificities
#' @return the list of number of cell type for each gene
cell.type.specificity <- function(object, gene.list, datasets, min.cells=10, min.fraction=.25)
{
    if(min.fraction >= 1 || min.fraction <= 0) stop("min.fraction reached limit, please use values > 0 and < 1.0.") else message("Calculating cell-types for each gene...")
    datasets <- if(missing(datasets)) object@datasets else select.datasets(object, datasets)
    if(missing(gene.list)) 
    {
        res <- object@index$geneSupportInCellTypes(object@index$genes(), datasets) 
    }
    else 
    {
        gene.list <- caseCorrect(object, gene.list)
        res <- object@index$geneSupportInCellTypes(gene.list, datasets)
    }
    
    res.tissue <- res
    names(res.tissue) <- gsub("\\.", "#", names(res.tissue))
    df <- cbind(stack(res), stack(unlist(res.tissue)))
    # df[,4] <- sub("^[^.]+\\.", "", df[,4])
    df[,1] <- object@index$getCellTypeSupport( sub("^[^.]+\\.", "", df[,4])) * min.fraction
    if(length(which(df[,1] < min.cells)) != 0) df[which(df[,1] < min.cells),1] <- min.cells
    if(nrow(df) != 0) df <- df[which(df[,3] > df[,1]),] else return(split(rep(0, length(gene.list)), gene.list))
    if(nrow(df) != 0) return(as.list(summary(df[,2], maxsum=nrow(df)))) else return(split(rep(0, length(gene.list)), gene.list))
}

#' @rdname findCellTypeSpecificities
#' @aliases findCellTypeSpecificities
setMethod("findCellTypeSpecificities", 
          signature(object = "SCFind"), 
          cell.type.specificity)


#' Find out how many tissues each gene is found
#' 
#' @param object the \code{SCFind} object
#' @param gene.list genes to be searched in the gene.index
#' @param min.cells threshold of cell hit of a tissue
#' 
#' @name findTissueSpecificities
#' @return the list of number of tissue for each gene
tissue.specificity <- function(object, gene.list, min.cells = 10)
{
    if(length(object@datasets) <= 1) stop("Index contains 1 dataset only.") else message("Calculating tissues for each gene...")
    if(missing(gene.list)) 
    {
        res  <- object@index$geneSupportInCellTypes(object@index$genes(), object@datasets)
    }
    else
    {
        gene.list <- caseCorrect(object, gene.list)
        res <- object@index$geneSupportInCellTypes(gene.list, object@datasets)
    }
    
    if(length(res) > 0) res.tissue <- res else return(split(rep(0, length(gene.list)), gene.list))
    names(res.tissue) <- gsub("\\.", "#", names(res.tissue))
    df <- cbind(stack(res), stack(unlist(res.tissue)))
    df[,5] <- gsub("^[^.]*\\.([^.]*)\\..*$","\\1",df[,4])
    df <- aggregate(df[,1], by=list(df[,5], df[,2]), FUN=sum)
    df <- df[which(df[,3] > min.cells),]
    
    if(nrow(df) != 0) return(as.list(summary(df[,2], maxsum=nrow(df)))) else return(split(rep(0, length(gene.list)), gene.list))
}

#' @rdname findTissueSpecificities
#' @aliases findTissueSpecificities
setMethod("findTissueSpecificities", 
          signature(object = "SCFind"), 
          tissue.specificity)

#' Find the set of genes that are ubiquitously expressed in a query of cell types
#' 
#' @param object the \code{SCFind} object
#' @param cell.types a list of cell types for the list to evaluated
#' @param min.recall threshold of minimun recall value
#' @param max.genes threshold of number of genes to be considered for each cell type
#' 
#' @name findHouseKeepingGenes
#' @return the list of gene that ubiquitously expressed in a query of cell types
#' 
house.keeping.genes <- function(object, cell.types, min.recall=.5, max.genes=1000) {
    if(min.recall >= 1 || min.recall <= 0) stop("min.recall reached limit, please use values > 0 and < 1.0.") 
    if(max.genes > length(object@index$genes())) stop(paste("max.genes exceeded limit, please use values > 0 and < ", length(object@index$genes()))) else message("Searching for house keeping genes...")
    df <- cellTypeMarkers(object, cell.types[1], top.k=max.genes, sort.field="recall", message=F)
    house.keeping.genes <- df$genes[which(df$recall>min.recall)]
    
    for (i in 2:length(cell.types)) {
        setTxtProgressBar(txtProgressBar(1, length(cell.types), style = 3), i) 
        df <- cellTypeMarkers(object, cell.types[i], top.k=max.genes, sort.field="recall", message=F)
        house.keeping.genes <- intersect(house.keeping.genes, df$genes[which(df$recall>min.recall)])
        if (length(house.keeping.genes)==0) { stop("No house keeping gene is found.") }
    }
    cat('\n')
    return( house.keeping.genes )
}


#' @rdname findHouseKeepingGenes
#' @aliases findHouseKeepingGenes
setMethod("findHouseKeepingGenes", 
          signature(object = "SCFind",
                    cell.types = "character"), 
          house.keeping.genes)

#'  Find the signature genes for a cell-type
#' 
#' @param object the \code{SCFind} object
#' @param cell.types a list of cell types for the list to evaluated
#' @param max.genes threshold of number of genes to be considered for each cell type
#' @param min.cells threshold of cell hit of a tissue
#' @param max.pval threshold of p-value
#' 
#' @name findGeneSignatures
#' @return the list of gene signatures in a query of cell types
#' 
gene.signatures <- function(object, cell.types, max.genes=1000, min.cells=10, max.pval=0) 
{
    message("Searching for gene signatures...")
    cell.types.all <- if(missing(cell.types)) object@index$getCellTypes() else cellTypeNames(object)[tolower(cellTypeNames(object)) %in% tolower(cell.types)]
    signatures <- list()
    if(length(cell.types.all) != 0)
    {
        for (i in 1:length(cell.types.all)) {
            if(i > 1) setTxtProgressBar(txtProgressBar(1, length(cell.types.all), style = 3), i)
            signatures[[cell.types.all[i]]] <- find.signature(object, cell.types.all[i], max.genes=max.genes, min.cells=min.cells, max.pval=max.pval)
        }
        cat('\n')
        return( signatures )
    } 
    else
    {
        return(message(paste0("Ignored ", toString(cell.types),". Cell type not found in index.")))
    }
}

#' @rdname findGeneSignatures
#' @aliases findGeneSignatures
setMethod("findGeneSignatures", 
          signature(object = "SCFind"), 
          gene.signatures)

#'  Look at all other genes and rank them based on the similarity of their expression pattern to the pattern defined by the gene query
#' 
#' @param object the \code{SCFind} object
#' @param gene.list genes to be searched in the gene.index
#' @param datasets the datasets that will be considered
#' @param top.k how many genes to retrieve
#' 
#' @name findSimilarGenes
#' @return the list of genes and their similarities presented in Jaccard indices
#' 
similar.genes <- function(object, gene.list, datasets, top.k=5) {
    message("Searching for genes with similar pattern...")
    datasets <- if(missing(datasets)) object@datasets else select.datasets(object, datasets)
    gene.list <- caseCorrect(object, gene.list)
    e <- object@index$findCellTypes(gene.list, datasets) #the cells expressing the genes in gene.list
    n.e <- length(unlist(e))
    if (n.e>0) {
        gene.names <- setdiff(object@index$genes(), gene.list)
        similarities <- rep(0, length(gene.names))
        ns <- rep(0, length(gene.names))
        ms <- rep(0, length(gene.names))
        for (i in 1:length(gene.names)) {
            setTxtProgressBar(txtProgressBar(1, length(gene.names), style = 3), i) 
            f <- object@index$findCellTypes(gene.names[i], datasets) #find expression pattern of other gene
            if (length(f)>0) {
                m <- rep(0, length(e))
                for (j in 1:length(names(e))) {
                    m[j] <- length(intersect(e[[j]], f[[names(e)[j]]]))
                }
                #calculate the Jaccard index for the similarity of the cells expressing the gene
                similarities[i] <- sum(m)/(n.e + length(unlist(f)) - sum(m))
                ns[i] <- length(unlist(f))
                ms[i] <- sum(m)            
            }
        }
        cat('\n')
        r <- sort(similarities, index.return=T)
        inds <- tail(r$ix, top.k)
        res <- data.frame("gene" = gene.names[inds], "Jaccard"=similarities[inds], "overlap"=ms[inds], "n"=ns[inds])
        return( res )
    }
    else
    {
        message(paste("Cannot find cell expressing", toString(gene.list), "in the index."))
        return( c() )
    }
}


#' @rdname findSimilarGenes
#' @aliases findSimilarGenes
setMethod("findSimilarGenes", 
          signature(object = "SCFind",
                    gene.list = "character"), 
          similar.genes)





