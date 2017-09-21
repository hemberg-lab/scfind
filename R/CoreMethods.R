#' Build gene index for a dataset
#' 
#' Calculates a fraction of expressed cells per gene per cell type
#'
#' @param object object of SingleCellExperiment class
#' @param cell_type_column column name in the colData slot of the object SingleCellExperiment 
#' containing the cell classification information
#' 
#' @name buildGeneIndex
#'
#' @return a `data.frame` containing calculated gene index
#'
#' @importFrom SummarizedExperiment rowData colData colData<-
#' @importFrom dplyr group_by summarise %>%
#' @importFrom reshape2 melt dcast
buildGeneIndex.SCESet <- function(object, cell_type_column) {
    if (is.null(object)) {
        stop("Please define a object using the `object` parameter!")
    }
    if (is.null(colData(object)[[cell_type_column]])) {
        stop("Please define a correct `cell_type_column` in the `colData` slot!")
    }
    gene <- cell_class <- exprs <- NULL
    object_local <- logcounts(object) > 0
    rownames(object_local) <- rowData(object)$feature_symbol
    colnames(object_local) <- colData(object)[[cell_type_column]]
    
    # calculate median feature expression in every cell class of object
    object_local <- reshape2::melt(object_local)
    colnames(object_local) <- c("gene", "cell_class", "exprs")
    object_local <- object_local %>% group_by(gene, cell_class) %>% summarise(gene_exprs_prob = sum(exprs)/length(exprs))
    object_local <- reshape2::dcast(object_local, gene ~ cell_class, value.var = "gene_exprs_prob")
    rownames(object_local) <- object_local$gene
    object_local <- object_local[, 2:ncol(object_local), drop = FALSE]
    return(object_local)
}

#' @rdname buildGeneIndex
#' @aliases buildGeneIndex
setMethod("buildGeneIndex", "SingleCellExperiment", buildGeneIndex.SCESet)

#' Find a cell type associated with a given gene list
#' 
#' Calculates p-values of a log-likelihood of a list of genes to be associated
#' with each cell type
#'
#' @param gene_index a data.frame with cell types in columns and genes in rows
#' @param gene_list genes that need to be searched in the gene_index
#' 
#' @name queryGeneList
#'
#' @return a named numeric vector containing p-values
#'
#' @importFrom stats pchisq
queryGeneList.data.frame <- function(gene_index, gene_list) {
    if (is.null(gene_index)) {
        stop("Please define a gene_index using the `gene_index` parameter!")
    }
    if (is.null(gene_list)) {
        stop("Please define a list of genes using the `gene_list` parameter!")
    }
    if (!"data.frame" %in% is(gene_index)) {
        stop("The gene_index must be a data.frame!")
    }
    if (!"character" %in% is(gene_list)) {
        stop("The gene_list must be a character vector!")
    }
    
    p0 <- colSums(gene_index) / nrow(gene_index)
    
    if (length(gene_list[!gene_list %in% rownames(gene_index)]) != 0) {
        warning(paste0("Genes: ", paste(gene_list[!gene_list %in% rownames(gene_index)], collapse = ", "), 
                       " were exluded from search since they are not present in the Gene Index!"))
        gene_list <- gene_list[gene_list %in% rownames(gene_index)]
    }
    
    if (length(gene_list) == 0) {
        stop("None of the genes in the gene_list are present in the gene_index!")
    }
    
    gene_index <- gene_index[gene_list, ]
    lambda <- 2 * log(apply(gene_index, 2, prod) / p0^(nrow(gene_index)))
    lambda[is.na(lambda)] <- NA
    lambda[is.infinite(lambda)] <- NA
    
    # compare the probability that each cell type expresses the genes of 
    # interest vs the baseline. Assuming that the genes are uncorrelated, 
    # we may calculate the probability of observing them all as a series of 
    # Bernoulli trials.
    # Return the p-value (not corrected for multiple testing)
    p_values <- pchisq(lambda, length(gene_list), lower.tail = F)
    return(p_values)
}

#' @rdname queryGeneList
#' @aliases queryGeneList
setMethod("queryGeneList", "data.frame", queryGeneList.data.frame)


#' @useDynLib scfind
#' @importFrom hash hash has.key
#' @importFrom bit as.bit
#' 
#' @export
parseDataset <- function(file.name, dataset.name=c()) {
    
    print(paste("Reading", file.name))
    d <- readRDS(file.name)
    size.orig <- object.size(d)
    if (length(dataset.name)==0) {
        #This bit needs to be made more robust to allow for more flexibility in how file-names are chosen
        tmp <- strsplit(strsplit(file.name, "\\.")[[1]], "/")[[1]]
        dataset.name <- tmp[length(tmp)]
    }
    genenames <- unique(rowData(d)$feature_symbol)
    inds.genes <- which(rowData(d)$feature_symbol %in% genenames)
    n.genes <- length(inds.genes)
    
    
    #Group cells by cell-type
    cell.types <- unique(colData(d)$cell_type1)
    #For keeping track of the parameters for the Elias-Fano code. Use l specific to each gene and cell-type
    m <- matrix(rep(0, n.genes*length(cell.types)), nrow=n.genes)
    #Keep track of how the cells were permuted
    permutation <- c()
    p0 <- c()
    #Keep track of the cluster for each cell. Use a numeric index to save space
    assigned.clusters <- rep(0, dim(d)[2])
    for (i in 1:dim(d)[2]) { assigned.clusters[i] <- which(cell.types==colData(d)$cell_type1[i]) }
    #Keep track of the number of cells in the cluster
    cluster.sizes <- rep(0, length(cell.types))
    cluster.names <- rep(0, length(cell.types))
    #Keep track of where the next cell-type starts for the H array for faster access
    H.start <- matrix(rep(1, n.genes*length(cell.types)), nrow=n.genes)
    #Use a hash to keep track of the codes from each gene
    H <- hash()
    L <- hash()
    exprs <- logcounts(d)
    for (i in 1:length(cell.types)) {
        inds.cell <- which(cell.types[i]==colData(d)$cell_type1)
        u <- length(inds.cell)
        cluster.sizes[i] <- u
        cluster.names[i] <- paste0(dataset.name, "_", cell.types[i])
        permutation <- c(permutation, inds.cell)
        #Calculate the baseline probability that a gene will be expressed in a cell
        p0 <- c(p0, sum(exprs[,inds.cell]>0)/(n.genes*u))
        for (j in 1:n.genes) {
            if (i>1) { H.start[j,i] <- length(H[[genenames[j]]]) + 1 }
            inds.exprs <- which(exprs[inds.genes[j],inds.cell]>0)
            m[j,i] <- length(inds.exprs)
            if (m[j,i]>0) {
                l <- floor(log2(u/m[j,i]))
                C <- eliasFanoCodingCpp(inds.exprs, l) #eliasFanoCoding(inds.exprs, l)
                if (has.key(genenames[j], H)) {
                    H[[genenames[j]]] <- c(H[[genenames[j]]], as.bit(C$H))
                    L[[genenames[j]]] <- c(L[[genenames[j]]], as.bit(C$L))
                } else {
                    H[[genenames[j]]] <- as.bit(C$H)
                    L[[genenames[j]]] <- as.bit(C$L)
                }
            }
        }
    }
    return( list( "experiment.name"=dataset.name, "assigned.clusters"=assigned.clusters, "m"=m, "permutation"=permutation, "genenames"=genenames, "cluster.sizes"=cluster.sizes, "p0"=p0, "n.clusters"=length(cluster.sizes), "n.cells"=length(permutation), "cluster.names"=cluster.names, "H.start"=H.start, "H"=H, "L"=L, "size.orig"=size.orig ) )
}

#' @importFrom hash hash has.key
#' 
#' @export
mergeDatasets <- function(efdb1, efdb2) {
    #For now, we are going to only use those genes that are present in both the old and the new
    efdb <- c()
    efdb$H <- hash()
    efdb$L <- hash()
    genenames.intersect <- intersect(efdb1$genenames, efdb2$genenames)
    #inds.efdb1 <- which(efdb1$genenames %in% genenames.intersect)
    #inds.efdb2 <- which(efdb2$genenames %in% genenames.intersect)
    #efdb$m <- cbind(efdb1$m[inds.efdb1,], efdb2$m[inds.efdb2,])
    efdb$cluster.sizes <- c(efdb1$cluster.sizes, efdb2$cluster.sizes)
    efdb$n.clusters <- c(efdb1$n.clusters, efdb2$n.clusters)
    efdb$n.cells <- c(efdb1$n.cells, efdb2$n.cells)
    efdb$H.start <- matrix(rep(0, length(genenames.intersect)*sum(efdb$n.clusters)), ncol=sum(efdb$n.clusters))
    efdb$m <- matrix(rep(0, length(genenames.intersect)*sum(efdb$n.clusters)), ncol=sum(efdb$n.clusters))
    efdb$assigned.clusters <- c(efdb1$assigned.clusters, sum(efdb1$n.clusters) + efdb2$assigned.clusters)
    efdb$cluster.names <- c(efdb1$cluster.names, efdb2$cluster.names)
    efdb$permutation <- c(efdb1$permutation, efdb2$permutation)
    efdb$p0 <- c(efdb1$p0, efdb2$p0)
    efdb$experiment.name <- c(efdb1$experiment.name, efdb2$experiment.name)
    efdb$size.orig <- c(efdb1$size.orig, efdb2$size.orig)
    efdb$genenames <- genenames.intersect
    for (i in 1:length(genenames.intersect)) {
        ind.efdb1 <- which(efdb1$genenames==genenames.intersect[i])
        ind.efdb2 <- which(efdb2$genenames==genenames.intersect[i])
        efdb$m[i,] <- c(efdb1$m[ind.efdb1,], efdb2$m[ind.efdb2,])
        if (has.key(genenames.intersect[i], efdb1$H) & has.key(genenames.intersect[i], efdb2$H)) {
            efdb$H[[genenames.intersect[i]]] <- c(efdb1$H[[genenames.intersect[i]]], efdb2$H[[genenames.intersect[i]]])
            efdb$L[[genenames.intersect[i]]] <- c(efdb1$L[[genenames.intersect[i]]], efdb2$L[[genenames.intersect[i]]])
            efdb$H.start[i,] <- c(efdb1$H.start[ind.efdb1,], length(efdb1$H[[genenames.intersect[i]]]) + efdb2$H.start[ind.efdb2,])
        }
        else if (has.key(genenames.intersect[i], efdb1$H) & !has.key(genenames.intersect[i], efdb2$H)) {
            efdb$H[[genenames.intersect[i]]] <- efdb1$H[[genenames.intersect[i]]]
            efdb$L[[genenames.intersect[i]]] <- efdb1$L[[genenames.intersect[i]]]
            efdb$H.start[i,] <- c(efdb1$H.start[ind.efdb1,], rep(length(efdb1$H[[genenames.intersect[i]]]), sum(efdb2$n.clusters)))
        }
        else if (!has.key(genenames.intersect[i], efdb1$H) & has.key(genenames.intersect[i], efdb2$H)) {
            efdb$H[[genenames.intersect[i]]] <- efdb2$H[[genenames.intersect[i]]]
            efdb$L[[genenames.intersect[i]]] <- efdb2$L[[genenames.intersect[i]]]
            efdb$H.start[i,] <- c(rep(1, sum(efdb1$n.clusters)), efdb2$H.start[ind.efdb2,])
        }
    }
    return( efdb )
}


