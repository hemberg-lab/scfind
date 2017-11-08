#' Build a cell type Index
#' 
#' Calculates a fraction of expressed cells per gene per cell type
#'
#' @param object object of SingleCellExperiment class
#' @param cell_type_column column name in the colData slot of the object SingleCellExperiment 
#' containing the cell classification information
#' 
#' @name buildCellTypeIndex
#'
#' @return a `data.frame` containing calculated gene index
#'
#' @importFrom SummarizedExperiment rowData colData colData<-
#' @importFrom SingleCellExperiment logcounts
#' @importFrom dplyr group_by summarise %>%
#' @importFrom reshape2 melt dcast
buildCellTypeIndex.SCESet <- function(object, cell_type_column) {
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

#' @rdname buildCellTypeIndex
#' @aliases buildCellTypeIndex
#' @importFrom SingleCellExperiment SingleCellExperiment
setMethod("buildCellTypeIndex", "SingleCellExperiment", buildCellTypeIndex.SCESet)

#' Find cell types associated with a given gene list
#' 
#' Calculates p-values of a log-likelihood of a list of genes to be associated
#' with each cell type. Log-likelihood is based on gene expression values.
#'
#' @param gene_index a data.frame with cell types in columns and genes in rows
#' @param gene_list genes that need to be searched in the gene_index
#' 
#' @name findCellType
#'
#' @return a named numeric vector containing p-values
#'
#' @importFrom stats pchisq
#' @importFrom methods is
findCellType.data.frame <- function(gene_index, gene_list) {
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
    
    p0 <- colSums(gene_index)/nrow(gene_index)
    
    if (length(gene_list[!gene_list %in% rownames(gene_index)]) != 0) {
        warning(paste0("Genes: ", paste(gene_list[!gene_list %in% rownames(gene_index)], collapse = ", "), 
            " were exluded from search since they are not present in the Gene Index!"))
        gene_list <- gene_list[gene_list %in% rownames(gene_index)]
    }
    
    if (length(gene_list) == 0) {
        stop("None of the genes in the gene_list are present in the gene_index!")
    }
    
    gene_index <- gene_index[gene_list, ]
    lambda <- 2 * log(apply(gene_index, 2, prod)/p0^(nrow(gene_index)))
    lambda[is.na(lambda)] <- NA
    lambda[is.infinite(lambda)] <- NA
    p_values <- pchisq(lambda, length(gene_list), lower.tail = FALSE)
    return(p_values)
}

#' @rdname findCellType
#' @aliases findCellType
setMethod("findCellType", "data.frame", findCellType.data.frame)

#' Build a cell Index
#' 
#' Creates a compressed cell Index
#'
#' @param object object of SingleCellExperiment class
#' containing the cell classification information
#' @param cell_type_column column name in the colData slot of the object SingleCellExperiment 
#' containing the cell classification information
#' 
#' @name buildCellIndex
#'
#' @return a `data.frame` containing calculated gene index
#' @useDynLib scfind
#' @importFrom hash hash
#' @importFrom bit as.bit
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment colData
#' @import Rcpp
buildCellIndex.SCESet <- function(object, cell_type_column) {
    if (is.null(object)) {
        stop("Please define a object using the `object` parameter!")
    }
    if (is.null(colData(object)[[cell_type_column]])) {
        stop("Please define a correct `cell_type_column` in the `colData` slot!")
    }
    gene_exprs <- logcounts(object) > 0
    l <- as.numeric(floor(log2(ncol(object)/rowSums(gene_exprs))))
    filter <- !is.infinite(l)
    l <- l[filter]
    gene_exprs <- gene_exprs[filter, ]
    inds <- lapply(apply(gene_exprs, 1, which), as.numeric)
    p0 <- vapply(unique(colData(object)[[cell_type_column]]), function(ct) {
        sum(gene_exprs[, colData(object)[[cell_type_column]] == ct])/
        (nrow(gene_exprs) * length(which(colData(object)[[cell_type_column]] == ct)))
    }, numeric(1))
    names(p0) <- unique(colData(object)[[cell_type_column]])
    f_symbs <- rowData(object)$feature_symbol[filter]
    rownames(gene_exprs) <- f_symbs
    codes <- eliasFanoCoding(inds, l)
    res <- Map(list, H = lapply(codes$H, as.bit), L = lapply(codes$L, as.bit), l = l)
    index <- hash(f_symbs, res)
    return(list(index = index, cell_types = colData(object)[[cell_type_column]], p0 = p0))
}

#' @rdname buildCellIndex
#' @aliases buildCellIndex
#' @importFrom SingleCellExperiment SingleCellExperiment
setMethod("buildCellIndex", "SingleCellExperiment", buildCellIndex.SCESet)

#' Find cells associated with a given gene list
#' 
#' Calculates p-values of a log-likelihood of a list of genes to be associated
#' with each cell type. Log-likelihood is based on gene expression values.
#'
#' @param input object of SingleCellExperiment class
#' @param genelist column name in the colData slot of the object SingleCellExperiment 
#' containing the cell classification information
#' @param statistics defines statistics to be used to calculate log-likelihood.
#' 'G' is the default. The second option is 'chisq'.
#' 
#' @name findCell
#'
#' @return a `list` containing calculated gene index
#' @useDynLib scfind
#' @import Rcpp
findCell.SCESet <- function(input, genelist, statistics) {
    if (is.null(input)) {
        stop("Please define an input parameter!")
    }
    if (is.null(genelist)) {
        stop("Please define a list of genes using the `genelist` parameter!")
    }
    if (!"list" %in% is(input)) {
        stop("The gene_index must be a list!")
    }
    if (!"character" %in% is(genelist)) {
        stop("The genelist must be a character vector!")
    }
    inds <- list()
    for (i in genelist) {
        tmp <- eliasFanoDecoding(as.numeric(input$index[[i]]$H), as.numeric(input$index[[i]]$L), 
            input$index[[i]]$l)
        inds[[i]] <- tmp
    }
    common_exprs_cells <- data.frame(cell_id = Reduce(intersect, inds), cell_type = input$cell_types[Reduce(intersect, 
        inds)])
    cell_types_p <- sapply(sapply(inds, function(x) {
        factor(input$cell_types[x], levels = unique(input$cell_types))
    }, simplify = FALSE), table)/as.vector(table(factor(input$cell_types, levels = unique(input$cell_types))))
    
    if(statistics == "G") {
        lambda <- 2 * apply(cell_types_p * log(cell_types_p / input$p0), 1, sum)
    } else {
        lambda <- 2 * apply(log(cell_types_p / input$p0), 1, sum)
    }
    lambda[is.na(lambda)] <- NA
    lambda[is.infinite(lambda)] <- NA
    p_values <- pchisq(lambda, length(genelist), lower.tail = FALSE)
    return(list(p_values = p_values, common_exprs_cells = common_exprs_cells))
}

#' @rdname findCell
#' @aliases findCell
setMethod("findCell", "list", findCell.SCESet)
