#' @export
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
#' # this is needed to calculate dropout rate for feature selection
#' # important: normcounts have the same zeros as raw counts (fpkm)
#' counts(sce) <- normcounts(sce)
#' logcounts(sce) <- log2(normcounts(sce) + 1)
#' # use gene names as feature symbols
#' rowData(sce)$feature_symbol <- rownames(sce)
#' isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' index <- buildCellTypeIndex(sce)
#' 
setGeneric("buildCellTypeIndex", signature = "object", function(object = NULL, cell_type_column = "cell_type1") {
    standardGeneric("buildCellTypeIndex")
})

#' @export
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
#' # this is needed to calculate dropout rate for feature selection
#' # important: normcounts have the same zeros as raw counts (fpkm)
#' counts(sce) <- normcounts(sce)
#' logcounts(sce) <- log2(normcounts(sce) + 1)
#' # use gene names as feature symbols
#' rowData(sce)$feature_symbol <- rownames(sce)
#' isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' index <- buildCellTypeIndex(sce)
#' res <- findCellType(index, gene_list = c('SOX6', 'SNAI3'))
#' 
setGeneric("findCellType", signature = "gene_index", function(gene_index = NULL, gene_list = NULL) {
    standardGeneric("findCellType")
})

#' @export
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
#' # this is needed to calculate dropout rate for feature selection
#' # important: normcounts have the same zeros as raw counts (fpkm)
#' counts(sce) <- normcounts(sce)
#' logcounts(sce) <- log2(normcounts(sce) + 1)
#' # use gene names as feature symbols
#' rowData(sce)$feature_symbol <- rownames(sce)
#' isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' index <- buildCellIndex(sce)
#' 
setGeneric("buildCellIndex", signature = "object", function(object = NULL, cell_type_column = "cell_type1") {
    standardGeneric("buildCellIndex")
})

#' @export
#' 
#' @examples
#' library(SingleCellExperiment)
#' sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
#' # this is needed to calculate dropout rate for feature selection
#' # important: normcounts have the same zeros as raw counts (fpkm)
#' counts(sce) <- normcounts(sce)
#' logcounts(sce) <- log2(normcounts(sce) + 1)
#' # use gene names as feature symbols
#' rowData(sce)$feature_symbol <- rownames(sce)
#' isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' index <- buildCellIndex(sce)
#' res <- findCell(index, genelist = c('SOX6', 'SNAI3'))
#' 
setGeneric("findCell", signature = "input", function(input = NULL, genelist = NULL,
                                                     statistics = "G") {
    standardGeneric("findCell")
})
