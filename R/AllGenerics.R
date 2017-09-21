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
#' isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' gene_index <- buildGeneIndex(sce)
#' 
setGeneric("buildGeneIndex", signature = "object", function(object = NULL, cell_type_column = "cell_type1") {
    standardGeneric("buildGeneIndex")
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
#' isSpike(sce, "ERCC") <- grepl("^ERCC-", rownames(sce))
#' # remove features with duplicated names
#' sce <- sce[!duplicated(rownames(sce)), ]
#' gene_index <- buildGeneIndex(sce)
#' res <- queryGeneList(gene_index, gene_list = c("ELMO2", "PNMA1"))
#' 
setGeneric("queryGeneList", signature = "gene_index", function(gene_index = NULL, gene_list = NULL) {
    standardGeneric("queryGeneList")
})
