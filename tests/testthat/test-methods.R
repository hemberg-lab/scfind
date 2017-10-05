test_that("testing buildCellTypeIndex()", {
    sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
    # this is needed to calculate dropout rate for feature selection
    # important: normcounts have the same zeros as raw counts (fpkm)
    counts(sce) <- normcounts(sce)
    logcounts(sce) <- log2(normcounts(sce) + 1)
    # use gene names as feature symbols
    rowData(sce)$feature_symbol <- rownames(sce)
    isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
    # remove features with duplicated names
    sce <- sce[!duplicated(rownames(sce)), ]
    index <- buildCellTypeIndex(sce)
    expect_is(index, "data.frame")
    expect_equal(ncol(index), 6)
    expect_equal(nrow(index), 20214)
    expect_error(index <- buildCellTypeIndex(yan))
})

test_that("testing findCellType()", {
    sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
    # this is needed to calculate dropout rate for feature selection
    # important: normcounts have the same zeros as raw counts (fpkm)
    counts(sce) <- normcounts(sce)
    logcounts(sce) <- log2(normcounts(sce) + 1)
    # use gene names as feature symbols
    rowData(sce)$feature_symbol <- rownames(sce)
    isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
    # remove features with duplicated names
    sce <- sce[!duplicated(rownames(sce)), ]
    index <- buildCellTypeIndex(sce)
    res <- findCellType(index, gene_list = c('SOX6', 'SNAI3'))
    expect_is(res, "numeric")
    expect_equal(length(res), 6)
    expect_error(res <- findCellType(sce, gene_list = c('SOX6', 'SNAI3')))
    expect_error(res <- findCellType(index, gene_list = c(1, 2)))
})

test_that("testing buildCellIndex()", {
    sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
    # this is needed to calculate dropout rate for feature selection
    # important: normcounts have the same zeros as raw counts (fpkm)
    counts(sce) <- normcounts(sce)
    logcounts(sce) <- log2(normcounts(sce) + 1)
    # use gene names as feature symbols
    rowData(sce)$feature_symbol <- rownames(sce)
    isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
    # remove features with duplicated names
    sce <- sce[!duplicated(rownames(sce)), ]
    index <- buildCellIndex(sce)
    expect_is(index, "list")
    expect_is(index$index, "hash")
    expect_is(index$cell_types, "factor")
    expect_is(index$p0, "numeric")
    expect_equal(length(index), 3)
    expect_equal(length(index$cell_types), 90)
    expect_equal(length(index$p0), 6)
    expect_error(index <- buildCellIndex(yan))
})

test_that("testing findCell()", {
    sce <- SingleCellExperiment(assays = list(normcounts = as.matrix(yan)), colData = ann)
    # this is needed to calculate dropout rate for feature selection
    # important: normcounts have the same zeros as raw counts (fpkm)
    counts(sce) <- normcounts(sce)
    logcounts(sce) <- log2(normcounts(sce) + 1)
    # use gene names as feature symbols
    rowData(sce)$feature_symbol <- rownames(sce)
    isSpike(sce, 'ERCC') <- grepl('^ERCC-', rownames(sce))
    # remove features with duplicated names
    sce <- sce[!duplicated(rownames(sce)), ]
    index <- buildCellIndex(sce)
    res <- findCell(index, genelist = c('SOX6', 'SNAI3'))
    expect_is(res, "list")
    expect_equal(length(res), 2)
    expect_equal(length(res$p_values), 6)
    expect_is(res$p_values, "numeric")
    expect_is(res$common_exprs_cells, "data.frame")
    expect_error(res <- findCell(sce, genelist = c('SOX6', 'SNAI3')))
    expect_error(res <- findCell(index, genelist = c(1, 2)))
})
