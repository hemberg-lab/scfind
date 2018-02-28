


plotCompressions <- function() {
    require("ggplot2")
    efdb.hs <- readRDS("efdb_hs.rds")
    efdb.mm <- readRDS("efdb_mm.rds")
    dat <- data.frame(compression=c(efdb.mm$size.orig/efdb.mm$size.compressed, efdb.hs$size.orig/efdb.hs$size.compressed), n=c(efdb.mm$n.cells, efdb.hs$n.cells))
    png(filename="Images/scfind_compression_scatter.png")
    gg <- ggplot(data=dat, aes(x=n, y=compression)) + geom_point() + scale_y_log10() + scale_x_log10() + ylab("Compression") + xlab("#cells") + theme_minimal() + theme(axis.title.y = element_text(size=20), axis.title.x = element_text(size=20), text = element_text(size=14))
    print(gg)
    dev.off()

}

plotSearchTimes <- function(max.genes=8, n.reps=10, rejection.sampling=T) {
    require("ggplot2")
    efdb.hs <- readRDS("efdb_hs.rds")
    efdb.mm <- readRDS("efdb_mm.rds")
    efdb.tenx <- readRDS("/lustre/scratch117/cellgen/team218/MH/efdb_tenx.rds")
    efdbs <- list(efdb.hs, efdb.mm, efdb.tenx)
    efdb.names <- c("Hs", "Mm", "1M")
    dat <- data.frame()
    ns <- c("efdb", "n.genes", "t.mean", "t.stderr")
    #names(dat) <- ns
    for (i in 1:length(efdbs)) {
        print(paste0("Considering dataset ", efdb.names[i]))
        for (j in 1:max.genes) {
            print(paste0("Searching for ", j))
            #Draw random genes, make sure that they have expression in at least one of the remaining cells
            ts <- testIndexSpeed(efdbs[[i]], j, n.reps, rejection.sampling)
            dat2 <- data.frame(efdb.names[i], j, mean(ts), sd(ts)/sqrt(n.reps))
            names(dat2) <- ns
            dat <- rbind(dat, dat2)
        }
    }
    png(filename="Images/scfind_search_times.png")
    gg <- ggplot(data=dat, aes(x=n.genes, y=t.mean, colour=efdb)) + geom_errorbar(aes(ymin=t.mean-t.stderr, ymax=t.mean+t.stderr), width=.1) + geom_line() + geom_point() + ylab("Time (s)") + xlab("#genes") + theme_minimal() + theme(axis.title.y = element_text(size=20), axis.title.x = element_text(size=20), text = element_text(size=14))
    print(gg)
    dev.off()
}

testIndexSpeed <- function(efdb, n.genes, n.mc, rejection.sampling) {
    genenames <- efdb$genenames
    ts <- rep(0, n.mc)
    for (i in 1:n.mc) {
        inds <- c()
        t <- 0
        while(length(inds)==0) {
            genes <- genenames[sample(length(genenames))[1:n.genes]]
            t <- Sys.time()
            inds <- findCellsWithGenes2(efdb, genes)
            t <- Sys.time() - t
            if (length(inds)==0 & rejection.sampling) { next }
        }
        ts[i] <- t
    }
    return( ts )
}

############################################

plotPancreasGeneSearch <- function(alpha=T) {
    efdb.muraro <- readRDS("efdb_muraro.rds")
    efdb.pancreas <- readRDS("efdb_pancreas.rds")
    efdb.hs <- readRDS("efdb_hs.rds")
    alpha.markers.muraro.authors <- c("GCG", "LOXL4", "PLCE1", "IRX2", "GC", "KLHL41", "CRYBA2", "TTR", "TM4SF4", "RGS4")

    findCellsWithGenes2(efdb.muraro
}

plotPancreasMarkerGenes <- function(max.genes=5, alpha=T) {
    efdb.muraro <- readRDS("efdb_muraro.rds")
    efdb.pancreas <- readRDS("efdb_pancreas.rds")
    efdb.hs <- readRDS("efdb_hs.rds")
    alpha.markers.muraro.authors <- c("GCG", "LOXL4", "PLCE1", "IRX2", "GC", "KLHL41", "CRYBA2", "TTR", "TM4SF4", "RGS4")
    beta.markers.segerstolpe.authors <- c("INS", "SCGN", "IAPP", "FXYD2", "RPL3", "G6PC2", "HSP90AB1", "PEBP1", "HSPA8", "EIR4A2", "ERO1LB", "SERINC1", "SLC30A8", "FAM159B", "SURF4", "NPTX2", "PFKFB2", "EDARADD", "MEG3", "HOPX", "LMO1", "PDX1", "PTEN", "SH3GL2", "ADCYAP1")
    cell.type <- "alpha"
    gene.list <- alpha.markers.muraro.authors
    if (!alpha) {
        cell.type <- "beta"
        gene.list <- beta.markers.segerstolpe.authors
    }
    #First find marker genes for beta cells in muraro only
    ret.muraro <- findMarkerGenes(efdb.muraro, max.genes, "*", cell.type)
    print(ret.muraro$genenames)
    #Now do beta cells in all pancreas datasets
    ret.pancreas <- findMarkerGenes(efdb.pancreas, max.genes, "*", cell.type)
    print(ret.pancreas$genenames)
    #Finally against the whole collection of Hs cells
    ret.hs <- findMarkerGenes(efdb.hs, max.genes, "*", cell.type)
    print(ret.hs$genenames)
    #######################################################
    #Plot precision-recall curves when using AND or OR for different number of 
    dat <- data.frame(precision=c(ret.muraro$precision[1,], ret.pancreas$precision[1,], ret.hs$precision[1,], ret.muraro$precision[2,], ret.pancreas$precision[2,], ret.hs$precision[2,]), recall=c(ret.muraro$recall[1,], ret.pancreas$recall[1,], ret.hs$recall[1,], ret.muraro$recall[2,], ret.pancreas$recall[2,], ret.hs$recall[2,]), logic=rep(c("AND", "OR"), each=3*max.genes), bg=rep(rep(c("Muraro", "pancreas", "Hs"), each=max.genes), 2), n.genes=rep(1:max.genes, 2*3))
    png(paste0(filename=paste0("Images/pancreas_",cell.type,"_markers.png")))
    gg <- ggplot(data=dat, aes(x=precision, y=recall, colour=bg, shape=logic, size=n.genes))+ geom_point() + xlim(0, 1) + ylim(0, 1) + xlab("Precision = TP/(TP+FP)") + ylab("Recall = TP/(TP+FN)") + theme_minimal() + theme(axis.title.y = element_text(size=20), axis.title.x = element_text(size=20), text = element_text(size=14), legend.text = element_text(size=14), legend.title=element_blank())
    print(gg)
    dev.off()
    #######################################################
    #Now we can evaluate the utility of other lists
    ret.muraro <- markerGeneScore(efdb.muraro, gene.list, "*", cell.type)
    print(ret.muraro$genenames)
    ret.pancreas <- markerGeneScore(efdb.pancreas, gene.list, "*", cell.type)
    print(ret.pancreas$genenames)
    ret.hs <- markerGeneScore(efdb.hs, gene.list, "*", cell.type)
    print(ret.hs$genenames)
    dat <- data.frame(precision=c(ret.muraro$precision[1,], ret.pancreas$precision[1,], ret.hs$precision[1,], ret.muraro$precision[2,], ret.pancreas$precision[2,], ret.hs$precision[2,]), recall=c(ret.muraro$recall[1,], ret.pancreas$recall[1,], ret.hs$recall[1,], ret.muraro$recall[2,], ret.pancreas$recall[2,], ret.hs$recall[2,]), logic=rep(c("AND", "OR"), each=3*max.genes), bg=rep(rep(c("Muraro", "pancreas", "Hs"), each=max.genes), 2), n.genes=rep(1:max.genes, 2*3))
    png(paste0(filename=paste0("Images/pancreas_",cell.type,"_markers_muraro.png")))
    gg <- ggplot(data=dat, aes(x=precision, y=recall, colour=bg, shape=logic, size=n.genes))+ geom_point() + xlim(0, 1) + ylim(0, 1) + xlab("Precision = TP/(TP+FP)") + ylab("Recall = TP/(TP+FN)") + theme_minimal() + theme(axis.title.y = element_text(size=20), axis.title.x = element_text(size=20), text = element_text(size=14), legend.text = element_text(size=14), legend.title=element_blank())
    print(gg)
    dev.off()
}

############################################

calculateQuantizationErrorsAndCompression <- function(save.file.name="quantization_compression_errors.rds") {
    #This takes several hours to run since a large number of cells must be indexed multiple times
    exclude.clusters <- c("not applicable", "alpha.contaminated", "beta.contaminated", "gamma.contaminated", "delta.contaminated", "miss")
    dir.name <- "/nfs/team218/sc-rna-seq-data/"
    #file.names <- c("segerstolpe", "muraro", "pollen")
    file.names <- c("baron-mouse", "fan", "deng-reads", "goolam") #, "biase", "klein", "chen", "campbell", "kolodziejczyk", "macosko", "shekhar", "manno_mouse", "marques", "romanov", "tasic-reads", "usoskin", "zeisel", "baron-human", "camp1", "darmanis", "manno_human", "segerstolpe", "muraro", "pollen", "xin", "yan")
    quantizations <- c(2, 4, 6)
    results <- list("file.names"=file.names, "gene.inds"=list(), "nrmse"=list(), "corrs"=list(), "compression"=list(), "p0"=c(), "quantizations"=quantizations, "n"=c())
    for (i in 1:length(file.names)) {
        #Load the file
        file.name <- paste0(dir.name, file.names[i], ".rds")
        d <- readRDS(file.name)
        e <- logcounts(d)
        #Figure out which columns to keep to avoid the excluded clusters
        cell.types.all <- colData(d)$cell_type1
        cell.types <- setdiff(unique(cell.types.all), exclude.clusters)
        inds.cell <- c()
        for (j in 1:length(cell.types)) {
            inds.cell <- c(inds.cell, which(cell.types[j]==cell.types.all))
        }
        e.excluded <- e[,sort(inds.cell)]
        genenames <- rowData(d)$feature_symbol
        nrmse <- matrix(rep(0, length(genenames)*length(quantizations)), ncol=length(quantizations))
        corrs <- matrix(rep(0, length(genenames)*length(quantizations)), ncol=length(quantizations))
        compressions <- matrix(rep(0, length(quantizations)))
        for (j in 1:length(quantizations)) {
        #Index
            efdb <- parseDataset(file.name, file.names[i], tfidf=quantizations[j], exclude.clusters=exclude.clusters)
            inds.keep <- c()
            for (k in 1:length(genenames)) {
                n <- length(which(e.excluded[k,]>0))
                if (n>1) {
                    ee <- extractExpressionLevels(efdb, genenames[k])
                    ek <- e[k, ee$inds]
                    cor.tmp <- cor(ee$exp.levels, ek)
                    if (!is.na(cor.tmp)) {
                        nrmse[k, j] <- sqrt(sum((ee$exp.levels - ek)^2)/n)/mean(ek)
                        corrs[k, j] <- cor.tmp
                        inds.keep <- c(inds.keep, k)
                    }
                }
            }
            results$gene.inds[[i]] <- inds.keep
            compressions[j] <- efdb$size.orig/object_size(efdb)
            results$p0[i] <- efdb$p0
            results$n[i] <- efdb$n.cells
        }
        results$nrmse[[i]] <- nrmse
        results$corrs[[i]] <- corrs
        results$compression[[i]] <- compressions
    }
    saveRDS(results, save.file.name)
}

plotQuantizationErrorsAndCompression <- function(load.file.name="quantization_compression_errors.rds") {
    require("ggplot2")
    results <- readRDS(load.file.name)
    #Put all of the results into a dataframe for plotting
    dat <- data.frame(compression=unlist(results$compression), quantization=as.character(rep(results$quantizations, length(results$file.names))), name=rep(results$file.names, each=length(results$quantizations)), n=rep(results$n, each=length(results$quantizations)))
    #Plot the compression
    png(filename="Images/scfind_compression_tfidf_scatter.png")
    gg <- ggplot(data=dat, aes(x=n, y=compression, colour=quantization)) + geom_point() + scale_y_log10() + scale_x_log10() + ylab("Compression") + xlab("#cells") + theme_minimal() + theme(axis.title.y = element_text(size=20), axis.title.x = element_text(size=20), text = element_text(size=14))
    print(gg)
    dev.off()
    #Plot the quantization errors
    nrmse <- c()
    for (i in 1:length(results$file.names)) {
        nrmse <- c(nrmse, colMeans(results$nrmse[[i]][results$gene.inds[[i]],]))
    }
    dat <- cbind(dat, nrmse)
    corrs <- c()
    for (i in 1:length(results$file.names)) {
        corrs <- c(corrs, colMeans(results$corrs[[i]][results$gene.inds[[i]],]))
    }
    dat <- cbind(dat, corrs)
    png(filename="Images/scfind_nrmse_scatter.png")
    gg <- ggplot(data=dat, aes(x=n, y=nrmse, colour=quantization)) + geom_point() + ylab("mean NRMSE") + xlab("#cells") + theme_minimal() + theme(axis.title.y = element_text(size=20), axis.title.x = element_text(size=20), text = element_text(size=14))
    print(gg)
    dev.off()
    png(filename="Images/scfind_corr_scatter.png")
    gg <- ggplot(data=dat, aes(x=n, y=corrs, colour=quantization)) + geom_point() + ylab("Pearson correlation") + xlab("#cells") + theme_minimal() + theme(axis.title.y = element_text(size=20), axis.title.x = element_text(size=20), text = element_text(size=14))
    print(gg)
    dev.off()
}

hypothalamusMarkerGenesSearch <- function() {
    #Use this method to
    #efdb.campbell <- parseDataset("/nfs/team218/sc-rna-seq-data/campbell.rds", "campbell", tfidf=4, exclude.clusters=c("miss"))
    #efdb.chen <- parseDataset("/nfs/team218/sc-rna-seq-data/chen.rds", "chen", tfidf=4)
    efdb.mm <- readRDS("efdb_mm.rds")
    
    bmi.gwas.n25.genes <- c("Bmp2", "Rarb", "Meox2", "Trp63", "Zfp366", "Tfap2b", "Gnat2", "Pde1c", "Nme5", "Dpf3")
    bmi.gwas.n32.genes <- c("Cck", "Lmo7", "Cps1", "Onecut1", "Plcd4", "4930452B06Rik", "4930453B24Rik", "Lrrc56", "Zfp366", "Stxbp6")
    obesity.n25.genes <- c("Sox8", "Fabp4", "Irs1", "Adra2a", "Npbwr1", "Ppargc1a", "Sgk1", "Apoc3", "Nr4a1", "Fosb", "Klf7", "Gfra2", "Rxrg", "Thrb")
    obesity.n32.genes <- c("Igf1", "Hrh1", "Igfbp6", "Rxrg", "Bbs1", "Cnr1", "Mc4r", "Bdnf", "Thrb", "Mmp2", "Ncoa3", "Lipa", "Ofd1", "Lpin1", "Flna")

    n.cells <- 100
    tfidf.bmi.n25 <- findBestCellTFIDF(efdb.campbell, bmi.gwas.n25.genes, n.cells=n.cells)
    tfidf.bmi.n32 <- findBestCellTFIDF(efdb.campbell, bmi.gwas.n32.genes, n.cells=n.cells)
    tfidf.obesity.n25 <- findBestCellTFIDF(efdb.campbell, bmi.obesity.n25.genes, n.cells=n.cells)
    tfidf.obesity.n32 <- findBestCellTFIDF(efdb.campbell, bmi.obesity.n32.genes, n.cells=n.cells)

    #Plot as heatmap to show enrichment of clusters
    
    #Find candidates from Chen et al as well and validate using scmap
    
}

demonstrateBatchCorrection <- function() {
    #Generate plots that shows how the fast batch-correction can improve the similarity of expression profiles
    genename.alpha <- "GCG"
    genename.beta <- "INS"
    efdb.hs <- readRDS("efdb_hs.rds")
}
