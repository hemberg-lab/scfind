findMarkerGenes <- function(efdb, n.markers, experiment.name, cluster.name, cross.val.fold=1) {
    #Use to identify the set of genes that will be best as marker genes. We use a greedy strategy where genes are first ranked based on 
    n.genes <- length(efdb$genenames)
    n.cells <- sum(efdb$cluster.sizes)
    inds <- list()
    #Find the indices of the cell-type that we are interested in
    tmp <- findClusterCellInds(efdb, experiment.name, cluster.name)
    cell.inds <- tmp$cell.inds
    cell.inds.test <- cell.inds
    cluster.inds <- tmp$cluster.inds
    #We identify the candidate genes based on the summary statistics based on the fraction of cells that they are expressed in
    print(paste("Calculating the overlap of each gene with the cell-type of interest"))
    support.gene <- rep(0, length(efdb$genenames))
    for (i in 1:length(efdb$genenames)) {
        #support.gene[i] <- -sum(efdb$m[i,]/efdb$cluster.sizes) + 2*sum(efdb$m[i,cluster.inds]/efdb$cluster.sizes[cluster.inds])
        support.gene[i] <- -sum(efdb$m[i,]) + 2*sum(efdb$m[i,cluster.inds])
    }
    #Create a mask for scoring these
    cell.mask <- rep(-1, n.cells)
    cell.mask[cell.inds] <- 1
    #Calculate the overlap of each gene with the cell.mask vector
    experiment.ends <- cumsum(efdb$n.cells)
    #If there are cells that we would like to ignore, then we should set the mask to 0
    tp.or <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    tn.or <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    fp.or <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    fn.or <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    tp.and <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    tn.and <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    fp.and <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    fn.and <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    if (cross.val.fold==1) {
        support.order <- sort(support.gene, decreasing=T, index.return=T)$ix[1:n.markers]
        #print(paste0("Candidate marker genes from greedy selection are: ", efdb$genenames[support.order]))
        ids.found.and <- 1:n.cells
        ids.found.or <- c()
        for (i in 1:n.markers) {
            inds[[i]] <- findCellsWithGenes2(efdb, efdb$genenames[support.order[i]])
            ids.found.and <- intersect(ids.found.and, inds[[i]])
            tp.and[i] <- length(which(cell.inds %in% ids.found.and))
            fn.and[i] <- length(cell.inds) - tp.and[i]
            fp.and[i] <- length(ids.found.and) - tp.and[i]
            tn.and[i] <- n.cells - tp.and[i] - fn.and[i] - fp.and[i]

            ids.found.or <- union(ids.found.or, inds[[i]])
            tp.or[i] <- length(which(cell.inds %in% ids.found.or))
            fn.or[i] <- length(cell.inds) - tp.or[i]
            fp.or[i] <- length(ids.found.or) - tp.or[i]
            tn.or[i] <- n.cells - tp.or[i] - fn.or[i] - fp.or[i]
        }
    }
    else {
        support.order <- sort(support.gene, decreasing=T, index.return=T)$ix[1:(2*n.markers)]
        print(efdb$genenames[support.order])
        gene.mask <- matrix( rep( -1, 2*n.markers*n.cells ), ncol=2*n.markers )
        for (i in 1:(2*n.markers)) {
            inds <- findCellsWithGenes2(efdb, efdb$genenames[support.order[i]])
            gene.mask[inds,i] <- 1
        }
        for (j in 1:cross.val.fold) {
            tmp <- sample(n.cells)
            cells.id.test <- tmp[1:round(n.cells/cross.val.fold)]
            cells.id.train <- tmp[(1+round(n.cells/cross.val.fold)):n.cells]
            cell.inds.test <- intersect(cells.id.test, cell.inds) #cells of interest in the test set
            support.gene.train <- matrix( rep(0, 2*n.markers), nrow=1)
            for (i in 1:2*n.markers) { support.gene.train[i] <- sum(cell.mask[cells.id.train]==gene.mask[cells.id.train,i]) }
            support.order.train <- support.order[sort(support.gene.train, decreasing=T, index.return=T)$ix[1:n.markers]]
            #print(paste0("Candidate marker genes from greedy selection are: ", efdb$genenames[support.order.train]))
            ids.found.and <- 1:n.cells
            ids.found.or <- c()
            for (i in 1:n.markers) {
                cell.inds <- intersect(cell.inds.test, which(gene.mask[support.order.train[i],support.order.train[i]]>0))
                ids.found.and <- intersect(ids.found.and, cell.inds)
                tp.and[j, i] <- length(intersect(cell.inds.test, ids.found.and))
                fn.and[j, i] <- length(cell.inds.test) - tp.and[j, i]
                fp.and[j, i] <- length(intersect(cells.id.test, ids.found.and)) - tp.and[j, i]
                tn.and[j, i] <- length(cells.id.test) - tp.and[j, i] - fn.and[j, i] - fp.and[j, i]

                ids.found.or <- union(ids.found.or, cell.inds)
                tp.or[j, i] <- length(intersect(cell.inds.test, ids.found.or))
                fn.or[j, i] <- length(cell.inds.test) - tp.or[j, i]
                fp.or[j, i] <- length(intersect(cells.id.test, ids.found.or)) - tp.or[j, i]
                tn.or[j, i] <- length(cells.id.test) - tp.or[j, i] - fn.or[j, i] - fp.or[j, i]
            }
        }
    }
    recall.or <- tp.or/length(cell.inds.test)
    precision.or <- tp.or/(tp.or+fp.or)
    fpr.or <- fp.or/(n.cells - length(cell.inds))
    f1.or <- 2/(recall.or^-1 + precision.or^-1)
    recall.and <- tp.and/length(cell.inds.test)
    precision.and <- tp.and/(tp.and+fp.and)
    fpr.and <- fp.and/(n.cells - length(cell.inds))
    f1.and <- 2/(recall.and^-1 + precision.and^-1)
    return( list( "inds"=inds, "cell.inds"=cell.inds, "f1" = rbind(f1.and, f1.or), "tp" = rbind(tp.and, tp.or), "fp" = rbind(fp.and, fp.or), "recall"=rbind(recall.and, recall.or), "precision"=rbind(precision.and, precision.or), "genenames" = efdb$genenames[support.order] ) ) 
}

extractExpressionLevels <- function(efdb, genename) {
    gene.ind <- which(efdb$genenames==genename)
    experiment.ends <- cumsum(efdb$n.cells)
    exp.levels <- c()
    inds <- c()
    if (length(gene.ind)==1) {
        tfs <- as.numeric(efdb$tf[[genename]])
        cell.inds.bycluster <- eliasFanoDecodingGeneByClusterCpp(as.numeric(efdb$H[[genename]]), as.numeric(efdb$L[[genename]]), efdb$m[gene.ind,], efdb$H.start[gene.ind,], efdb$cluster.sizes, efdb$permutation)
        #Since the quantized values were stored in order of the clusters, that is how they need to be retrieved
        start <- 0
        experiment.ind <- 1
        clusters.cumsum <- cumsum(efdb$n.clusters)
        q <- efdb$tf.quantization[experiment.ind]
        for (i in 1:length(efdb$cluster.sizes)) {
            if (length(cell.inds.bycluster[[i]])>0) {
                inds <- c(inds, cell.inds.bycluster[[i]])
                n.quant.bits <- q*length(cell.inds.bycluster[[i]])
                tmp <- dequantizeCpp(tfs[start+(1:n.quant.bits)], q)/2^q + .5/2^q
                exp.levels <- c(exp.levels, 10^-qnorm(tmp, efdb$tf.mean[gene.ind, experiment.ind], efdb$tf.sd[gene.ind,experiment.ind])/efdb$norm.factor[cell.inds.bycluster[[i]]])
                start <- start + n.quant.bits
            }
            if (i>clusters.cumsum[experiment.ind]) {
                experiment.ind <- experiment.ind + 1
                q <- efdb$tf.quantization[experiment.ind]
            }
        }
    }
    return( list("exp.levels"=exp.levels, "inds"=inds) )
}

extractExpressionLevelsByCluster <- function(efdb, genename) {
    gene.ind <- which(efdb$genenames==genename)
    experiment.ends <- cumsum(efdb$n.cells)
    exp.levels <- list()
    inds <- list()
    cluster.inds <- c()
    if (length(gene.ind)==1) {
        tfs <- as.numeric(efdb$tf[[genename]])
        cell.inds.bycluster <- eliasFanoDecodingGeneByClusterCpp(as.numeric(efdb$H[[genename]]), as.numeric(efdb$L[[genename]]), efdb$m[gene.ind,], efdb$H.start[gene.ind,], efdb$cluster.sizes, efdb$permutation)
        #Since the quantized values were stored in order of the clusters, that is how they need to be retrieved
        start <- 0
        experiment.ind <- 1
        clusters.cumsum <- cumsum(efdb$n.clusters)
        q <- efdb$tf.quantization[experiment.ind]
        for (i in 1:length(efdb$cluster.sizes)) {
            if (length(cell.inds.bycluster[[i]])>0) {
                inds[[length(inds) + 1]] <- cell.inds.bycluster[[i]]
                n.quant.bits <- q*length(cell.inds.bycluster[[i]])
                tmp <- dequantizeCpp(tfs[start+(1:n.quant.bits)], q)/2^q + .5/2^q
                exp.levels[[length(exp.levels) + 1]] <- 10^-qnorm(tmp, efdb$tf.mean[gene.ind, experiment.ind], efdb$tf.sd[gene.ind,experiment.ind])/efdb$norm.factor[cell.inds.bycluster[[i]]]
                start <- start + n.quant.bits
                cluster.inds[length(cluster.inds) + 1] <- i
            }
            if (i>clusters.cumsum[experiment.ind]) {
                experiment.ind <- experiment.ind + 1
                q <- efdb$tf.quantization[experiment.ind]
            }
        }
    }
    return( list("exp.levels"=exp.levels, "inds"=inds, "cluster.inds"=cluster.inds) )
}


fastBatchCorrection <- function(efdb, genenames, experiment.names, cluster.names) {
    #Use the simplified version of Combat proposed by Tung et al (2017) to carry out simplified batch correct the quantized expression levels
    clust.starts <- c(1, cumsum(efdb$cluster.sizes)+1)
    #Find all of the cells across experiments that correspond to the clusters of interest
    cell.inds <- list()
    cluster.inds <- list()
    lm.res <- list()
    for (i in 1:length(cluster.names)) {
        tmp <- findClusterCellInds(efdb, experiment.names, cluster.names[i])
        cell.inds[[i]] <- tmp$cell.inds
        cluster.inds[[i]] <- tmp$cluster.inds
    }
    n.cells <- sum(efdb$cluster.sizes[unlist(cluster.inds)])
    select.cluster.starts <- c(1, 1+cumsum(clust.starts[unlist(cluster.inds)+1] - clust.starts[unlist(cluster.inds)]))[1:length(unlist(cluster.inds))]
    #Extract the approximate expression levels of the genes of interest
    all.exp.levels <- matrix(rep(0, n.cells*length(genenames)), ncol=length(genenames))
    all.exp.levels.estimate <- matrix(rep(0, n.cells*length(genenames)), ncol=length(genenames))
    for (j in 1:length(genenames)) {
    #The cells that have expression of a given gene. First index is for cluster and second for gene
        genename.cell.inds <- list()
        exp.levels <- list()
        inds <- list()
        ee <- list()
        for (i in 1:length(cluster.names)) {
            genename.cell.inds[[i]] <- list()
            exp.levels[[i]] <- list()
            ee[[i]] <- extractExpressionLevelsByCluster(efdb, genenames[j])
            inds[[i]] <- which(ee[[i]]$cluster.inds %in% cluster.inds[[i]])
            for (k in 1:length(inds[[i]])) {
                genename.cell.inds[[i]][[k]] <- ee[[i]]$inds[[inds[[i]][k]]]
            }
        }
    #Carry out the approximate batch effects correction using a linear model with coefficients for cluster membership
        #
        exp.levels <- c()
        cell.type <- c()
        exp.inds <- list()
        tmp.inds <- rep(0, n.cells)
        for (i in 1:length(cluster.names)) {
            for (k in 1:length(inds[[i]])) {
                tmp <- matrix(rep(0, efdb$cluster.sizes[cluster.inds[[i]][k]]))
                #tmp[which(ee[[i]]$inds[[k]] %in% efdb$permutation[clust.starts[inds[[i]][k]]:(clust.starts[inds[[i]][k]+1]-1)])] <- ee[[i]]$exp.levels[[inds[[i]][k]]]
                which(ee[[i]]$inds[[k]] %in% efdb$permutation[clust.starts[inds[[i]][k]]:(clust.starts[inds[[i]][k]+1]-1)])
                #print(length(tmp))
                exp.levels <- c(exp.levels, ee[[i]]$exp.levels[[inds[[i]][k]]])
                cell.type <- c(cell.type, rep(efdb$cluster.names[cluster.inds[[i]][k]], length(ee[[i]]$exp.levels[[inds[[i]][k]]])))
            }
        }
        lm.res[[j]] <- lm(exp.levels ~ cell.type)
        all.exp.levels.estimate[tmp.inds,j] <- exp.levels - lm.res[[j]]$fitted.values
        all.exp.levels[tmp.inds,j] <- exp.levels
    }
}

findClusterCellInds <- function(efdb, experiment.name, cluster.name) {
    cell.inds <- c()
    cluster.starts <- c(1, cumsum(efdb$cluster.sizes))
    experiment.starts <- cumsum(efdb$n.cells)
    experiment.cluster.starts <- cumsum(efdb$n.clusters)
    cluster.ind <- c()
    experiment.ind <- c()
    if (experiment.name=="*") {
        for (i in 1:length(efdb$experiment.name)) {
            tmp <- which(efdb$cluster.names==paste0(efdb$experiment.name[i], "_", cluster.name))
            if (length(tmp)>0) {
                cluster.ind <- c(cluster.ind, tmp)
                experiment.ind <- c(experiment.ind, min(which(tmp<experiment.cluster.starts)))
            }
        }
    }
    else {
        cluster.ind <- which(efdb$cluster.names==paste0(experiment.name, "_", cluster.name))
        experiment.ind <- min(which(cluster.ind<experiment.cluster.starts))
    }
    for (i in 1:length(cluster.ind)) {
        cell.inds <- c(cell.inds, cluster.starts[cluster.ind[i]]:cluster.starts[cluster.ind[i]+1])
    }
    #This is an ugly hack that needs to be fixed when the index is created!
    if (length(experiment.starts)==1) { cell.inds <- efdb$permutation[cell.inds] }
    return( list("cell.inds"=cell.inds, "cluster.inds"=cluster.ind) )
}


generateMarkerGenesPlotsPancreas <- function(method="f1") {
    #Use this method to generate plots for marker gene identification in the human pancreas
    require("ggplot2")
    alpha.markers.muraro.authors <- c("GCG", "LOXL4", "PLCE1", "IRX2", "GC", "KLHL41", "CRYBA2", "TTR", "TM4SF4", "RGS4")
    #Find the marker genes from just looking at the Muraro dataset
    ret.muraro <- findMarkerGenes(efdb.muraro, 10, "*", "alpha", 5)
    #Find the marker genes from looking at all pancreas datasets
    ret.pancreas <- findMarkerGenes(efdb.pancreas, 10, "*", "alpha", 5)
    #Find marker genes from looking at all Hs datasets
    ret.hs <- findMarkerGenes(efdb.hs, 10, "*", "alpha", 5)

    #Compare how good the different sets of genes are across other backgrounds
    #alpha.markers.score <- evaluateMarkerGenes
    dat <- data.frame(precision=c(colSums(ret.muraro$precision)/cross.val, colSums(ret.pancreas$precision)/cross.val, colSums(ret.hs$precision)/cross.val), recall=c(colSums(ret.muraro$recall)/cross.val, colSums(ret.pancreas$recall)/cross.val, colSums(ret.hs$recall)/cross.val), bg=rep(c("Muraro", "pancreas", "Hs"), each=10))
    png(paste0(filename="Images/pancreas_alpha_markers_",method,".png"))
    gg <- ggplot(data=dat, aes(x=precision, y=recall, colour=bg))+ geom_point() + xlim(0, 1) + ylim(0, 1) + ylab("Precision = TP/(TP+FP)") + xlab("Recall = TP/(TP+FN)") + theme_minimal() + theme(axis.title.y = element_text(size=20), axis.title.x = element_text(size=20), text = element_text(size=14), legend.text = element_text(size=14), legend.title=element_blank())
    print(gg)
    dev.off()
}


markerGeneScore <- function(efdb, marker.genenames, experiment.name, cluster.name, n.markers=5, method="f1", cross.val.fold=1) {
    #Find the ids of the cells of interests
    n.genes <- length(efdb$genenames)
    n.cells <- sum(efdb$cluster.sizes)
    tmp <- findClusterCellInds(efdb, experiment.name, cluster.name)
    cell.inds <- tmp$cell.inds
    cluster.inds <- tmp$cluster.inds
    support.gene <- rep(0, length(marker.genenames))
    for (i in 1:length(marker.genenames)) {
        gene.ind[i] <- which(efdb$genenames==marker.genenames[i])
        support.gene[i] <- -sum(efdb$m[gene.ind[i],]) + 2*sum(efdb$m[gene.ind[i],cluster.inds])
    }
    support.order <- gene.ind[sort(support.gene, decreasing=T, index.return=T)$ix[1:n.markers]]
    ids.found.and <- 1:n.cells
    ids.found.or <- c()
    inds <- list()
    tp.or <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    tn.or <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    fp.or <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    fn.or <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    tp.and <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    tn.and <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    fp.and <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    fn.and <- matrix( rep(0, n.markers*cross.val.fold), nrow=cross.val.fold)
    print(efdb$genenames[support.order])
    for (i in 1:n.markers) {
        inds[[i]] <- findCellsWithGenes2(efdb, efdb$genenames[support.order[i]])
        ids.found.and <- intersect(ids.found.and, inds[[i]])
        tp.and[i] <- length(which(cell.inds %in% ids.found.and))
        fn.and[i] <- length(cell.inds) - tp.and[i]
        fp.and[i] <- length(ids.found.and) - tp.and[i]
        tn.and[i] <- n.cells - tp.and[i] - fn.and[i] - fp.and[i]
        
        ids.found.or <- union(ids.found.or, inds[[i]])
        tp.or[i] <- length(which(cell.inds %in% ids.found.or))
        fn.or[i] <- length(cell.inds) - tp.or[i]
        fp.or[i] <- length(ids.found.or) - tp.or[i]
        tn.or[i] <- n.cells - tp.or[i] - fn.or[i] - fp.or[i]
    }
    recall.or <- tp.or/length(cell.inds)
    precision.or <- tp.or/(tp.or+fp.or)
    fpr.or <- fp.or/(n.cells - length(cell.inds))
    f1.or <- 2/(recall.or^-1 + precision.or^-1)
    recall.and <- tp.and/length(cell.inds)
    precision.and <- tp.and/(tp.and+fp.and)
    fpr.and <- fp.and/(n.cells - length(cell.inds))
    f1.and <- 2/(recall.and^-1 + precision.and^-1)
    return( list( "inds"=inds, "cell.inds"=cell.inds, "f1" = rbind(f1.and, f1.or), "tp" = rbind(tp.and, tp.or), "fp" = rbind(fp.and, fp.or), "recall"=rbind(recall.and, recall.or), "precision"=rbind(precision.and, precision.or), "genenames" = efdb$genenames[support.order[1:n.markers]] ) ) 
}


calculateInformationGain <- function(n.gene.cell, n.other.cell, n.cell, n.other) {
    information.total <- -n.cell*log2(n.cell/(n.cell+n.other))/(n.cell+n.other) - n.other*log2(n.other/(n.cell+n.other))/(n.cell+n.other)
    information.cell <- -n.gene.cell*log2(n.gene.cell/(n.gene.cell+n.other.cell))/(n.gene.cell+n.other.cell) - n.other.cell*log2(n.other.cell/(n.gene.cell+n.other.cell))/(n.gene.cell+n.other.cell)
    information.not.cell <- -(n.cell-n.gene.cell)*log2((n.cell-n.gene.cell)/(n.cell-n.gene.cell+n.other-n.other.cell))/(n.cell-n.gene.cell+n.other-n.other.cell) - (n.other-n.other.cell)*log2((n.other-n.other.cell)/(n.cell-n.gene.cell+n.other-n.other.cell))/(n.cell-n.gene.cell+n.other-n.other.cell)
    E <- (n.gene.cell+n.other.cell)*information.cell/(n.cell+n.other) + (n.cell-n.gene.cell+n.other-n.other.cell)*information.not.cell/(n.cell+n.other)
    information.gain <- information.total - E
    Iv <- -(n.gene.cell+n.other.cell)*log2((n.gene.cell+n.other.cell)/(n.cell+n.other))-(n.cell-n.gene.cell+n.other-n.other.cell)*log2((n.cell-n.gene.cell+n.other-n.other.cell)/(n.cell+n.other))
    return( list( "information.gain"=information.gain, "Iv"=Iv ) )
}

findCellsWithGenes2 <- function(efdb, genelist) {
  #Use this method to run the AND query on the big database
  #find the number of cells for each gene, start by searching the smallest ones first
  require('bit')
  require('hash')
  require('Rcpp')
  sourceCpp('eliasFanoCoding.cpp')
  n <- rep(0, length(genelist))
  for (i in 1:length(genelist)) { n[i] <- numberOfCells(efdb, genelist[i]) }
  order <- sort.int(n, index.return=T)$ix
  cell.inds.byexperiment <- vector("list", length(genelist))
  inds <- vector("list", length(genelist))
  experiment.ends <- cumsum(efdb$n.cells)
  ids.unique <- c()
  n.genes.not.found <- 0
  clusters.found <- 1:length(efdb$cluster.sizes)
  gene.ind <- rep(1, length(genenames))
  if (length(order)>=1) {
      for (i in 1:length(order)) {
          tmp <- which(genelist[order[i]]==efdb$genenames)
          if (length(tmp)>0) {
              gene.ind[i] <- tmp[1]  
              cell.inds.byexperiment[[i]] <- eliasFanoDecodingGeneByExperimentCpp(as.numeric(efdb$H[[genelist[order[i]]]]), as.numeric(efdb$L[[genelist[order[i]]]]), efdb$m[gene.ind[i],], efdb$H.start[gene.ind[i],], efdb$cluster.sizes, efdb$permutation, cumsum(efdb$n.cells))
              if (length(experiment.ends)>1) {
                  inds[[i]] <- which(efdb$permutation[1:experiment.ends[1]] %in% cell.inds.byexperiment[[i]][[1]])
                  for (j in 2:length(experiment.ends)) {
                      inds[[i]] <- c(inds[[i]], experiment.ends[j-1] + which(efdb$permutation[(experiment.ends[j-1]+1):experiment.ends[j]] %in% cell.inds.byexperiment[[i]][[j]]))
                  }
              }
              else {
                  inds[[i]] <- cell.inds.byexperiment[[i]][[1]]
              }
          }
          if (i==1) {
              ids.unique <- inds[[i]]
          }
          else {
              ids.unique <- intersect(inds[[i]], ids.unique)
          }
      }
  }
  return( ids.unique )
}

calculatePvaluesBinarySearch <- function(efdb, ids.unique, n.genes.eff) {
    #Use this method to calculate the p-values using different methods for the basic search function
    clusters <- efdb$assigned.clusters[ids.unique]
    #Calculate the significance of the enrichment for each cell-type. For a single gene, this corresponds to a hypergeometric test with x=number of cells with gene, n=number of cells in cluster, m=number of cells expressing gene and N=total number of cells. This can be generalized for larger number of genes. Note that this tests for the significance of the found cells being in the same cluster, NOT that we would find that many cells expressing the set of genes
    tmp <- rle(sort(clusters))
    clusters.found <- tmp$values
    cluster.counts <- tmp$lengths
    p.vals.hyper <- rep(1, length(clusters.found))
    for (i in 1:length(clusters.found)) {
      p.vals.hyper[i] <- 1 - phyper(cluster.counts[i], efdb$cluster.sizes[clusters.found[i]], sum(efdb$cluster.sizes) - efdb$cluster.sizes[clusters.found[i]], length(ids.unique))
    }
    #As an option, use instead the likelihood ratio test that was used for gcompass. The advantage of using this method is that it accounts for the overall expression in each cluster. The problem with this method is for small gene-lists where it is very hard to achieve significance due to how the test is constructed
    p1 <- cluster.counts/efdb$cluster.sizes[clusters.found]
    p0 <- efdb$p0[clusters.found]
    p.vals.ratio <- 1 - pchisq(-2*log(p1) + (n.genes.eff*log(p0)), 1, lower.tail=FALSE)
    #For short lists, it is useful to consider the binomial test
    p.vals.binom <- rep(0, length(cluster.counts))
    for (i in 1:length(clusters.found)) { p.vals.binom[i] <- binom.test(cluster.counts[i], efdb$cluster.sizes[clusters.found[i]], p0[i]^n.genes.eff, "greater")$p.value }
    #Calculate the summary statistics 
    return( list( "ids.unique"=ids.unique, "clusters"=clusters, "clusters.found"=clusters.found, "cluster.counts"=cluster.counts, "cluster.names"=efdb$cluster.names[clusters.found], "p.vals.ratio"=p.vals.ratio, "p.vals.binom"=p.vals.binom, "p.vals.hyper"=p.vals.hyper, "p1"=p1, "p0"=p0 ) )
}


findAssociationRulesCellsFPgrowth <- function(efdb, min.support=.5, min.cells=5, min.genes=10, transpose=T) {
    require("Rcpp")
    sourceCpp("eliasFanoCoding.cpp")
    #Calculate the support, ie the fraction of cells where each gene appears
    n.cells <- sum(efdb$cluster.sizes)
    genenames <- efdb$genenames
    support <- rowSums(efdb$m)/n.cells
    ret.sort <- sort.int(support, decreasing=T, index.return=T)
    support.order <- ret.sort$ix[which(ret.sort$x>min.support)]
    print(paste0("Found ",length(support.order)," candidate genes."))
    Hss <- vector("list", length(support.order))
    Lss <- vector("list", length(support.order))
    for (i in 1:length(support.order)) {
        Hss[[i]] <- as.numeric(efdb$H[[genenames[support.order[i]]]])
        Lss[[i]] <- as.numeric(efdb$L[[genenames[support.order[i]]]])
    }
    n.cells.cumsum <- cumsum(efdb$n.cells)
    tmp <- findAssociationRulesCellsFPgrowthCpp(min.genes, min.cells, Hss, Lss, efdb$m[support.order,], efdb$H.start[support.order,], efdb$cluster.sizes, efdb$permutation, n.cells.cumsum, transpose)
    ret <- vector("list", length(tmp$potential.itemsets))
    for (i in 1:length(ret)) {
        ret[[i]]$genenames <- genenames[support.order[tmp$potential.itemsets[[i]]$genes]]
        ret[[i]]$cells <- tmp$potential.itemsets[[i]]$cells
        #Find the breakdown by experiment
        ret[[i]]$experiment <- rep(0, length(efdb$n.cells))
        for (j in 2:length(n.cells.cumsum)) {
            ret[[i]]$experiment[j-1] <- length(which(tmp$potential.itemsets[[i]]$cells>n.cells.cumsum[j-1] & tmp$potential.itemsets[[i]]$cells<=n.cells.cumsum[j]))
        }
    }
    return( ret )
}

bitToInt <- function(x) {
    ret <- 0
    len <- length(x)
    for (i in 1:len) {
        if (x[i]) { ret <- ret + 2^(len-i) }
    }
    return( ret )
}

findBestCellTFIDF <- function(efdb, genenames, n.cells=10, score.method="sum") {
    #Use this method to search for cells using the TF-IDF strategy. Will return the n.cells best matches
    require("bit")
    N <- sum(efdb$cluster.sizes)
    experiment.ends <- cumsum(efdb$n.cells)
    cell.inds.byexperiment <- vector("list", length(genenames))
    tf <- vector("list", length(genenames))
    idf <- rep(0, length(genenames))
    gene.ind <- rep(1, length(genenames))
    #Find all of the genes that have a non-zero value, this is one list per experiment, containing inds of all cells that express at least one of the genes
    allcells.bygenename.indexes <- vector("list", length(experiment.ends))
    for (i in 1:length(genenames)) { #Obtain the scores for each of the genes
        tmp <- which(genenames[i]==efdb$genenames)
        tf[[i]] <- vector("list", length(experiment.ends))
        if (length(tmp)>0) {
            gene.ind[i] <- tmp[1]
            #This is a list of lists where the 1st dimension is genename and the second is experiments
            cell.inds.byexperiment[[i]] <- eliasFanoDecodingGeneByExperimentCpp(as.numeric(efdb$H[[genenames[i]]]), as.numeric(efdb$L[[genenames[i]]]), efdb$m[gene.ind[i],], efdb$H.start[gene.ind[i],], efdb$cluster.sizes, efdb$permutation, cumsum(efdb$n.cells))
            tfs <- as.numeric(efdb$tf[[genenames[i]]])
            #decode the quantized expression levels
            inds <- which(cell.inds.byexperiment[[i]][[1]] %in% efdb$permutation[1:experiment.ends[1]])
            tmp <- dequantizeCpp(tfs[1:(efdb$tf.quantization[1]*length(inds))], efdb$tf.quantization[1])/2^efdb$tf.quantization[1] + .5/2^efdb$tf.quantization[1]
            tf[[i]][[1]] <- 10^-qnorm(tmp, efdb$tf.mean[gene.ind[i],1], efdb$tf.sd[gene.ind[i],1])
            allcells.bygenename.indexes[[1]] <- union(allcells.bygenename.indexes[[1]], cell.inds.byexperiment[[i]][[1]])
            first.ind <- 1 + efdb$tf.quantization[1]*length(inds)
            if (length(experiment.ends)>1) {
                for (j in 2:length(experiment.ends)) {
                    inds <- experiment.ends[j-1] + which(cell.inds.byexperiment[[i]][[j]] %in% efdb$permutation[(experiment.ends[j-1]+1):experiment.ends[j]])
                    tmp <- dequantizeCpp(tfs[first.ind:(first.ind + efdb$tf.quantization[j]*length(inds))], efdb$tf.quantization[j])/2^efdb$tf.quantization[j] + .5/2^efdb$tf.quantization[j]
                    tf[[i]][[j]] <- 10^-qnorm(tmp, efdb$tf.mean[gene.ind[i],j], efdb$tf.sd[gene.ind[i],j])
                    allcells.bygenename.indexes[[j]] <- union(allcells.bygenename.indexes[[j]], cell.inds.byexperiment[[i]][[j]])
                    first.ind <- efdb$tf.quantization[j]*length(inds) + first.ind
                }
            }
            #calculate the inverse frequency
            idf[i] <- 1 + log(N/(1+length(unlist(cell.inds.byexperiment[[i]]))))
        }
    }
    #Now we can calculate the cosine similarity with the query
    cell.score <- c() #get a list with indexes for each experiment
    if (score.method=="cosine") {
        cell.score <- calculateCosineSimilarity(genenames, allcells.bygenename.indexes, tf, idf, cell.inds.byexperiment)
    }
    else if (score.method=="sum") {
        cell.score <- calculateTFIDFSum(genenames, allcells.bygenename.indexes, tf, idf, cell.inds.byexperiment)
    }
    #Find the top n.cells and figure out their cell-types
    cluster.starts <- cumsum(efdb$cluster.sizes)
    merged.cell.score <- unlist(cell.score)
    score.cutoff <- sort(merged.cell.score, decreasing=T)[n.cells]
    inds <- which(efdb$permutation[1:experiment.ends[1]] %in% allcells.bygenename.indexes[[1]][which(cell.score[[1]]>=score.cutoff)])
    cluster <- c()
    cluster.counts <- rep(0, length(efdb$cluster.names))
    cluster.ind <- c()
    if (length(inds)>0) {
        for (i in 1:length(inds)) {
            cluster.ind <- c(cluster.ind, min(which(inds[i]<cluster.starts)))
            cluster <- c(cluster, efdb$cluster.names[tail(cluster.ind,1)])
            cluster.counts[tail(cluster.ind,1)] <- cluster.counts[tail(cluster.ind,1)] + 1
        }
    }
    if (length(experiment.ends)>1) {
        for (j in 2:length(experiment.ends)) {
            inds <- which(efdb$permutation[(experiment.ends[j-1]+1):experiment.ends[j]] %in% allcells.bygenename.indexes[[j]][which(cell.score[[j]]>=score.cutoff)])
            if (length(inds)>0) {
                for (i in 1:length(inds)) {
                    cluster.ind <- c(cluster.ind, min(which(experiment.ends[j-1]+inds[i]<cluster.starts)))
                    cluster <- c(cluster, efdb$cluster.names[tail(cluster.ind,1)])
                    cluster.counts[tail(cluster.ind,1)] <- cluster.counts[tail(cluster.ind,1)] + 1
                }
            }
        }
    }
    #calculate the p-value for the different clusters using a hypergeometric test
    tmp <- rle(cluster)
    clusters.found <- tmp$values
    cluster.counts <- tmp$lengths
    p.vals.hyper <- rep(1, length(clusters.found))
    cluster.sizes <- rep(1, length(clusters.found))
    for (i in 1:length(p.vals.hyper)) {
        cluster.sizes[i] <- efdb$cluster.sizes[which(efdb$cluster.names==clusters.found[i])[1]]
      p.vals.hyper[i] <- 1 - phyper(cluster.counts[i], cluster.sizes[i], N - cluster.sizes[i], length(cluster))
    }
    #store summary in a data frame
    cluster.data <- data.frame(cluster.names=clusters.found, count=cluster.counts, size=cluster.sizes, p.vals=p.vals.hyper)
    return( list( "cluster"=cluster, "inds"=inds, "cell.score"=cell.score, "tf"=tf, "idf"=idf, "allcells.indexes"=allcells.bygenename.indexes, "cell.inds"=cell.inds.byexperiment, "cluster.data"=cluster.data ) )
}

findBestCellTFIDF2 <- function(efdb, genenames, n.cells=10, score.method="sum") {
    #Use this method to search for cells using the TF-IDF strategy. Will return the n.cells best matches. Same as before, but more compact since it uses the extractExpressionLevels function
    require("bit")
    N <- sum(efdb$cluster.sizes)
    experiment.ends <- cumsum(efdb$n.cells)
    cell.inds.byexperiment <- vector("list", length(genenames))
    tf <- vector("list", length(genenames))
    idf <- rep(0, length(genenames))
    gene.ind <- rep(1, length(genenames))
    #Find all of the genes that have a non-zero value, this is one list per experiment, containing inds of all cells that express at least one of the genes
    allcells.bygenename.indexes <- vector("list", length(experiment.ends))
    for (i in 1:length(genenames)) { #Obtain the scores for each of the genes
        tmp <- which(genenames[i]==efdb$genenames)
        tf[[i]] <- vector("list", length(experiment.ends))
        if (length(tmp)>0) {
            ee <- extractExpressionLevels(efdb, genenames[i], no.norm=T)
            idf[i] <- 1 + log(N/(1+length(ee$inds)))
            allcells.bygenename.indexes[[i]] <- ee$inds
            tf[[i]] <- ee$exp.levels
            cell.inds.byexperiment[[i]] <- eliasFanoDecodingGeneByExperimentCpp(as.numeric(efdb$H[[genenames[i]]]), as.numeric(efdb$L[[genenames[i]]]), efdb$m[gene.ind[i],], efdb$H.start[gene.ind[i],], efdb$cluster.sizes, efdb$permutation, cumsum(efdb$n.cells))
        }
    }
    #Now we can calculate the cosine similarity with the query
    cell.score <- c() #get a list with indexes for each experiment
    if (score.method=="cosine") {
        cell.score <- calculateCosineSimilarity(genenames, allcells.bygenename.indexes, tf, idf, cell.inds.byexperiment)
    }
    else if (score.method=="sum") {
        cell.score <- calculateTFIDFSum(genenames, allcells.bygenename.indexes, tf, idf, cell.inds.byexperiment)
    }
    #Find the top n.cells and figure out their cell-types
    cluster.starts <- cumsum(efdb$cluster.sizes)
    merged.cell.score <- unlist(cell.score)
    score.cutoff <- sort(merged.cell.score, decreasing=T)[n.cells]
    inds <- which(efdb$permutation[1:experiment.ends[1]] %in% allcells.bygenename.indexes[[1]][which(cell.score[[1]]>=score.cutoff)])
    cluster <- c()
    cluster.counts <- rep(0, length(efdb$cluster.names))
    cluster.ind <- c()
    if (length(inds)>0) {
        for (i in 1:length(inds)) {
            cluster.ind <- c(cluster.ind, min(which(inds[i]<cluster.starts)))
            cluster <- c(cluster, efdb$cluster.names[tail(cluster.ind,1)])
            cluster.counts[tail(cluster.ind,1)] <- cluster.counts[tail(cluster.ind,1)] + 1
        }
    }
    if (length(experiment.ends)>1) {
        for (j in 2:length(experiment.ends)) {
            inds <- which(efdb$permutation[(experiment.ends[j-1]+1):experiment.ends[j]] %in% allcells.bygenename.indexes[[j]][which(cell.score[[j]]>=score.cutoff)])
            if (length(inds)>0) {
                for (i in 1:length(inds)) {
                    cluster.ind <- c(cluster.ind, min(which(experiment.ends[j-1]+inds[i]<cluster.starts)))
                    cluster <- c(cluster, efdb$cluster.names[tail(cluster.ind,1)])
                    cluster.counts[tail(cluster.ind,1)] <- cluster.counts[tail(cluster.ind,1)] + 1
                }
            }
        }
    }
    #calculate the p-value for the different clusters using a hypergeometric test
    tmp <- rle(cluster)
    clusters.found <- tmp$values
    cluster.counts <- tmp$lengths
    p.vals.hyper <- rep(1, length(clusters.found))
    cluster.sizes <- rep(1, length(clusters.found))
    for (i in 1:length(p.vals.hyper)) {
        cluster.sizes[i] <- efdb$cluster.sizes[which(efdb$cluster.names==clusters.found[i])[1]]
      p.vals.hyper[i] <- 1 - phyper(cluster.counts[i], cluster.sizes[i], N - cluster.sizes[i], length(cluster))
    }
    #store summary in a data frame
    cluster.data <- data.frame(cluster.names=clusters.found, count=cluster.counts, size=cluster.sizes, p.vals=p.vals.hyper)
    return( list( "cluster"=cluster, "inds"=inds, "cell.score"=cell.score, "tf"=tf, "idf"=idf, "allcells.indexes"=allcells.bygenename.indexes, "cell.inds"=cell.inds.byexperiment, "cluster.data"=cluster.data ) )
}

getClusterSummary <- function(efdb, inds) {
    #Use this method to get information about cluster representation amongst a set of cells found
    cluster.starts <- cumsum(efdb$cluster.sizes)
    cluster <- c()
    cluster.counts <- rep(0, length(efdb$cluster.names))
    cluster.ind <- c()
    for (i in 1:length(inds)) {
        cluster.ind <- c(cluster.ind, min(which(inds[i]<cluster.starts)))
        cluster <- c(cluster, efdb$cluster.names[tail(cluster.ind,1)])
        cluster.counts[tail(cluster.ind,1)] <- cluster.counts[tail(cluster.ind,1)] + 1
    }
    #calculate the p-value for the different clusters using a hypergeometric test
    tmp <- rle(cluster)
    clusters.found <- tmp$values
    cluster.counts <- tmp$lengths
    p.vals.hyper <- rep(1, length(clusters.found))
    cluster.sizes <- rep(1, length(clusters.found))
    for (i in 1:length(p.vals.hyper)) {
        cluster.sizes[i] <- efdb$cluster.sizes[which(efdb$cluster.names==clusters.found[i])[1]]
      p.vals.hyper[i] <- 1 - phyper(cluster.counts[i], cluster.sizes[i], N - cluster.sizes[i], length(cluster))
    }
    #store summary in a data frame
    cluster.data <- data.frame(cluster.names=clusters.found, count=cluster.counts, size=cluster.sizes, p.vals=p.vals.hyper)
    return ( cluster.data )
}

calculateTFIDFSum <- function(genenames, allcells.bygenename.indexes, tf, idf, cell.inds.byexperiment) {
    cell.score <- list()
    for (i in 1:length(allcells.bygenename.indexes)) { #Keep track of the scores for each experiment separately
        cell.score[[i]] <- rep(0, length(allcells.bygenename.indexes[[i]]))
    }
    for (i in 1:length(genenames)) {
        if (idf[i]>0) {
            for (j in 1:length(cell.inds.byexperiment[[i]])) { 
                #Find the indexes for the cells where this gene is expressed
                tmp.inds <- which(allcells.bygenename.indexes[[j]] %in% cell.inds.byexperiment[[i]][[j]])
                cell.score[[j]][tmp.inds] <- cell.score[[j]][tmp.inds] + tf[[i]][[j]]
            }
        }
    }
    return( cell.score )
}

calculateCosineSimilarity <- function(genenames, allcells.indexes, tf, idf, cell.inds) {
    tf.query <- log(1+rep(1/length(genenames), length(genenames)))
    idf.query <- rep(1, length(genenames))
    tfidf.query <- tf.query*idf.query
    cell.score <- rep(0, length(allcells.indexes))
    denominator <- rep(0, length(allcells.indexes))
    denominator.query <- rep(0, length(allcells.indexes))
    for (i in 1:length(genenames)) {
        if (idf[i]>0) {
            tmp.query <- rep(tfidf.query[i], length(tf[[i]]))
            tmp <- tf[[i]]*idf[i]
        #Find the indexes for the cells where this gene is expressed
            tmp.inds <- which(allcells.indexes %in% cell.inds[[i]])
            cell.score[tmp.inds] <- cell.score[tmp.inds] + tmp*tmp.query
            denominator[tmp.inds] <- denominator[tmp.inds] + tmp^2
            denominator.query <- denominator.query + tfidf.query[i]^2
        }
    }
    return( cell.score/(sqrt(denominator)*sqrt(denominator.query)) )
}

findCellsWithGOtermTFIDF <- function(efdb, goterm, goterm2genename=NULL, n.cells=10, score.method="sum") {
  if (is.null(goterm2genename)) { goterm2genename <- parseGOterms("genenames_go_terms.tsv") }
  #Find the cells associated with the goterm
  return( findBestCellTFIDF(efdb, goterm2genename[[goterm]], n.cells=n.cells, score.method=score.method) )
}

findGOtermsTFIDF <- function(efdb, goterm2genename=NULL, score.method="sum") {
    #Rank the GO terms for each cell
}
