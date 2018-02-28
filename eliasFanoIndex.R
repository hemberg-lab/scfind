#Use these methods to build an index which is compressed using the Elias-Fano strategy to search for cells expressed in specfic genes

buildMmIndex <- function(save.file.name="efdb_mm.rds") {
  #Use this method to build a big database
  efdb.fan <- parseDataset("/nfs/team218/sc-rna-seq-data/fan.rds", "fan")
  efdb.deng <- parseDataset("/nfs/team218/sc-rna-seq-data/deng-reads.rds", "deng")
  efdb.goolam <- parseDataset("/nfs/team218/sc-rna-seq-data/goolam.rds", "goolam")
  efdb.biase <- parseDataset("/nfs/team218/sc-rna-seq-data/biase.rds", "biase")
  efdb.klein <- parseDataset("/nfs/team218/sc-rna-seq-data/klein.rds", "klein")
  efdb.chen <- parseDataset("/nfs/team218/sc-rna-seq-data/chen.rds", "chen")
  efdb.campbell <- parseDataset("/nfs/team218/sc-rna-seq-data/campbell.rds", "campbell", exclude.clusters="miss")
  efdb.kolodziejczyk <- parseDataset("/nfs/team218/sc-rna-seq-data/kolodziejczyk.rds", "kolodziejczyk")
  efdb.macosko <- parseDataset("/nfs/team218/sc-rna-seq-data/macosko.rds", "macosko")
  efdb.shekhar <- parseDataset("/nfs/team218/sc-rna-seq-data/shekhar.rds", "shekhar")

  #efdb.camp <- parseDataset("/nfs/team218/sc-rna-seq-data/camp2.rds", "camp")
  efdb.baron <- parseDataset("/nfs/team218/sc-rna-seq-data/baron-mouse.rds", "baron")

  efdb.manno <- parseDataset("/nfs/team218/sc-rna-seq-data/manno_mouse.rds", "manno")
  efdb.marques <- parseDataset("/nfs/team218/sc-rna-seq-data/marques.rds", "marques")
  efdb.romanov <- parseDataset("/nfs/team218/sc-rna-seq-data/romanov.rds", "romanov")
  efdb.tasic <- parseDataset("/nfs/team218/sc-rna-seq-data/tasic-reads.rds", "tasic")
  efdb.usoskin <- parseDataset("/nfs/team218/sc-rna-seq-data/usoskin.rds", "usoskin")
  efdb.zeisel <- parseDataset("/nfs/team218/sc-rna-seq-data/zeisel.rds", "zeisel")

  size.compressed <- c(object_size(efdb.fan), object_size(efdb.deng), object_size(efdb.goolam), object_size(efdb.biase), object_size(efdb.klein), object_size(efdb.chen), object_size(efdb.campbell), object_size(efdb.kolodziejczyk), object_size(efdb.macosko), object_size(efdb.shekhar), object_size(efdb.baron), object_size(efdb.manno), object_size(efdb.marques), object_size(efdb.romanov), object_size(efdb.tasic), object_size(efdb.usoskin), object_size(efdb.zeisel))
  
  #Do the actual merging
  print("Merging datasets")
  efdb.mm <- mergeDatasets(efdb.fan, efdb.deng)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.goolam)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.biase)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.klein)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.chen)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.campbell)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.kolodziejczyk)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.macosko)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.shekhar)
  #efdb.mm <- mergeDatasets(efdb.mm, efdb.camp)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.baron)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.manno)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.marques)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.romanov)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.tasic)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.usoskin)
  efdb.mm <- mergeDatasets(efdb.mm, efdb.zeisel)
  efdb.mm$size.compressed <- size.compressed
  print(paste("Saving merged database to", save.file.name))
  saveRDS(efdb.mm, save.file.name)
  return( efdb.mm )
}

buildTabulaMurisIndex <- function(save.file.name="efdb_tm.rds", tfidf=4) {
  #Use this method to build a big database
    dir.path <- "/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TabulaMuris/"
    tissues <- c("Bladder_10X", "Bladder_FACS", "Brain_Microglia_FACS", "Brain_Neurons_FACS", "Colon_FACS", "Fat_FACS", "Heart_10X", "Heart_FACS", "Kidney_10X", "Kidney_FACS", "Liver_10X", "Liver_FACS", "Lung_10X", "Lung_FACS", "Mammary_10X", "Mammary_FACS", "Marrow_10X", "Marrow_FACS", "Muscle_10X", "Muscle_FACS", "Pancreas_FACS", "Skin_FACS", "Spleen_10X", "Spleen_FACS", "Thymus_10X", "Thymus_FACS", "Tongue_10X", "Tongue_FACS", "Trachea_10X", "Trachea_FACS")

    efdb <- parseDataset(paste0(dir.path, tissues[1], ".rds"), tabula.muris=T, dataset.name=tissues[1], tfidf=tfidf)
    size.compressed <- object_size(efdb)
    for (i in 2:length(tissues)) {
        efdb2 <- parseDataset(paste0(dir.path, tissues[i], ".rds"), tabula.muris=T, dataset.name=tissues[i], tfidf=tfidf)
        size.compressed <- c(size.compressed, object_size(efdb2))
        efdb <- mergeDatasets(efdb, efdb2)
    }
    efdbsize.compressed <- size.compressed
    print(paste("Saving Tabula Muris database to", save.file.name))
    saveRDS(efdb, save.file.name)
    return( efdb )
}

buildHsIndex <- function(save.file.name="efdb_hs.rds", tfidf=4) {
  #Use this method to build a big database for the human data stored on the farm
  efdb.camp_hs <- parseDataset("/nfs/team218/sc-rna-seq-data/camp1.rds", "camp_hs", tfidf=tfidf)
  efdb.darmanis <- parseDataset("/nfs/team218/sc-rna-seq-data/darmanis.rds", "darmanis", tfidf=tfidf)
  #Exclude since this is single nucleus sequencing
  #efdb.lake <- parseDataset("/nfs/team218/sc-rna-seq-data/lake.rds", "lake", tfidf=tfidf)
  efdb.manno_hs <- parseDataset("/nfs/team218/sc-rna-seq-data/manno_human.rds", "manno_hs", tfidf=tfidf)

  efdb.baron_hs <- parseDataset("/nfs/team218/sc-rna-seq-data/baron-human.rds", "baron_hs", tfidf=tfidf)
  efdb.muraro <- parseDataset("/nfs/team218/sc-rna-seq-data/muraro.rds", "muraro", tfidf=tfidf)
  efdb.segerstolpe <- parseDataset("/nfs/team218/sc-rna-seq-data/segerstolpe.rds", "segerstolpe", tfidf=tfidf, exclude.clusters="not applicable")
  efdb.xin <- parseDataset("/nfs/team218/sc-rna-seq-data/xin.rds", "xin", tfidf=tfidf, exclude.clusters=c("alpha.contaminated", "beta.contaminated", "gamma.contaminated", "delta.contaminated"))

  efdb.pollen <- parseDataset("/nfs/team218/sc-rna-seq-data/pollen.rds", "pollen", tfidf=tfidf)
  efdb.yan <- parseDataset("/nfs/team218/sc-rna-seq-data/yan.rds", "yan", tfidf=tfidf)

  size.compressed <- c(object_size(efdb.camp_hs), object_size(efdb.darmanis), object_size(efdb.manno_hs), object_size(efdb.baron_hs), object_size(efdb.muraro), object_size(efdb.segerstolpe), object_size(efdb.xin), object_size(efdb.pollen), object_size(efdb.yan))

  #Do the actual merging
  print("Merging datasets")
  efdb.hs <- mergeDatasets(efdb.camp_hs, efdb.darmanis)
  #efdb.hs <- mergeDatasets(efdb.hs, efdb.lake)
  efdb.hs <- mergeDatasets(efdb.hs, efdb.manno_hs)
  efdb.hs <- mergeDatasets(efdb.hs, efdb.baron_hs)
  efdb.hs <- mergeDatasets(efdb.hs, efdb.muraro)
  efdb.hs <- mergeDatasets(efdb.hs, efdb.segerstolpe)
  efdb.hs <- mergeDatasets(efdb.hs, efdb.xin)
  efdb.hs <- mergeDatasets(efdb.hs, efdb.pollen)
  efdb.hs <- mergeDatasets(efdb.hs, efdb.yan)
  efdb.hs$size.compressed <- size.compressed
  print(paste("Saving merged database to", save.file.name))
  saveRDS(efdb.hs, save.file.name)
  return( efdb.hs )
}

parseh5Dataset <- function(file.name="/nfs/team218/VK/10x/1M_neurons/data/1M_neurons_filtered_gene_bc_matrices_h5.h5", name="mm10") {
  require("rhdf5")
  H <- hash()
  L <- hash()

  h <- h5read(file.name, name)
  genenames <- h[["gene_names"]]
  n.genes <- length(genenames)
  inds.cell.starts <- h[["indptr"]]
  gene.start <- 1
  m <- rep(0, n.genes)
  H.start <- rep(1, n.genes)
  for (i in 1:length(genenames)) {
    inds.exprs <- c()
    inds.cell <- which(h[["indices"]]==(gene.start-1+i-1))
    prev <- 0
    ind <- 1
    for (j in 1:length(inds.cell)) {
      while((ind<=length(inds.cell.starts)) & (inds.cell.starts[ind]<inds.cell[j])) { ind <- ind + 1 }
      if (ind<=length(inds.cell.starts)) { inds.exprs <- c(inds.exprs, ind - 1) }
    }
    #Now we can store the indices using the Elias-Fano code. Note that we don't have access to cell-types for the 10X dataset, so we just store everything without the same types of pointers that we use for parseDataset below
    u <- length(inds.cell)
    m[i] <- length(inds.exprs)
    if (m[i]>0) {
      l <- floor(log2(u/m[i]))
      C <- eliasFanoCodingCpp(inds.exprs, l) #eliasFanoCoding(inds.exprs, l)
      H[[genenames[i]]] <- as.bit(C$H)
      L[[genenames[i]]] <- as.bit(C$L)
    }
  }
  return( list( "experiment.name"=dataset.name, "assigned.clusters"=assigned.clusters, "m"=m, "permutation"=permutation, "genenames"=genenames, "cluster.sizes"=cluster.sizes, "p0"=p0, "n.clusters"=1, "n.cells"=length(permutation), "cluster.names"=cluster.names, "H.start"=H.start, "H"=H, "L"=L, "size.orig"=0 ) )
}

parseh5DatasetJulia <- function(file.name.h5="/nfs/team218/VK/10x/1M_neurons/data/1M_neurons_filtered_gene_bc_matrices_h5.h5", file.start="/lustre/scratch117/cellgen/team218/MH/1Mh5/1M_", n.files=280, genes.per.file=100) {
  #Use this method to parse the h5 file read by Julia code in rice_inverse_index.jl
    require("rhdf5")
  require('SingleCellExperiment')
  require('bit')
  require("pryr")
  require('hash')
  require('Rcpp')
  sourceCpp('eliasFanoCoding.cpp')
    genenames <- h5read(file.name.h5, "mm10/gene_names")
  H <- hash()
  L <- hash()
  u <- h5read(file.name.h5, "mm10/shape")[2]
    m <- matrix(rep(0, length(genenames)))
  for (i in 1:n.files) {
      print(i)
      file.name <- paste0(file.start, i, ".h5")
    for (j in 1:genes.per.file) {
      ind <- (i-1)*100+j
      #since the H is unary encoding, it is possible to infer m, the number of expressed genes and from that l which is required to figure out the number of low bits
      if (file.exists(file.name)) {
          if (length(which(paste0(genenames[ind],"_H")==h5ls(file.name)[,2]))>0) {
              h <- h5read(file.name, paste0(genenames[ind], "_H"))
              m[ind] <- length(which(h>0))
              H[[genenames[ind]]] <- as.bit(h)
              if (length(h5read(file.name, paste0(genenames[ind], "_L")))>0) {
                  L[[genenames[ind]]] <- as.bit(h5read(file.name, paste0(genenames[ind], "_L")))
              }
          }
      }
    }
    H5close()
  }
    efdb.tenx <- list( "experiment.name"="!0Xneurons", "assigned.clusters"=rep(1, u), "m"=m, "permutation"=1:u, "genenames"=genenames, "cluster.sizes"=c(u), "p0"=(sum(m)/u)/length(m), "n.clusters"=1, "n.cells"=u, "cluster.names"=c("10Xneurons"), "H.start"=matrix(rep(1, length(m))), "H"=H, "L"=L, "size.orig"=0 )
    saveRDS(efdb.tenx, "efdb_tenx.rds")
  return( efdb.tenx )
}

parseDataset <- function(file.name, dataset.name=c(), cell.type.2=F, tabula.muris=F, MCA=F, tfidf=4, exclude.clusters=c()) {
  #This function parses a scater object and encodes it using the Elias-Fano code. Several different pointers are also stored to allow us to skip to the cluster that we are intersted in only. This will allow for faster access at a relatively modest increase of disk use
  require('SingleCellExperiment')
  require('bit')
  require("pryr")
  require('hash')
  require('Rcpp')
  sourceCpp('eliasFanoCoding.cpp')
  
  print(paste("Reading", file.name))
  d <- readRDS(file.name)
  size.orig <- tryCatch( object_size(logcounts(d)), error=function(e) object_size(log2(1+counts(d))))
  if (length(dataset.name)==0) {
  #This bit needs to be made more robust to allow for more flexibility in how file-names are chosen
    tmp <- strsplit(strsplit(file.name, "\\.")[[1]], "/")[[1]]
    dataset.name <- tmp[length(tmp)]
  }
  genenames <- unique(rowData(d)$feature_symbol)
  inds.genes <- which(rowData(d)$feature_symbol %in% genenames)
  if (tabula.muris | MCA) {
      genenames <- unique(rownames(d))
      inds.genes <- which(rownames(d) %in% genenames)
  }
  n.genes <- length(inds.genes)
  n.cells <- dim(d)[2]
  exprs.norm <- c()
  tf.mean <- matrix(rep(0, n.genes), nrow=n.genes)
  tf.sd <- matrix(rep(0, n.genes), nrow=n.genes)
  #store the normalization factor so that we can retrieve the original values
  norm.factor <- matrix(rep(0, n.cells), nrow=n.cells)
  if (tfidf>0) { #calculate the normalized frequency of each gene for each cell, a positive number indicates the number of bins that we use for storing the expression
      if (tabula.muris | MCA) { exprs.norm <- log2(1+counts(d)) }
      else {       exprs.norm <- logcounts(d) }
      for (i in 1:dim(d)[2]) {
          #normalize for the sequencing depth as well by multiplying by nGenes/10000
          norm.factor[i] <- 10^(length(which(exprs.norm[,i]>0))/10000)*10000/sum(exprs.norm[,i])
          exprs.norm[,i] <- exprs.norm[,i]*norm.factor[i]
      }
      #For storing the expression levels, we need to carry out a transformation so that we can make the distribution more uniform and thereby maximizing the dynamic range and the information content for the encoding.
      exprs.norm.log.inv <- log10(1/(exprs.norm))
  }
  #Group cells by cell-type
  cell.types.all <- colData(d)$cell_type1
  if (cell.type.2) { cell.types.all <- unique(colData(d)$cell_type2) }
  if (tabula.muris) { cell.types.all <- colData(d)$type }
  if (MCA) { cell.types.all <- colData(d)$ClusterID }
  cell.types <- setdiff(unique(cell.types.all), exclude.clusters)
  print(paste0("Found ", length(cell.types), " clusters."))
  #For keeping track of the parameters for the Elias-Fano code. Use l specific to each gene and cell-type
  m <- matrix(rep(0, n.genes*length(cell.types)), nrow=n.genes)
  #Keep track of how the cells were permuted
  permutation <- c()
  p0 <- c()
  #Keep track of the cluster for each cell. Use a numeric index to save space
  assigned.clusters <- rep(0, dim(d)[2])
  for (i in 1:dim(d)[2]) {
      tmp <- which(cell.types==cell.types.all[i])
      if (length(tmp)>0) { assigned.clusters[i] <- tmp }
  }
  #Keep track of the number of cells in the cluster
  cluster.sizes <- rep(0, length(cell.types))
  cluster.names <- rep(0, length(cell.types))
  #Keep track of where the next cell-type starts for the H array for faster access
  H.start <- matrix(rep(1, n.genes*length(cell.types)), nrow=n.genes)
  #Use a hash to keep track of the codes from each gene
  H <- hash()
  L <- hash()
  tf <- hash()
  for (i in 1:length(cell.types)) {
      print(paste("Indexing", cell.types[i]))
      inds.cell <- which(cell.types[i]==cell.types.all)
      u <- length(inds.cell)
      cluster.sizes[i] <- u
      cluster.names[i] <- paste0(dataset.name, "_", cell.types[i])
      permutation <- c(permutation, inds.cell)
      #Calculate the baseline probability that a gene will be expressed in a cell
      exprs <- c()
      if (tabula.muris | MCA) { exprs <- log2(1+counts(d)) }
      else {       exprs <- logcounts(d) }
      p0 <- c(p0, sum(exprs[,inds.cell]>0)/(n.genes*u))
      for (j in 1:n.genes) {
          if (i>1) { H.start[j,i] <- length(H[[genenames[j]]]) + 1 }
          inds.logcounts <- which(exprs[inds.genes[j],inds.cell]>0)
          m[j,i] <- length(inds.logcounts)
          if (m[j,i]>0) {
              l <- floor(log2(u/m[j,i]))
              C <- eliasFanoCodingCpp(inds.logcounts, l) #eliasFanoCoding(inds.logcounts, l)
              if (has.key(genenames[j], H)) {
                  H[[genenames[j]]] <- c(H[[genenames[j]]], as.bit(C$H))
                  L[[genenames[j]]] <- c(L[[genenames[j]]], as.bit(C$L))
              } else {
                  H[[genenames[j]]] <- as.bit(C$H)
                  L[[genenames[j]]] <- as.bit(C$L)
              }
              if (tfidf>0) {
                  tf.mean[j] <- mean(exprs.norm.log.inv[inds.genes[j], which(exprs.norm[inds.genes[j],]>0)])
                  tf.sd[j] <- sd(exprs.norm.log.inv[inds.genes[j], which(exprs.norm[inds.genes[j],]>0)])
                  #Store the quantized gene expression in the same order as in the H/L arrays
                  tmp <- as.bit(quantizeCpp(pnorm(exprs.norm.log.inv[inds.genes[j],inds.cell[inds.logcounts]], tf.mean[j], tf.sd[j]), tfidf))
                  if (has.key(genenames[j], tf)) { tf[[genenames[j]]] <- c(tf[[genenames[j]]], tmp) }
                  else { tf[[genenames[j]]] <- tmp }
              }
              #else { tf[[genenames[j]]] <- c() }
          }
      }
  }
  return( list( "experiment.name"=dataset.name, "assigned.clusters"=assigned.clusters, "m"=m, "permutation"=permutation, "genenames"=genenames, "cluster.sizes"=cluster.sizes, "p0"=p0, "n.clusters"=length(cluster.sizes), "n.cells"=length(permutation), "cluster.names"=cluster.names, "H.start"=H.start, "H"=H, "L"=L, "size.orig"=size.orig, "tf"=tf, "tf.quantization"=tfidf, "tf.mean"=tf.mean, "tf.sd"=tf.sd, "norm.factor"=norm.factor ) )
}

quantize <- function(x, n) {
    tmp <- as.bit(intToBitVect(floor(x*2^n-1e-10)))
    #Pad so that we use n bits
    return( c(as.bit(rep(FALSE, n-length(tmp))), tmp) )
}

mergeDatasets <- function(efdb1, efdb2) {
    #For now, we are going to only use those genes that are present in both the old and the new
  efdb <- c()
  efdb$H <- hash()
  efdb$L <- hash()
  efdb$tf <- hash()
  efdb$tf.quantization <- c(efdb1$tf.quantization, efdb2$tf.quantization)
  genenames.intersect <- intersect(efdb1$genenames, efdb2$genenames)
  #inds.efdb1 <- which(efdb1$genenames %in% genenames.intersect)
  #inds.efdb2 <- which(efdb2$genenames %in% genenames.intersect)
  #efdb$m <- cbind(efdb1$m[inds.efdb1,], efdb2$m[inds.efdb2,])
  efdb$experiment.name <- c(efdb1$experiment.name, efdb2$experiment.name)
  efdb$cluster.sizes <- c(efdb1$cluster.sizes, efdb2$cluster.sizes)
  efdb$n.clusters <- c(efdb1$n.clusters, efdb2$n.clusters)
  efdb$n.cells <- c(efdb1$n.cells, efdb2$n.cells)
  efdb$H.start <- matrix(rep(0, length(genenames.intersect)*sum(efdb$n.clusters)), ncol=sum(efdb$n.clusters))
  efdb$m <- matrix(rep(0, length(genenames.intersect)*sum(efdb$n.clusters)), ncol=sum(efdb$n.clusters))
  efdb$tf.mean <- matrix(rep(0, length(genenames.intersect)*length(efdb$experiment.name)), ncol=length(efdb$experiment.name))
  efdb$tf.sd <- matrix(rep(0, length(genenames.intersect)*length(efdb$experiment.name)), ncol=length(efdb$experiment.name))
  efdb$assigned.clusters <- c(efdb1$assigned.clusters, sum(efdb1$n.clusters) + efdb2$assigned.clusters)
  efdb$cluster.names <- c(efdb1$cluster.names, efdb2$cluster.names)
  efdb$permutation <- c(efdb1$permutation, efdb2$permutation)
  efdb$p0 <- c(efdb1$p0, efdb2$p0)
  efdb$norm.factor <- c(efdb1$norm.factor, efdb2$norm.factor)
  efdb$size.orig <- c(efdb1$size.orig, efdb2$size.orig)
  efdb$genenames <- genenames.intersect
  for (i in 1:length(genenames.intersect)) {
    ind.efdb1 <- which(efdb1$genenames==genenames.intersect[i])
    ind.efdb2 <- which(efdb2$genenames==genenames.intersect[i])
    efdb$m[i,] <- c(efdb1$m[ind.efdb1,], efdb2$m[ind.efdb2,])
    efdb$tf.mean[i,] <- c(efdb1$tf.mean[ind.efdb1,], efdb2$tf.mean[ind.efdb2,])
    efdb$tf.sd[i,] <- c(efdb1$tf.sd[ind.efdb1,], efdb2$tf.sd[ind.efdb2,])
    if (has.key(genenames.intersect[i], efdb1$H) & has.key(genenames.intersect[i], efdb2$H)) {
      efdb$H[[genenames.intersect[i]]] <- c(efdb1$H[[genenames.intersect[i]]], efdb2$H[[genenames.intersect[i]]])
      efdb$L[[genenames.intersect[i]]] <- c(efdb1$L[[genenames.intersect[i]]], efdb2$L[[genenames.intersect[i]]])
      efdb$H.start[i,] <- c(efdb1$H.start[ind.efdb1,], length(efdb1$H[[genenames.intersect[i]]]) + efdb2$H.start[ind.efdb2,])
      efdb$tf[[genenames.intersect[i]]] <- c(efdb1$tf[[genenames.intersect[i]]], efdb2$tf[[genenames.intersect[i]]])
  }
    else if (has.key(genenames.intersect[i], efdb1$H) & !has.key(genenames.intersect[i], efdb2$H)) {
      efdb$H[[genenames.intersect[i]]] <- efdb1$H[[genenames.intersect[i]]]
      efdb$L[[genenames.intersect[i]]] <- efdb1$L[[genenames.intersect[i]]]
      efdb$H.start[i,] <- c(efdb1$H.start[ind.efdb1,], rep(length(efdb1$H[[genenames.intersect[i]]]), sum(efdb2$n.clusters)))
      efdb$tf[[genenames.intersect[i]]] <- efdb1$tf[[genenames.intersect[i]]]
    }
    else if (!has.key(genenames.intersect[i], efdb1$H) & has.key(genenames.intersect[i], efdb2$H)) {
      efdb$H[[genenames.intersect[i]]] <- efdb2$H[[genenames.intersect[i]]]
      efdb$L[[genenames.intersect[i]]] <- efdb2$L[[genenames.intersect[i]]]
      efdb$H.start[i,] <- c(rep(1, sum(efdb1$n.clusters)), efdb2$H.start[ind.efdb2,])
      efdb$tf[[genenames.intersect[i]]] <- efdb2$tf[[genenames.intersect[i]]]
    }
  }
  return( efdb )
}

parseGOterms <- function(filename) {
  #Use this method to parse the list of GO terms so that they can be used for searches
  d <- read.table(filename)
  goterm2genename <- hash()
  for (i in 1:dim(d)[1]) {
    if (!has.key(as.character(d[i,2]), goterm2genename)) { goterm2genename[[d[i,2]]] <- c() }
    goterm2genename[[as.character(d[i,2])]] <- c(goterm2genename[[as.character(d[i,2])]], as.character(d[i,1]))
  }
  return( goterm2genename )
}

##################################################################

numberOfCells <- function(efdb, genename) {
  ind <- which(efdb$genenames==genename)
  return( sum(efdb$m[ind,]) )
}

findCellsWithGenes <- function(efdb, genelist, experiments.allowed=c(), no.p.vals=F) {
  #Use this method to run the AND query on the big database
  #find the number of cells for each gene, start by searching the smallest ones first
  require('bit')
  require('hash')
  require('Rcpp')
  sourceCpp('eliasFanoCoding.cpp')
  n <- rep(0, length(genelist))
  for (i in 1:length(genelist)) { n[i] <- numberOfCells(efdb, genelist[i]) }
  order <- sort.int(n, index.return=T)$ix
  ids.unique <- NULL
  ids.original <- c()
  n.genes.not.found <- 0
  clusters.found <- 1:length(efdb$cluster.sizes)
  if (length(order)>=1) {
    for (i in 1:length(order)) {
      tmp <- eliasFanoDecodingGene(efdb, genelist[order[i]], experiments.allowed, clusters.found)
      if (!is.null(tmp)) {
        if (is.null(ids.unique)) {
          inds.keep <- tmp$ids.unique
          ids.unique <- tmp$ids.unique
          ids.original <- tmp$ids.original
        }
        else {
          inds.keep <- which(ids.unique %in% tmp$ids.unique)
          ids.unique <- ids.unique[inds.keep] #intersect(ids.unique, tmp$unique) #
          ids.original <- ids.original[inds.keep]
        }
        clusters.found <- tmp$clusters.found
        if (length(ids.unique)==0) { break }
      }
      else { n.genes.not.found <- n.genes.not.found + 1 }
    }
  }
  if (no.p.vals) { return( ids.unique ) }
  return( calculatePValues(efdb, ids.unique, length(genelist) - n.genes.not.found) )
}

calculatePValues <- function(efdb, ids.unique, n.genes.eff) {
  if (length(ids.unique)>0) { #Now we need to find out which clusters these cells come from
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
  return( list( "ids.unique"=ids.unique, "clusters"=c(), "clusters.found"=c(), "cluster.counts"=c(), "cluster.names"=c(), "p.vals"=c() ) )
}

findExperiments <- function(ids, cluster.sizes, n.clusters, experiment.name) {
  #Use this method to find out which experiments are represented amongst ids. This information can be used to speed up conjunctive searches since we do not even have to decode the indices corresponding to those cell-types
  experiments.found <- c()
  experiment.ind <- 1
  n.clusters.cumsum <- cumsum(n.clusters)
  for (i in 1:length(ids)) {
    if (length(which(experiments.found==experiment.name[experiment.ind]))==0) { experiments.found <- c(experiments.found, experiment.name[experiment.ind]) }
    while (ids[i]>n.clusters.cumsum[experiment.ind]) {
      experiment.ind <- experiment.ind + 1
    }        
  }
  return( experiments.found )
}

findCellsWithGOterm <- function(efdb, goterm, goterm2genename=NULL, mismatches.allowed=0) {
  if (is.null(goterm2genename)) { goterm2genename <- parseGOterms("genenames_go_terms.tsv") }
  #Find the cells associated with the goterm
  if (mismatches.allowed==0) {  return( findCellsWithGenes(efdb, goterm2genename[[goterm]]) ) }
  return( findCellsWithGenesMismatch(efdb, goterm2genename[[goterm]], mismatches.allowed) )
}

findCellsWithGenesMismatch <- function(efdb, genelist, mismatches.allowed, experiments.allowed=c()) {
  #Same as before, except that we now allow for up to n.mismatches.allowed genes from genelist to be missing for each cell. This method is very slow which is most likely due to the two sets of for-loops. If this is to be used in practice it needs to be sped up dramatically.
  cells2genes <- rep(0, sum(efdb$n.cells)) #Keep track of the number of genes found for each cell
  n.genes.not.found <- 0
  for (i in 1:length(genelist)) {
    f <- eliasFanoDecodingGene(efdb, genelist[i], experiments.allowed, 1:length(efdb$cluster.sizes))
    if (is.null(f)) {
        n.genes.not.found <- n.genes.not.found + 1
        print(paste0("Could not find ", genelist[i]))
    }
    else {
      if (length(f$ids.unique)>0) { #update the array that keeps track of how many genes were present in each of the cells so far
        cells2genes[f$ids.unique] <- cells2genes[f$ids.unique] + 1
      }
    }
  }
  #Find the cells that have the minimal number of genes required
  min.genes <- 0
  if (mismatches.allowed<1) {
      min.genes <- round(mismatches.allowed*(length(genelist) - n.genes.not.found))
  }
  else { min.genes <- length(genelist) - mismatches.allowed - n.genes.not.found }
  #print(paste("Requiring", as.character(min.genes), "for cell to be included."))
  ids.unique <- which(cells2genes>min.genes)
  return( calculatePValues(efdb, ids.unique, min.genes) )
}

findCellTypeMatchingGOterm <- function(efdb, goterms, goterm2genename=NULL, mismatches.allowed=0) {
  if (is.null(goterm2genename)) { goterm2genename <- parseGOterms("genenames_go_terms.tsv") }
  #Find the fraction of terms for which each type is considered significant
  cluster.terms.count <- rep(0, length(efdb$cluster.sizes))
  cluster.terms.pvals <- c()
  cluster.terms.fractions <- c()
  for (i in 1:length(goterms)) {
      ret <- findCellsWithGenesMismatch(efdb, goterms[i], goterm2genename, mismatches.allowed)
      inds <- which(ret$p.vals.binom<.05/length(efdb$cluster.sizes))
      if (length(inds)>0) {
          cluster.terms.count[inds] <- cluster.terms.count[inds] + 1
          for (j in 1:length(inds)) {
              cluster.terms.pvals[[inds[j]]] <- c(cluster.terms.pvals[[inds[j]]], ret$p.vals.binom[inds[j]])
          }
      }
  }
  #return( )
}

##################################################################


findAssociationRulesCellsFPgrowthOld <- function(efdb, min.support=.5) {
    #For this method we should also add the option to filter for the experiments and clusters that one is interested in so that one does not need to run on the whole database.
    #This method would strongly benefit from Rcpp treatment since it is quite slow
    require("data.tree")
    sourceCpp("eliasFanoCoding.cpp")
    #Calculate the support, ie the fraction of cells where each gene appears
    nCells <- sum(efdb$cluster.sizes)
    genenames <- efdb$genenames
    support <- rowSums(efdb$m)/nCells
    ret <- sort.int(support, decreasing=T, index.return=T)
    support.order <- ret$ix[which(ret$x>min.support)]
    #Now we need to build a transposed datastructure to be able to search by cell instead. Since we have been able to remove a lot of genes that are cell-type specific, the matrix should not be too large. However, if nCells is very large then we may need to store it in a more compact manner
    print(paste0("Extracting frequency matrix for ", as.character(length(support.order)), " genes."))
    mat.fpgrowth <- extractFrequencyMatrix(efdb, nCells, genenames[support.order]);
    print("Building FP tree")
    fp <- buildFPTree(mat.fpgrowth, round(nCells*min.support))
    #Now, we can build a tree in R and use the proper gene names
    #fpgrowth.root <- Node$new("fpgrowth.root")
    #For the FP growth algorith, use a tree structure
    #fpgrowth.root <- Node$new("fpgrowth.root")
    ##Go through the cells and add genes in the order dictated by the support
    #for(i in 1:nCells) {
    #    tmp.n <- fpgrowth.root
    #    for (j in 1:length(support.order)) {
    #        if (mat[i,j]>0) {
    #            childInd <- hasChild(tmp.n$children, genenames[support.order[j]])
    #            if (childInd>0) {
    #                tmp.n$children[[childInd]]$freq <- tmp.n$children[[childInd]]$freq + 1
    #                tmp.n <- tmp.n$children[[childInd]]
    #            }
    #            else {
    #                tmp.n <- tmp.n$AddChild(genenames[support.order[j]], freq=1)
    #            }
    #        }
    #    }
    #}
    ##Now we can prune the tree by removing the nodes where the frequency is too low
    #Prune(fpgrowth.root, function(x) x$freq > min.support*nCells)
    return( list( "fp" = fp, "genenames" = genenames[support.order] ) ) #fpgrowth.root )
}

findMarkerGenes <- function(efdb, cluster.ind, max.genes=20, cross.val.fold=5) {
    #Use to identify the set of genes that will be best as marker genes
    sourceCpp("eliasFanoCoding.cpp")
    nGenes <- length(efdb$genenames)
    nCells <- sum(efdb$cluster.sizes)
    #Find the indices of the cell-type that we are interested in
    cell.inds <- efdb$permutation[1:efdb$cluster.sizes[1]]
    if (cluster.ind>1) {
        cell.inds <- efdb$permutation[cumsum((1+efdb$cluster.sizes)[cluster.ind]):efdb$cluster.sizes[cluster.ind-1]]
    }
    cell.mask <- rep(-1, nCells)
    cell.mask[cell.inds] <- 1
    #Calculate the overlap of each gene with the cell.mask vector
    print(paste("Calculating the overlap of each gene with the cell-type of interest"))
    support.gene <- rep(0, length(efdb$genenames))
    gene.mask <- matrix(rep(-1, nCells*nGenes), nCells)
    mask <- matrix(rep(-1, nCells*nGenes), nCells)
    for (i in 1:length(efdb$genenames)) {
        gene.mask[findCellsWithGenes(efdb, efdb$genenames[i], no.p.vals=T),i] <- 1
        support.gene[i] <- sum(cell.mask==gene.mask[,i])
        mask[cell.mask==gene.mask[,i], i] <- 1
    }
    #If there are cells that we would like to ignore, then we should set the mask to 0
    tp <- matrix( rep(0, max.genes*cross.val.fold), nrow=cross.val.fold)
    tn <- matrix( rep(0, max.genes*cross.val.fold), nrow=cross.val.fold)
    fp <- matrix( rep(0, max.genes*cross.val.fold), nrow=cross.val.fold)
    fn <- matrix( rep(0, max.genes*cross.val.fold), nrow=cross.val.fold)
    if (cross.val.fold==1) {
        for (i in 1:length(efdb$genenames)) { support.gene[i] <- sum(cell.mask==gene.mask[,i]) }
        ret <- sort.int(support.gene, decreasing=T, index.return=T)
        support.order <- ret$ix[which(ret$x>min.support*nCells)]
        for (i in 1:max.genes) {
            ids.found <- findCellsWithGenesMismatch(efdb, efdb$genenames[support.order[1:i]], i/max.genes)$ids.unique
            tp[i] <- length(which(cell.inds %in% ids.found))
            fn[i] <- length(cell.inds) - tp[i]
            fp[i] <- length(ids.found) - tp[i]
            tn[i] <- nCells - tp[i] - fn[i] - fp[i]
        }
    }
    else {
        for (j in 1:cross.val.fold) {
            tmp <- sample(nCells)
            cells.id.test <- tmp[1:round(nCells/cross.val.fold)]
            cells.id.train <- tmp[(1+round(nCells/cross.val.fold)):nCells]
            cell.inds.test <- intersect(cells.id.test, cell.inds) #cells of interest in the test set
            for (i in 1:length(efdb$genenames)) { support.gene[i] <- sum(cell.mask[cells.id.train]==gene.mask[cells.id.train,i]) }
            ret <- sort.int(support.gene, decreasing=T, index.return=T)
            support.order <- ret$ix[which(ret$x>min.support*length(cells.id.train))]
            for (i in 1:max.genes) {
                ids.found <- findCellsWithGenesMismatch(efdb, efdb$genenames[support.order[1:i]], i/max.genes)$ids.unique
                tp[j, i] <- length(intersect(cell.inds.test, ids.found))
                fn[j, i] <- length(cell.inds.test) - tp[j, i]
                fp[j, i] <- length(intersect(cells.id.test, ids.found)) - tp[j, i]
                tn[j, i] <- length(cells.id.test) - tp[j, i] - fn[j, i] - fp[j, i]
            }
        }
    }
    recall <- tp/length(cell.inds.test)
    precision <- tp/(tp+fp)
    fpr <- fp/(nCells - length(cell.inds))
    f1 <- 2/(recall^-1 + precision^-1)
    return( list( "fp" = fp, "genenames" = genenames[support.order] ) ) #fpgrowth.root )
}


extractFrequencyMatrix <- function(efdb, nCells, genenames) {
    mat <- matrix(rep(0, nCells*length(genenames)), nrow=nCells) #use bit class instead to save memory
    for (i in 1:length(genenames)) { mat[findCellsWithGenes(efdb, genenames[i], no.p.vals=T),i] <- 1  }
    return( mat );
}

hasChild <- function(children, genename) {
    if( is.null(children) ) { return( -1 ) }
    for (i in 1:length(children)) { if (children[[i]]$name==genename) { return( i ) } }
    return( 0 )
}

findAssociationRulesChi2 <- function(efdb, min.support=.5, alpha=.05, min.chi2=100) {
    #Use the Brin et al (1997) method for identifying correlated pairs of positive and negative relationships
    #Calculate the support, ie the fraction of cells where each gene appears
    nCells <- sum(efdb$cluster.sizes)
    genenames <- efdb$genenames
    support <- rowSums(efdb$m)
    support.inds <- which(support>min.support*nCells)
    mat <- extractFrequencyMatrix(efdb, nCells, efdb$genenames[support.inds])
    significant.pairs <- c()
    for (i in 1:length(support.inds)) {
        for (j in (i+1):length(support.inds)) {
            #construct a contingency table
            n <- sum(mat[,i]*mat[,j])
            t <- rbind(c(n, support[support.inds[i]] - n), c(support[support.inds[j]] - n, nCells - support[support.inds[i]] - support[support.inds[j]] + n))
            if (length(which(t>nCells*min.support))>=1) {
                chi2 <- chisq.test(t)
                if (!is.nan(chi2$p.value) & chi2$p.value<alpha*2/(length(support.inds)*(length(support.inds)-1)) & chi2$statistic>min.chi2) {
                    significant.pairs <- rbind(significant.pairs, c(genenames[support.inds[c(i, j)]], chi2$statistic))
                }
            }
        }
    }
    return( significant.pairs )
}

#######################################################

findBiclustersFPgrowth <- function(efdb, minGenes=10, minCells=100) {
    #This method is inspired by the method by Buerher et al for finding biclusters based on the FP growth algorithm. The nice feature of this method is that it is linear in the size of the database and that it only needs to make one pass over the collection (provided that we have pre-computed the marginals). The approximation can also be used for boolean matrix factorization
    require("Rcpp")
    mat <- t(extractFrequencyMatrix(efdb, sum(efdb$cluster.sizes), efdb$genenames))
    genes.order <- sort.int(rowSums(mat), decreasing=T, index.return=T)$ix
    cells.order <- sort.int(colSums(mat), decreasing=T, index.return=T)$ix
    m <- mat[genes.order, cells.order]
    interesting.patterns <- buildFPTreeT(m, minGenes, minCells)
    return( interesting.patterns )
}

booleanMatrixFactorizationQuality <- function(m, interesting.patterns) {
    p <- length(interesting.patterns)
    n1 <- matrix(rep(0, p*dim(m)[1]), nrow=dim(m)[1])
    n2 <- matrix(rep(0, p*dim(m)[2]), nrow=p)
    for (i in 1:p) {
        n1[interesting.patterns[[i]]$genes,i] <- 1
        n2[i,interesting.patterns[[i]]$cells] <- 1
    }
    n <- 1*((n1%*%n2)>0)
    print(sum(n==m)/prod(dim(m)))
    return( n )
}
    
findDecisionTreeID3 <- function(efdb, celltype) {
    #Use the ID3 algorithm to find a decision tree that can be used to identify celltype

    #Find which cellIDs correspond to the celltype of interest

    
}

informationGain <- function(efdb, genename, celltype, cell.ids) {
    #get the entropy of all the cells when considering the binary distinction of celltype vs not celltype
    p.celltype

}
    
#######################################################

greedyPack <- function(efdb) {
    #Implement the GreedyPack method by Tatti and Vreeken (2008) which can identify representative patterns based on both presence and absence

    #Generate the trivial trees for each gene

    #The graph that ensures that we will have proper decision trees
    
    while (T) {
        for (i in 1:length(trees)) {
            cand[[i]] <- trees[[i]]
            leaves <- leaves(cand)
        }
    }
    for (i in 1:length(trees)) {
        costDiff[i] <- treeCost(cand[[i]]) - treeCost(trees[[i]])
    }
    m <- min(costDiff)
    k <- which(costDiff==m)
    if (m<0) {

    }
}

treeCost <- function() {
    #Use this method to compute the cost of encoding a tree t
    if (isLeaf()) { return( 1 + 0 ) }
    else {
    }
    return( c )
}

costMDL <- function(M) {
    c <- rep(0, M+1)
    for (k in 0:M) {
        c[k+1] <- choose(M, k)*(k/M)^k*((M-k)/M)^(M-k)
    }
    return( log2(sum(c)) )
}

bmfasso <- function(efdb, k, tau) {
    #Use the Asso algorithm to carry out boolean matrix factorization
    nCells <- sum(efdb$cluster.sizes)
    nGenes <- length(efdb$genenames)
    B <- matrix(rep(0, nGenes*k), nrow=nGenes)
    C <- matrix(rep(0, nCells*k), nrow=k)
    Y <- matrix(rep(0, nGenes^2), nrow=nGenes)
    for (i in 1:nGenes) {
        inds1 <- findCellsWithGenes(efdb, efdb$genenames[i], no.p.vals=T)
        for (j in 1:nGenes) {
            inds2 <- findCellsWithGenes(efdb, efdb$genenames[c(i,j)], no.p.vals=T)
            Y[j,i] <- length(inds2)/length(inds1)
        }
    }
    X <- 1*(Y>=tau)
    for (l in 1:k) {
        
    }
}

cover <- function(efdb, B, C) {
    #Calculate the number of correct elements covered in the approximation A = B*C minus the number of incorrect elements
    for (i in 1:nGenes) {

    }
}

##################################################################

intToBitVect <- function(x){
  tmp <- rev(as.bit(as.numeric(intToBits(x))))
  id <- seq_len(match(1,tmp,length(tmp))-1)
  tmp[-id]
}

eliasFanoCoding <- function(ids, l) {
  require('bit')
  L <- bit(0)
  H <- bit(0)
  ids <- c(0, ids);
  for (i in 2:length(ids)) {
    c <- intToBitVect(ids[i]) 
    if (length(c)<l) { c <- c(bit(l-length(c)), c) }
    L <- c(L, c[(length(c)-l+1):length(c)])
    #c <- intToBin(ids[i])
    #if (nchar(c)<l) { c <- c(bit(l-length(c)), c) }
    #now use a unary code for the high bits, ie floor(ids[i]/2^l) - floor(ids[i-]/2^l)
    m = floor(ids[i]/2^l) - floor(ids[i-1]/2^l)
    h <- bit(m+1)
    h[m+1] <- T
    H <- c(H, h)
  }
  return( list( "H"=H, "L"=as.bit(L) ) )
}

#######################################################
  
eliasFanoDecodingGene <- function(efdb, genename, experiments.allowed=c(), clusters.found) {
  ind <- which(efdb$genenames==genename)
  if (length(ind)==0) {
    print(paste("Warning, could not find gene", genename))
    return( NULL )
  }
  ids <- eliasFanoDecodingClusters(efdb$H[[genename]], efdb$L[[genename]], efdb$m[ind,], efdb$H.start[ind,], efdb$cluster.sizes, efdb$n.clusters, efdb$permutation, efdb$experiment.name, experiments.allowed, clusters.found)
  return( ids )
}

eliasFanoDecodingClusters <- function(Hs, Ls, ms, i.h, ns, n.clusters, permutation, experiment.name, experiments.allowed=c(), clusters.found) {
  #Need to figure out how to decode H and L since the parameter l varies. We know that there are ls[i]*ms[i] bits used in L for the ith entry. Similarly, we know that we need to read until we have encountered ms[i] 1s in H.
  ids.original <- c()
  ids.unique <- c()
  i.l <- 1
  ids.offset <- 0
  ids.offset.dataset <- 0 #To make sure that the ids are unique globally
  #Keep track of which experiment we are evaluating. This is useful since we may want to check if the experiment is on the list of the ones that we have restricted ourselves to.
  experiment.ind <- 1
  n.clusters.cumsum <- cumsum(n.clusters)
  experiment <- experiment.name[experiment.ind]
  for (i in 1:length(ns)) {
    if (ms[i]>0 & length(which(clusters.found==i))>0) {
      #Check if the experiment is on the white-list
      if ((length(experiments.allowed)==0) | (length(which(experiments.allowed==experiment))>0)) {
        l <- floor(log2(ns[i]/ms[i]))
        #Now we need to look in the permutation vector to figure out the original ID for the cell
        #tmp <- permutation[ids.offset + eliasFanoDecoding(Hs, Ls, ms[i], l, i.h[i], i.l)]
        tmp <- eliasFanoDecodingCpp(as.numeric(Hs), as.numeric(Ls), ms[i], l, i.h[i]-1, i.l-1)
        #if (i==1) { print(tmp) ; print(permutation[tmp]) }
        ids.original <- c(ids.original, permutation[ids.offset + tmp])
        ids.unique <- c(ids.unique, ids.offset.dataset + permutation[ids.offset + tmp])
        #Figure out the index where we should start in the L vector
        i.l <- i.l + l*ms[i]
      }
    }
    else { clusters.found <- setdiff(clusters.found, i) }
    ids.offset <- ids.offset + ns[i]
    if (i>=n.clusters.cumsum[experiment.ind]) {
      experiment.ind <- experiment.ind + 1
      experiment <- experiment.name[experiment.ind]
      ids.offset.dataset <- ids.offset
    }
  }
  return( list( "ids.unique"=ids.unique, "ids.original"=ids.original, "clusters.found"=clusters.found ) )
}
  
#This function should also be re-written using Rcpp
eliasFanoDecoding <- function(H, L, m, l, i.h=1, i.l=1) {
  ids <- rep(0, m)
  nZeros <- 0
  n.cells.found <- 0
  j <- 0
  i <- i.h
  prevH <- 0
  while (n.cells.found<m) {
    if (!H[i]) { nZeros <- nZeros + 1 }
    else {
      if (l>0) { #Calculate the value stored in the low bits
        for (k in 1:l) {
          if (L[i.l-1+j*l+k]) { ids[n.cells.found+1] <- ids[n.cells.found+1] + 2^(l-k) }
        }
      }
      #Calculate the value stored in the high bits
      h = intToBitVect(nZeros + prevH)
      for (k in 1:length(h)) {
        if (h[k]) { ids[n.cells.found+1] <- ids[n.cells.found+1] + 2^(length(h)-k+l) }
      }
      prevH <- floor(ids[n.cells.found+1]/2^l)
      j <- j + 1
      nZeros <- 0
      n.cells.found <- n.cells.found + 1
    }
    i <- i + 1
  }
  return( ids )
}
  
  eliasFanoDecoding2 <- function(H, L, l) {
    n <- length(L)/l
    ids <- rep(0, n)
  nZeros <- 0
  j <- 1
  i <- 1
  prevH <- 0
  while (i<=length(H)) {
    if (!H[i]) { nZeros <- nZeros + 1 }
    else {
      #Calculate the value stored in the low bits
      for (k in 1:l) {
        if (L[(j-1)*l+k]) { ids[j] <- ids[j] + 2^(l-k) }
      }
      #Calculate the value stored in the high bits
      h = intToBitVect(nZeros + prevH)
      for (k in 1:length(h)) {
        if (h[k]) { ids[j] <- ids[j] + 2^(length(h)-k+l) }
      }
      prevH <- floor(ids[j]/2^l)
      j <- j + 1
      nZeros <- 0
    }
    i <- i + 1
  }
  return( ids )
}
