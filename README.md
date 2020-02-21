<img src=https://genat.uk/img/scfind2_colour.png height="200">

## scfind - Fast searches of large collections of single cell data

Single cell technologies have made it possible to profile millions of cells, but for these resources to be useful they must be easy to query and access. To facilitate interactive and intuitive access to single cell data we have developed scfind (source available at https://github.com/hemberg-lab/scfind), a __search engine__ for cell atlases. Scfind can be used to evaluate marker genes, to perform in silico gating, and to identify both cell-type specific and housekeeping genes. An interactive interface website with 9 single cell datasets is available at https://scfind.sanger.ac.uk. 

__Q__: What is this?

__A__: __scfind__ is a search engine that makes single cell data accessible to a wide range of users by enabling sophisticated queries for large datasets through an interface which is both very fast and familiar to users from any background.

__Q__: How to install/run __scfind__?

__A__: If you would like to install the latest development version of scfind please install it from the GitHub repository:

```
# Linux and Mac users, run this in your R session:
install.packages("devtools")
devtools::install_github("hemberg-lab/scfind")

library("scfind")

# For Windows users:
# Please install the latest version of Rtools at https://cran.r-project.org/bin/windows/Rtools/ prior to installation of scfind
```

__Update__ The latest version (`3.5.1`) of __scfind__ released on 3rd October 2019 has provided 2 datasets and 2 pre-processed __scfind__ indexes as example. To update the latest version:

```
install.packages("devtools")
devtools::install_github("hemberg-lab/scfind", force = TRUE)
```

__Q__: Where can I find the `scfind` example datasets and indexes?

__A__: The latest version of the package has provided a list of example `SingleCellExperiment` objects and __scfind__ indexes created from the [The Tabula Muris Consortium](https://doi.org/10.1038/s41586-018-0590-4) for your first __scfind__ experience:

```
library("scfind")

# List of `Tabula Muris (FACS)` `SingleCellExperiment` objects
data(tmfacs)

# List of `Tabula Muris (10X)` `SingleCellExperiment` objects
data(tm10x)
```

The detail of building __scfind__ index from `SingleCellExperiment` object is described in [this page](https://github.com/hemberg-lab/scfind/blob/master/Vignettes/scfind.Rmd). 

```
library("scfind")
library("SingleCellExperiment")

# To build the `Bladder` index
sce.bladder <- readRDS(url(tmfacs[1]))
scfind.index <-  buildCellTypeIndex(sce = sce.bladder, 
                             cell.type.label = "cell_type1",
                             dataset.name = "Bladder", 
                             assay.name = "counts")
```

You can use the `mergeDataset` function to combine more than one dataset into one super index. The function `saveObject` allows you to save your index for future use.

To Quick Start __scfind__ with pre-computed indexes:

```
# `scfind` index of the `Tabula Muris (FACS)` dataset
data(ExampleIndex_TabulaMurisFACS)

scfind.index.tmfacs <- loadObject(file = url(ExampleIndex_TabulaMurisFACS))

# `scfind` index of the `Tabula Muris (10X)` dataset
data(ExampleIndex_TabulaMuris10X)

scfind.index.tm10x <- loadObject(file = url(ExampleIndex_TabulaMuris10X))
```

__Q__: Where can I report bugs, comments, issues or suggestions?

__A__: Please use [this page](https://github.com/hemberg-lab/scfind/issues).

__Q__: Is __scfind__ published?  

__A__: Not yet, but a copy of __scmap__ manuscript is available on [bioRxiv](https://doi.org/10.1101/788596).

__Q__: What is __scfind__ licence?

__A__: GPL-3
