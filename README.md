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

# For Windows users:
# Please install the latest version of Rtools at https://cran.r-project.org/bin/windows/Rtools/ prior to installation of scfind
```

__Q__: Where can I report bugs, comments, issues or suggestions?
__A__: Please use [this page](https://github.com/hemberg-lab/scfind/issues).

__Q__: Is __scfind__ published?  
__A__: Not yet

__Q__: What is __scfind__ licence?  
__A__: GPL-3
