% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/CoreMethods.R
\docType{methods}
\name{markerGenes}
\alias{markerGenes}
\alias{find.marker.genes}
\alias{markerGenes,SCFind,character-method}
\title{Performs query optimization and return the best candidate gene sets}
\usage{
markerGenes(object, gene.list, datasets, exhaustive = FALSE, support.cutoff = -1)

find.marker.genes(object, gene.list, datasets, exhaustive, support.cutoff)

\S4method{markerGenes}{SCFind,character}(object, gene.list, datasets, exhaustive = FALSE, support.cutoff = -1)
}
\arguments{
\item{object}{SCFind object}

\item{gene.list}{A list of nGenes existing in the database}

\item{datasets}{the datasets of the objects to be considered}

\item{exhaustive}{use bootstrap method instead of the default frequent pattern growth}

\item{support.cutoff}{define result cutoff}
}
\value{
hierarchical list of queries and their respective scores
}
\description{
This function can be used with quite long gene lists
that otherwise would have no cell hits in the database
}
