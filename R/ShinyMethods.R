#' Opens \code{scfind} index in an interactive session in a web browser.
#'
#' Runs interactive \code{shiny} session of \code{scfind} based on the indexed project.
#'
#' @name scfindShiny
#' 
#' @param object an object of \code{SCFind} class
#' @param dictionary a word2vec model
#'
#' @return Opens a browser window with an interactive \code{shiny} app and visualize queries on the dataset.
#' 
#'
#' @importFrom shiny shinyApp
#' @importFrom graphics plot
scfind.interactive <- function(object, dictionary) {
    
    if (is.null(object@index) || length(object@datasets) == 0) {
        warning("You should consider the option of indexing a dataset")
        return()
    }
    
    if(!missing(dictionary) && length(names(dictionary)) == 0){
        warning("Please download a valid dictionary at `https://github.com/hemberg-lab/scfind`")
        return()
    }
    
    shinyApp(
        ui = ui.scfind(),
        server = if(missing(dictionary) || all(startsWith(tolower(object@index$genes()), "chr"))) server.scfind(object) else server.scfind.w2v(object, dictionary)
    )
}

#' @rdname scfindShiny
#' @aliases scfindShiny
setMethod("scfindShiny", signature(object = "SCFind"), scfind.interactive)
