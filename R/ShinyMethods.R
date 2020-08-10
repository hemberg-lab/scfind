#' Opens \code{scfind} index in an interactive session in a web browser.
#'
#' Runs interactive \code{shiny} session of \code{scfind} based on the indexed project.
#'
#' @name scfindShiny
#' 
#' @param object an object of \code{SCFind} class
#'
#' @return Opens a browser window with an interactive \code{shiny} app and visualize queries on the dataset.
#' 
#'
#' @importFrom shiny shinyApp
#' @importFrom graphics plot
scfind.interactive <- function(object) {
    
    if (is.null(object@index) || length(object@datasets) == 0) {
        warning("You should consider the option of indexing a dataset")
        return()
    }
    
    shinyApp(
        ui = ui.scfind(),
        server = server.scfind(object)
    )
}

#' @rdname scfindShiny
#' @aliases scfindShiny
setMethod("scfindShiny", signature(object = "SCFind"), scfind.interactive)
