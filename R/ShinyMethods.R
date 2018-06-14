#' UI handler for the shiny front end of scfind
#'
#' @param object SCFind object
#'
#' @name ui.scfind
#' @aliases ui.scfind
#'
#' @importFrom shiny titlePanel sidebarLayout textInput mainPanel plotOutput fluidPage h1 h2 h3 h4 sidebarPanel uiOutput checkboxGroupInput plotOutput
#' 
ui.scfind <- function(object)
{
    return(
        fluidPage(
            titlePanel("Scfind index"),
            sidebarLayout(
                sidebarPanel(
                    textInput("geneList", label = h3("Gene List"), value = ""),
                    uiOutput("geneCheckbox"),
                    checkboxGroupInput("datasetCheckbox", h3("Datasets"), choices = object@datasets, selected = object@datasets)
                ),
                mainPanel(
                    plotOutput("cellTypesHisto")
                )
            )
        )
    )
}



#' Server handler for scfind
#'
#' @param input handler item for the ShinyApp
#' @param output handel item for the ShinyApp
#' @param session handler item for the stateful ShinyApp
#'
#' @name server.scfind
#' @aliases server.scfind
#'
#' @importFrom shiny renderPlot stopApp checkboxGroupInput
#' @importFrom ggplot2 ggplot geom_bar ggtitle xlab ylab aes coord_flip theme_minimal
server.scfind <- function(object)
{

    return(
        function(input,output,session)
        {
            output$cellTypesHisto <- renderPlot({
                ## Render a barplot
                print(length(input$geneCheckbox))
                print(input$geneCheckbox)
                if (length(input$geneCheckbox) != 0)
                {
                    result <- findCellTypes(object, input$geneCheckbox)
                    print(result)
                    result <- setNames(unlist(result, use.names=F),rep(names(result), lengths(result)))
                    df <- data.frame(cell_type = names(result), cell_id = result)
                    g <- ggplot(df, aes(x=cell_type)) +
                        xlab("Cell Type") +
                        ylab("Cells") +
                        geom_bar(color = "blue") +
                        ggtitle(paste0(input$geneCheckbox,collapse=",")) +
                        coord_flip() +
                        theme_minimal()
                    g
                }
                else
                {
                    barplot(table(mtcars$gear), 
                        ## main=input$region,
                        ylab="Number of Cells",
                        xlab="Cell Type")
                }
            })

            output$geneCheckbox <- renderUI({
                genes <-  unlist(strsplit(gsub("\\s", "", input$geneList), ","))
                checkboxGroupInput("geneCheckbox","Select Genes", choices = genes, selected = genes, inline = T)

            })

            session$onSessionEnded(function() {
                stopApp()
            })
        })
}


#' Opens \code{scfind} index in an interactive session in a web browser.
#'
#' Runs interactive \code{shiny} session of \code{scfind} based on the indexed project.
#'
#' @param object an object of \code{SCFind} class
#'
#' @return Opens a browser window with an interactive \code{shiny} app and visualize queries on the dataset.
#' 
#' @name scfind.interactive
#' @aliases scfind.interactive
#'
#' @importFrom shiny shinyApp
#' @importFrom graphics plot
scfind.interactive <- function(object) {

    if (is.null(object@index) || length(object@datasets) == 0) {
        warning("You should consider the option of indexing a dataset")
        return()
    }
    
    
    shinyApp(
        ui = ui.scfind(object),
        server = server.scfind(object),
        options = list(launch.browser = TRUE)
    )
}

#' @rdname scfind.interactive
#' @aliases scfind.interactive
setMethod("scfind_interactive", signature(object = "SCFind"), scfind.interactive)
