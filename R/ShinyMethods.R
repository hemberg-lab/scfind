#' UI handler for the shiny front end of scfind
#'
#' @param object SCFind object
#'
#' @name ui.scfind
#' @aliases ui.scfind
#'
#' @importFrom shiny titlePanel sidebarLayout textInput mainPanel plotOutput fluidPage h1 h2 h3 h4 sidebarPanel uiOutput checkboxGroupInput plotOutput
#' @importFrom DT dataTableOutput
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
                    checkboxGroupInput("datasetCheckbox",
                                       h3("Datasets"),
                                       choices = object@datasets,
                                       selected = object@datasets,
                                       inline = T),
                     dataTableOutput("queryOptimizer")
                ),
                   
                mainPanel(
                    plotOutput("cellTypesHisto", height = 800),
                    dataTableOutput("cellTypesData")
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
#' @importFrom DT renderDataTable datatable
#' @importFrom data.table as.data.table data.table
#' @importFrom ggplot2 ggplot geom_bar ggtitle xlab ylab aes coord_flip theme_minimal
server.scfind <- function(object)
{

    return(
        function(input,output,session)
        {
            

            gene.list <- reactive({
                genes <-  unlist(strsplit(gsub("\\s", "", input$geneList), ","))
                genes
            })

            checkbox.selection <- reactive({
                selected.index <- input$queryOptimizer_rows_selected
                if (!is.null(selected.index))
                {
                    available.queries <-  recommended.queries()
                    selected.query <- available.queries[selected.index, 'query']
                    genes <-  unlist(strsplit(gsub("\\s", "", selected.query), ","))
                }
                else
                {
                    if (is.null(input$geneCheckbox))
                    {
                        genes <- gene.list()
                    }
                    else
                    {
                        genes <- input$geneCheckbox
                    }
                }
                genes
                ## updateSelectInput(session, "geneCheckbox", selected  = genes)
            })

            
            output$geneCheckbox <- renderUI({
                
                ## Select genes
                checkboxGroupInput("geneCheckbox", h4("Select Genes"), choices = gene.list(), selected = checkbox.selection(), inline = T)
            })
           
            
            
            
            recommended.queries <- reactive({
                
                selected.genes <- gene.list()
                selected.datasets <- input$datasetCheckbox
                if (length(selected.genes) != 0)
                {
                    available.queries <-  markerGenes(object, selected.genes, selected.datasets)
                }
                else
                {
                    available.queries <- c()
                }
                available.queries <- as.data.table(available.queries)
                available.queries
            })


            output$queryOptimizer <- renderDataTable({
                datatable(recommended.queries(), selection = 'single')
            })

            
            
            cell.types <- reactive({
                selection <- input$geneCheckbox
                if (length(input$geneCheckbox) != 0){
                    print(selection)
                    result <- findCellTypes(object, selection, input$datasetCheckbox)
                    ## print(result)
                    result <- setNames(unlist(result, use.names=F), rep(names(result), lengths(result)))
                    df <- data.frame(cell_type = names(result), cell_id = result)
                    df
                }
                else
                {
                    data.frame(cell_type = c(), cell_id = c())
                }
            })
            
            output$cellTypesData <- renderDataTable({
                df <- cell.types()
                cell.types.df <- aggregate(cell_id ~ cell_type, df, FUN = length)
                datatable(cell.types.df, selection = 'single')
                
            })
            
            output$cellTypesHisto <- renderPlot({
                ## Render a barplot
                ## print(length(input$geneCheckbox))
                ## print(input$geneCheckbox)
                df <- cell.types()
                if (nrow(df) != 0)
                {
                    g <- ggplot(df, aes(x=cell_type)) +
                        xlab("Cell Type") +
                        ylab("Cells") +
                        geom_bar(color = "blue") +
                        ggtitle(paste0(input$geneCheckbox, collapse = ",")) +
                        coord_flip() +
                        theme_minimal()
                }
                else
                {
                    g <- plot(0,type='n',axes=FALSE,ann=FALSE)
                }
                g
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
