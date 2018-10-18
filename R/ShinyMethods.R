#' UI handler for the shiny front end of scfind
#'
#' @name ui.scfind
#' @aliases ui.scfind
#'
#' @importFrom shiny titlePanel sidebarLayout textInput mainPanel plotOutput fluidPage h1 h2 h3 h4 sidebarPanel uiOutput checkboxGroupInput plotOutput
#' @importFrom DT dataTableOutput
ui.scfind <- function()
{
    fluidPage(
        titlePanel("Scfind index"),
        sidebarLayout(
            sidebarPanel(
                textInput("geneList", label = h3("Gene List"), value = ""),
                uiOutput("geneCheckbox"),
                uiOutput("datasetCheckbox"),
                plotOutput("geneSupportHisto", height = 400),
                dataTableOutput("queryOptimizer")
            ),
            mainPanel(
                dataTableOutput("cellTypesData")
            )
        )
    )
}


#' Server handler for scfind
#'
#' @param object An SCFind object that the shiny app will be deployed upon
#'
#' @name server.scfind
#' @aliases server.scfind
#'
#' @importFrom shiny renderPlot stopApp checkboxGroupInput observeEvent reactiveVal
#' @importFrom DT renderDataTable datatable
#' @importFrom data.table as.data.table data.table
#' @importFrom ggplot2 ggplot geom_bar geom_col ggtitle xlab ylab aes coord_flip theme_minimal
server.scfind <- function(object)
{

    return(
        function(input, output, session)
        {

            last.query.state <- reactiveVal("genelist")
            initial.datasets <- "initial"
            gene.list <- reactiveVal(c())

            
            observeEvent(
                input$geneCheckbox,
                {
                    last.query.state("checkbox")
                })
            
            observeEvent(
                input$queryOptimizer_rows_selected,
                {   
                    last.query.state("query_optimizer")
                })
           
            observeEvent(input$geneList,{
                text <- gsub("\\s", "", input$geneList)
                gene.list.input <- unlist(strsplit(text, ","))
                last.query.state("genelist")
                print(paste("GeneList",gene.list.input))
                gene.list(gene.list.input)
            })

            recommended.queries <- reactive({
                selected.genes <- gene.list()
                selected.datasets <- input$datasetCheckbox
                if (length(selected.genes) != 0)
                {
                    ## print(paste("QO gene:",selected.genes))
                    ## print(paste("QO selected:",selected.datasets))
                    available.queries <-  markerGenes(object, selected.genes, selected.datasets)
                }
                else
                {
                    available.queries <- c()
                }
                available.queries <- as.data.table(available.queries)
                available.queries
            })

            qo.output <- reactive({
                selected.index <- input$queryOptimizer_rows_selected
                available.queries <-  recommended.queries()
                selected.query <- available.queries[selected.index, 'Query']
                ## print(paste0('selected query', selected.query))
                unlist(strsplit(gsub("\\s", "", selected.query), ","))
            })

            
            
            output$geneCheckbox <- renderUI({
                if(last.query.state() == "query_optimizer")
                {
                    checkboxGroupInput("geneCheckbox", h4("Select Genes"), choices = gene.list(), selected = qo.output(), inline = T)
                }
                else if (last.query.state() == "checkbox")
                {
                    checkboxGroupInput("geneCheckbox", h4("Select Genes"), choices = gene.list(), selected = input$geneCheckbox, inline = T)                    
                }
                else
                {
                    checkboxGroupInput("geneCheckbox", h4("Select Genes"), choices = gene.list(), selected = gene.list(), inline = T)
                }
                
            })


            
            output$datasetCheckbox <-  renderUI({
                datasets <- object@datasets
                box.selection <-  object@datasets
                if (initial.datasets == "initial")
                {
                    ## set the flag down
                    initial.datasets <- "initialized"
                }
                else
                {
                    box.selection <-  input@datasetCheckbox
                }
                checkboxGroupInput("datasetCheckbox", h3(paste("Datasets", length(box.selection), "/", length(datasets))), choices = datasets, selected = box.selection, inline = T)
            })
           
            

            
            
            output$queryOptimizer <- renderDataTable({
                datatable(recommended.queries(), selection = 'single')
            })

            
            cell.types <- reactive({
                selection <- input$geneCheckbox
                
                if (length(selection) != 0){
                    df <- query.result.as.dataframe(findCellTypes(object, selection, input$datasetCheckbox))
                    df
                }
                else
                {
                    data.frame(cell_type = c(), cell_id = c())
                }
            })

            gene.support <- reactive({
                gene.selection <- gene.list()
                dataset.selection <- input$datasetCheckbox
                gene.support <- as.data.frame(object@index$genesSupport(gene.selection, dataset.selection))
                dimnames(gene.support)[[2]] <- 'support'
                gene.support$genes <- rownames(gene.support)
                gene.support
            })
            
            
            output$cellTypesData <- renderDataTable({       
                df <- cell.types()
                datatable(phyper.test(object, df, input$datasetCheckbox), selection = 'single')
            })
            
            output$geneSupportHisto <- renderPlot({
                ## Render a barplot
                ## print(length(input$geneCheckbox))
                ## print(input$geneCheckbox)
                df <- gene.support()
                ## print(df)
                if (nrow(df) != 0)
                {
                    g <- ggplot(df, aes(x=genes, y= support)) +
                        xlab("Gene") +
                        ylab("Cells") +
                        geom_col(color = "blue") +
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

#' @rdname scfindShinyServer
#' @aliases scfindShinyServer
setMethod("scfindShinyServer", signature(object = "SCFind"), server.scfind)



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
        ui = ui.scfind(),
        server = server.scfind(object),
        options = list(launch.browser = FALSE)
    )
}

#' @rdname scfind.interactive
#' @aliases scfind.interactive
setMethod("scfindShiny", signature(object = "SCFind"), scfind.interactive)




scfind.get.genes.in.db <- function(object){
    
    return(object@index$genes())

}


#' @rdname scfind.get.genes.in.db
#' @aliases scfind.get.genes.in.db
setMethod("scfindGenes", signature(object = "SCFind"), scfind.get.genes.in.db)
