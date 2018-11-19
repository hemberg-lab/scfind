#' UI handler for the shiny front end of scfind
#'
#' @name ui.scfind
#' @aliases ui.scfind
#'
#' @importFrom shiny  sidebarLayout actionButton textInput navbarPage navbarMenu plotOutput fluidPage fluidRow column h1 h2 h3 h4 a tags$h4 tags$head tags$style HTML sidebarPanel uiOutput checkboxGroupInput actionLink plotOutput verbatimTextOutput
#' @importFrom DT dataTableOutput
ui.scfind <- function()
{
    fluidPage(
          title = "SCfind",
          
          
          shiny::tags$head(
            shiny::tags$style(shiny::HTML("
                            body {
          
                            background-image: url('https://scfind.sanger.ac.uk/img/scfind.png');
                            background-size: 200px;
                            background-attachment: fixed;
                            background-repeat: no-repeat;
                            background-position: center;
                            }
                            
                            #homeBtn {
                            position: fixed;
                            top: 20px;
                            right: 20px;
                            height: 80px;
                            width: 80px;
                            font-size: 14px;
                            z-index: 3;
                            }

                            .row {
                            width: 100%;
                            }
                            
                            .navbar-brand {
                            font-weight: bold;
                            }
                            
                            #search {
                            
                            position: fixed;
                            top: 8%;
                            padding-left: 20%;
                            padding-right:20%;
                            background-color: rgba(256,256,256, 0.6);
                            z-index: 2;
                            }
                            
                            #search h3 {
                            color: rgb(0, 180, 204);
                            font-size: 20px;
                            }
                            
                            #geneList {
                            height: 35px;
                            border: 2px solid #00B4CC;
                            border-radius: 5px;
                            padding: 5px;
                            color: #9DBFAF;
                            outline: none;
                            }
                            
                            #geneList:focus {
                            color: #00B4CC;
                            }
                            
                            #main {
                            position: absolute;
                            top: 20%;
                            }
                            
                            
                            #main h3 {
                            margin-top: 0;
                            text-align: center;
                            color: rgb(0, 180, 204);
                            }
                            
                            #geneCheckbox{
                            position: relative;
                            z-index: 1;
                            width: 300px;
                            }
                            
                            #geneCheckbox .shiny-options-group {
                            height: 440px;
                            border:0;
                            padding:0;
                            margin:0;
                            display: flex;
                            flex-direction: column;
                            justify-content: space-between;
                            }
                            
                            h4 {
                            font-size: 20px;
                            color: #455254;
                            
                            }
                            
                            .checkbox {
                            width: 30%;
                            padding-left: 20px;
                            padding-top: 5px;
                            padding-bottom: 5px;
                            background-color: white;
                            
                            }
                            .checkbox span {
                            font-size: 14px;
                            font-weight: bold;
                            
                            }
                            
                            #geneSupportHisto {
                            border-left: 80px sold white;
                            padding-top: 30px;
                            background-color: transparent;
                            }
                            
                            .datasetsBox {
                            position: fixed; 
                            bottom: -200px; 
                            left: 30%;
                            right:30%;
                            height: 250px;
                            width: auto;
                            border: 2px solid #00B4CC;
                            border-radius: 5px 5px 0 0; 
                            padding: 15px; 
                            font-size: 16px; 
                            font-weight: bolder;
                            color: #455254; 
                            background-color: rgba(0, 180, 204, 0.3); 
                            transition: 0.8s; 
                            overflow: scroll;
                            z-index: 3;
                            }
                            
                            
                            .datasetsBox:hover{
                            bottom: 0;
                            }

                            .titleDataset{
                            position: fixed;
                            z-index:4;
                            }

                            #datasetCheckbox{
                            margin-top: 40px;
                            text-align: right;
                            }
                            
                            #selectGenes {
                            position: absolute;
                            }
                            
                            #query {
                            position: absolute;
                            right: 35%;
                            width: 37%; 
                            
                            }
                            
                            #celltype {
                            position: absolute;
                            right: 0;
                            width: 35%;
                            }
                            
                            table {
                            table-layout: fixed;
                            }
                            
                            table td{
                            word-wrap: break-word;
                            }
                            
                            #ctData {
                            position: fixed;
                            bottom: 5px;
                            left: -50px;
                            height: 20%;
                            border: 2px solid #00B4CC;
                            padding-left: 50px;
                            border-radius: 5px 5px 0 0; 
                            text-align: left;
                            color: #455254; 
                            background-color: rgba(0, 180, 204, 0.3); 
                            overflow: scroll;
                            z-index: 3;
                            }

                            
                            "))
            ),
                      shiny::actionButton(inputId = 'homeBtn', label = "Home", shiny::icon('th'), onclick = "window.open('https://scfind.sanger.ac.uk', '_self')"),
                      shiny::fluidRow(id = "search",
                              shiny::column(12, align="center",
                                            
                                             textInput("geneList", label = h3(paste("What's your gene list today?")), value = "")
                              )
                      ),
                    shiny::fluidRow(id = "main",
                              shiny::column(3,
                                        shiny::fluidRow(id = "selectGenes",
                                                  shiny::column(1,
                                                            uiOutput("geneCheckbox")
                                                  ),
                                                  shiny::column(2,
                                                            plotOutput("geneSupportHisto", width = 350, height = 500),
                                                            shiny::verbatimTextOutput("ctData")
                                                  )
                                        )
                              ),
                              shiny::column(4, id = "query",
                                        shiny::tags$div(class = "datasetsBox",
                                            shiny::tags$span(class="titleDataset",                                            
                                            shiny::tags$h3(uiOutput("selectedDataset")),
                                            shiny::actionLink("selectall","Select/Deselect All")),
                                            uiOutput("datasetCheckbox")),
                                        shiny::tags$h4(uiOutput("suggestHyper")),
                                        dataTableOutput("queryOptimizer")
                              ),
                              shiny::column(5, id = "celltype",
                                        shiny::tags$h4(uiOutput("selectedQuery")),
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
#' @importFrom shiny renderPlot stopApp checkboxGroupInput observeEvent observe reactiveVal renderText renderPrint
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
                if (length(selected.genes) > 1 && length(selected.datasets) != 0)
                {
                    #print(paste("QO gene:",selected.genes))
                    #print(paste("QO selected:",selected.datasets))
                    available.queries <-  markerGenes(object, selected.genes, selected.datasets)
                }
                else
                {
                    available.queries <- c()
                }
                available.queries <- as.data.table(available.queries)
                if(length(available.queries != 0)){
                    available.queries$tfidf <- signif(available.queries$tfidf, digits = 4) 
                    available.queries
                } else {
                    available.queries  
                }
            })
            

            qo.output <- reactive({
                selected.index <- input$queryOptimizer_rows_selected
                available.queries <-  recommended.queries()
                selected.query <- available.queries[selected.index, 'Query']
                unlist(strsplit(gsub("\\s", "", selected.query), ","))
            })

            output$geneCheckbox <- renderUI({
                genes.in.list <- as.character(gene.list())
                ## else fallback
                selection <- gene.list()
                
                if(last.query.state() == "query_optimizer")
                {
                    selection <- qo.output()
                }
                else if (last.query.state() == "checkbox")
                {
                    selection <-  input$geneCheckbox
                }
                if(length(gene.list()) != 0 && length(input$datasetCheckbox) != 0){
                          checkboxGroupInput("geneCheckbox", h4(paste("Select Genes", length(selection), "/", length(genes.in.list))), choices = genes.in.list , selected = selection, inline = F)
                } else {
                          checkboxGroupInput("geneCheckbox", label='', choices = c() , selected = c(), inline = F)
                }
            })


            output$selectedDataset <- renderText({
                datasetName <- gsub("/", "", session$clientData$url_pathname)
                datasets <- object@datasets
                box.selection <-  input$datasetCheckbox
                
                if(length(datasetName) != 0 && datasetName != '#'){
                    datasetName <- gsub("mca", "Mouse Cell Atlas ", datasetName)
                    datasetName <- gsub("tm-10X", "Tabula Muris - 10X ", datasetName)
                    datasetName <- gsub("tm-facs", "Tabula Muris - FACS ", datasetName)
                    datasetName <- gsub("brain", "Mouse Brain Atlas ", datasetName)
                    datasetName <- gsub("malaria", "Malaria Cell Atlas ", datasetName)
                    datasetName <- gsub("liver", "Human Liver Atlas ", datasetName)
                    datasetName <- gsub("atacseq", "Mouse sci-ATAC-seq Atlas ", datasetName)
                } else {
                    datasetName <- ''
                }
                if(length(box.selection) == 0){
                    paste("Please select at least 1 dataset below:")
                }
                else if(length(box.selection) < 2){
                    paste(datasetName, length(box.selection), "/", length(datasets), "dataset")
                } else {
                    paste(datasetName, length(box.selection), "/", length(datasets), "datasets")
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

                checkboxGroupInput("datasetCheckbox", label = '', choices = datasets, selected = box.selection, inline = T)
                     
                
            })
            
            observe({
                exampleCT <- sample(cellTypeNames(object), 1)
                exampleQuery <- cellTypeMarkers(object, exampleCT)$genes
                updateTextInput(session, "geneList", value = "", placeholder = exampleQuery)
            })
            
            observe({
                datasets <- object@datasets
                box.selection <- object@datasets
                
                if(input$selectall == 0) return(NULL) 
                else if (input$selectall%%2 == 0)
                {
                    updateCheckboxGroupInput(session, "datasetCheckbox", label = '', choices = datasets, selected = box.selection, inline = T)
                }
                else
                {
                    updateCheckboxGroupInput(session, "datasetCheckbox", label = '', choices = datasets, inline = T)
                    
                }
            })

            
            output$queryOptimizer <- renderDataTable({
                datatable(recommended.queries(), selection = 'single', options = list(autoWidth = TRUE))
            })

            cell.types <- reactive({

                selection <- input$geneCheckbox

                if (length(selection) != 0 && length(input$datasetCheckbox) != 0){
                    df <- query.result.as.dataframe(findCellTypes(object, selection, input$datasetCheckbox))
                    #print("yo")
                    #print(df)
                }
                else
                {
                    df <- data.frame(cell_type = c(), cell_id = c())
                }
                df
            })

            gene.support <- reactive({
                      gene.selection <- gene.list()
                      dataset.selection <- input$datasetCheckbox
                      if(length(dataset.selection) != 0){
                          gene.support <- as.data.frame(object@index$genesSupport(gene.selection, dataset.selection))
                          dimnames(gene.support)[[2]] <- 'support'
                          gene.support$genes <- rownames(gene.support)
                      } else {
                          gene.support <- data.frame(support = c(), genes = c())
                          print('Initiating datasets or dataset not defined')
                      }
                      gene.support
            })
            
            
            output$cellTypesData <- renderDataTable({
                
                selection <- input$geneCheckbox
                df <- cell.types()
                df
                if (nrow(df) != 0 && length(input$datasetCheckbox) != 0)
                {
                    rdt <- phyper.test(object, df, input$datasetCheckbox)
                    rdt <- data.frame(rdt)
                }
                else
                {
                  rdt <- data.table()
                }

            }, server = TRUE)
            
            output$ctData = renderPrint({
                      s = input$cellTypesData_rows_selected
                      selection <- input$geneCheckbox
                      df <- cell.types()
                      df
                      if (nrow(df) != 0 && length(input$datasetCheckbox) != 0)
                      {
                                rdt <- phyper.test(object, df, input$datasetCheckbox)
                                mge <- evaluateMarkers(object, selection, rdt$cell_type)
                                mge <- mge[order(mge$genes),]
                                mge <- mge[order(mge[,1]),]

                                mge$precision <- round((mge$precision), digits = 4)
                                mge$recall <- round((mge$recall), digits = 4)
                                mge$f1 <- round((mge$f1), digits = 4)
                                mge <- aggregate(.~cellType, mge[,c("cellType", "genes", "precision", "recall", "f1")], FUN=toString)

                                
                                for(i in 1: nrow(mge)){
                                          for(j in 1: length(selection)){
                                                    mge$genes[i]<- gsub(j, selection[j], mge$genes[i])
                                          }
                                }
                                
                                ct <- c()
                                for(i in 1: nrow(mge)){
                                          
                                          cellName <- as.character(rdt$cell_type[i])
                                          cellInfo <- paste(cellName, 
                                                            "\n| genes: ", 
                                                            mge$genes[i],
                                                            "\n| precision:", 
                                                            mge$precision[i], 
                                                            "\n| recall:", 
                                                            mge$recall[i],
                                                            "\n| f1:", 
                                                            mge$f1[i],"\n\n")
                                          ct <- c(ct, cellInfo)
                                          
                                          
                                }
                                if(length(s)!=0) {
                                          cat(ct[s])
                                          
                                } else {
                                          print("Select cell types to view marker evaluation")
                                }
                      } else {
                                ct <- c()
                      }
            })
            
            
            
            output$selectedQuery <- renderText({
                      
                      if(nrow(cell.types()) != 0){
                                selectedList <- paste(input$geneCheckbox, collapse=", ")
                              paste("Scfind found ", length(unique(cell.types()$cell_type)), " cell types for: ", selectedList)
                      } else if(nrow(cell.types()) == 0 && length(gene.list()) != 0 && length(input$datasetCheckbox) != 0) {
                                paste("Expecting specific cell types?", "Select our recommended query!", sep="<br>")
            } else {
                                ""
                      }

            })
            
            output$suggestHyper <- renderText({
                if(length(input$datasetCheckbox) == 0 && length(gene.list()) != 0){
                    "Oops! dataset is missing."
                }
                else if(length(gene.list()) != 0){
                    "Recommended query for the initial genesets: "
                } else {
                    ""
                }
            })
          
            output$geneSupportHisto <- renderPlot({
                      df <- gene.support()
                      if (nrow(df) != 0)
                      {
                                g <- ggplot(df, aes(x=genes, y= support)) +
                                          xlab("") + 
                                          ylab("Cells") +
                                          geom_col(color = "blue") +
                                          coord_flip() +
                                          ggplot2::theme() +
                                          ggplot2::scale_x_discrete(limits = rev(df$genes))
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
