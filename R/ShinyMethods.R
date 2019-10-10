#' UI handler for the shiny front end of scfind
#'
#' @name ui.scfind
#' @aliases ui.scfind
#'
#' @importFrom shiny  sidebarLayout actionButton textInput navbarPage navbarMenu plotOutput fluidPage fluidRow column h1 h2 h3 h4 a  HTML sidebarPanel uiOutput checkboxGroupInput actionLink plotOutput verbatimTextOutput
#' @importFrom DT dataTableOutput
#' @export
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
                            
                            .shiny-notification {
                            position:fixed;
                            height: 60px;
                            width: 400px;
                            top: 12%;
                            left: calc(50% - 200px);;
                            }

                            .progress-text {
                            text-align: center;
                            font-size: 18px;
                            font-style: italic; 
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
                            background-color: rgba(256,256,256, 0.6);
                            z-index: 2;
                            }
                            
                            #search h3 {
                            color: rgb(0, 180, 204);
                            font-size: 20px
                            }
                            
                            #geneList {
                            left: calc(50% - 150px);;
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
                            
                            #geneCheckbox {
                            position: relative;
                            width: 300px;
                            z-index: 0;
                            }

                            #geneCheckbox:hover {
                            z-index: 1;
                            }
                            
                            #geneCheckbox .shiny-options-group {
                            left: 0;
                            border: 0;
                            padding: 0;
                            margin: 0;
                            display: flex;
                            flex-direction: column;
                            justify-content: space-between;
                            opacity: 0.6;
                            }

                            #geneCheckbox .shiny-options-group:hover {
                            opacity: 1;
                            }
                            
                            #geneCheckbox .checkbox {
                            width: 150px;
                            height: 50px;
                            overflow-wrap: break-word;
                            }

                            h4 {
                            font-size: 18px;
                            color: #455254;
                            background-color: rgba(256,256,256, 0.6);
                            }
                            
                            .checkbox {
                            width: 30%;
                            padding-left: 20px;
                            padding-top: 5px;
                            padding-bottom: 5px;
                            background-color: rgba(255,255,255,0.8);
                            
                            }
                            .checkbox span {
                            font-size: 14px;
                            font-weight: bold;
                            
                            }
                            
                            #lociRegions {
                            width: 500px;
                            }

                            .ucsc {
                            color: red;
                            }
                            
                            .border-l {
                            border-left: 1px solid;
                            }

                            .border-lb {
                            border-left: 1px solid;
                            border-bottom: 1px solid;
                            }                            

                            #geneSupportHisto {
                            padding-top: 30px;
                            margin-left: 80px;
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

                            .datasetsBox:hover {
                            bottom: 0;
                            }

                            .titleDataset {
                            position: fixed;
                            z-index:4;
                            }

                            #datasetCheckbox {
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
                            
                            table td {
                            word-wrap: break-word;
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
                                            uiOutput("lociRegions"),
                                        shiny::fluidRow(id = "selectGenes",
                                                  shiny::column(1,
                                                            uiOutput("geneCheckbox")
                                                  ),
                                                  shiny::column(2,
                                                            uiOutput("geneHisto")
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
                                        dataTableOutput("queryOptimizer"),
                                        shiny::tags$h4(uiOutput("selectedQuery")),
                                        dataTableOutput("cellTypesData")
                                        
                              ),
                              shiny::column(5, id = "celltype",
                                            shiny::tags$h4(uiOutput("ctData")),
                                            uiOutput("evaluateSum")
                              )
                    )
          )
}


#' Server handler for scfind
#'
#' @param object An SCFind object that the shiny app will be deployed upon
#'
#' @name scfindShinyServer
#' @aliases scfindShinyServer
#'
#' @importFrom shiny reactive updateCheckboxGroupInput renderPlot stopApp checkboxGroupInput observeEvent observe reactiveVal renderText renderPrint renderUI Progress updateTextInput updateSelectInput
#' @importFrom DT renderDataTable datatable
#' @importFrom data.table as.data.table data.table
#' @importFrom grDevices rainbow
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
                
# This is necessary, because users are allowed to input chromosome as the format chrX:XXXXX-YYYYY instead of chrX_XXXXX_YYYYY
                        text <- gsub("\\s|,", ",", input$geneList)
                        gene.list.input <- unlist(strsplit(text, ","))
                        chr.list.input <- c()

                            if(all(startsWith(object@index$genes(), "chr") == T)){
                                gene.names.all <- object@index$genes()

                                p.pos <- list(p.chr = c(), p.start = 0, p.end = 0)
                                p.pos$p.chr <- gsub("_.*$", "", gene.names.all)
                                p.pos$p.start <- as.numeric(gsub(".*_(.*)\\_.*", "\\1", gene.names.all))
                                p.pos$p.end <- as.numeric(gsub(".*_", "", gene.names.all))
                                p.pos <- data.frame(p.pos)

                                for(i in 1: length(gene.list.input)){
                                    chr.list.input <- c(chr.list.input, gene.names.all[which(p.pos$p.chr == tolower(gsub(":.*$", "", gene.list.input[i])) & p.pos$p.start >= as.numeric(gsub(".*:(.*)\\-.*", "\\1", gene.list.input[i])) & p.pos$p.end <= as.numeric(gsub(".*-", "", gene.list.input[i])))])
                                }
                            }

                            gene.list.input <- if(length(chr.list.input) != 0) gsub(':|-','_',chr.list.input) else gene.list.input
                            last.query.state("genelist")
                            # print(paste("GeneList",gene.list.input))
                            gene.list.input <- caseCorrect(object, gene.list.input) 
                            gene.list(gene.list.input)
                            
            })

            recommended.queries <- reactive({
                
                selected.genes <- gene.list()
                selected.datasets <- input$datasetCheckbox
                
                if (length(selected.genes) > 1 && length(grep(TRUE, (caseCorrect(object, selected.genes) %in% object@index$genes()))) != 0 && length(selected.datasets) != 0)
                { 
                    #print(paste("QO gene:",selected.genes))
                    #print(paste("QO selected:",selected.datasets))
                    available.queries <-  as.data.table(markerGenes(object, selected.genes, selected.datasets))
                }
                else
                {
                    available.queries <- data.table()
                }
                
                if(length(available.queries) != 0){
                    available.queries$tfidf <- signif(available.queries$tfidf, digits = 5) 
                    available.queries
                } else {
                    available.queries  
                }
                
            })
            

            qo.output <- reactive({
                selected.index <- input$queryOptimizer_rows_selected
                available.queries <-  recommended.queries()
                
                selected.query <- if(all(startsWith(object@index$genes(), "chr") == T)) gsub(":|-", "_", available.queries$Query[selected.index]) else available.queries$Query[selected.index]
                unlist(strsplit(gsub("\\s", "", selected.query), ","))
            })

            output$geneCheckbox <- renderUI({
                genes.in.list <- caseCorrect(object, as.character(gene.list()))
                ## else fallback

                selection <- caseCorrect(object, gene.list())
                
                list.indx <- intersect(grep('chr', genes.in.list), grep('_', genes.in.list))
                if(length(list.indx) != 0) {
                    genes.in.list <- gsub('_', '-', sub('_', ':', genes.in.list[list.indx]))
                    selection <- gsub('_', '-', sub('_', ':', selection[list.indx]))
                }
                
                if(last.query.state() == "query_optimizer")
                {
                    selection <- qo.output()
                }
                else if (last.query.state() == "checkbox")
                {
                    selection <-  input$geneCheckbox
                }
                
                if(length(list.indx) != 0) selection <- gsub('_', '-', sub('_', ':', selection[list.indx]))
                
                if(length(gene.list()) != 0  && sum(gene.support()$support) != 0 && length(input$datasetCheckbox) != 0) {
                          checkboxGroupInput("geneCheckbox", h4(paste("Select", if(length(list.indx) != 0) "Peaks" else "Genes", length(selection), "/", length(genes.in.list))), choices = genes.in.list , selected = selection, inline = F)
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
                    datasetName <- gsub("spinalcord", "Mouse Spinal Cord Atlas ", datasetName)
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
                progress <- Progress$new(session, min=0)
                on.exit(progress$close())
                progress$set(message = 'Initialising index...',
                             detail = 'This may take a few seconds')
                
                exampleCT <- sample(cellTypeNames(object), 1)
                exampleQuery <- cellTypeMarkers(object, exampleCT)$genes
                
                if(all(startsWith(object@index$genes(), "chr") == T)) exampleQuery <- gsub("_", "-", sub("_", ":", exampleQuery))
                
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
            
             output$evaluateSum <- renderUI ({
                 s = input$cellTypesData_rows_selected # return row number, NULL
                 
                 if(!is.null(s) && length(object@metadata) != 0){
                     shiny::tabsetPanel(type = "tabs",
                                        shiny::tabPanel("UMAP", 
                                                        shiny::selectInput('subdataset', 'Datasets', selected = NULL, choices = NULL, selectize=TRUE),
                                                        plotOutput("cellUMAP", width = 400, height = 500)),
                                        shiny::tabPanel("Evaluate Markers",  
                                                        dataTableOutput("evaluateCtMarkers"))
                                        
                     )
                 } else {
                     ""
                 }
             })

            
            observe({
                s = input$cellTypesData_rows_selected # return row number, NULL
                df <- cell.types()
                if(!is.null(s) && length(object@metadata) != 0 && nrow(df) != 0){ #!#
                    
                    rdt <- phyper.test(object, df, input$datasetCheckbox)
                    selectedCellTypes <- as.character(rdt$cell_type[s])
                    subdataset <- unique(sub("\\..*", "", selectedCellTypes)) # get datasetName. from datasetName.cellType
                    # print(subdataset)
                    if(length(subdataset) > 1){
                        updateSelectInput(session, 'subdataset', label = 'Datasets', choices = subdataset,  selected = subdataset[length(subdataset)])
                    } else {
                        updateSelectInput(session, 'subdataset', label = 'Datasets', choices = subdataset,  selected = subdataset)
                    }
                        
                }
            })
            
            output$queryOptimizer <- renderDataTable({
                if(!is.null(recommended.queries())){
                    col <- if(all(startsWith(object@index$genes(), "chr") == T)) "Peaks" else "Genes"
                    
                    datatable(recommended.queries(), selection = 'single',
                              options = list(columnDefs = list(list(width = '70px', targets = c(2, 3, 4)), list(width = '10px', targets = c(0))), pageLength = 5,
                                             autoWidth = TRUE,
                                             dom = 'Bfrtip',
                                             buttons = c('copy', 'csv', 'excel')),
                              extensions = 'Scroller', colnames = c(col , 'Sub-queries', 'TF-IDF', 'Cell No.', 'Cell types'),
                              rownames = F)
                } else {
                    data.table()
                }
            })

            cell.types <- reactive({
                selection <- if(all(startsWith(object@index$genes(), "chr") == T)) input$geneCheckbox else caseCorrect(object, input$geneCheckbox)
                
                if (length(selection) != 0 && length(input$datasetCheckbox) != 0 ){ 
                    if(all(startsWith(object@index$genes(), "chr") == T)) {
                        selection <- gsub(":|-", "_", selection)
                    }
                    df <- query.result.as.dataframe(findCellTypes.geneList(object, selection, input$datasetCheckbox))
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
                gene.selection <- if(all(startsWith(object@index$genes(), "chr") == T)) gene.list() else caseCorrect(object,gene.list())
                      dataset.selection <- input$datasetCheckbox
                      if(length(dataset.selection) != 0 && length(gene.selection) != 0){
                          gene.support <- as.data.frame(object@index$genesSupport(gene.selection, dataset.selection))
                          dimnames(gene.support)[[2]] <- 'support'
                          gene.support$genes <- rownames(gene.support)
                      } else {

                          gene.support <- data.frame()
                          gene.support$support <- c()
                          gene.support$genes <-  c()
                          print('dataset not defined')

                      }
                      gene.support
            })
            
            
            output$cellTypesData <- renderDataTable({
                df <- cell.types()
                #df
                if (nrow(df) != 0 && length(input$datasetCheckbox) != 0) 
                {
                    rdt <- phyper.test(object, df, input$datasetCheckbox)
                    rdt$pval <- as.numeric(signif(rdt$pval, digits = 6))
                    rdt <- data.frame(rdt)
                }
                else
                {
                  rdt <- data.table()
                }

            }, server = TRUE, extensions = c('Scroller', 'Buttons')
            , options = list(columnDefs = list(list(width = '60px', targets = c(1, 2)), list(width = '70px', targets = c(3))),
                deferRender = TRUE,
                scrollY = 200,
                scroller = TRUE, 
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
            ), colnames = c('Cell types', 'Cell hits', 'Cell No.', 'p-value'), rownames = F)
            
            
            output$evaluateCtMarkers <- renderDataTable({

                s = input$cellTypesData_rows_selected
                
                selection <- input$geneCheckbox
                if(all(startsWith(object@index$genes(), "chr") == T)) {
                    selection <- gsub(":|-", "_", selection)
                }
                
                df <- cell.types()

                if(!is.null(df$cell_type) && nrow(df) != 0) { #!#
                    rdt <- phyper.test(object, df, input$datasetCheckbox)
                    mge <- if(length(unique(df$cell_type)) < 2) evaluateMarkers(object, selection, as.character(unique(df$cell_type))) else evaluateMarkers(object, selection, rdt$cell_type[s])
                    if(nrow(mge) != 0){
                        mge <- mge[order(mge$genes),]
                        mge <- mge[order(mge[,1]),]
                        mge$genes <- if(all(startsWith(object@index$genes(), "chr") == T)) gsub("_", "-", sub("_", ":", mge$genes)) else mge$genes
                        mge$precision <- paste(as.character(mge$genes), ": ", round((mge$precision), digits = 3))
                        mge$recall <- paste(mge$genes, ": ", round((mge$recall), digits = 3))
                        mge$f1 <- paste(mge$genes, ": ", round((mge$f1), digits = 3))
                        mge <- aggregate(.~cellType, mge[,c("cellType", "precision", "recall", "f1")], FUN=toString)
                    }
                } else {
                    data.table()
                }

            }, extensions = 'Buttons' , options = list(columnDefs = list(list(width = '90px', targets = c(1, 2, 3))), dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')), 
            colnames = c('Cell types', 'Precision', 'Recall', 'F1 score'), 
            rownames = F)
            
            output$ctData <- renderText({
                s = input$cellTypesData_rows_selected
                selection <- toString(input$geneCheckbox)
                
                df <- cell.types()
                if (nrow(df) != 0 && length(input$datasetCheckbox) != 0)
                {
                    rdt <- phyper.test(object, df, input$datasetCheckbox)
                    
                    if(nrow(rdt) > 0){
                        if(nrow(rdt) == 1){
                            paste0("Marker evaluation of ", selection, " in ", rdt$cell_type, ": ")
                        } else {
                            if(!is.null(selection) && !is.null(s)){
                                paste0("Marker evaluation of ", selection, " in ", length(s), if(length(s) < 2) " cell type:" else " cell types:")
                            } else {
                                paste0("<span style = 'color: rgb(0, 180, 204);'>Select cell types from cell type table to start evaluation of ", selection, "</span>")
                            }
                        }
                        
                    } else {
                        ""
                    }
                }
       
            })

            
            output$selectedQuery <- renderText({
                selectedList = input$queryOptimizer_rows_selected
                df <- cell.types()
                
                if(!is.null(df$cell_type)) {
                    list.indx <- intersect(grep('chr', input$geneCheckbox), grep('_', input$geneCheckbox))
                    
                    paste("Scfind found ", length(unique(cell.types()$cell_type)), if ((length(unique(cell.types()$cell_type))) < 2) " cell type" else " cell types for the <I>selected query</I>:")
                } else {
                    if (nrow(recommended.queries()) != 0 && length(gene.list()) != 0 && length(input$datasetCheckbox) != 0 && nrow(gene.support()) != 0)  paste("<span style = 'color: rgb(0, 180, 204);'>Expecting specific cell types?", "Select above recommended query!</span>", sep="<br>") else ""
                }
            })
            
            output$suggestHyper <- renderText({
                if(length(input$datasetCheckbox) == 0 && input$geneList != ""){
                    "Oops! dataset is missing."
                } else if(input$geneList != "" && sum(gene.support()$support) == 0){
                    
                    if(all(startsWith(object@index$genes(), "chr") == T)){
                        paste("Please follow the format below for your query: ", 
                              "",
                              "chr13:56515712-56516189 , chr5:136247792-136248036", 
                              "",
                              "for search of", 
                              "'chromosome 13, Start at 56515712, End at 56516189", 
                              "and chromosome 5, Start at 136247792, End at 136248036'", 
                              sep="<br>")
                    } else {
                        paste("Please follow the format below for your query: ",
                              "",
                              paste(sample(object@index$genes(), 5), collapse = ' , '),
                              sep="<br>")
                    }

                } else if(length(gene.list()) != 0){
                    if(sum(gene.support()$support) != 0 && length(recommended.queries()) != 0){ 
                        "Recommended queries for the initial gene set: "
                    } 
                } else {
                    ""
                }
            })
            
            output$lociRegions <- renderText({
                if(all(startsWith(object@index$genes(), "chr") == T)){
                    text <- gsub("\\s|,", ",", input$geneList)
                    gene.list.input <- unlist(strsplit(text, ","))
                    gene.list.input <- unique(grep("chr",gsub("_", "-", sub("_", ":", tolower(gene.list.input))), value = T))
                    gene.list.input <- grep(":", gene.list.input, value = T)
                    gene.list.input <- gene.list.input[which(as.numeric(gsub(".*:(.*)\\-.*", "\\1", gene.list.input)) < as.numeric(gsub(".*-", "", gene.list.input)))]
                    
                    
                    if(length(gene.list.input) != 0){
                        peak.indx <- list()
                        gene.names.all <- object@index$genes()
                        
                        p.pos <- list(p.chr = c(), p.start = 0, p.end = 0)
                        p.pos$p.chr <- gsub("_.*$", "", gene.names.all)
                        p.pos$p.start <- as.numeric(gsub(".*_(.*)\\_.*", "\\1", gene.names.all))
                        p.pos$p.end <- as.numeric(gsub(".*_", "", gene.names.all))
                        p.pos <- data.frame(p.pos)
                        
                        for(i in 1: length(gene.list.input)){
                            peak.tmp <- which(p.pos$p.chr == tolower(gsub(":.*$", "", gene.list.input[i])) & p.pos$p.start >= as.numeric(gsub(".*:(.*)\\-.*", "\\1", gene.list.input[i])) & p.pos$p.end <= as.numeric(gsub(".*-", "", gene.list.input[i])))
                            peak.indx$index <- c(peak.indx$index, peak.tmp)
                            for(j in 1:length(peak.tmp)){
                                peak.indx$input <- c(peak.indx$input, gene.list.input[i])
                                peak.indx$loc_1 <- c(peak.indx$loc_1, as.numeric(gsub(".*:(.*)\\-.*", "\\1", gene.list.input[i])))
                                peak.indx$loc_4 <- c(peak.indx$loc_4, as.numeric(gsub(".*-", "", gene.list.input[i])))
                                
                            }
                        }
                        
                        peak.indx <- data.frame(peak.indx)
                        peak.indx$chromosome <- p.pos$p.chr[peak.indx$index]
                        peak.indx$loc_2 <- as.numeric(p.pos$p.start[peak.indx$index])
                        peak.indx$loc_3 <- as.numeric(p.pos$p.end[peak.indx$index])
                        peak.indx <- peak.indx[c(3,6,7,4,5,2)]
                        
                        ucsc <- "<a href='https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&position="
                        
                        chr.chart <- paste0("<br><h4>Available Peak Locations</h4><p>(Click locus to view in UCSC Genome Browser)</p><table>")
                        
                        for(i in 1: nrow(peak.indx)){
                            
                            locus <- peak.indx[i,]
                            
                            ucsc_end <- paste0("' class='ucsc' target='_blank' title='Click to view Mouse mm9: ", locus$input, " UCSC Genome Browser'>", "[")
                            ucsc_a <- paste(ucsc, locus$chromosome, "%3A", locus[2], "%2D", locus[3], ucsc_end, sep="")

                            gap <- locus[2] - locus[1]
                            
                            chromosome <- if(gap < 200) c() else paste0(locus[1], (if(gap<=400) paste(rep("&#8594;", (gap/200))) else paste("&#8594;", "...", "&#8594;")))
    
                            gap <- locus[3] - locus[2]        
                            chromosome <- paste0(chromosome, ucsc_a, locus[2], paste(rep("&#187;", (gap/100)),collapse = ""), locus[3], "]", "</a>" )
                            
                            gap <- locus[4] - locus[3]
                            
                            chromosome <- if(gap < 200) chromosome else paste0(chromosome, (if(gap<=400) paste(rep("&#8594;", (gap/200))) else paste("&#8594;", "...", "&#8594;")), locus[4]) 
                            chromosome <- paste("<tr><td>", locus$chromosome, "</td><td", if(i == nrow(peak.indx)) " class='border-lb'>" else " class='border-l'>", "...", chromosome, "...", "</td></tr>")
                            chr.chart <- paste(chr.chart, chromosome, (if(i == nrow(peak.indx))"<tr><td></td><td style = 'text-align: center; font-size: 10px;'>base (Mouse mm9)</td></tr></table>" else ""))
                        }
                        
                        chr.chart
                    }
                } else {
                    ''
                }
                
                
            })
            output$geneHisto <- renderUI({
                
                number.of.choices <- length(caseCorrect(object, as.character(gene.list())))
                if(number.of.choices < 4){
                    histoHeight = 60*4
                } else {
                    histoHeight = 250 + (number.of.choices-4) * 60
                }
                plotOutput("geneSupportHisto", width = 300, height = histoHeight)
            })
            
            
            output$geneSupportHisto <- renderPlot({
                      df <- gene.support()
                      
                      if (nrow(df) != 0 && length(grep(TRUE, (df$genes %in% object@index$genes()))) != 0)
                      {
                                list.indx <- intersect(grep('chr', df$genes), grep('_', df$genes))
                                if(length(list.indx) != 0) df$genes <- gsub('_', '-', sub('_', ':', df$genes[list.indx]))

                                max.axis <- max(df$support)
                                g <- ggplot(df, aes(x = genes, y = support), colour = support) +
                                          xlab("") + 
                                          ylab("Cells") +
                                          ggplot2::ylim(0, max.axis) +
                                          ggplot2::geom_text(aes(label=df$genes), vjust = -3, hjust = "inward", size = 4) +
                                          geom_col(fill = rainbow(length(df$genes)), position = "dodge", width = .5, alpha = .5, colour = "black") + 
                                          theme_minimal() + 
                                          coord_flip() +
                                          ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(),
                                                         panel.background = ggplot2::element_blank(), axis.line = ggplot2::element_line(colour = "black")) +
                                          ggplot2::scale_x_discrete(limits = rev(df$genes))
                      }
                      else
                      {
                                g <- plot(0,type='n',axes=FALSE,ann=FALSE)
                      }
                      g
                      
            })
            
            
            output$cellUMAP <- renderPlot({
                s = input$cellTypesData_rows_selected # return row number, NULL
                selection <- if(all(startsWith(object@index$genes(), "chr") == T)) gsub(":|-", "_", input$geneCheckbox) else caseCorrect(object, input$geneCheckbox)
                highlightCells <- c()
                hits.summary <- c()
                
                df <- cell.types()
                if(!is.null(s) && length(object@metadata) != 0 && nrow(df) != 0){
                      
                      rdt <- phyper.test(object, df, input$datasetCheckbox)
                      
                      umapChoice = input$subdataset[length(input$subdataset)]
                      getDatasetUmap <- object@metadata[[object@metadata[[1]]$umap[which(object@metadata[[1]]$dataset == umapChoice)]]]
                      getDatasetUmap <- getDatasetUmap[which(!is.na(rownames(getDatasetUmap))),]

                      #umapChoice = paste0(umapChoice, ".")
                      
                      selectedCellTypes <- as.character(rdt$cell_type[s])
                      subCellTypes <- grep(paste0(umapChoice, "."), selectedCellTypes, value = T)
                      subCellTypes <- sub(paste0(umapChoice, "."), "", subCellTypes)

                    # get coordinations of selected celltypes
                    for(j in 1: length(subCellTypes)){
                        highlightCells <- c(highlightCells, grep(paste0("^",subCellTypes[j],"$"), rownames(getDatasetUmap)))
                    }

                    # get hits
                    if (length(selection) != 0 && length(subCellTypes) != 0){ 
                        
                        true.hits <- object@index$findCellTypes(selection, umapChoice)
                        hit.ct <- paste0(umapChoice,".",subCellTypes)

                        for(i in 1: length(hit.ct)){
                            subct.id <- true.hits[[grep(paste0("^",hit.ct[i],"$"), names(true.hits))]]
                            hits <- rep(F, length(grep(paste0("^",subCellTypes[i],"$"), rownames(getDatasetUmap))))
                            hits[subct.id] <- T
                            hits.summary <- c(hits.summary, hits)
                        }
                    } 
                      
                    
                      
                    umap_plot <- data.frame(x = getDatasetUmap[,1], y = getDatasetUmap[,2], col = rownames(getDatasetUmap))
                    umap_highlight <- data.frame(x = getDatasetUmap[highlightCells,1], y = getDatasetUmap[highlightCells,2], col = rownames(getDatasetUmap)[highlightCells], shape = hits.summary)
                    umap_truehits <- data.frame(x = umap_highlight[grep(T, hits.summary),1], y = umap_highlight[grep(T, hits.summary),2], col =umap_highlight[grep(T, hits.summary),3])
                    # umap_falsehits <- data.frame(x = umap_highlight[grep(F, hits.summary),1], y = umap_highlight[grep(F, hits.summary),2], col =umap_highlight[grep(F, hits.summary),3])
                    g <- ggplot() + ggplot2::geom_point(aes(x, y, group = col), data = umap_plot, colour = ggplot2::alpha("grey", .3)) + 
                        ggplot2::geom_point(aes(x, y, colour = col), data = umap_highlight, alpha = .5)+
                        ggplot2::geom_point(aes(x, y), data = umap_truehits, shape = 21, colour = "black")+
                        ggplot2::labs(x = "UMAP1", y = "UMAP2", color = "") + 
                        ggplot2::theme_bw() + ggplot2::theme(axis.line = ggplot2::element_line(colour = "black"),
                                                                                      panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=2),
                                                                                      plot.background = ggplot2::element_blank(),
                                                                                      legend.position = "bottom",
                                                                                      axis.text.x = ggplot2::element_blank(),
                                                                                      axis.text.y = ggplot2::element_blank(),
                                                                                      axis.ticks = ggplot2::element_blank(),
                                                                                      panel.grid.major = ggplot2::element_blank(),
                                                                                      panel.grid.minor = ggplot2::element_blank(),
                                                                                      text = ggplot2::element_text(size=20),
                                                                                      legend.text = ggplot2::element_text(size=14),
                                                                                      aspect.ratio = 1) +
                        ggplot2::guides(colour = ggplot2::guide_legend(nrow = 4,byrow = TRUE))
            
                    g
                } else {
                    g <- plot(0,type='n',axes=FALSE,ann=FALSE)
                }
                
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


