#' Server handler for scfind
#'
#' @param object An SCFind object that the shiny app will be deployed upon
#'
#' @name scfindShinyServer
#' @aliases scfindShinyServer
#'
#' @importFrom shiny reactive updateCheckboxGroupInput renderPlot stopApp checkboxGroupInput observeEvent observe reactiveVal renderText renderPrint
#' @importFrom DT renderDataTable datatable
#' @importFrom data.table as.data.table data.table
#' @importFrom ggplot2 ggplot geom_bar geom_col ggtitle xlab ylab aes coord_flip theme_minimal
server.scfind <- function(object)
{
    message(getwd())
    return(
        function(input, output, session)
        {
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
            
            output$title <- renderText({
                if(all(startsWith(object@index$genes(), "chr") == T)) "What's your peak list today ?" else "What's your gene list today ?"
            })
            
            observeEvent(input$keyPressed,{
                TxtInput(input$geneList)
                
                # This is necessary, because users are allowed to input chromosome as the format chrX:XXXXX-YYYYY instead of chrX_XXXXX_YYYYY
                text <- gsub("\\s|,", ",", TxtInput())
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
                
                gene.list.input <- caseCorrect(object, gene.list.input)
                gene.list(gene.list.input)
                
            })
            
            recommended.queries <- reactive({
                progress <- Progress$new(session, min=0)
                on.exit(progress$close())
                progress$set(message = 'Dream, dream, dream...',
                             detail = 'Genes come true!')
                selected.genes <- gene.list()
                selected.datasets <- input$datasetCheckbox

                if (length(selected.genes) > 1 && any(caseCorrect(object, selected.genes) %in% object@index$genes()) && !is.null(selected.datasets))
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
                
                
                
                if(last.query.state() == "query_optimizer" && !is.null(input$queryOptimizer_rows_selected))
                {
                    selection <- qo.output()
                }
                else if (last.query.state() == "checkbox" )
                {
                    selection <-  input$geneCheckbox
                }
                
                if(length(list.indx) != 0) selection <- gsub('_', '-', sub('_', ':', selection[list.indx]))
                
                if(length(gene.list()) != 0  && sum(gene.support()$support) != 0 && !is.null(input$datasetCheckbox)) {
                    checkboxGroupInput("geneCheckbox", label='', choices = genes.in.list , selected = selection, inline = F)
                }
                # else
                # {
                #           checkboxGroupInput("geneCheckbox", label='', choices = c() , selected = c(), inline = F)
                # }
            })
            
            
            output$selectedDataset <- renderText({
                datasetName <- gsub("/", "", session$clientData$url_pathname)
                datasets <- object@datasets
                box.selection <-  input$datasetCheckbox

                datasetName <- if(length(datasetName) != 0 && datasetName != '#' && datasetName != '') dataset_def_names[[datasetName]] else ''
                
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
                             detail = 'This may take a few minutes')
                
                exampleCT <- sample(cellTypeNames(object), 1)
                exampleQuery <- cellTypeMarkers(object, exampleCT)$genes
                
                if(all(startsWith(object@index$genes(), "chr") == T)) exampleQuery <- gsub("_", "-", sub("_", ":", exampleQuery))
                
                updateTextInput(session, "geneList", value = "", placeholder = exampleQuery)
            })
            
            observe({
                datasets <- object@datasets
                box.selection <- object@datasets
                
                if(input$selectall == 0) {
                    return(NULL)
                } else if (input$selectall%%2 == 0)  {
                    updateCheckboxGroupInput(session, "datasetCheckbox", label = '', choices = datasets, selected = box.selection, inline = T)
                    
                }  else {
                    updateCheckboxGroupInput(session, "datasetCheckbox", label = '', choices = datasets, inline = T)
                    
                }
            })
            
            output$logicPanel <- renderUI ({
                df <- gene.support()
                
                number.of.choices <- length(gene.list())
                if(number.of.choices < 4){
                    histoHeight = 52.5*4
                } else {
                    histoHeight = 100 + (number.of.choices-4) * 30
                }
                
                if (nrow(df) != 0 && any(df$genes %in% object@index$genes()))
                {
                    tabsetPanel(type = "tabs",
                                tabPanel(paste(if(all(startsWith(tolower(df$genes), "chr"))) "Peaks" else "Genes", "Summary"),
                                         uiOutput("geneCheckbox"),
                                         plotOutput("geneSupportHisto", width = 250, height = histoHeight)
                                )
                    )
                } else {
                    ""
                }
            })
            
            output$evaluateSum <- renderUI ({
                s = input$cellTypesData_rows_selected # return row number, NULL
                
                if(!is.null(s) && length(object@metadata) != 0){
                    tabsetPanel(type = "tabs",
                                tabPanel("UMAP",
                                         selectInput('subdataset', 'Datasets', selected = NULL, choices = NULL, selectize=TRUE),
                                         plotOutput("cellUMAP", width = 400, height = 500)),
                                tabPanel("Evaluate Markers",
                                         dataTableOutput("evaluateCtMarkers"))
                                
                    )
                } else {
                    ""
                }
            })
            
            
            observe({
                s = input$cellTypesData_rows_selected # return row number, NULL
                df <- cell.types()
                if(!is.null(s) && length(object@metadata) != 0 && nrow(df) != 0 && !is.null(input$datasetCheckbox)){ 
                    
                    rdt <- phyper.test(object, df, input$datasetCheckbox)
                    
                    selectedCellTypes <- as.character(rdt$cell_type[s])
                    subdataset <- unique(sub("\\..*", "", selectedCellTypes)) # get datasetName. from datasetName.cellType
                    
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
                              options = list(columnDefs = list(list(width = '70px', targets = c(2, 3)), list(width = '10px', targets = c(0))), pageLength = 5,
                                             autoWidth = TRUE,
                                             dom = 'Bfrtip',
                                             buttons = c('copy', 'csv', 'excel')),
                              extensions = 'Scroller', colnames = c(col , 'Sub-queries', 'TF-IDF', 'Cell No.'),
                              rownames = F)
                } else {
                    data.table()
                }
            })
            
            cell.types <- reactive({
                selection <- if(all(startsWith(object@index$genes(), "chr") == T)) input$geneCheckbox else caseCorrect(object, input$geneCheckbox)
                
                if (length(selection) != 0 && !is.null(input$datasetCheckbox)){
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
                if(!is.null(dataset.selection) && length(gene.selection) != 0){
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
                
                if (nrow(df) != 0 && !is.null(input$datasetCheckbox) && TxtInput() != ''  && nrow(gene.support()) != 0)
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
                
                if(!is.null(df$cell_type) && nrow(df) != 0 && !is.null(input$datasetCheckbox)) { #!#
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
                if (nrow(df) != 0 && !is.null(input$datasetCheckbox) && TxtInput() != '' && nrow(gene.support()) != 0)
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
                
                if(!is.null(df$cell_type) && TxtInput() != '' && nrow(gene.support()) != 0) {
                    list.indx <- intersect(grep('chr', input$geneCheckbox), grep('_', input$geneCheckbox))
                    
                    paste("Scfind found ", length(unique(cell.types()$cell_type)), if ((length(unique(cell.types()$cell_type))) < 2) " cell type" else " cell types for the <I>selected query</I>:")
                } else {
                    if (nrow(recommended.queries()) != 0 && length(gene.list()) != 0 && !is.null(input$datasetCheckbox) && nrow(gene.support()) != 0)  paste("<span style = 'color: rgb(0, 180, 204);'>Expecting specific cell types?", "Select above recommended query!</span>", sep="<br>") else ""
                }
            })
            
            output$suggestHyper <- renderText({
                if(!is.null(TxtInput())) {
                    if(is.null(input$datasetCheckbox) && TxtInput() != '' ){
                        "Oops! dataset is missing."
                    } else if(TxtInput() != '' && sum(gene.support()$support) == 0){
                        
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
                        else
                        {
                            paste0(
                                "<p>Hmm,<br>", "There's no cell type representing your query in<br>",
                                toString(input$datasetCheckbox), ".<br></p>",
                                "<p>Please try other datasets or gene combinations</p>"
                            )
                        }
                    } else {
                        ""
                    }
                } else {
                    ""
                }
                    
            })
            
            
            output$lociRegions <- renderText({
                if(all(startsWith(object@index$genes(), "chr") == T)){
                    text <- gsub("\\s|,", ",", TxtInput())
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
                            
                            ucsc_end <- paste0("' class='ucsc' target='_blank' title='Click to view Mouse mm9: ", locus$input, " UCSC Genome Browser'>", "【")
                            ucsc_a <- paste(ucsc, locus$chromosome, "%3A", locus[2], "%2D", locus[3], ucsc_end, sep="")
                            
                            gap <- locus[2] - locus[1]
                            
                            chromosome <- if(gap < 200) c() else paste0(locus[1], (if(gap<=400) paste(rep("&#8594;", (gap/200))) else paste("&#8594;", "...", "&#8594;")))
                            
                            gap <- locus[3] - locus[2]
                            chromosome <- paste0(chromosome, ucsc_a, locus[2], paste(rep("&#187;", (gap/100)),collapse = ""), locus[3], "】", "</a>" )
                            
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
                
                number.of.choices <- length(gene.list())
                if(number.of.choices < 4){
                    histoHeight = 52.5*4
                } else {
                    histoHeight = 100 + (number.of.choices-4) * 30
                }
                plotOutput("geneSupportHisto", width = 250, height = histoHeight)
            })
            
            
            output$geneSupportHisto <- renderPlot({
                df <- gene.support()
                
                if (nrow(df) != 0 && length(grep(TRUE, (df$genes %in% object@index$genes()))) != 0)
                {
                    max.axis <- max(df$support)
                    g <- ggplot(df, aes(x = genes, y = support), colour = support) +
                        xlab("") +
                        ylab("Cells") +
                        ylim(0, max.axis) +
                        geom_text(aes(label=df$genes), vjust = "inward", hjust = "inward", size = 4) +
                        geom_col(fill = rainbow(length(df$genes)), position = "dodge", width = 1, alpha = .5, colour = "black") +
                        theme_minimal() +
                        coord_flip() +
                        theme(panel.grid.minor = element_blank(), axis.text.y = element_blank(),
                              panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                        scale_x_discrete(limits = rev(df$genes))
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
                
                umap_highlight <- data.frame()
                
                
                df <- cell.types()
                
                if(!is.null(s) && length(object@metadata) != 0 && nrow(df) != 0 && !is.null(input$datasetCheckbox)){
                    
                    rdt <- phyper.test(object, df, input$datasetCheckbox)
                    
                    
                    umapChoice = input$subdataset[length(input$subdataset)]
                    
                    getDatasetUmap <- object@metadata[[paste0("umap.", umapChoice)]]
                    
                    selectedCellTypes <- as.character(rdt$cell_type[s])
                    
                    subCellTypes <- grep(paste0(umapChoice, "."), selectedCellTypes, value = T)
                    subCellTypes <- sub(paste0(umapChoice, "."), "", subCellTypes)
                    # get hits
                    if (length(selection) != 0 && length(subCellTypes) != 0){ 
                        true.hits <- findCellTypes.geneList(object, selection, umapChoice)
                        
                        for( j in subCellTypes)
                        {
                            getCluster <- subset(getDatasetUmap, rownames(getDatasetUmap) == j )
                            getCluster <- data.frame(x = getCluster[,1], y = getCluster[,2], col = rownames(getCluster), shape = F)
                            getCluster$shape[true.hits[[paste0(umapChoice, ".", j)]]] <- T
                            umap_highlight <- data.frame(rbind(umap_highlight, getCluster))
                        }
                    } 
                    
                    umap_plot <- data.frame(x = getDatasetUmap[,1], y = getDatasetUmap[,2], col = rownames(getDatasetUmap))
                    umap_truehits <- data.frame(subset(umap_highlight, shape == T))
                    
                    
                    g <- ggplot() + geom_point(aes(x, y, group = col), data = umap_plot, colour = alpha("grey", .3)) +
                        geom_point(aes(x, y, colour = col), data = umap_highlight, alpha = .5, shape = 21)+
                        geom_point(aes(x, y), data = umap_truehits, shape = 21, colour = "black")+
                        labs(x = "UMAP1", y = "UMAP2", color = "") +
                        theme_bw() + 
                        theme(axis.line = element_line(colour = "black"),
                              panel.border = element_rect(colour = "black", fill=NA, size=2),
                              plot.background = element_blank(),
                              legend.position = "bottom",
                              axis.text.x = element_blank(),
                              axis.text.y = element_blank(),
                              axis.ticks = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              text = element_text(size=20),
                              legend.text = element_text(size=14),
                              aspect.ratio = 1) +
                        guides(colour = guide_legend(nrow = 4,byrow = TRUE))
                    
                    g
                } else {
                    g <- plot(0,type='n',axes=FALSE,ann=FALSE)
                }
                
            })
            
            session$onSessionEnded(function() {
                stopApp()
            })
        }
    )
}




#' @rdname scfindShinyServer
#' @aliases scfindShinyServer
setMethod("scfindShinyServer", signature(object = "SCFind"), server.scfind)


#' Server handler for scfind free text search
#'
#' @param object An SCFind object that the shiny app will be deployed upon
#'
#' @name scfindShinyW2VServer
#' @aliases scfindShinyServer
#'
#' @importFrom shiny reactive updateCheckboxGroupInput renderPlot stopApp checkboxGroupInput observeEvent observe reactiveVal renderText renderPrint
#' @importFrom DT renderDataTable datatable
#' @importFrom data.table as.data.table data.table
#' @importFrom ggplot2 ggplot geom_bar geom_col ggtitle xlab ylab aes coord_flip theme_minimal
server.scfind.w2v <- function(object, dictionary)
{
    
    return(
        function(input, output, session)
        {
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
            
            output$title <- renderText({
                "Search by genes or natural language<br><span style='font-size: 15px;'>(version 1.0)</span>"
            })
            
            observeEvent(input$keyPressed,{
                TxtInput(input$geneList)
            })
            
            observe({
                if(!is.null(TxtInput()))
                {
                    progress <- Progress$new(session, min=0)
                    on.exit(progress$close())
                    progress$set(message = 'Dream, dream, dream...',
                                 detail = 'Genes come true!')
                    gene.list.input <- suppressMessages(query2genes(object = object, dictionary = dictionary, query = TxtInput(), greedy = greedy(), strict = T))
                    
                    
                    # Take care of duplicated genes
                    input.and <- setdiff(gene.list.input$and, c(gene.list.input$or, gene.list.input$not, gene.list.input$ornot))
                    input.or <- setdiff(gene.list.input$or, c(gene.list.input$not, gene.list.input$ornot))
                    input.ex <- setdiff(gene.list.input$not, gene.list.input$ornot)
                    
                    gene.list(unique(input.and))
                    gene.list.or(unique(input.or))
                    gene.list.ex(unique(input.ex))
                    gene.list.ex.or(unique(gene.list.input$ornot))
                }

            })
            
            observeEvent(input$greedy,{
                if(gsub("\\s", "", input$geneList) != "")
                {
                    updateTextInput(session, "geneList", value = paste0(TxtInput(), " "), placeholder = sample(exampleQuery,1))
                }
                else
                {
                    updateTextInput(session, "geneList", value = "", placeholder = sample(exampleQuery,1))
                }
                last.query.state("genelist")
                greedy(input$greedy)
            })
            
            recommended.queries <- reactive({
                selected.genes <- gene.list()
                selected.datasets <- input$datasetCheckbox

                if (length(selected.genes) > 1 && any(selected.genes %in% object@index$genes()) && !is.null(selected.datasets))
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
                    available.queries <- available.queries[order(available.queries$tfidf, decreasing = T),]
                    available.queries$tfidf <- signif(available.queries$tfidf, digits = 5) 
                    available.queries
                } else {
                    available.queries  
                }
                
            })
            
            
            qo.output <- reactive({
                selected.index <- input$queryOptimizer_rows_selected
                available.queries <-  recommended.queries()
                
                selected.query <- available.queries$Query[selected.index]
                unlist(strsplit(gsub("\\s", "", selected.query), ","))
            })
            
            output$geneCheckbox <- renderUI({
                genes.in.list <- gene.list()
                ## else fallback
                selection <- gene.list()
                
                if(last.query.state() == "query_optimizer" && !is.null(input$queryOptimizer_rows_selected))
                {
                    selection <- qo.output()
                }
                else if (last.query.state() == "checkbox")
                {
                    selection <-  input$geneCheckbox
                }
                
                
                if(length(gene.list()) != 0  && sum(gene.support()$support) != 0 && !is.null(input$datasetCheckbox)) {
                    checkboxGroupInput("geneCheckbox", label = '', choices = genes.in.list , selected = selection, inline = F)
                } 
                # else 
                # {
                #     checkboxGroupInput("geneCheckbox", label='', choices = c() , selected = c(), inline = F)
                # }
            })
            
            
            observe({
                if(length(gene.list()) != 0)
                {
                    if(is.null(input$selectall.and) ) return (NULL)
                    if(input$selectall.and%%2 == 0 || last.query.state() == "genelist")
                    {
                        updateCheckboxGroupInput(session, "geneCheckbox", label = '', choices = gene.list(), selected = gene.list(), inline = F)
                    }
                    else
                    {
                        updateCheckboxGroupInput(session, "geneCheckbox", label = '', choices = gene.list(), inline = F)
                        
                    }
                }
            })
            
            observe({
                if(length(gene.list.or()) != 0)
                {
                    if(is.null(input$selectall.or) ) return (NULL)
                    if(input$selectall.or%%2 == 0 || last.query.state() == "genelist" )
                    {
                        updateCheckboxGroupInput(session, "geneOrCheckbox", label = '', choices = gene.list.or(), selected = gene.list.or(), inline = F)
                    }
                    else
                    {
                        updateCheckboxGroupInput(session, "geneOrCheckbox", label = '', choices = gene.list.or(), inline = F)
                        
                    }
                }
            })
            
            observe({
                if(length(gene.list.ex()) != 0)
                {
                    if(is.null(input$selectall.not) ) return (NULL)
                    if(input$selectall.not%%2 == 0 || last.query.state() == "genelist")
                    {
                        updateCheckboxGroupInput(session, "geneExCheckbox", label = '', choices = gene.list.ex(), selected = gene.list.ex(), inline = F)
                    }
                    else
                    {
                        updateCheckboxGroupInput(session, "geneExCheckbox", label = '', choices = gene.list.ex(), inline = F)
                        
                    }
                }
            })
            
            observe({
                if(length(gene.list.ex.or()) != 0)
                {
                    if(is.null(input$selectall.ornot) ) return (NULL)
                    if(input$selectall.ornot%%2 == 0 || last.query.state() == "genelist" )
                    {
                        updateCheckboxGroupInput(session, "geneExOrCheckbox", label = '', choices = gene.list.ex.or(), selected = gene.list.ex.or(), inline = F)
                    }
                    else
                    {
                        updateCheckboxGroupInput(session, "geneExOrCheckbox", label = '', choices = gene.list.ex.or(), inline = F)
                        
                    }
                }
            })
            
            
            observe({
                if(length((gene.list.or()) != 0 && !is.null(input$datasetCheckbox)))
                {
                    updateCheckboxGroupInput(session, "geneOrCheckbox", label = '', choices = gene.list.or(), selected = gene.list.or(), inline = F)
                }
                else
                {
                    updateCheckboxGroupInput(session, "geneOrCheckbox", label = '', choices = gene.list.or(), inline = F)
                    
                }
                
                if(length(gene.list.ex()) != 0 && !is.null(input$datasetCheckbox))
                {
                    updateCheckboxGroupInput(session, "geneExCheckbox", label = '', choices = gene.list.ex(), selected = gene.list.ex(), inline = F)
                }
                else
                {
                    updateCheckboxGroupInput(session, "geneExCheckbox", label = '', choices = gene.list.ex(), inline = F)
                    
                }
                
                if(length(gene.list.ex.or()) != 0 && !is.null(input$datasetCheckbox))
                {
                    updateCheckboxGroupInput(session, "geneExOrCheckbox", label = '', choices = gene.list.ex.or(), selected = gene.list.ex.or(), inline = F)
                }
                else
                {
                    updateCheckboxGroupInput(session, "geneExOrCheckbox", label = '', choices = gene.list.ex.or(), inline = F)
                    
                }
            })
            
            ### OR GENES
            
            output$geneOrCheckbox <-  renderUI({
                genes.in.list <- gene.list.or()
                box.selection <-  gene.list.or()
                
                if (initial.OR == "initial")
                {
                    ## set the flag down
                    initial.OR <- "initialized"
                }
                else
                {
                    box.selection <-  input@geneOrCheckbox
                    
                }
                
                checkboxGroupInput("geneOrCheckbox", label = '', choices = genes.in.list, selected = box.selection, inline = F)
                
                
            })
            ### EX GENES
            output$geneExCheckbox <-  renderUI({
                genes.in.list <- gene.list.ex()
                box.selection <-  gene.list.ex()
                
                if (initial.EX == "initial")
                {
                    ## set the flag down
                    initial.EX <- "initialized"
                }
                else
                {
                    box.selection <-  input@geneExCheckbox
                    
                }
                
                checkboxGroupInput("geneExCheckbox", label = '', choices = genes.in.list, selected = box.selection, inline = F)
                
                
            })
            
            ### EX-OR GENES
            output$geneExOrCheckbox <-  renderUI({
                genes.in.list <- gene.list.ex.or()
                box.selection <-  gene.list.ex.or()
                
                if (initial.EXOR == "initial")
                {
                    ## set the flag down
                    initial.EXOR <- "initialized"
                }
                else
                {
                    box.selection <-  input@geneExOrCheckbox
                    
                }
                
                checkboxGroupInput("geneExOrCheckbox", label = '', choices = genes.in.list, selected = box.selection, inline = F)
                
                
            })
            
            ###
            
            
            output$selectedDataset <- renderText({
                datasetName <- gsub("/", "", session$clientData$url_pathname)
                datasets <- object@datasets
                box.selection <-  input$datasetCheckbox
                
                datasetName <- if(length(datasetName) != 0 && datasetName != '#' && datasetName != '') dataset_def_names[[datasetName]] else ''
                
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
            
            
            
            # Progress bar while initializing index
            observe({
                progress <- Progress$new(session, min=0)
                on.exit(progress$close())
                progress$set(message = 'Initialising index...',
                             detail = 'This may take a few minutes')
                
                updateTextInput(session, "geneList", value = "", placeholder = sample(exampleQuery,1))
            })
            
            observe({
                datasets <- object@datasets
                box.selection <- object@datasets
                
                if(input$selectall == 0) {
                    return(NULL) 
                } else if (input$selectall%%2 == 0) {
                    updateCheckboxGroupInput(session, "datasetCheckbox", label = '', choices = datasets, selected = box.selection, inline = T)
                    
                    last.query.state("query_optimizer")

                                    } else {
                    updateCheckboxGroupInput(session, "datasetCheckbox", label = '', choices = datasets, inline = T)
                    print('dataset not defined')
                }
                
            })
            
            # UI of the evaluation panels
            output$paraPanel <- renderUI ({
                fluidRow(
                    tags$h3(renderText("Parameters")),
                    sliderInput("greedy", label = "Greedy:", min = 0, max = 1, value = .6, step = .1)
                )
            })
            
            output$logicPanel <- renderUI ({
                df <- gene.support()
                
                
                if (nrow(df) != 0 && any(df$genes %in% object@index$genes()))
                {
                    tabsetPanel(type = "tabs",
                                tabPanel("Query",
                                         if(!is.null(operatorfootprint()))
                                         {
                                             selectInput('operating', '', selected = operatorfootprint()[1],  choices = operatorfootprint(), selectize = T)
                                         }
                                         else 
                                         {
                                             selectInput('operating', '', selected = NULL, choices = NULL, selectize=TRUE)    
                                         },
                                         
                                         uiOutput("operatorUI")
                                ),
                                tabPanel("Summary",
                                         renderUI({
                                             fluidRow(
                                                 tags$h4(renderText(paste0("Greedy: ", greedy()))),
                                                 if(!is.null(input$geneCheckbox))
                                                 {
                                                     if(last.query.state() == "query_optimizer")
                                                     {
                                                         tags$h4(renderText(paste0("AND genes: ", toString(qo.output()))))
                                                     }
                                                     else if (last.query.state() == "checkbox" )
                                                     {
                                                         tags$h4(renderText(paste0("AND genes: ", toString(input$geneCheckbox))))
                                                     }
                                                     
                                                 },
                                                 
                                                 tags$h4(renderText(paste0(if(length(gene.list.or()) != 0) "OR genes: " else "", if(is.null(input$geneOrCheckbox)) toString(gene.list.or()) else toString(input$geneOrCheckbox)))),
                                                 tags$h4(renderText(paste0(if(length(gene.list.ex()) != 0) "NOT genes: " else "",  if(is.null(input$geneExCheckbox)) toString(gene.list.ex()) else toString(input$geneExCheckbox)))),
                                                 tags$h4(renderText(paste0(if(length(gene.list.ex.or()) != 0) "ORNOT genes: " else "",  if(is.null(input$geneExOrCheckbox)) toString(gene.list.ex.or()) else toString(input$geneExOrCheckbox))))
                                             )
                                         })
                                )
                    )
                } else {
                    ""
                }
            })
            
            
            output$evaluateSum <- renderUI ({
                s = input$cellTypesData_rows_selected # return row number, NULL
                
                if(!is.null(s) && length(object@metadata) != 0 ){
                    tabsetPanel(type = "tabs",
                                tabPanel("UMAP", 
                                         selectInput('subdataset', 'Datasets', selected = NULL, choices = NULL, selectize=TRUE),
                                         plotOutput("cellUMAP", width = 400, height = 500)),
                                tabPanel("Evaluate Markers",  
                                         dataTableOutput("evaluateCtMarkers")),
                                tabPanel("WordCloud",
                                         plotOutput("wordcloud"))
                                
                    )
                } else {
                    ""
                }
            })
            
            
            
            output$operatorUI <- renderUI({
                if(input$operating == LogicChoices[['and']] || last.query.state() == "query_optimizer")
                {
                    fluidRow(
                        tags$h3(actionLink("selectall.and","Select/Deselect All")),
                        tags$h4(renderText(paste0("Selected genes: ", length(input$geneCheckbox), "/", length(gene.list())))),
                        uiOutput("geneCheckbox"),
                        uiOutput("geneHisto")
                    )
                }
                else if(input$operating == LogicChoices[['or']])
                {
                    fluidRow(
                        tags$h3(actionLink("selectall.or","Select/Deselect All")),
                        tags$h4(renderText(paste0("Selected genes: ", length(input$geneOrCheckbox), "/", length(gene.list.or())))),
                        uiOutput("geneOrCheckbox"),
                        uiOutput("geneHisto")
                    )
                }
                else if(input$operating == LogicChoices[['not']])
                {
                    fluidRow(
                        tags$h3(actionLink("selectall.not","Select/Deselect All")),
                        tags$h4(renderText(paste0("Selected genes: ", length(input$geneExCheckbox), "/", length(gene.list.ex())))),
                        uiOutput("geneExCheckbox"),
                        uiOutput("geneHisto")
                    )
                }
                else if(input$operating == LogicChoices[['ornot']])
                {
                    fluidRow(
                        tags$h3(actionLink("selectall.ornot","Select/Deselect All")),
                        tags$h4(renderText(paste0("Selected genes: ", length(input$geneExOrCheckbox), "/", length(gene.list.ex.or())))),
                        uiOutput("geneExOrCheckbox"),
                        uiOutput("geneHisto")
                    )
                }
                
                
                
            })
            
            observe({
                operators <- c(
                    if(length(gene.list())!= 0) LogicChoices[['and']],
                    if(length(gene.list.or())!= 0) LogicChoices[['or']],
                    if(length(gene.list.ex())!= 0) LogicChoices[['not']],
                    if(length(gene.list.ex.or())!= 0) LogicChoices[['ornot']]
                )
                
                if(length(operators) != 0)
                {
                    updateSelectInput(session, 'operating', label = '', choices = operators)
                }
                
                operatorfootprint(operators)
                
            })
            
            observe({
                s = input$cellTypesData_rows_selected # return row number, NULL
                df <- cell.types()
                if(!is.null(s) && length(object@metadata) != 0 && nrow(df) != 0 && !is.null(input$datasetCheckbox)){ #!#
                    
                    rdt <- phyper.test(object, df, input$datasetCheckbox)
                    rdt <- data.frame(rdt[order(rdt$pval, decreasing = F), ])
                    
                    selectedCellTypes <- as.character(rdt$cell_type[s])
                    subdataset <- unique(sub("\\..*", "", selectedCellTypes)) # get datasetName. from datasetName.cellType
                    
                    if(length(subdataset) > 1){
                        updateSelectInput(session, 'subdataset', label = 'Datasets', choices = subdataset,  selected = subdataset[length(subdataset)])
                    } else {
                        updateSelectInput(session, 'subdataset', label = 'Datasets', choices = subdataset,  selected = subdataset)
                    }
                    
                }
                
                
            })
            
            output$queryOptimizer <- renderDataTable({
                if(!is.null(recommended.queries())){
                    
                    
                    datatable(recommended.queries(), selection = 'single',
                              options = list(columnDefs = list(list(width = '70px', targets = c(2, 3)), list(width = '10px', targets = c(0))), pageLength = 5,
                                             autoWidth = TRUE,
                                             dom = 'Bfrtip',
                                             buttons = c('copy', 'csv', 'excel')),
                              extensions = 'Scroller', colnames = c("Genes" , 'Sub-queries', 'TF-IDF', 'Cell No.'),
                              rownames = F)
                } else {
                    data.table()
                }
            })
            
            cell.types <- reactive({
                selection <- c(input$geneCheckbox,
                               if(length(gene.list.or()) != 0) { paste0("*", if(is.null(input$geneOrCheckbox)) gene.list.or() else input$geneOrCheckbox ) },
                               if(length(gene.list.ex()) != 0) { paste0("-", if(is.null(input$geneExCheckbox)) gene.list.ex() else input$geneExCheckbox ) },
                               if(length(gene.list.ex.or()) != 0) { paste0("*-", if(is.null(input$geneExOrCheckbox)) gene.list.ex.or() else input$geneExOrCheckbox )})
                
                if (length(selection) != 0 && !is.null(input$datasetCheckbox) ){ 
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
                gene.selection <- gene.list()
                dataset.selection <- input$datasetCheckbox
                if(!is.null(dataset.selection) && length(gene.selection) != 0){
                    gene.support <- as.data.frame(object@index$genesSupport(gene.selection, dataset.selection))
                    dimnames(gene.support)[[2]] <- 'support'
                    gene.support$genes <- rownames(gene.support)
                } else {
                    gene.support <- data.frame()
                    gene.support$support <- c()
                    gene.support$genes <-  c()
                }
                gene.support
            })
            
            
            output$cellTypesData <- renderDataTable({
                df <- cell.types()
                selection <- c(input$geneCheckbox,
                               if(is.null(input$geneOrCheckbox)) gene.list.or() else input$geneOrCheckbox,
                               if(is.null(input$geneExCheckbox)) gene.list.ex() else input$geneExCheckbox,
                               if(is.null(input$geneExOrCheckbox)) gene.list.ex.or() else input$geneExOrCheckbox)
                if (nrow(df) != 0 && !is.null(input$datasetCheckbox) && length(c(gene.list(), gene.list.or(), gene.list.ex(), gene.list.ex.or())) != 0 && length(selection) != 0 && TxtInput() != ''  && nrow(gene.support()) != 0) 
                {
                    rdt <- phyper.test(object, df, input$datasetCheckbox)
                    rdt$pval <- as.numeric(signif(rdt$pval, digits = 6))
                    rdt <- data.frame(rdt[order(rdt$pval, decreasing = F), ])
                    
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
                
                df <- cell.types()
                
                if(!is.null(df$cell_type) && nrow(df) != 0 && length(gene.list) != 0 && !is.null(input$datasetCheckbox)) { #!#
                    rdt <- phyper.test(object, df, input$datasetCheckbox)
                    rdt <- data.frame(rdt[order(rdt$pval, decreasing = F), ])
                    
                    mge <- if(length(unique(df$cell_type)) < 2) evaluateMarkers(object, selection, as.character(unique(df$cell_type))) else evaluateMarkers(object, selection, rdt$cell_type[s])
                    if(nrow(mge) != 0){
                        mge <- mge[order(mge$genes),]
                        mge <- mge[order(mge[,1]),]
                        mge$genes <- mge$genes
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
                selection <- c(input$geneCheckbox,
                               if(is.null(input$geneOrCheckbox)) gene.list.or() else input$geneOrCheckbox,
                               if(is.null(input$geneExCheckbox)) gene.list.ex() else input$geneExCheckbox,
                               if(is.null(input$geneExOrCheckbox)) gene.list.ex.or() else input$geneExOrCheckbox)
                selection.text <- toString(input$geneCheckbox)
                
                df <- cell.types()
                if (nrow(df) != 0 && !is.null(input$datasetCheckbox) && TxtInput() != '' && nrow(gene.support()) != 0)
                {
                    rdt <- phyper.test(object, df, input$datasetCheckbox)
                    rdt <- data.frame(rdt[order(rdt$pval, decreasing = F), ])
                    
                    if(nrow(rdt) > 0 && length(c(gene.list(), gene.list.or(), gene.list.ex(), gene.list.ex.or())) != 0 && length(selection) != 0){
                        if(nrow(rdt) == 1){
                            if(is.null(s)) '' else paste0("Marker evaluation of ", selection.text, " in ", rdt$cell_type, ": ")
                        } else {
                            if(!is.null(selection.text) && length(s) != 0){
                                paste0("Marker evaluation of ", selection.text, " in ", length(s), if(length(s) < 2) " cell type:" else " cell types:")
                            } else {
                                paste0("<span style = 'color: rgb(0, 180, 204);'>Select cell types from cell type table to start evaluation of your query.")
                            }
                        }
                        
                    } else {
                        ""
                    }
                }
                
            })
            
            
            output$selectedQuery <- renderText({
                selectedList = input$queryOptimizer_rows_selected
                selection <- c(input$geneCheckbox,
                               if(is.null(input$geneOrCheckbox)) gene.list.or() else input$geneOrCheckbox,
                               if(is.null(input$geneExCheckbox)) gene.list.ex() else input$geneExCheckbox,
                               if(is.null(input$geneExOrCheckbox)) gene.list.ex.or() else input$geneExOrCheckbox)
                df <- cell.types()
                if(length(c(gene.list(), gene.list.or(), gene.list.ex(), gene.list.ex.or())) != 0 && length(selection) != 0)
                {
                    if(!is.null(df$cell_type) && TxtInput() != '' && nrow(gene.support()) != 0) {        
                        paste(if( TxtInput() != '' && sum(gene.support()$support) == 0 ) "But",
                              "Scfind found ", length(unique(cell.types()$cell_type)), if ((length(unique(cell.types()$cell_type))) < 2) " cell type" else " cell types for your query:")
                    } else {
                        if (nrow(recommended.queries()) != 0 && length(gene.list()) != 0 && !is.null(input$datasetCheckbox) && nrow(gene.support()) != 0)  paste("<span style = 'color: rgb(0, 180, 204);'>Expecting specific cell types?", "Select above recommended query!</span>", sep="<br>") else ""
                    }
                }
                
            })
            
            output$suggestHyper <- renderText({
                if(!is.null(TxtInput())) {
                    if(is.null(input$datasetCheckbox) && TxtInput() != '' ){
                        "Oops! dataset is missing."
                    } 
                    else if(TxtInput() != '' && sum(gene.support()$support) == 0 && length(gene.list.ex()) != 0)
                    {
                        # Add did you mean here ?
                        paste("Oops! No best match gene set for your query.")
                        
                    } 
                    else if(length(gene.list()) != 0)
                    {
                        if(sum(gene.support()$support) != 0 && length(recommended.queries()) != 0){ 
                            "Optimized queries for AND genes: "
                        } 
                        else
                        {
                            paste0(
                                "<p>Hmm,<br>", "There's no cell type representing your query in<br>",
                                toString(input$datasetCheckbox), ".<br></p>",
                                "<p>Please try other datasets or gene combinations</p>"
                            )
                        }
                    } 
                    else 
                    {
                        '<p>Have you tried to use operators "NOT", "OR" & "ORNOT" in your query?<br>By the way, `scfind` also supports query with RS and MeSH IDs.</p>Have fun with your single cell data analysis!<br>(Note: Query without operators will be treated as "AND")'
                    }
                } else {
                    ""
                }
            })
            
            
            # histogram on left hand side
            output$geneHisto <- renderUI({
                
                number.of.choices <- if(input$operating == LogicChoices[['and']]) length(gene.list()) else if (input$operating == LogicChoices[['or']]) length(gene.list.or()) else if (input$operating == LogicChoices[['not']]) length(gene.list.ex()) else if (input$operating == LogicChoices[['ornot']]) length(gene.list.ex.or())
                if(number.of.choices < 4){
                    histoHeight = 52.5*4
                } else {
                    histoHeight = 100 + (number.of.choices-4) * 30
                }
                plotOutput("geneSupportHisto", width = 250, height = histoHeight)
            })
            
            
            
            output$geneSupportHisto <- renderPlot({
                
                dataset.selection <- input$datasetCheckbox
                if(input$operating == LogicChoices[['and']] )
                {
                    plotCellNumberHisto(gene.list(), dataset.selection, object)
                }
                else if (input$operating == LogicChoices[['or']])
                {
                    plotCellNumberHisto(gene.list.or(), dataset.selection, object)
                }
                else if (input$operating == LogicChoices[['not']])
                {
                    plotCellNumberHisto(gene.list.ex(), dataset.selection, object)
                }
                else if (input$operating == LogicChoices[['ornot']])
                {
                    plotCellNumberHisto(gene.list.ex.or(), dataset.selection, object)
                }

                
            })
            
            
            output$cellUMAP <- renderPlot({
                s = input$cellTypesData_rows_selected 
                
                selection <- c(input$geneCheckbox,
                               if(length(gene.list.or()) != 0) { paste0("*", if(is.null(input$geneOrCheckbox)) gene.list.or() else input$geneOrCheckbox ) },
                               if(length(gene.list.ex()) != 0) { paste0("-", if(is.null(input$geneExCheckbox)) gene.list.ex() else input$geneExCheckbox ) },
                               if(length(gene.list.ex.or()) != 0) { paste0("*-", if(is.null(input$geneExOrCheckbox)) gene.list.ex.or() else input$geneExOrCheckbox )})
                
                umap_highlight <- data.frame()
                
                df <- cell.types()
                
                if(!is.null(s) && length(object@metadata) != 0 && nrow(df) != 0 && !is.null(input$datasetCheckbox)){
                    
                    rdt <- phyper.test(object, df, input$datasetCheckbox)
                    rdt <- data.frame(rdt[order(rdt$pval, decreasing = F), ])
                    
                    umapChoice = input$subdataset[length(input$subdataset)]
                    
                    getDatasetUmap <- object@metadata[[paste0("umap.", umapChoice)]]
                    
                    selectedCellTypes <- as.character(rdt$cell_type[s])
                    
                    subCellTypes <- grep(paste0(umapChoice, "."), selectedCellTypes, value = T)
                    subCellTypes <- sub(paste0(umapChoice, "."), "", subCellTypes)
                    
                    # get hits
                    if (length(selection) != 0 && length(subCellTypes) != 0){ 
                        true.hits <- findCellTypes.geneList(object, selection, umapChoice)
                        
                        for( j in subCellTypes)
                        {
                            getCluster <- subset(getDatasetUmap, rownames(getDatasetUmap) == j )
                            getCluster <- data.frame(x = getCluster[,1], y = getCluster[,2], col = rownames(getCluster), shape = F)
                            getCluster$shape[true.hits[[paste0(umapChoice, ".", j)]]] <- T
                            umap_highlight <- data.frame(rbind(umap_highlight, getCluster))
                        }
                    } 
                    
                    
                    umap_plot <- data.frame(x = getDatasetUmap[,1], y = getDatasetUmap[,2], col = rownames(getDatasetUmap))
                    umap_truehits <- data.frame(subset(umap_highlight, shape == T))
                    
                    g <- ggplot() + geom_point(aes(x, y, group = col), data = umap_plot, colour = alpha("grey", .3)) +
                        geom_point(aes(x, y, colour = col), data = umap_highlight, alpha = .5, shape = 21)+
                        geom_point(aes(x, y), data = umap_truehits, shape = 21, colour = "black")+
                        labs(x = "UMAP1", y = "UMAP2", color = "") +
                        theme_bw() + 
                        theme(axis.line = element_line(colour = "black"),
                              panel.border = element_rect(colour = "black", fill=NA, size=2),
                              plot.background = element_blank(),
                              legend.position = "bottom",
                              axis.text.x = element_blank(),
                              axis.text.y = element_blank(),
                              axis.ticks = element_blank(),
                              panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              text = element_text(size=20),
                              legend.text = element_text(size=14),
                              aspect.ratio = 1) +
                        guides(colour = guide_legend(nrow = 4,byrow = TRUE))
                    
                    g
                } else {
                    g <- plot(0,type='n',axes=FALSE,ann=FALSE)
                }
                
            })
            
            output$wordcloud <- renderPlot({
                selection <- c(input$geneCheckbox, if(length(gene.list.or()) != 0) {if(is.null(input$geneOrCheckbox)) gene.list.or() else input$geneOrCheckbox})
                s = input$cellTypesData_rows_selected
                
                if(length(selection) != 0 && !is.null(s))
                {
                    
                    selection <- if(length(selection) > 15) sample(selection,10) else selection
                    
                    withProgress(message = "Hold on tight, we're landing on your wordcloud...", value = 0, {
                        
                        suppressMessages(gene2wordcloud(dictionary = dictionary, gene.list = selection))
                    })
                    
                }
                else
                {
                    plot(0,type='n',axes=FALSE,ann=FALSE)
                }
            })
            
            
            session$onSessionEnded(function() {
                stopApp()
            })
        })
}




#' @rdname scfindShinyW2VServer
#' @aliases scfindShinyW2VServer
setMethod("scfindShinyW2VServer", signature(object = "SCFind"), server.scfind.w2v)
