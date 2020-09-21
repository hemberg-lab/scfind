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
                                      tags$h3(uiOutput("title")),
                                      tags$script(js),
                                      tags$script(inactivity),
                                      textInput("geneList", label = "", value = "")
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
