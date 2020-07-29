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
        
        
        tags$head(
            tags$style(HTML("
                            body {
                            background-image: url('https://scfind.sanger.ac.uk/img/scfind.png');
                            background-size: 200px;
                            background-attachment: fixed;
                            background-repeat: no-repeat;
                            background-position: center;
                            }
                            
                            .shiny-notification {
                            position:fixed;
                            height: 100px;
                            width: 400px;
                            top: 12%;
                            left: calc(50% - 200px);
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
                            font-size: 20px;
                            margin-bottom:0;
                            }
                            
                            #geneList {
                            position: fixed; 
                            margin-top:0;
                            left: 35%;
                            right:35%;
                            height: 35px;
                            width:30%;
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
                            /*top: 20%;*/
                            top:25%;
                            }

                            #main h3 {
                            margin-top: 0;
                            text-align: center;
                            color: rgb(0, 180, 204);
                            }

                            #logicPanel h3 {
                            font-size: 14px;
                            }                               

                            #logicPanel h4 {
                            font-size: 14px;
                            padding-left: 5%;
                            }                          

                            #geneCheckbox,
                            #geneOrCheckbox,
                            #geneExCheckbox,
                            #geneExOrCheckbox {
                            float:left;  
                            width:25%; 
                            z-index: 0;
                            }

                            #geneCheckbox:hover,
                            #geneOrCheckbox:hover,
                            #geneExCheckbox:hover,
                            #geneExOrCheckbox:hover {
                            z-index: 1;
                            }
                            
                            #geneCheckbox .shiny-options-group,
                            #geneOrCheckbox .shiny-options-group,
                            #geneExCheckbox .shiny-options-group,
                            #geneExOrCheckbox .shiny-options-group
                            {
                            left: 0;
                            border: 0;
                            padding: 0;
                            margin: 0;
                            display: flex;
                            flex-direction: column;
                            justify-content: space-between;
                            opacity: 0.6;
                            }

                            #geneCheckbox .shiny-options-group:hover,
                            #geneOrCheckbox .shiny-options-group:hover,
                            #geneExCheckbox .shiny-options-group:hover,
                            #geneExOrCheckbox .shiny-options-group:hover
                            {
                            opacity: 1;
                            }
                            
                            #geneCheckbox .checkbox,
                            #geneOrCheckbox .checkbox,
                            #geneExCheckbox .checkbox,
                            #geneExOrCheckbox .checkbox
                            {
                            width: 150px;
                            height: 20px;
                            overflow-wrap: break-word;
                            }

                            h4 {
                            font-size: 18px;
                            color: #455254;
                            background-color: rgba(256,256,256, 0.6);
                            }
                            
                            .checkbox {
                            width: 20%;
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
                            float:left; 
                            width:50%; 
                            padding-top:30px;
                            /*padding-top: 30px;
                            margin-top: 20px;
                            margin-left: 80px;
                            background-color: transparent;*/
                            }

                            #wordcloud img {
                            float:left;
                            height:100%;
                            width:100%;
                            padding-top:30px;
                            }
                            
                            .datasetsBox {
                            position: fixed; 
                            bottom: -200px; 
                            left: 15%;
                            right:15%;
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

                            .parameters {
                            float: right;
                            padding-left:5%;
                            width:30%;
                            }

                            
                            "))
        ),
        actionButton(inputId = 'homeBtn', label = "Home", icon('th'), onclick = "window.open('https://scfind.sanger.ac.uk', '_self')"),
        fluidRow(id = "search",
                 column(12, align="center",
                        tags$h3(uiOutput("title")),
                        tags$script(js),
                        textInput("geneList", label = "", value = "", placeholder = "Initialising. This may take a few minutes...")
                 )
        ),
        fluidRow(id = "main",
                 column(3,
                        uiOutput("lociRegions"),
                        fluidRow(id = "selectGenes",
                                 uiOutput("logicPanel")
                        )
                 ),
                 column(4, id = "query",
                        tags$div(class = "datasetsBox",
                                 tags$div( class = "parameters", uiOutput("paraPanel")),
                                 tags$span(class="titleDataset",                                            
                                           tags$h3(uiOutput("selectedDataset")),
                                           actionLink("selectall","Select/Deselect All")),
                                 uiOutput("datasetCheckbox")),
                        tags$h4(uiOutput("suggestHyper")),
                        dataTableOutput("queryOptimizer"),
                        tags$h4(uiOutput("selectedQuery")),
                        dataTableOutput("cellTypesData")
                        
                 ),
                 column(5, id = "celltype",
                        tags$h4(uiOutput("ctData")),
                        uiOutput("evaluateSum")
                 )
        )
    )
}
