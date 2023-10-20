## Global variables and functions for the shiny app of scfind

utils::globalVariables(c("support", "Progress", "renderUI", "updateTextInput", "updateSelectInput", "shape", "x", "y"))

dataset_def_names <- setNames(c("Mouse Cell Atlas ", "Tabula Muris - 10X ", "Tabula Muris - FACS ", "Mouse Brain Atlas ", "Malaria Cell Atlas ", "Human Liver Atlas ", "Mouse Spinal Cord Atlas ", "Mouse sci-ATAC-seq Atlas ", "Human Kidney Atlas "), c("mca", "tm-10X", "tm-facs", "brain", "malaria", "liver", "spinalcord", "atacseq", "kidney"))

# To prevent query optimization runs automatically when the query is incomplete
js <- '
$(document).on("keyup", function(e) {
  if(e.keyCode == 13){
    Shiny.onInputChange("keyPressed", Math.random());
  }
});
'
TxtInput <- reactiveVal()
last.query.state <- reactiveVal("genelist")
initial.datasets <- "initial"
gene.list <- reactiveVal(c())
plotCellNumberHisto <- function(gene.list, dataset, object)
{
    
    if(!is.null(dataset) && length(gene.list) != 0){
        df <- as.data.frame(object@index$genesSupport(gene.list, dataset))
        dimnames(df)[[2]] <- 'support'
        df$`genes` <- rownames(df)
        max.axis <- max(df$support)
        return(ggplot(df, aes(x = genes, y = support), colour = support) +
                   ylab("Cells") + 
                   ylim(0, max.axis) +
                   geom_text(aes(label=df$`genes`), vjust = "inward", hjust = "inward", size = 4) +
                   geom_col(fill = rainbow(length(df$`genes`)), position = "dodge", width = 1, alpha = .5, colour = "black") + 
                   theme_minimal() + 
                   coord_flip() +
                   theme(panel.grid.minor = element_blank(), axis.text.y = element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black")) +
                   scale_x_discrete(limits = rev(df$`genes`), breaks = NULL))
        
    } else {
        return(plot(0,type='n',axes=FALSE,ann=FALSE))
    }
    
}

# Reset timeout to prevent shiny app quit when idle, 12 hour.
inactivity <- "function idleTimer() {
  var t = setTimeout(logout, 4.32e7);
  window.onmousemove = resetTimer; // catches mouse movements
  window.onmousedown = resetTimer; // catches mouse movements
  window.onclick = resetTimer;     // catches mouse clicks
  window.onscroll = resetTimer;    // catches scrolling
  window.onkeypress = resetTimer;  //catches keyboard actions

  function logout() {
    window.close();  //close the window
  }

  function resetTimer() {
    clearTimeout(t);
    t = setTimeout(logout, 4.32e7);  // time is in milliseconds (1000 is 1 second)
  }
}
idleTimer();"
