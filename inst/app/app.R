############# app setting #############
# number of cores used in tryptic digestion
ncores  <- 8
# root directory
root <- "~/projects"

##########################

library(pepset)
library(plotly)
library(ggplot2)
library(ggdendro)
library(shiny)
library(DT)
library(shinyFiles)

ui <- fluidPage(
  fluidRow(
    column(12, h3('PepSet - Checking fasta files in folder: ')),
    column(12, verbatimTextOutput("path")),
    column(2, plotlyOutput("p1"), offset = 0, style='padding-right:0px;'),
    column(6, plotlyOutput("p2"), offset = 0, style='padding-left:0px;'),
    column(4, plotlyOutput("p3")),
    column(12, hr()),
    column(6, wellPanel(DTOutput("tab1"))),
    column(6, wellPanel(DTOutput("tab2")))
  )
)

server <- function(input, output, session) {
  
  ############## folder - need R & W permission ###########
  fastaPath <- pepset:::input_popup(id = "pop0", dir = root)
  
  obj <- reactive({
    req(fastaPath())
    showModal(modalDialog(
      "If this is the first time you selec this folder, 
      it may take several minutes to in-silico digest
      all the fasta files and stats. 
      Please do not close this window.\n",
      "Next time when you open this folder, you do not need to reprocess the fasta files and
      will be faster. ",
      title = "Load data ...",
      footer = NULL,
      size = "m",
      easyClose = FALSE,
      fade = TRUE,
      style = "z-index: 9999"
    ))
    ll <- pepset:::loadObj(f = fastaPath(), ncores = ncores)
    removeModal()
    ################ dendrogram #############
    h <- ll$hcl
    dhc <- as.dendrogram(h)
    data <- dendro_data(dhc, type = "rectangle")
    fs <- length(unique(ll$bar$fasta))
    p1 <- ggplot(segment(data)) +
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
      coord_flip() +
      scale_y_reverse() +
      coord_flip(xlim = c(0.5, fs+0.5), ylim = c(1.01, 0), expand=FALSE) +
      theme_minimal() +
      theme(
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
      ) +
      labs(y = "Jaccard distance")
    
    ############ bar chart ############
    
    k <- ll$bar
    p2 <- ggplot(k, aes(fill=Find_in, y=N, x=fasta)) +
      geom_bar(stat="identity") +
      coord_flip(expand = 0) +
      theme_minimal() +
      labs(y = "Number of peptides") +
      theme(
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y = element_text(hjust=0)
      )
    
    ############ pie chart ############
    lx <- ll$binmat[which(ll$binmat$n == 1), -(1:4)]
    lx <- colSums(lx)
    data <- data.frame(
      fasta = c(names(lx), "Shared"),
      peptide_number = c(lx, sum(ll$binmat$n > 1, na.rm = TRUE))
    )
    p3 <- plotly::plot_ly(
      data=data,labels=~fasta, values=~peptide_number, type="pie"
    )
    c(ll, list(p1 = p1, p2 = p2, p3 = p3))
  })
  # 
  
  # # #### output
  
  output$path <- renderText(fastaPath())
  
  output$p1 <- renderPlotly(ggplotly(obj()$p1))
  output$p2 <- renderPlotly(ggplotly(obj()$p2))
  output$p3 <- renderPlotly(ggplotly(obj()$p3))
  
  output$tab1 <- renderDT({
    req(obj()$seqmat)
    pepset:::formatTab(obj()$seqmat)
  }, server = TRUE)
  output$tab2 <- renderDT({
    req(obj()$binmat)
    pepset:::formatTab (obj()$binmat)
  }, server = TRUE)
}

shinyApp(ui, server)


