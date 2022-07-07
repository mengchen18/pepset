#' check the folder, if no object, then process. If exist, then load directly
#' @param f a folder containing the fasta files. Need read and write permission
#' @param ncores n cores passed to mclapply
#' @return an object to be visualized
#' 
loadObj <- function(f, ncores) {
  if (file.exists(ff <- file.path(f, "zzz_pepset_sumstats.RDS"))) {
    ll <- readRDS(ff)
  } else {
    fs <- list.files(f, full.names = TRUE, pattern = ".fasta|.fa")
    ll <- peptideSets(fs, mc.cores = ncores)
    sta <- peptideStats(ll$binmat)
    ll <- c(ll, sta)
    saveRDS(ll,  file = file.path(f, "zzz_pepset_sumstats.RDS"))
  }
  ll
}


#' Format DT 
#' @param tab a table-like object, such as data.frame
#' 
formatTab <- function(tab) {    
  ci <- unname(which(vapply(tab, inherits, c('factor', "character"), FUN.VALUE = logical(1))))
  if (length(ci) > 0)
    tab[ci] <- lapply(tab[ci], function(x) {
      x[is.na(x)] <- ""
      x
    })
  dt <- DT::datatable( 
    tab,
    selection = "none", #list(mode = c("single", "multiple")[as.integer(sel)+1], selected = tab_status()$rows_selected, target = "row"),
    rownames = FALSE,
    filter = "top",
    class="table-bordered compact nowrap",
    options = list(
      scrollX = TRUE, pageLength = 10, dom = 'tip',
      columnDefs = list(list(
        targets = ci-1,
        render = DT::JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 50 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
          "}")
      )),
      stateSave = TRUE,  stateDuration = -1
    )
  )
  DT::formatStyle(dt, columns = seq_len(ncol(tab)), fontSize = '90%')
}

#' pop up window to choose fasta folder
#' @param id id
#' @param dir directory to monitor
#' @import shinyFiles
input_popup <- function(id, dir) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns
      rt <- structure(normalizePath(dir), names = "Current root folder")
      
      showModal(modalDialog(
        shinyDirButton(
          id = ns('fastaDir'), label = "Select folder containing the fasta files", 
          title = "Acceptable file/format: .fa/.fasta", multiple = FALSE),
        title = "Select a folder containing all the fasta files ...",
        footer = NULL,#modalButton("Open"),
        size = "m",
        easyClose = FALSE,
        fade = TRUE,
        style = "z-index: 9999"
      ))
      
      shinyDirChoose(
        input = input, 
        id = 'fastaDir', 
        roots = rt,
        defaultRoot = names(rt)[1], 
        session = session, 
        filetypes=c('', 'fasta', "fa")
      )
      
      observeEvent( input$fastaDir, {
        req(input$fastaDir)
        req(!inherits(input$fastaDir, "integer"))
        removeModal()
      })
      
      reactive({
        req(input$fastaDir)
        req(!inherits(input$fastaDir, "integer"))
        v <- do.call(file.path, c(
           rt[[input$fastaDir$root]], input$fastaDir$path
        ))
        normalizePath(v)
      })
    }
  )}
