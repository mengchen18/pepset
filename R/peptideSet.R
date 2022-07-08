#' Get peptides set
#' @param files a named character vector indicates the path of input fasta files. 
#'  The names of the vector will be used in the figure. If the vector does not 
#'  have name, name will be created from the file name. 
#' @param mc.cores number of cores used to digest peptides
#' @param subset either a character vector of peptides or a \code{data.frame} of two columns. 
#'   if a character vector is given, then the digested peptides will only be kept
#'   if the peptides are present in the input vector. If a \code{data.frame} is given,
#'   the first column should be a vector vector (see above) and the second column 
#'   should be a numeric vector, which usually indicates the intensity of peptides
#'   in an experiment. Then, the statistic is calculated not only from the presence 
#'   or absence of peptide, also the intensity information will be used. 
#' @param ... other parameters passed to \code{trypsinDigestList}
#' @importFrom parallel mclapply
#' @importFrom seqinr read.fasta
#' @import stringr
#' @import reshape2
#' @import philentropy
#' @rawNamespace import(plotly, except = last_plot)
#' @import ggplot2
#' @import ggdendro
#' @import DT
#' @rawNamespace import(shiny, except = c(dataTableOutput,renderDataTable))
#' @importFrom stats as.dist hclust na.omit


peptideSets <- function(files, mc.cores = 1, subset = NULL, ...) {
  
  if (!all(i <- file.exists(files))) 
    stop("Following file(s) not exist:", files[!i])
  
  weight <- NULL
  if (is.null(subset)) {
    weight <- 1
  } else {
    if (is.character(subset)) {
      weight <- 1
    } else if (is.data.frame(subset)) {
      if (ncol(subset) == 2) {
        if (is.character(subset[[1]]) && is.numeric(subset[[2]])) {
          weight <- tapply(subset[[2]], subset[[1]], sum, na.rm = TRUE)
          subset <- subset[[1]]
        }
      }
    } 
  } 
  if (is.null(weight))
    stop(
      "subset should be either a character vector or a data.frame of two columns 
         (first column is a character vector and second column is a numeric vector)!"
    )
  
  if (is.null(names(files)))
    names(files) <- sub(".fasta$|.fa$", "", basename(files))
  
  label <- names(files)
  
  files <- lapply(files, read.fasta, seqtype = "AA", as.string = TRUE)
  
  cat("Digesting protein ... \n")
  fseq <- lapply(files, trypsinDigestList, mc.cores = mc.cores)
  names(fseq) <- label
  cat("Digesting protein done!")
  
  fseq <- lapply(names(fseq), function(x) {
    d <- fseq[[x]]
    d$fasta <- x
    unique(d)
  })
  
  seqmat <- do.call(rbind, fseq)
  
  as <- unique(seqmat$seq)
  binmat <- sapply(fseq, function(x) {
    as %in% x$seq
  })
  rownames(binmat) <- as
  colnames(binmat) <- label
  
  we <- 1
  if (is.character(subset)) {
    if ( length(weight) == 1 )
      we <- as.integer(as %in% subset) else {
       we <- weight[as]
       we[is.na(we)] <- 0
      }
    } 
  
  dff <- data.frame(
    seq = as,
    length = nchar(as),
    n = rowSums(binmat), 
    weight = we, 
    binmat,
    stringsAsFactors = FALSE
  )
  
  list(seqmat = seqmat, binmat = dff)
}

#' Stats peptides set
#' @param bmat binary matrix returned by \code{peptideSet} function
#' @param filter whether to filter out peptide whose weight < 0
#' @importFrom philentropy distance
peptideStats <- function(bmat, filter = FALSE) {
  
  if (filter)
    bmat <- bmat[bmat$weight > 0, ]
  
  # clustering
  x <- as.matrix(bmat[, -(1:4)])
  if (ncol(x) > 2) {
    d <- distance(t(x), method = "jaccard")
    rownames(d) <- colnames(d) <- colnames(x)
    d <- as.dist(d)
    h <- hclust(d)
  } else 
    h <- list(order = seq_len(ncol(x)))
  
  # barplot obj
  an <- as.character(sort(unique(bmat$n)))
  k <- sapply(colnames(bmat)[-(1:4)], function(i) {
    i <- bmat[, i]
    tapply(bmat$weight[i], bmat$n[i], sum)[an]
  })
  k[is.na(k)] <- 0
  k <- melt(k)
  k$Var2 <- factor(k$Var2, levels = colnames(bmat)[-(1:4)][h$order])
  k$Var1 <- factor(k$Var1, levels = rev(unique(k$Var1)))
  colnames(k) <- c("Find_in", "fasta", "N")
  
  list(hcl = h, bar = k)
}