# Heatmap of ranks
rankheatmap <- function(rankmat, colpal=NULL, yl=NULL, xl=NULL, 
  htitle=NULL,  tcol="white", dendro="row", cellnotes=FALSE){
  
  if(!is.matrix(rankmat)) rankmat <- as.matrix(rankmat)
  rankmat <- rankmat[, !apply(rankmat, 2, function(col) any(is.nan(col)))] 
  
  p <- ncol(rankmat)
  n <- nrow(rankmat)
  
  if (n < 2 || p < 2) {
      stop("Invalid matrix. At least 2 rows and 2 columns are required.")
  }

  colpal1 <- colorRampPalette(c("red","yellow","springgreen","dodgerblue"))
  colpal2 <- colorRampPalette(c("yellow", "green","dodgerblue"))
  colpal3 <- colorRampPalette(c("skyblue1","deeppink4"))
  colpal4 <- monochromeR::generate_palette("darkgreen", modification = "go_lighter",
      n_colours = ncol(rankmat), view_palette = FALSE, view_labels = TRUE)
  colpal5 <- monochromeR::generate_palette("darkgreen", modification = "go_lighter",
      n_colours = ncol(rankmat), view_palette = FALSE, view_labels = TRUE)
  colpal6 <- monochromeR::generate_palette(c(255, 0, 0), modification = "go_darker",
     n_colours = ncol(rankmat), view_palette = FALSE, view_labels = TRUE)
  colpald <- monochromeR::generate_palette("skyblue", modification = "go_darker",
       n_colours = ncol(rankmat), view_palette = FALSE, view_labels = TRUE)
       
  if(is.null(colpal)){
    colpal <- colpal4
  }
  
  if(is.numeric(colpal)){
    if(colpal==1){
      colpal <- colpal1
    }else if(colpal==2){
      colpal <- colpal2
    }else if(colpal==3){
      colpal <- colpal3
    }else if(colpal==4){
      colpal <- colpal4     
    }else if(colpal==5){
      colpal <- colpal5
    }else if(colpal==6){
      colpal <- colpal6
    }else{
      colpal <- colpald
    }
  }
  
  if(is.null(xl)) xl <- ""
  if(is.null(yl)) yl <- ""
  if(is.null(htitle)) htitle <- ""
  
  if(cellnotes){
    cn <- round(rankmat, 3) 
  } else{
    cn <- matrix(NA, nrow=n, ncol=p)
  }

  gplots::heatmap.2(rankmat, dendrogram=dendro, col=colpal, 
    key=TRUE, keysize=0.75, key.par = list(cex=0.5), 
    trace="none", density.info="none", 
    cellnote=cn, notecol=tcol, notecex=1.1,
    colRow=rep(1, n), colCol=rep(1, p),
    cexRow=0.8, cexCol=0.9,
    ylab=yl, xlab=xl, main=htitle
  )
}