boxplotmcda <- function(datamatrix, labels=NULL, colors=NULL,  ax=NULL, mt=NULL) {
    if (is.null(ax)) {
        ax <- dev.cur()
    }
    
    datamatrix <- as.matrix(datamatrix)
    
    if (is.null(colors)) {
       colors <- rainbow(ncol(datamatrix))
    }
    
    do.call("boxplot", c(list(datamatrix, col=colors, xaxt="n", main=mt)))
    if (is.null(labels)) {
       labels <- paste0("C", seq_len(ncol(datamatrix)))
    }
    
    axis(1, at=1:ncol(datamatrix), labels=labels, las=2)
    grid(lty="dotted", col="grey")
    invisible(ax)
}
