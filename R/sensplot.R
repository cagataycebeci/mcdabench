sensplot <- function(senstable, topn = 10, 
   colpal = NULL, type = "bar", mtitle = NULL) {
   
   tab <- senstable[order(-senstable$Percent), ]
   tab$Pattern <- factor(tab$Pattern, levels = tab$Pattern)

    if (nrow(tab) > topn) {
     tab <- tab[1:topn, ]
   }

   if (is.null(colpal)) {
     colpal <- rainbow(nrow(tab))
   }
   
   if (is.null(mtitle)) {
     mtitle <- "Sensitivity Analysis"
   }

   if (type == "bar") {
     ymax <- max(tab$Percent) * 1.15
     bp <- barplot(tab$Percent,
                   names.arg = tab$Pattern,
                   las = 2,
                   col = colpal,
                   main = mtitle,
                   ylab = "Percentage (%)",
                   cex.names = 0.8,
                   ylim = c(0, ymax))
     text(x = bp, y = tab$Percent, labels = paste0(tab$Percent, "%"), 
          pos = 3, cex = 0.7, col = "black")   
   } else if (type == "pie") {
      pie(tab$Percent,
         labels = paste0(tab$Pattern, " (", tab$Percent, "%)"),
         col = colpal,
         main = mtitle)
   } else {
     warning("Invalid plot type. Please use 'bar' or 'pie'.")
   }
}
