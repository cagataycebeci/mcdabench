parcorplot <- function(x, xl=NULL, yl=NULL, lt=NULL, colpal=NULL) {
  if(is.null(xl)) xl <- "x"
  if(is.null(yl)) xl <- "y"
  if(is.null(lt)) xl <- "x"

  n <- nrow(x)
  m <- ncol(x)
  
  par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE) 
  plot(1, type = "n",  xlab=xl, ylab=yl,
       xlim = c(1, m), ylim = range(x), xaxt = "n")
  axis(1, at = 1:m, labels = colnames(x))

  precolors <-  c("#FF0000", "#00FF00", "#0000FF", "pink", "darkgray",
       "#6A5ACD", "#FFA500", "#FF22CC", "#0000CD", "#DEB887", 
       "#2F4F4F", "#556B2F", "#8B4513", "#6B8E23", "#800000",
       "#A0FFAC", "#006400", "#708090", "#483D8B", "#3CB371", 
       "#FFFF00", "#008080", "#B8860B", "#4682B4", "#D2691E", 
       "#CD5C5C", "#00008B", "#32CD32", "#7F007F", "#8FBC8F",
       "#B03060", "#9932CC", "#FF0000", "#FF8C00", "#FFD700",
       "#00FF7F", "#E9967A", "#DC143C", "#00FFFF", "#00BFFF",
       "#BC8F8F", "#0F0F0F", "#A020F0", "#ADFF2F", "#FF6347", 
       "#FF00FF", "#F0E68C", "#6495ED", "#DDA0DD", "#90EE90",
       "#87CEEB", "#FF1493", "#7FFFD4", "#FF69B4", "#9ACD32",
       "#DA70D6")
  if(is.null(colpal)){
     if(n > length(precolors)){
       colpal <- sample(topo.colors(n))
     } else {
       colpal <- precolors[1:n]
     }
   }
   
  for (i in 1:n) {
    lines(1:m, as.numeric(x[i, ]), col = colpal[i], lwd = 3)
  }
    
  legend("topright", 
         inset = c(-0.2, 0),
         legend = row.names(x), 
         col = colpal, 
         lwd = 3, 
         ncol = 1, 
         bty = "n", 
         title = lt)
}
