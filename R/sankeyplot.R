sankeyplot <- function(rankmat, fs=12, nw=30){
   methods <- rownames(rankmat)
   alternatives <- colnames(rankmat)
   k <- length(methods)
   
   longdf <- data.frame(
     Method = rep(methods, each = length(alternatives)),
     Alternative = rep(alternatives, times = length(methods)),
     Rank = as.vector(t(rankmat)),
     stringsAsFactors = FALSE
   )

   nodes <- data.frame(name = c(methods, alternatives), stringsAsFactors = FALSE)

   longdf$source <- match(longdf$Method, nodes$name) - 1
   longdf$target <- match(longdf$Alternative, nodes$name) - 1

   longdf$value <- k - longdf$Rank 

   networkD3::sankeyNetwork(Links = longdf[, c("source", "target", "value")],
      Nodes = nodes,
      Source = "source",
      Target = "target",
      Value = "value",
      NodeID = "name",
      fontSize = fs,
      nodeWidth = nw)
 }



