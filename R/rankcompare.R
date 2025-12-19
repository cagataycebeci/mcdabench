rankcompare <- function(rankmat, nperms = 1000, nboot=1000, entropyopt = "jsd", 
   j=1, k=3, alpha = 0.05, padjmethod = "fdr", biplot=TRUE){

res1 <- rankspearman(rankmat)$displaymat   
res2 <- rankentperm(rankmat, nperms = nperms, entropyopt = entropyopt, padjmethod=padjmethod)$displaymat
res3 <- rankentboot(rankmat, nboot=nboot, entropyopt = entropyopt, padjmethod=padjmethod)$displaymat
res4 <- rankwssim(rankmat)$displaymat
res5 <- rankmia(rankmat, j = j)$displaymat
res6 <- rankwilcox(rankmat, alpha = alpha, padjmethod=padjmethod)$displaymat
res7 <- rankrangesim(rankmat, k=k, plotsim=FALSE)$percentage
res8 <- rankgamma(rankmat)$gammamat
res9 <- rankpca(rankmat, biplot = biplot)
results <- list(src=res1, entperm=res2, entboot=res3, wsrs=res4, mia=res5, wilcox=res6, rangesim=res7, gamma=res8, geladipca=res9)
return(results)
}
