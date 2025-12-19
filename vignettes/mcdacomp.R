options("width" = 260)
knitr::opts_chunk$set(fig.width=11, fig.height=9)

# install.packages("mcdabench", dep=TRUE)

library(mcdabench)

 # Load the data set
 data(egrids)

 # Extract the decision matrix, benefit-cost vector and weights
 dmat <- egrids$dmat
 bc <- egrids$bcvec
 userwei <- egrids$weights
 print(egrids)

nmatrix1 <- calcnormal(dmat, bcvec=bc, type="maxmin")
nmatrix2 <- calcnormal(dmat, bcvec=bc, type="sum")
nmatrix3 <- calcnormal(dmat, bcvec=bc, type="vector")
nmatrix4 <- calcnormal(dmat, bcvec=bc, type="zavadskas")

round(nmatrix1, 3) # MaxMin normalized matrix
corplot(nmatrix1, xlab="Alternative", ylab="Criterion", title="MaxMin Normalized Matrix")

opar <- par(mfrow=c(2,2))
boxplotmcda(nmatrix1, mt = "MaxMin")
boxplotmcda(nmatrix2, mt = "Sum")
boxplotmcda(nmatrix3, mt = "Vector")
boxplotmcda(nmatrix4, mt = "Zavadskas-Turskis")
par(opar)


critwei <- calcweights(nmatrix3, bcvec=bc, type="critic")
entwei <- calcweights(nmatrix3, bcvec=bc, type="entropy")
equwei <- calcweights(nmatrix3, bcvec=bc, type="equal")
giniwei <- calcweights(nmatrix3, bcvec=bc, type="gini")
sdevwei <- calcweights(nmatrix3, bcvec=bc, type="sdev")
merecwei <- calcweights(nmatrix3, bcvec=bc, type="merec")
mpsiwei <- calcweights(nmatrix3, bcvec=bc, type="mpsi")
geomwei <- calcweights(nmatrix3, bcvec=bc, type="geom")
rocwei <- calcweights(nmatrix3, bcvec=bc, type="roc")
rswei <- calcweights(nmatrix3, bcvec=bc, type="rs")
wmatrix <- cbind(Equal=equwei, Merec=merecwei, Geometric=geomwei, Mpsi=mpsiwei,
           Gini=giniwei, Critic=critwei, Entropy=entwei, StdDev=sdevwei, Rs=rswei, Roc=rocwei)
print(round(wmatrix,3))

parcorplot(wmatrix, xl="Weighting Methods", yl="Weight", lt="Criteria")
corplot(wmatrix, xlab="Weighting Methods", ylab="Weight", title="Weights", 
   colpal=c("gray","dodgerblue", "orange"))

    paramlist <- list(
      aras      = list(),
      aroman    = list(lambda = 0.5, beta = 0.5),
      codas     = list(thr = 0.1),
      cocoso    = list(lambda = 0.5),
      electre4  = list(p = 0.6, q = 0.4, v = 0.1),
      fuca      = list(),
      gra       = list(idesol = NULL, grdmethod = "sum", rho = 0.5),
      mabac     = list(),
      macont6   = list(p = 0.5, q = 0.5, delta = 0.5, theta = 0.5),
      marcos    = list(),
      mairca    = list(),
      maut      = list(utilfuncs = NULL, normutil = TRUE, ss = 1),
      mavt      = list(valfuncs = NULL, normvals = TRUE, ss = 1),
      megan     = list(normethod = "maxmin", thr = 0, tht = "sdev"),
      megan2    = list(normethod = "ratio", thr = NULL, tht = "p25"),
      moora     = list(),
      ocra      = list(),
      oreste    = list(domplot = FALSE),
      promethee1 = list(),
      promethee2 = list(),
      promethee3 = list(strict = FALSE),
      promethee4 = list(alpha = 0.2),
      promethee5 = list(g = 0, l = 100),
      promethee6 = list(varmethod = "abs_sum"),
      ram       = list(normethod = "sum"),
      rov       = list(normethod = "maxmin"),
      smart     = list(),
      topsis    = list(normethod = "maxmin"),
      vikor     = list(normethod = "maxmin", v = 0.5),
      waspas    = list(normethod = "linear", v = 0.5),
      wpm       = list(normethod = "vector")
    )

methodlist <- c("aras", "edas", "elect4", "fuca", "gra", "mabac", "codas", "marcos", "megan", 
                "moora", "promt2", "smart",  "topsis", "vikor", "waspas")

equwei <- calcweights(dmat, bcvec = bc, type = "equal")

resmcda <- methodbench(dmatrix = dmat, bcvec = bc, weights = equwei, 
     mcdm = methodlist, params = paramlist)

str(resmcda)  # Structure of benchmarking object

rankmat <- resmcda$rankmat  # Ranking matrix
print(rankmat)

rankheatmap(rankmat, colpal=1, cellnotes=TRUE, tcol="black")

corplot(rankmat, xlab="MCDM Methods", ylab="Alternatives", title="MCDA Methods", colpal=c("gray","green","dodgerblue"))
parcorplot(rankmat, xl="Alternatives", yl="Ranks", lt="MCDA Methods")

rescomp <- rankcompare(rankmat, nperms = 100, nboot=100, entropyopt = "jsd", 
   alpha = 0.05, padjmethod = "fdr", biplot=FALSE)
print(rescomp$src) # Spearman rank correlations matrix
print(rescomp$wsrs) # WS similarity matrix
print(rescomp$rangesim) # Rank range similarity matrix
print(rescomp$wilcox) # Wilcoxon rank sum test matrix
print(rescomp$entper) # Rank entropy matrix with permutations
print(rescomp$entboot) # Rank entropy matrix with bootstrap
rankheatmap(rescomp$rangesim, colpal=1, cellnotes=TRUE, tcol="black")
sccmatrix <- rankspearman(rankmat)$cormat
corplot(sccmatrix, xlab="MCDM Methods", ylab="MCDM Methods", title="Spearman Correlation Matrix")

ressens <- sensana(rankmat)
print(ressens$stabtable) # Stability
print(ressens$sensscores) # Sensitivity score

respref <- rankaggregate(rankmat, topk=3)
print(respref$preference_ranking)
print(respref$preference_table)

fod <- respref$preference_ranking
flowplot(fod["BORDACNT",], colpal=terrain.colors(12), txtcol = "black", orientation = "vertical") 
