knitr::opts_chunk$set(fig.width=11, fig.height=9)

# install.packages("mcdabench", dep=TRUE)

library("mcdabench")

 data(egrids)
 # Extract the decision matrix, benefit-cost vector and weights

 dmatrix <- egrids$dmat
 bcvec <- egrids$bcvec
 weights <- egrids$weights
 egrids

nmatrix1 <- calcnormal(dmatrix, bcvec=bcvec, type="maxmin")
nmatrix2 <- calcnormal(dmatrix, bcvec=bcvec, type="sum")
nmatrix3 <- calcnormal(dmatrix, bcvec=bcvec, type="vector")
nmatrix4 <- calcnormal(dmatrix, bcvec=bcvec, type="zavadskas")
nmatrix5 <- calcnormal(dmatrix, bcvec=bcvec, type="ratio")
nmatrix6 <- calcnormal(dmatrix, bcvec=bcvec, type="enacc")

opar <- par(mfrow=c(3,2))
boxplotmcda(nmatrix1, mt="MaxMin")
boxplotmcda(nmatrix2, mt="Sum")
boxplotmcda(nmatrix3, mt="Vector")
boxplotmcda(nmatrix4, mt="Zavadskas")
boxplotmcda(nmatrix5, mt="Ratio")
boxplotmcda(nmatrix6, mt="Enhanced Accuracy")
par(opar)

nmatrix <- calcnormal(dmatrix, bcvec, type="vector")
critwei <- calcweights(nmatrix, bcvec=bcvec, type="critic")
entwei <- calcweights(nmatrix, bcvec=bcvec, type="entropy")
equwei <- calcweights(nmatrix, bcvec=bcvec, type="equal")
giniwei <- calcweights(nmatrix, bcvec=bcvec, type="gini")
sdevwei <- calcweights(nmatrix, bcvec=bcvec, type="sdev")
merecwei <- calcweights(nmatrix, bcvec=bcvec, type="merec")
mpsiwei <- calcweights(nmatrix, bcvec=bcvec, type="mpsi")
geomwei <- calcweights(nmatrix, bcvec=bcvec, type="geom")
rocwei <- calcweights(nmatrix, bcvec=bcvec, type="roc")
rswei <- calcweights(nmatrix, bcvec=bcvec, type="rs")

wmatrix <- cbind(Equal=equwei, Merec=merecwei, Geometric=geomwei, Gini=giniwei, 
           Critic=critwei, Mpsi=mpsiwei, StdDev=sdevwei, Rs=rswei, Roc=rocwei, Entropy=entwei)
print(round(wmatrix, 3))

corplot(wmatrix, colpal=c("orange", "dodgerblue","red"), xlab="Criterion", ylab="Weighting Method", title="Weights by Criteria", shape="square")
parcorplot(wmatrix, xl="Weighting Methods", yl="Weight", lt="Criteria")

resmegan_1 <- megan(dmatrix, bcvec=bcvec, weights="gini", normethod="maxmin", thr=0, tht=NULL)
str(resmegan_1)

# Superiority scores
 superiority <- resmegan_1$superiority
 print(superiority)
# Ranking of Alternatives
 rank <- resmegan_1$rank
 print(rank)

ratepar <- seq(0.01, 0.5, 0.05)
giniwei <- calcweights(nmatrix, bcvec=bcvec, type="gini")
resmegan_2 <- megan2(dmatrix, bcvec=bcvec, weights=giniwei, normethod="maxmin", 
  weimethod="gradual", rp=ratepar, thr=0)
  
 str(resmegan_2)
 rank <- resmegan_2$rank
 print(rank)

rankt0 <- megan(dmatrix, bcvec=bcvec, weights=giniwei, thr=0)$rank
rankt10 <- megan(dmatrix, bcvec=bcvec, weights=giniwei, thr=10)$rank
rankp5 <- megan(dmatrix, bcvec=bcvec, weights=giniwei, thr=NULL, tht="p5")$rank
rankp25 <- megan(dmatrix, bcvec=bcvec, weights=giniwei, thr=NULL, tht="p25")$rank
ranksdev <- megan(dmatrix, bcvec=bcvec, weights=giniwei, thr=NULL, tht="sdev")$rank

meganranks_2 <- cbind(T0=rankt0, T10=rankt10, P5=rankp5, P25=rankp25, SD=ranksdev)
meganranks_2 <- t(meganranks_2)
print(meganranks_2)

rankheatmap(meganranks_2, colpal=1, cellnotes=TRUE, tcol="black", dendro="row")

thr0 <- megan(dmatrix, bcvec=bcvec, weights=giniwei, thr=0)$thresholdval
thr10 <- megan(dmatrix, bcvec=bcvec, weights=giniwei, thr=10)$thresholdval
thrp5 <- megan(dmatrix, bcvec=bcvec, weights=giniwei, thr=NULL, tht="p5")$thresholdval
thrp25 <- megan(dmatrix, bcvec=bcvec, weights=giniwei, thr=NULL, tht="p25")$thresholdval
thrsdev <- megan(dmatrix, bcvec=bcvec, weights=giniwei, thr=NULL, tht="sdev")$thresholdval
thrp5 <- unname(thrp5); thrp25 <- unname(thrp25)

meganthresholds <- c(T0=thr0, T10=thr10, P5=thrp5, P25=thrp25, SD=thrsdev)
print(meganthresholds)

thrval <- 0L
thtmet <- "p5"
critrank <- megan(dmatrix, bcvec=bcvec, weights=critwei, thr=thrval, tht=thtmet)$rank
entrank <- megan(dmatrix, bcvec=bcvec, weights=entwei, thr=thrval, tht=thtmet)$rank
equrank <- megan(dmatrix, bcvec=bcvec, weights=equwei, thr=thrval, tht=thtmet)$rank
ginirank <- megan(dmatrix, bcvec=bcvec, weights=giniwei, thr=thrval, tht=thtmet)$rank
mpsirank <- megan(dmatrix, bcvec=bcvec, weights=mpsiwei, thr=thrval, tht=thtmet)$rank
sdevrank <- megan(dmatrix, bcvec=bcvec, weights=sdevwei, thr=thrval, tht=thtmet)$rank
geomrank <- megan(dmatrix, bcvec=bcvec, weights=geomwei, thr=thrval, tht=thtmet)$rank
merecrank <- megan(dmatrix, bcvec=bcvec, weights=merecwei, thr=thrval, tht=thtmet)$rank

meganranks_3 <- cbind(Equal=equrank, Critic=critrank, StdDev=sdevrank, Gini=ginirank,
  Geometric=geomrank, Merec=merecrank, Mpsi=mpsirank, Entropy=entrank)
rownames(meganranks_3) <- rownames(dmatrix)
print(meganranks_3)

respref <- rankaggregate(t(meganranks_3), topk=1)

print(respref$preference_table)  # Preference table
print(respref$preference_ranking) # Preference ranking

rankheatmap(meganranks_3, colpal=1, cellnotes=TRUE, tcol="black", dendro="column")

a <- rankpca(meganranks_3, biplot=TRUE, reverse=FALSE)

thrval <- NULL
thtmet <- "p5"
enhrank <- megan(dmatrix, bcvec=bcvec, weights=giniwei, normethod="enhanced", thr=thrval, tht=thtmet)$rank
ratiorank <- megan(dmatrix, bcvec=bcvec, weights=giniwei, normethod="ratio", thr=thrval, tht=thtmet)$rank
maxminrank <- megan(dmatrix, bcvec=bcvec, weights=giniwei, normethod="maxmin", thr=thrval, tht=thtmet)$rank
sumrank <- megan(dmatrix, bcvec=bcvec, weights=giniwei, normethod="sum", thr=thrval, tht=thtmet)$rank
vecrank <- megan(dmatrix, bcvec=bcvec, weights=giniwei, normethod="vector", thr=thrval, tht=thtmet)$rank
zavrank <- megan(dmatrix, bcvec=bcvec, weights=giniwei, normethod="zavadskas", thr=thrval, tht=thtmet)$rank
meganranks_4 <- cbind(Enhanced=enhrank, MaxMin=maxminrank, 
   Sum=sumrank, Ratio=ratiorank, Vector=vecrank, Zavadskas=zavrank)
rownames(meganranks_4) <- rownames(dmatrix)
print(meganranks_4)

rankheatmap(meganranks_4, colpal=1, cellnotes=TRUE, tcol="black", dendro="both")

a <- rankpca(meganranks_4, biplot=TRUE, reverse=FALSE)

respref <- rankaggregate(t(meganranks_4), topk=1)
print(respref$preference_table)  # Outranking table
print(respref$preference_ranking)  # Outranking flow

mp <- list(thr=0)
wp <- list(rp = c(0.01, 0.05, 0.10, -0.01, -0.05, -0.10))
megangrawei <- weisana(dmatrix = dmatrix, bcvec = bcvec, weights=critwei,
     weimethod = "gradual", weipars = wp,
     mcdamethod = megan, methodpars = mp,
     sensplot=FALSE)
rankmat <- megangrawei$ranking_matrix    # Ranking matrix
senstable <- megangrawei$sensitivity_table # Different rankings summary table
head(rankmat)  # Rank matrix
print(senstable)


rankheatmap(rankmat, colpal=4, cellnotes=FALSE, tcol="black", dendro="both")
a <- rankpca(rankmat, biplot=TRUE, reverse=FALSE)
