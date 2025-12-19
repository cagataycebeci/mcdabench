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

nmat <- calcnormal(dmat, bc, type="vector")
critwei <- calcweights(dmatrix=nmat, bcvec=bc, type="critic")
entwei <- calcweights(dmatrix=nmat, bcvec=bc, type="entropy")
equwei <- calcweights(dmatrix=nmat, bcvec=bc, type="equal")
giniwei <- calcweights(dmatrix=nmat,bcvec=bc, type="gini")
sdevwei <- calcweights(dmatrix=nmat, bcvec=bc, type="sdev")
merecwei <- calcweights(dmatrix=nmat, bcvec=bc, type="merec")
mpsiwei <- calcweights(dmatrix=nmat, bcvec=bc, type="mpsi")
geomwei <- calcweights(dmatrix=nmat, bcvec=bc, type="geom")
rocwei <- calcweights(dmatrix=nmat, bcvec=bc, type="roc")
rswei <- calcweights(dmatrix=nmat, bcvec=bc, type="rs")
wmatrix <- cbind(Equal=equwei, Merec=merecwei, Geometric=geomwei, Gini=giniwei, 
    Critic=critwei, Mpsi=mpsiwei, Entropy=entwei, StdDev=sdevwei, Rs=rswei, Roc=rocwei )
print(round(wmatrix,3))

parcorplot(wmatrix, xl="Weighting Methods", yl="Weight", lt="Criteria")

critrank <- topsis(dmatrix=dmat, bcvec=bc, weights=critwei)$rank
entrank <- topsis(dmatrix=dmat, bcvec=bc, weights=entwei)$rank
equrank <- topsis(dmatrix=dmat, bcvec=bc, weights=equwei)$rank
ginirank <- topsis(dmatrix=dmat, bcvec=bc, weights=giniwei)$rank
sdevrank <- topsis(dmatrix=dmat, bcvec=bc, weights=sdevwei)$rank
geomrank <- topsis(dmatrix=dmat, bcvec=bc, weights=geomwei)$rank
merecrank <- topsis(dmatrix=dmat, bcvec=bc, weights=merecwei)$rank
mpsirank <- topsis(dmatrix=dmat, bcvec=bc, weights=mpsiwei)$rank
rocrank <- topsis(dmatrix=dmat, bcvec=bc, weights=rocwei)$rank
rsrank <- topsis(dmatrix=dmat, bcvec=bc, weights=rswei)$rank

topsisranks <- cbind(Equal=equrank, Critic=critrank, StdDev=sdevrank, Gini=ginirank, Geometric=geomrank, 
    Merec=merecrank, Mpsi=mpsirank, Entropy=entrank, Roc=rocrank, Rs=rsrank)

rownames(topsisranks) <- rownames(dmat)
print(t(topsisranks))

rankheatmap(t(topsisranks), colpal=1, cellnotes=TRUE, tcol="black", dendro="row")

pca <- rankpca(t(topsisranks), biplot=TRUE)

ressens1 <- sensana(t(topsisranks))
print(ressens1)

rescomp <- rankcompare(t(topsisranks), biplot=FALSE, nperms = 100, nboot=100, 
   entropyopt = "jsd", alpha = 0.05, padjmethod = "bonferroni")


print(rescomp$src) # Spearman rank correlations matrix

print(rescomp$wsrs) # WS similarity matrix

print(rescomp$wilcox) # Wilcox test matrix

print(rescomp$entper) # Rank entropy matrix with permutations

print(rescomp$entboot) # Rank entropy matrix with bootstrap

mp <- list()
wp <- list(rp = seq(0.01, 0.5, 0.05))
topsisgrawei <- weisana(dmatrix = dmat, bcvec = bc, weights=critwei,
     weimethod = "gradual", weipars = wp,
     mcdamethod = topsis, methodpars = mp)
str(topsisgrawei)

gradweimat <- topsisgrawei$weights_matrix # Modified weights matrix
rankmat <- topsisgrawei$ranking_matrix    # Ranking matrix
senstable <- topsisgrawei$sensitivity_table # Different rankings summary table
head(round(gradweimat, 3)) 
head(rankmat)  # Rank matrix

print(senstable)


pca <- rankpca(rankmat, biplot=TRUE)

ressens2 <- sensana(rankmat)
print(ressens2$stabtable)

respref <- rankaggregate(t(topsisranks), topk=1, tiesmethod="average")
preftable <- respref$preference_table
prefrank <- respref$preference_ranking
print(prefrank)
print(preftable)
