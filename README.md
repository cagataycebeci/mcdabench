---
title: "mcdabench: Benchmarking for Multi-Criteria Decision Analysis"
author: "Cagatay Cebeci"
date: "2025-06-22"
#output: html_document
---

# Introduction

The `mcdabench` package provides implementations of various Multi-Criteria Decision Analysis (MCDA) methods in R. These methods are designed to evaluate and rank alternatives based on 
multiple criteria, applying normalization, weighting, and aggregation techniques. The package includes popular decision-making methods such as  ARAS, AROMAN, COCOSO, CODAS, COPRAS, EDAS, 
ELECTRE family (I-IV), FUCA, GRA, MABAC, MAIRCA, MARCOS, MAUT, MAVT, MEGAN, MOORA, OCRA, ORETES, PROMETHEE family (I - VI), RAM, ROV, SMART, TOPSIS, VIKOR, WASPAS, WPM, WSM and many more, 
facilitating flexible and efficient analyses for multi-criteria problems.

# Install mcdabench

You can easily install the `mcdabench` package from CRAN:

```r
install.packages("mcdabench", dep=TRUE)
```

After installing the package, load it into your R session by running:
library(mcdabench)

```r
library(mcdabench)
```

# Quickly mcdabench
## Load data
```r
data(egrids)
dmat < egrids$dmat # Decision matrix
bc <- egrids$bcvec # Benefit-cost vector
uw <- egrids$weights # Criteria weights 

## Normalization
```r
nmatrix1 <- normalize(dmatrix, bcvec=bc, type="maxmin")
nmatrix1

nmatrix2 <- normalize(dmatrix, bcvec=bc, type="sum")
nmatrix2
```

## Calculate weights
```r
equwei <- calcweights(nmatrix1, bcvec=bc, type="equal")
equwei

entwei <- calcweights(nmatrix1, bcvec=bc, type="entropy")
entwei
```

## Apply MCDA Methods
### Run MEGAN
```r
resmegan <- megan(dmatrix, bcvec=bc, weights=uw, normethod="maxmin")
print(resmegan)$rank
```

### Run TOPSIS
```r
restopsis <- topsis(dmatrix, bcvec=bc, weights=uw, normethod="maxmin")
print(restopsis)
```

### Run VIKOR
```r
resvikor <- vikor(dmatrix, bcvec=bc, weights=uw, v=0.8)
print(resvikor$rank)
```

## Weight Sensitivity Analsyis
### Gradual weight modification for Testing VIKOR
```r
mp <- list(v=0.5)
wp <- list(rp = seq(0.01, 0.5, 0.05))

vikorgrawei <- weisana(dmatrix = dmat, bcvec = bc,
     weimethod = "gradual", weipars = wp,
     mcdamethod = vikor, methodpars = mp, sensplot=FALSE)
print(vikorgrawei)
sensplot(vikorgrawei$sensitivity_table, 
   mtitle="Weight Sensivity Analysis for VIKOR", colpal=terrain.colors(10))
```

### Test WASPAS with gradual and random weighting
```r
mp <- list(v=0.5, normethod="linear", tiesmethod="average")
wp <- list(rp = seq(0.01, 0.6, 0.01))

waspasgrawei <- weisana(dmatrix = dmat, bcvec = bc,
     weimethod = "gradual", weipars = wp,
     mcdamethod = waspas, methodpars = mp)
print(waspasgrawei)
sensplot(waspasgrawei$sensitivity_table, 
   mtitle="Weight Sensivity Analysis for WASPAS", colpal=terrain.colors(10))
   
waspasrandwei <- weisana(dmatrix = dmat, bcvec = bc,
     weimethod = "random", weipars = list(ss=0.05, niters=50),
     mcdamethod = waspas, methodpars = mp)
print(waspasrandwei)
sensplot(waspasrandwei$sensitivity_table, 
   mtitle="Weight Sensivity Analysis (random) for WASPAS", colpal=terrain.colors(10))
```

# Bechmarking the Methods
```r
# Sample decision matrix
dm <- matrix(c(
  10, 20, 30, 1.5, 102, 55,
  15, 25, 35, 1.6, 90, 60,
  12, 22, 32, 1.7, 100, 58,
  13, 24, 33, 1.8, 95, 57,
  14, 26, 37, 1.9, 98, 59,
  11, 23, 31, 1.65, 101, 56,
  16, 27, 36, 1.55, 97, 61,
  17, 28, 38, 1.7, 99, 63,
  18, 29, 39, 1.8, 94, 62,
  19, 30, 40, 1.75, 96, 64
), nrow = 10, byrow = TRUE)
colnames(dm) <- paste0("C", 1:ncol(dm))
rownames(dm) <- paste0("ALT", 1:nrow(dm))

# Benefit-Cost vector
bc <- c(1, -1, 1, -1, 1, 1)

# User-defined weights
userwei <- c(0.3, 0.1, 0.2, 0.1, 0.2, 0.1)

prmlist <- list(
      aras      = list(),
      aroman    = list(lambda = 0.5, beta = 0.5),
      cocoso    = list(lambda = 0.5),
      codas     = list(thr = 0.1),  
      smart     = list(),
      topsis    = list(normethod = "maxmin"),
      vikor     = list(normethod = "maxmin", v = 0.8),
      waspas    = list(normethod = "linear", v = 0.5),
      wpm       = list(normethod = "vector"),
      wsm       = list()
)

# Compare selected methods with 'gini' weights
giniwei <- calcweights(dm, bcvec=bc, type="gini")
rescomp_3 <- methodbench(dmatrix = dm, bcvec = bc, weights = giniwei,
   mcdm = methodlist, params=prmlist)
print(rescomp_3$rankmat)
rankheatmap(rescomp_3$rankmat, colpal=1, cellnotes=TRUE, tcol="white")

# Overall ranks and outranking
resoverall <- rankaggregate(rescomp_3$rankmat, tiesmethod="average", topk = 3,
   damp = 0.5, niters = 200, tol = 1e-4)
print(resoverall)
```

# Citation

To cite the mcdabench package in publications, please run the following code in R.

```r 
citation("mcdabench")
```

