methodbench <- function(dmatrix, bcvec, weights, normethod = "maxmin",
   mcdm      = "all", params = NULL, tiesmethod = "average") {

  built_in <- c("aras","aroman","cocoso","codas","copras","edas", 
                 "elect1", "elect2", "elect3", "elect4",
                "fuca","gra","mabac","macont6","marcos","mairca","maut","mavt",
                "megan","megan2","moora","ocra","oretes",
                "promt1","promt2","promt3","promt4","promt5","promt6",
                "ram","rov","smart","topsis","vikor","waspas","wpm","wsm")

  alias <- c(
      aromn="aroman", 
      macont="macont6", 
      mega="megan",
      elect1="electre1",
      elect2="electre2",
      elect3="electre3",
      elect4="electre4",
      promt1="promethee1", promet1="promethee1", 
      promt2="promethee2", promet2="promethee2", 
      promt3="promethee3", promet3="promethee3",
      promt4="promethee4", promet4="promethee4", 
      promt5="promethee5", promet5="promethee5", 
      promt6="promethee6", promet6="promethee6"
  )

  if (is.null(params)) {
    params <- list(
      aras      = list(),
      aroman    = list(lambda = 0.5, beta = 0.5),
      cocoso    = list(lambda = 0.5),
      codas     = list(thr = 0.1),
      copras    = list(),
      electre1  = list(),
      electre2  = list(),
      electre3  = list(),
      electre4  = list(),
      fuca      = list(),
      gra       = list(idesol = NULL, grdmethod = "sum", rho = 0.5),
      mabac     = list(),
      macont6   = list(p = 0.5, q = 0.5, delta = 0.5, theta = 0.5),
      marcos    = list(),
      mairca    = list(),
      maut      = list(utilfuncs = NULL, normutil = TRUE, ss = 1),
      mavt      = list(valfuncs = NULL, normvals = TRUE, ss = 1),
      megan     = list(normethod = "ratio", thr = 0, tht = "sdev"),
      megan2    = list(normethod = "ratio", thr = 0, tht = "sdev"),
      moora     = list(),
      ocra      = list(),
      oretes    = list(domplot = FALSE),
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
      wpm       = list(normethod = "vector"),
      wsm       = list()
    )
  }

  pick <- tolower(if (length(mcdm) == 1 && mcdm == "all") built_in else mcdm)
  methods_run <- unique(c(intersect(pick, built_in), setdiff(pick, built_in)))

  is_fun <- function(f) exists(f, mode = "function")

  trim <- function(fun, a) {
    fm <- names(formals(fun))
    if ("..." %in% fm) a else a[names(a) %in% fm]
  }
  mergeArgs <- function(core, extra) {
    if (is.null(extra)) extra <- list()
    modifyList(core, extra, keep.null = TRUE)
  }
  grab_rank <- function(res) {
    if (is.list(res))
      for (tg in c("rank","ranking","Rank","R","rank_score"))
        if (!is.null(res[[tg]])) return(res[[tg]])
    if (is.atomic(res)) return(res)
    NULL
  }

  results <- list()
  failed  <- character()
  n_alt   <- nrow(dmatrix)
  rnames  <- rownames(dmatrix)

  methods_run <- sort(methods_run)
  
  for (m in methods_run) {

    fun <- unname(if (m %in% names(alias)) alias[m] else m)
    if (!is_fun(fun)) { failed <- c(failed, m); next }

    base <- list(dmatrix = dmatrix, bcvec = bcvec, weights = weights,
            normethod = normethod, tiesmethod = tiesmethod)
            
    call <- trim(fun, mergeArgs(base, params[[fun]]))

    res <- tryCatch(do.call(fun, call),
           error = function(e) { failed <<- c(failed, m); NULL })
    if (is.null(res)) next

    rvec <- grab_rank(res)
    if (is.null(rvec) || length(rvec) != n_alt) { failed <- c(failed, m); next }

    if (is.data.frame(rvec)) rvec <- as.matrix(rvec)
    if (is.matrix(rvec))     rvec <- rvec[, 1]

    results[[toupper(m)]] <- as.numeric(rvec)
  }

  rankmat <- if (length(results)) {
               mat <- do.call(rbind, results)
               rownames(mat) <- names(results); colnames(mat) <- rnames
               mat
             } else matrix(numeric(0), nrow = 0, ncol = n_alt)

if (length(failed))
  message("The failed methods with problems or non-exist were excluded: ",
          paste(toupper(failed), collapse = ", "))
          
  list(dmatrix = dmatrix, weights = weights, rankmat = rankmat)
}
