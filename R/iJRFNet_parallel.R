#' Joint Random Forest for the simultaneous estimation of multiple related networks
#'
#' MAIN FUNCTION -- > iJRFNet


importance <- function(x,  scale=TRUE) {
  # --- Function importance is a modified version of function importance from R package randomForest
  
  type=NULL;
  class=NULL;
  if (!inherits(x, "randomForest"))
    stop("x is not of class randomForest")
  classRF <- x$type != "regression"
  hasImp <- !is.null(dim(x$importance)) || ncol(x$importance) == 1
  hasType <- !is.null(type)
  if (hasType && type == 1 && !hasImp)
    stop("That measure has not been computed")
  allImp <- is.null(type) && hasImp
  if (hasType) {
    if (!(type %in% 1:2)) stop("Wrong type specified")
    if (type == 2 && !is.null(class))
      stop("No class-specific measure for that type")
  }
  
  imp <- x$importance
  if (hasType && type == 2) {
    if (hasImp) imp <- imp[, ncol(imp), drop=FALSE]
  } else {
    if (scale) {
      SD <- x$importanceSD
      imp[, -ncol(imp)] <-
        imp[, -ncol(imp), drop=FALSE] /
        ifelse(SD < .Machine$double.eps, 1, SD)
    }
    if (!allImp) {
      if (is.null(class)) {
        ## The average decrease in accuracy measure:
        imp <- imp[, ncol(imp) - 1, drop=FALSE]
      } else {
        whichCol <- if (classRF) match(class, colnames(imp)) else 1
        if (is.na(whichCol)) stop(paste("Class", class, "not found."))
        imp <- imp[, whichCol, drop=FALSE]
      }
    }
  }
  imp<-imp[,2]
  imp
}


# --- Functions called by iJRFNet
"iJRF_onetarget" <-
  function(x, y=NULL,  xtest=NULL, ytest=NULL, ntree,
           sampsize,              
           totsize = if (replace) ncol(x) else ceiling(.632*ncol(x)),
           mtry=if (!is.null(y) && !is.factor(y))
             max(floor(nrow(x)/3), 1) else floor(sqrt(nrow(x))),
           replace=TRUE, classwt=NULL, cutoff, strata,
           nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
           maxnodes=NULL,
           importance=FALSE, localImp=FALSE, nPerm=1,
           proximity, oob.prox=proximity,
           norm.votes=TRUE, do.trace=FALSE,
           keep.forest=!is.null(y) && is.null(xtest), corr.bias=FALSE,
           keep.inbag=FALSE, nclasses, sw,...) {
    
    ww=1/sampsize;
    nclass=mylevels=ipi=NULL
    addclass <- is.null(y)
    classRF <- addclass || is.factor(y)
    if (!classRF && length(unique(y)) <= 5) {
      warning("The response has five or fewer unique values.  Are you sure you want to do regression?")
    }
    if (classRF && !addclass && length(unique(y)) < 2)
      stop("Need at least two classes to do classification.")
    
    n <- ncol(y)           # number of samples
    p <- nrow(x)/nclasses  # number of variables
    
    if (n == 0) stop("data (x) has 0 rows")
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)
    
    keep.forest=!is.null(y) 
    xtest=NULL; ytest=NULL
    testdat <- !is.null(xtest)
    if (testdat) {
      if (ncol(x) != ncol(xtest))
        stop("x and xtest must have same number of columns")
      ntest <- nrow(xtest)
      xts.row.names <- rownames(xtest)
    }
    
    prox <- proxts <- double(1)
    ## Check for NAs.
    if (any(is.na(x))) stop("NA not permitted in predictors")
    if (testdat && any(is.na(xtest))) stop("NA not permitted in xtest")
    if (any(is.na(y))) stop("NA not permitted in response")
    if (!is.null(ytest) && any(is.na(ytest))) stop("NA not permitted in ytest")
    
    if (is.data.frame(x)) {
      xlevels <- lapply(x, mylevels)
      ncat <- sapply(xlevels, length)
      ## Treat ordered factors as numerics.
      ncat <- ifelse(sapply(x, is.ordered), 1, ncat)
      x <- data.matrix(x)
      if(testdat) {
        if(!is.data.frame(xtest))
          stop("xtest must be data frame if x is")
        xfactor <- which(sapply(xtest, is.factor))
        if (length(xfactor) > 0) {
          for (i in xfactor) {
            if (any(! levels(xtest[[i]]) %in% xlevels[[i]]))
              stop("New factor levels in xtest not present in x")
            xtest[[i]] <-
              factor(xlevels[[i]][match(xtest[[i]], xlevels[[i]])],
                     levels=xlevels[[i]])
          }
        }
        xtest <- data.matrix(xtest)
      }
    } else {
      ncat <- rep(1, p)
      xlevels <- as.list(rep(0, p))
    }
    maxcat <- max(ncat)
    if (maxcat > 32)
      stop("Can not handle categorical predictors with more than 32 categories.")
    
    addclass <- FALSE
    
    proximity <- addclass
    
    
    
    impout <- matrix(0.0, p*nclasses, 2)
    impSD <- matrix(0.0, p*nclasses, 1)
    #  names(impSD) <- x.col.names
    
    
    
    nsample <- if (addclass) 2 * n else n
    Stratify <- length(n) > 1
    
    nodesize=5;
    nrnodes <- 2 * trunc(n/max(1, nodesize - 4)) + 1
    
    maxnodes=NULL
    if (!is.null(maxnodes)) {
      ## convert # of terminal nodes to total # of nodes
      maxnodes <- 2 * maxnodes - 1
      if (maxnodes > nrnodes) warning("maxnodes exceeds its max value.")
      nrnodes <- min(c(nrnodes, max(c(maxnodes, 1))))
    }
    
    
    ## Compiled code expects variables in rows and observations in columns.
    # x <- t(x)
    storage.mode(x) <- "double"
    
    xtest <- double(1)
    ytest <- double(1)
    ntest <- 1
    labelts <- FALSE
    nt <- if (keep.forest) ntree else 1
    
    
    nPerm=1
    do.trace=F; oob.prox=F
    corr.bias=FALSE
    keep.inbag=FALSE
    impmat <- double(1)
    replace=T
    
    
    rfout <- .C("iJRF_regRF",
                x,
                y, ww,
                as.integer(c(totsize, p)),
                sampsize=as.integer(sampsize), as.integer(totsize),
                as.integer(nodesize),
                as.integer(nrnodes),
                as.integer(ntree),
                as.integer(mtry),
                as.integer(c(importance, localImp, nPerm)),
                as.integer(ncat),
                as.integer(maxcat),
                as.integer(do.trace),
                as.integer(proximity),
                as.integer(oob.prox),
                as.integer(corr.bias),
                ypred = double(n * nclasses),
                impout = impout,
                impmat = impmat,
                impSD = impSD,
                prox = prox,
                ndbigtree = integer(ntree),
                nodestatus = matrix(integer(nrnodes * nt * nclasses), ncol=nt),
                leftDaughter = matrix(integer(nrnodes * nt * nclasses), ncol=nt),
                rightDaughter = matrix(integer(nrnodes * nt * nclasses), ncol=nt),
                nodepred = matrix(double(nrnodes * nt * nclasses), ncol=nt),
                bestvar = matrix(integer(nrnodes * nt * nclasses), ncol=nt),
                xbestsplit = matrix(double(nrnodes * nt * nclasses), ncol=nt),
                mse = double(ntree * nclasses),
                keep = as.integer(c(keep.forest, keep.inbag)),
                replace = as.integer(replace),
                testdat = as.integer(testdat),
                xts = xtest,
                ntest = as.integer(ntest),
                yts = as.double(ytest),
                labelts = as.integer(labelts),
                ytestpred = double(ntest),
                proxts = proxts,
                msets = double(if (labelts) ntree else 1),
                coef = double(2),
                oob.times = integer(n),
                inbag = if (keep.inbag)
                  matrix(integer(n * ntree), n) else integer(1), as.integer(nclasses), 
                sw = as.double(sw))[c(16:28, 36:41)]
    #         ## Format the forest component, if present.
    if (keep.forest) {
      max.nodes <- max(rfout$ndbigtree)
      rfout$nodestatus <-
        rfout$nodestatus[1:max.nodes, , drop=FALSE]
      rfout$bestvar <-
        rfout$bestvar[1:max.nodes, , drop=FALSE]
      rfout$nodepred <-
        rfout$nodepred[1:max.nodes, , drop=FALSE]
      rfout$xbestsplit <-
        rfout$xbestsplit[1:max.nodes, , drop=FALSE]
      rfout$leftDaughter <-
        rfout$leftDaughter[1:max.nodes, , drop=FALSE]
      rfout$rightDaughter <-
        rfout$rightDaughter[1:max.nodes, , drop=FALSE]
    }
    cl <- match.call()
    cl[[1]] <- as.name("randomForest")
    #         ## Make sure those obs. that have not been OOB get NA as prediction.
    ypred <- rfout$ypred
    if (any(rfout$oob.times < 1)) {
      ypred[rfout$oob.times == 0] <- NA
    }
    out <- list(call = cl,
                type = "regression",
                predicted =0,
                mse = rfout$mse,
                rsq = 1 - rfout$mse / (var(y[1,]) * (n-1) / n),
                oob.times = rfout$oob.times,
                importance = if (importance) matrix(rfout$impout, p * nclasses, 2) else
                  matrix(rfout$impout, ncol=1),
                importanceSD=if (importance) rfout$impSD else NULL,
                localImportance = if (localImp)
                  matrix(rfout$impmat, p, n, dimnames=list(x.col.names,
                                                           x.row.names)) else NULL,
                proximity = if (proximity) matrix(rfout$prox, n, n,
                                                  dimnames = list(x.row.names, x.row.names)) else NULL,
                ntree = ntree,
                mtry = mtry,
                forest = if (keep.forest)
                  c(rfout[c("ndbigtree", "nodestatus", "leftDaughter",
                            "rightDaughter", "nodepred", "bestvar",
                            "xbestsplit")],
                    list(ncat = ncat), list(nrnodes=max.nodes),
                    list(ntree=ntree), list(xlevels=xlevels)) else NULL,
                coefs = if (corr.bias) rfout$coef else NULL,
                y = y,
                test = if(testdat) {
                  list(predicted = structure(rfout$ytestpred,
                                             names=xts.row.names),
                       mse = if(labelts) rfout$msets else NULL,
                       rsq = if(labelts) 1 - rfout$msets /
                         (var(ytest) * (n-1) / n) else NULL,
                       proximity = if (proximity)
                         matrix(rfout$proxts / ntree, nrow = ntest,
                                dimnames = list(xts.row.names,
                                                c(xts.row.names,
                                                  x.row.names))) else NULL)
                } else NULL,
                inbag = if (keep.inbag)
                  matrix(rfout$inbag, nrow(rfout$inbag),
                         dimnames=list(x.row.names, NULL)) else NULL)
    
    
    #      print(rfout$mse)
    class(out) <- "randomForest"
    return(out)
    
  }

"irafnet_onetarget" <-
  function(x, y=NULL,  xtest=NULL, ytest=NULL, ntree,
           mtry=if (!is.null(y) && !is.factor(y))
             max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))),
           replace=TRUE, classwt=NULL, cutoff, strata,
           sampsize = if (replace) nrow(x) else ceiling(.632*nrow(x)),
           nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
           maxnodes=NULL,
           importance=FALSE, localImp=FALSE, nPerm=1,
           proximity, oob.prox=proximity,
           norm.votes=TRUE, do.trace=FALSE,
           keep.forest=!is.null(y) && is.null(xtest), corr.bias=FALSE,
           keep.inbag=FALSE, sw) {
    addclass <- is.null(y)
    classRF <- addclass || is.factor(y)
    if (!classRF && length(unique(y)) <= 5) {
      warning("The response has five or fewer unique values.  Are you sure you want to do regression?")
    }
    if (classRF && !addclass && length(unique(y)) < 2)
      stop("Need at least two classes to do classification.")
    n <- nrow(x)
    p <- ncol(x)
    if (n == 0) stop("data (x) has 0 rows")
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)
    
    ## overcome R's lazy evaluation:
    keep.forest <- keep.forest
    
    testdat <- !is.null(xtest)
    if (testdat) {
      if (ncol(x) != ncol(xtest))
        stop("x and xtest must have same number of columns")
      ntest <- nrow(xtest)
      xts.row.names <- rownames(xtest)
    }
    
    ## Make sure mtry is in reasonable range.
    if (mtry < 1 || mtry > p)
      warning("invalid mtry: reset to within valid range")
    mtry <- max(1, min(p, round(mtry)))
    if (!is.null(y)) {
      if (length(y) != n) stop("length of response must be the same as predictors")
      addclass <- FALSE
    } else {
      if (!addclass) addclass <- TRUE
      y <- factor(c(rep(1, n), rep(2, n)))
      x <- rbind(x, x)
    }
    
    ## Check for NAs.
    if (any(is.na(x))) stop("NA not permitted in predictors")
    if (testdat && any(is.na(xtest))) stop("NA not permitted in xtest")
    if (any(is.na(y))) stop("NA not permitted in response")
    if (!is.null(ytest) && any(is.na(ytest))) stop("NA not permitted in ytest")
    
    ncat <- rep(1, p)
    xlevels <- as.list(rep(0, p))
    maxcat <- max(ncat)
    if (maxcat > 32)
      stop("Can not handle categorical predictors with more than 32 categories.")
    
    if (classRF) {
      nclass <- length(levels(y))
      ## Check for empty classes:
      if (any(table(y) == 0)) stop("Can't have empty classes in y.")
      if (!is.null(ytest)) {
        if (!is.factor(ytest)) stop("ytest must be a factor")
        if (!all(levels(y) == levels(ytest)))
          stop("y and ytest must have the same levels")
      }
      if (missing(cutoff)) {
        cutoff <- rep(1 / nclass, nclass)
      } else {
        if (sum(cutoff) > 1 || sum(cutoff) < 0 || !all(cutoff > 0) ||
              length(cutoff) != nclass) {
          stop("Incorrect cutoff specified.")
        }
        if (!is.null(names(cutoff))) {
          if (!all(names(cutoff) %in% levels(y))) {
            stop("Wrong name(s) for cutoff")
          }
          cutoff <- cutoff[levels(y)]
        }
      }
      if (!is.null(classwt)) {
        if (length(classwt) != nclass)
          stop("length of classwt not equal to number of classes")
        ## If classwt has names, match to class labels.
        if (!is.null(names(classwt))) {
          if (!all(names(classwt) %in% levels(y))) {
            stop("Wrong name(s) for classwt")
          }
          classwt <- classwt[levels(y)]
        }
        if (any(classwt <= 0)) stop("classwt must be positive")
        ipi <- 1
      } else {
        classwt <- rep(1, nclass)
        ipi <- 0
      }
    } else addclass <- FALSE
    
    if (missing(proximity)) proximity <- addclass
    if (proximity) {
      prox <- matrix(0.0, n, n)
      proxts <- if (testdat) matrix(0, ntest, ntest + n) else double(1)
    } else {
      prox <- proxts <- double(1)
    }
    
    if (localImp) {
      importance <- TRUE
      impmat <- matrix(0, p, n)
    } else impmat <- double(1)
    
    if (importance) {
      if (nPerm < 1) nPerm <- as.integer(1) else nPerm <- as.integer(nPerm)
      if (classRF) {
        impout <- matrix(0.0, p, nclass + 2)
        impSD <- matrix(0.0, p, nclass + 1)
      } else {
        impout <- matrix(0.0, p, 2)
        impSD <- double(p)
        names(impSD) <- x.col.names
      }
    } else {
      impout <- double(p)
      impSD <- double(1)
    }
    
    nsample <- if (addclass) 2 * n else n
    Stratify <- length(sampsize) > 1
    if ((!Stratify) && sampsize > nrow(x)) stop("sampsize too large")
    if (Stratify && (!classRF)) stop("sampsize should be of length one")
    if (classRF) {
      if (Stratify) {
        if (missing(strata)) strata <- y
        if (!is.factor(strata)) strata <- as.factor(strata)
        nsum <- sum(sampsize)
        if (length(sampsize) > nlevels(strata))
          stop("sampsize has too many elements.")
        if (any(sampsize <= 0) || nsum == 0)
          stop("Bad sampsize specification")
        ## If sampsize has names, match to class labels.
        if (!is.null(names(sampsize))) {
          sampsize <- sampsize[levels(strata)]
        }
        if (any(sampsize > table(strata)))
          stop("sampsize can not be larger than class frequency")
      } else {
        nsum <- sampsize
      }
      nrnodes <- 2 * trunc(nsum / nodesize) + 1
    } else {
      ## For regression trees, need to do this to get maximal trees.
      nrnodes <- 2 * trunc(sampsize/max(1, nodesize - 4)) + 1
    }
    if (!is.null(maxnodes)) {
      ## convert # of terminal nodes to total # of nodes
      maxnodes <- 2 * maxnodes - 1
      if (maxnodes > nrnodes) warning("maxnodes exceeds its max value.")
      nrnodes <- min(c(nrnodes, max(c(maxnodes, 1))))
    }
    ## Compiled code expects variables in rows and observations in columns.
    x <- t(x)
    storage.mode(x) <- "double"
    if (testdat) {
      xtest <- t(xtest)
      storage.mode(xtest) <- "double"
      if (is.null(ytest)) {
        ytest <- labelts <- 0
      } else {
        labelts <- TRUE
      }
    } else {
      xtest <- double(1)
      ytest <- double(1)
      ntest <- 1
      labelts <- FALSE
    }
    nt <- if (keep.forest) ntree else 1
    
    rfout <- .C("iRafNet_regRF",
                x,
                as.double(y),
                as.integer(c(n, p)),
                as.integer(sampsize),
                as.integer(nodesize),
                as.integer(nrnodes),
                as.integer(ntree),
                as.integer(mtry),
                as.integer(c(importance, localImp, nPerm)),
                as.integer(ncat),
                as.integer(maxcat),
                as.integer(do.trace),
                as.integer(proximity),
                as.integer(oob.prox),
                as.integer(corr.bias),
                ypred = double(n),
                impout = impout,
                impmat = impmat,
                impSD = impSD,
                prox = prox,
                ndbigtree = integer(ntree),
                nodestatus = matrix(integer(nrnodes * nt), ncol=nt),
                leftDaughter = matrix(integer(nrnodes * nt), ncol=nt),
                rightDaughter = matrix(integer(nrnodes * nt), ncol=nt),
                nodepred = matrix(double(nrnodes * nt), ncol=nt),
                bestvar = matrix(integer(nrnodes * nt), ncol=nt),
                xbestsplit = matrix(double(nrnodes * nt), ncol=nt),
                mse = double(ntree),
                keep = as.integer(c(keep.forest, keep.inbag)),
                replace = as.integer(replace),
                testdat = as.integer(testdat),
                xts = xtest,
                ntest = as.integer(ntest),
                yts = as.double(ytest),
                labelts = as.integer(labelts),
                ytestpred = double(ntest),
                proxts = proxts,
                msets = double(if (labelts) ntree else 1),
                coef = double(2),
                oob.times = integer(n),
                inbag = if (keep.inbag)
                  matrix(integer(n * ntree), n) else integer(1), sw = as.double(sw))[c(16:28, 36:41)]
    ## Format the forest component, if present.
    if (keep.forest) {
      max.nodes <- max(rfout$ndbigtree)
      rfout$nodestatus <-
        rfout$nodestatus[1:max.nodes, , drop=FALSE]
      rfout$bestvar <-
        rfout$bestvar[1:max.nodes, , drop=FALSE]
      rfout$nodepred <-
        rfout$nodepred[1:max.nodes, , drop=FALSE]
      rfout$xbestsplit <-
        rfout$xbestsplit[1:max.nodes, , drop=FALSE]
      rfout$leftDaughter <-
        rfout$leftDaughter[1:max.nodes, , drop=FALSE]
      rfout$rightDaughter <-
        rfout$rightDaughter[1:max.nodes, , drop=FALSE]
    }
    cl <- match.call()
    cl[[1]] <- as.name("randomForest")
    ## Make sure those obs. that have not been OOB get NA as prediction.
    ypred <- rfout$ypred
    if (any(rfout$oob.times < 1)) {
      ypred[rfout$oob.times == 0] <- NA
    }
    out <- list(call = cl,
                type = "regression",
                predicted = structure(ypred, names=x.row.names),
                mse = rfout$mse,
                rsq = 1 - rfout$mse / (var(y) * (n-1) / n),
                oob.times = rfout$oob.times,
                importance = if (importance) matrix(rfout$impout, p, 2,
                                                    dimnames=list(x.col.names,
                                                                  c("%IncMSE","IncNodePurity"))) else
                                                                    matrix(rfout$impout, ncol=1,
                                                                           dimnames=list(x.col.names, "IncNodePurity")),
                importanceSD=if (importance) rfout$impSD else NULL,
                localImportance = if (localImp)
                  matrix(rfout$impmat, p, n, dimnames=list(x.col.names,
                                                           x.row.names)) else NULL,
                proximity = if (proximity) matrix(rfout$prox, n, n,
                                                  dimnames = list(x.row.names, x.row.names)) else NULL,
                ntree = ntree,
                mtry = mtry,
                forest = if (keep.forest)
                  c(rfout[c("ndbigtree", "nodestatus", "leftDaughter",
                            "rightDaughter", "nodepred", "bestvar",
                            "xbestsplit")],
                    list(ncat = ncat), list(nrnodes=max.nodes),
                    list(ntree=ntree), list(xlevels=xlevels)) else NULL,
                coefs = if (corr.bias) rfout$coef else NULL,
                y = y,
                test = if(testdat) {
                  list(predicted = structure(rfout$ytestpred,
                                             names=xts.row.names),
                       mse = if(labelts) rfout$msets else NULL,
                       rsq = if(labelts) 1 - rfout$msets /
                         (var(ytest) * (n-1) / n) else NULL,
                       proximity = if (proximity)
                         matrix(rfout$proxts / ntree, nrow = ntest,
                                dimnames = list(xts.row.names,
                                                c(xts.row.names,
                                                  x.row.names))) else NULL)
                } else NULL,
                inbag = if (keep.inbag)
                  matrix(rfout$inbag, nrow(rfout$inbag),
                         dimnames=list(x.row.names, NULL)) else NULL)
    
    class(out) <- "randomForest"
    return(out)
  }

"JRF_onetarget" <-
  function(x, y=NULL,  xtest=NULL, ytest=NULL, ntree,
           sampsize,              
           totsize = if (replace) ncol(x) else ceiling(.632*ncol(x)),
           mtry=if (!is.null(y) && !is.factor(y))
             max(floor(nrow(x)/3), 1) else floor(sqrt(nrow(x))),
           replace=TRUE, classwt=NULL, cutoff, strata,
           nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
           maxnodes=NULL,
           importance=FALSE, localImp=FALSE, nPerm=1,
           proximity, oob.prox=proximity,
           norm.votes=TRUE, do.trace=FALSE,
           keep.forest=!is.null(y) && is.null(xtest), corr.bias=FALSE,
           keep.inbag=FALSE, nclasses, ...) {
    
    ww=1/sampsize;
    nclass=mylevels=ipi=sw=NULL
    addclass <- is.null(y)
    classRF <- addclass || is.factor(y)
    if (!classRF && length(unique(y)) <= 5) {
      warning("The response has five or fewer unique values.  Are you sure you want to do regression?")
    }
    if (classRF && !addclass && length(unique(y)) < 2)
      stop("Need at least two classes to do classification.")
    
    n <- ncol(y)           # number of samples
    p <- nrow(x)/nclasses  # number of variables
    
    if (n == 0) stop("data (x) has 0 rows")
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)
    
    keep.forest=!is.null(y) 
    xtest=NULL; ytest=NULL
    testdat <- !is.null(xtest)
    if (testdat) {
      if (ncol(x) != ncol(xtest))
        stop("x and xtest must have same number of columns")
      ntest <- nrow(xtest)
      xts.row.names <- rownames(xtest)
    }
    
    prox <- proxts <- double(1)
    ## Check for NAs.
    if (any(is.na(x))) stop("NA not permitted in predictors")
    if (testdat && any(is.na(xtest))) stop("NA not permitted in xtest")
    if (any(is.na(y))) stop("NA not permitted in response")
    if (!is.null(ytest) && any(is.na(ytest))) stop("NA not permitted in ytest")
    
    if (is.data.frame(x)) {
      xlevels <- lapply(x, mylevels)
      ncat <- sapply(xlevels, length)
      ## Treat ordered factors as numerics.
      ncat <- ifelse(sapply(x, is.ordered), 1, ncat)
      x <- data.matrix(x)
      if(testdat) {
        if(!is.data.frame(xtest))
          stop("xtest must be data frame if x is")
        xfactor <- which(sapply(xtest, is.factor))
        if (length(xfactor) > 0) {
          for (i in xfactor) {
            if (any(! levels(xtest[[i]]) %in% xlevels[[i]]))
              stop("New factor levels in xtest not present in x")
            xtest[[i]] <-
              factor(xlevels[[i]][match(xtest[[i]], xlevels[[i]])],
                     levels=xlevels[[i]])
          }
        }
        xtest <- data.matrix(xtest)
      }
    } else {
      ncat <- rep(1, p)
      xlevels <- as.list(rep(0, p))
    }
    maxcat <- max(ncat)
    if (maxcat > 32)
      stop("Can not handle categorical predictors with more than 32 categories.")
    
    addclass <- FALSE
    
    proximity <- addclass
    
    
    
    impout <- matrix(0.0, p*nclasses, 2)
    impSD <- matrix(0.0, p*nclasses, 1)
    #  names(impSD) <- x.col.names
    
    
    
    nsample <- if (addclass) 2 * n else n
    Stratify <- length(n) > 1
    
    nodesize=5;
    nrnodes <- 2 * trunc(n/max(1, nodesize - 4)) + 1
    
    maxnodes=NULL
    if (!is.null(maxnodes)) {
      ## convert # of terminal nodes to total # of nodes
      maxnodes <- 2 * maxnodes - 1
      if (maxnodes > nrnodes) warning("maxnodes exceeds its max value.")
      nrnodes <- min(c(nrnodes, max(c(maxnodes, 1))))
    }
    
    
    ## Compiled code expects variables in rows and observations in columns.
    # x <- t(x)
    storage.mode(x) <- "double"
    
    xtest <- double(1)
    ytest <- double(1)
    ntest <- 1
    labelts <- FALSE
    nt <- if (keep.forest) ntree else 1
    
    
    nPerm=1
    do.trace=F; oob.prox=F
    corr.bias=FALSE
    keep.inbag=FALSE
    impmat <- double(1)
    replace=T
    
    rfout <- .C("JRF_regRF",
                x,
                y, ww,
                as.integer(c(totsize, p)),
                sampsize=as.integer(sampsize), as.integer(totsize),
                as.integer(nodesize),
                as.integer(nrnodes),
                as.integer(ntree),
                as.integer(mtry),
                as.integer(c(importance, localImp, nPerm)),
                as.integer(ncat),
                as.integer(maxcat),
                as.integer(do.trace),
                as.integer(proximity),
                as.integer(oob.prox),
                as.integer(corr.bias),
                ypred = double(n * nclasses),
                impout = impout,
                impmat = impmat,
                impSD = impSD,
                prox = prox,
                ndbigtree = integer(ntree),
                nodestatus = matrix(integer(nrnodes * nt * nclasses), ncol=nt),
                leftDaughter = matrix(integer(nrnodes * nt * nclasses), ncol=nt),
                rightDaughter = matrix(integer(nrnodes * nt * nclasses), ncol=nt),
                nodepred = matrix(double(nrnodes * nt * nclasses), ncol=nt),
                bestvar = matrix(integer(nrnodes * nt * nclasses), ncol=nt),
                xbestsplit = matrix(double(nrnodes * nt * nclasses), ncol=nt),
                mse = double(ntree * nclasses),
                keep = as.integer(c(keep.forest, keep.inbag)),
                replace = as.integer(replace),
                testdat = as.integer(testdat),
                xts = xtest,
                ntest = as.integer(ntest),
                yts = as.double(ytest),
                labelts = as.integer(labelts),
                ytestpred = double(ntest),
                proxts = proxts,
                msets = double(if (labelts) ntree else 1),
                coef = double(2),
                oob.times = integer(n),
                inbag = if (keep.inbag)
                  matrix(integer(n * ntree), n) else integer(1), as.integer(nclasses))[c(16:28, 36:41)]
    #         ## Format the forest component, if present.
    if (keep.forest) {
      max.nodes <- max(rfout$ndbigtree)
      rfout$nodestatus <-
        rfout$nodestatus[1:max.nodes, , drop=FALSE]
      rfout$bestvar <-
        rfout$bestvar[1:max.nodes, , drop=FALSE]
      rfout$nodepred <-
        rfout$nodepred[1:max.nodes, , drop=FALSE]
      rfout$xbestsplit <-
        rfout$xbestsplit[1:max.nodes, , drop=FALSE]
      rfout$leftDaughter <-
        rfout$leftDaughter[1:max.nodes, , drop=FALSE]
      rfout$rightDaughter <-
        rfout$rightDaughter[1:max.nodes, , drop=FALSE]
    }
    cl <- match.call()
    cl[[1]] <- as.name("randomForest")
    #         ## Make sure those obs. that have not been OOB get NA as prediction.
    ypred <- rfout$ypred
    if (any(rfout$oob.times < 1)) {
      ypred[rfout$oob.times == 0] <- NA
    }
    out <- list(call = cl,
                type = "regression",
                predicted =0,
                mse = rfout$mse,
                rsq = 1 - rfout$mse / (var(y[1,]) * (n-1) / n),
                oob.times = rfout$oob.times,
                importance = if (importance) matrix(rfout$impout, p * nclasses, 2) else
                  matrix(rfout$impout, ncol=1),
                importanceSD=if (importance) rfout$impSD else NULL,
                localImportance = if (localImp)
                  matrix(rfout$impmat, p, n, dimnames=list(x.col.names,
                                                           x.row.names)) else NULL,
                proximity = if (proximity) matrix(rfout$prox, n, n,
                                                  dimnames = list(x.row.names, x.row.names)) else NULL,
                ntree = ntree,
                mtry = mtry,
                forest = if (keep.forest)
                  c(rfout[c("ndbigtree", "nodestatus", "leftDaughter",
                            "rightDaughter", "nodepred", "bestvar",
                            "xbestsplit")],
                    list(ncat = ncat), list(nrnodes=max.nodes),
                    list(ntree=ntree), list(xlevels=xlevels)) else NULL,
                coefs = if (corr.bias) rfout$coef else NULL,
                y = y,
                test = if(testdat) {
                  list(predicted = structure(rfout$ytestpred,
                                             names=xts.row.names),
                       mse = if(labelts) rfout$msets else NULL,
                       rsq = if(labelts) 1 - rfout$msets /
                         (var(ytest) * (n-1) / n) else NULL,
                       proximity = if (proximity)
                         matrix(rfout$proxts / ntree, nrow = ntest,
                                dimnames = list(xts.row.names,
                                                c(xts.row.names,
                                                  x.row.names))) else NULL)
                } else NULL,
                inbag = if (keep.inbag)
                  matrix(rfout$inbag, nrow(rfout$inbag),
                         dimnames=list(x.row.names, NULL)) else NULL)
    
    
    #      print(rfout$mse)
    class(out) <- "randomForest"
    return(out)
    
  }

"iJRF_parallel" <-
  function(X, W=NULL,ntree=NULL,mtry=NULL,genes.name=NULL,parallel) {
    
    p<-dim(X[[1]])[1];
    num.par<-parallel[1]; num.targets<-parallel[2]
    if (is.null(mtry)) mtry=sqrt(p)
    if (is.null(ntree)) ntree=1000
    if (is.null(genes.name)) genes.name=paste("G",seq(1,p),sep="")
    
    nclasses<-length(X)
    sampsize<-rep(0,nclasses)
    
    index.target<-c(round(p/num.targets)*(num.par-1)+1,round(p/num.targets)*num.par);
    if (index.target[2]>p) index.target[2]<-p
    if (index.target[1]>p) {print("Error: Maximum number of parallelizations excedeed")} else {
      
      imp<-array(0,c(p,index.target[2]-index.target[1]+1,nclasses))
      index<-seq(1,p)
      
      for (j in 1:nclasses) {  X[[j]] <- t(apply(X[[j]], 1, function(x) { (x - mean(x)) / sd(x) } ))    
                               sampsize[j]<-dim(X[[j]])[2] }
      tot<-max(sampsize);
      print(is.null(W))
      
      
      if (is.null(W)) { # -- implement standard JRF
        j.s=0
        for (j in index.target[1]:index.target[2]){
          
          j.s=j.s+1;
          covar<-matrix(0,(p-1)*nclasses,tot)             
          y<-matrix(0,nclasses,tot)             
          
          for (c in 1:nclasses)  {
            y[c,seq(1,sampsize[c])]<-as.matrix(X[[c]][j,])
            covar[seq((c-1)*(p-1)+1,c*(p-1)),seq(1,sampsize[c])]<-X[[c]][-j,]
          }
          jrf.out<-JRF_onetarget(x=covar,y=y,mtry=mtry,importance=TRUE,sampsize=sampsize,nclasses=nclasses,ntree=ntree)
          for (s in 1:nclasses) imp[-j,j.s,s]<-importance(jrf.out,scale=FALSE)[seq((p-1)*(s-1)+1,(p-1)*(s-1)+p-1)]  #- save importance score for net1
          
        }
      } else { # -- implement iJRF (integrative JRF)
        
        j.s=0
        for (j in index.target[1]:index.target[2]){
          j.s<-j.s+1
          weights.rf<-as.matrix(W[,j]); 
          weights.rf[j]<-0
          weights.rf<-weights.rf/sum(weights.rf);
          
          w.sorted<-sort(weights.rf,decreasing = FALSE,index.return=T)
          index<-w.sorted$ix
          w.sorted<-w.sorted$x
          
          covar<-matrix(0,p*nclasses,tot)             
          y<-matrix(0,nclasses,tot)             
          for (c in 1:nclasses)  {
            y[c,seq(1,sampsize[c])]<-X[[c]][j,]
            covar[seq((c-1)*(p)+1,c*p),seq(1,sampsize[c])]<-X[[c]][index,]
          }
          
          jrf.out<-iJRF_onetarget(x=covar,y=y,mtry=mtry,importance=TRUE,sampsize=sampsize,
                                  nclasses=nclasses,ntree=ntree,sw=as.double(w.sorted))
          
          for (s in 1:nclasses) imp[index,j.s,s]<-importance(jrf.out,scale=FALSE)[seq(p*(s-1)+1,p*s)]  #- save importance score for net1
          
        }
      }
      
      out<-list(num.par=num.par,model="iJRF",importance=imp,genes.name=genes.name)
      return(out)
      
    }    
  }


"iRafNet_parallel" <- function(X,W,ntree=NULL,mtry=NULL,genes.name,parallel) {
  
  X<-t(X[[1]])
  p<-dim(X)[2]
  if (is.null(ntree)) ntree=1000
  if (is.null(genes.name)) genes.name=paste("G",seq(1,p),sep="")
  if (is.null(W)) W=matrix(1,p,p)
  if (is.null(mtry)) mtry=sqrt(p)
  
  
  num.par<-parallel[1]; num.targets<-parallel[2]
  index.target<-c(round(p/num.targets)*(num.par-1)+1,round(p/num.targets)*num.par);
  if (index.target[2]>p) index.target[2]<-p
  if (index.target[1]>p) {print("Error: Maximum number of parallelizations excedeed")} else {
    
    imp<-matrix(0,p,index.target[2]-index.target[1]+1)
    X <- (apply(X, 2, function(x) { (x - mean(x)) / sd(x) } ))    
  
    j.s=0;
    
    for (j in index.target[1]:index.target[2]){
      j.s=j.s+1;
      y<-X[,j];
      
      weights.rf<-as.matrix(W[,j]); 
      weights.rf[j]<-0
      weights.rf<-weights.rf/sum(weights.rf);
      
      w.sorted<-sort(weights.rf,decreasing = FALSE,index.return=T)
      index<-w.sorted$ix
      x.sorted<-X[,index]
      w.sorted<-w.sorted$x
      
      rout<-irafnet_onetarget(x=x.sorted,y=as.double(y),importance=TRUE,mtry=round(sqrt(p-1)),ntree=1000,
                              sw=as.double(w.sorted))
      
      imp[index,j.s]<-c(importance(rout))
    }
    out<-list(num.par=num.par,model="iRafNet",importance=imp,genes.name=genes.name)
  }
  return(out)
}

"ptmJRF_parallel" <-
  function(X, ntree,mtry,genes.name,ptm.name,parallel) {
    
    p<-length(genes.name); ptm.p<-length(ptm.name)
    if (is.null(mtry)) mtry=sqrt(p)
    num.par<-parallel[1]; num.targets<-parallel[2]
    nclasses<-length(X)
    sampsize<-rep(0,nclasses)
    for (j in 1:nclasses) {  X[[j]] <- t(apply(X[[j]], 1, function(x) { (x - mean(x)) / sd(x) } ))    
                             sampsize[j]<-dim(X[[j]])[2] }
    
    index.target<-c(round(p/num.targets)*(num.par-1)+1,round(p/num.targets)*num.par);
    if (index.target[2]>p) index.target[2]<-p
    if (index.target[1]>p) {print("Error: Maximum number of parallelizations excedeed")} else {
      
      # --- reorder rows in PTM object X[[1]]
      X.ptm<-X[[1]]; s=0 ; locptm<-numptm<-rep(0,p)
      ptm.new<-ptm.name
      for (j in 1:p){ 
        ptm.j<-X[[1]][ptm.name==genes.name[j],]
        n.j<-sum(ptm.name==genes.name[j])
        X.ptm[seq(s+1,s+n.j),]<-ptm.j
        locptm[j]<-(s+1)
        numptm[j]<-n.j
        ptm.new[seq(s+1,s+n.j)]<-rep(genes.name[j],n.j)
        s<-s+n.j
      }
      X[[1]]<-X.ptm
      ptm.name<-ptm.new
      index<-seq(1,p)
      
      imp<-array(0,c(p,index.target[2]-index.target[1]+1,nclasses))
      
      # -- for loop over target proteins (1)
      j.s=0; 
      for (jj in index.target[1]:index.target[2]){ 
        j.s=j.s+1;
        index.ptm<-seq(1,ptm.p); 
        index.ptm<-index.ptm[ptm.name==genes.name[jj]] # - ptm corresponding to target protein
        imp.momentary<-array(0,c(p,numptm[jj],nclasses)); 
        k.s<-0
        
        # --- extract response and predictors for non-ptm data
        covar.nptm<-matrix(0,(p-1)*(nclasses-1),max(sampsize))             
        y.nptm<-matrix(0,nclasses-1,max(sampsize))             
        for (c in 2:nclasses) {
          y.nptm[c-1,seq(1,sampsize[c])]<-as.matrix(X[[c]][jj,])
          covar.nptm[seq((p-1)*(c-2)+1,(p-1)*(c-1)),seq(1,sampsize[c])]<-as.matrix(X[[c]][-jj,])
        }
        
        # -- for loop over ptm corresponding to target protein (2)
        for (k in 1:numptm[jj]){
          k.s<-k.s+1
          y.ptm<-matrix(0,1,max(sampsize)); covar.ptm<-matrix(0,ptm.p-numptm[jj],max(sampsize))  
          y.ptm[,seq(1,sampsize[1])]<-as.matrix(X[[1]][index.ptm[k],])
          covar.ptm[,seq(1,sampsize[1])]<-X[[1]][ptm.name!=genes.name[jj],]
          covar<-rbind(covar.ptm,covar.nptm)
          y<-rbind(y.ptm,y.nptm)
          
          numptm.j<-numptm[-jj]
          locptm.j<-locptm; 
          if (jj < p) locptm.j[seq(jj+1,p)]<-locptm.j[seq(jj+1,p)]-numptm[jj]
          locptm.j<-locptm.j[-jj]
          
          rfout<-ptmJRF_onetarget(x=covar,y=y,p=(p-1),mptm=ptm.p-numptm[jj],
                                  mtry=sqrt(p-1),importance=TRUE,sampsize=sampsize,
                                  nclasses=nclasses,
                                  ntree=ntree,numptm=numptm.j,locptm=locptm.j)
          
          imp.rfout<-importance(rfout)
          for (s in 1:nclasses) imp.momentary[-jj,k.s,s]<-imp.rfout[seq((p-1)*(s-1)+1,(p-1)*s)] 
        } # -- End for loop (2)
        
        if (numptm[jj]>1)  {
          for (c in 1:nclasses) imp[,j.s,c]<-apply(imp.momentary[,,c], 1, function(x) { mean(x) } )
        } else {
          for (c in 1:nclasses) imp[,j.s,c]<-imp.momentary[,1,c]        
        }
        
      } # -- End for loop (1)
      
      out<-list(num.par=num.par,model="ptmJRF",importance=imp,genes.name=genes.name,ptm.name=ptm.name)
      return(out)
    }
  }


# --- MAIN function
"iJRFNet_parallel" <-
  function(X, W=NULL,ntree=NULL,mtry=NULL,model=NULL,genes.name,ptm.name=NULL,parallel) {
    
     if (is.null(model)) {print("Error: Specify Model")} else { 
      if (is.null(ntree)) ntree=1000
      
      if (model=="iJRF")  out<-iJRF_parallel(X=X,W=W,ntree=ntree,mtry=mtry,genes.name=genes.name,parallel=parallel)
      if (model=="iRafNet") out<-iRafNet_parallel(X=X, W=W,ntree=ntree,mtry=mtry,genes.name=genes.name,parallel=parallel)
      if (model=="ptmJRF") out<-ptmJRF_parallel(X=X,ntree=ntree,mtry=mtry,genes.name=genes.name,ptm.name=ptm.name,parallel=parallel)
      return(out)
    }
  }