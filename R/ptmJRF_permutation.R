#' Joint Random Forest for the simultaneous estimation of multiple related networks
#'
#' MAIN FUNCTION -- > JRF
#' 
#' INPUT
#' 
#' X            list object containing data for each class
#' ntree        number of trees
#' mtry         number of variables to be sampled at each node
#' genes.name   list of gene names 
#' 
#' OUTPUT: importance score of interactions
#'
#'
#' OTHER FUNCTIONS -- > importance  and  JRF_onetarget
#' 
#' importance     compute importance score for an object of class JRF 
#' (this file is a modified version of file importance contained in package randomForest, A. Liaw and M. Wiener (2002))
#' 
#' JRF_onetarget  for each class, model the expression of a target gene as a function of the expression of other genes via random forest. 
#'                class specific tree ensemble are designed to borrow information across them. 
#' (this file is a modified version of file randomForest contained in package randomForest, A. Liaw and M. Wiener (2002))
#'   
#"JRF" <-  function(X, ...)UseMethod("JRF")


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


# --- Function JRF_onetarget is a modified version of function randomForest from R package randomForest

"ptmJRF_onetarget" <-
  function(x, y=NULL, p,mptm, xtest=NULL, ytest=NULL, ntree,
           sampsize,mtry,
           replace=TRUE, classwt=NULL, cutoff, strata,
           nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
           maxnodes=NULL,
           importance=FALSE, localImp=FALSE, nPerm=1,
           proximity, oob.prox=proximity,
           norm.votes=TRUE, do.trace=FALSE,
           keep.forest=!is.null(y) && is.null(xtest), corr.bias=FALSE,
           keep.inbag=FALSE, nclasses, numptm, locptm, ...) {
    #--- p is the number of genes    
    #--- mptm number of ptm sites
    
    ww=1/sampsize;
    totsize<-max(sampsize)
    nclass=mylevels=ipi=sw=NULL
    addclass <- is.null(y)
    classRF <- addclass || is.factor(y)
    if (!classRF && length(unique(y)) <= 5) {
      warning("The response has five or fewer unique values.  Are you sure you want to do regression?")
    }
    if (classRF && !addclass && length(unique(y)) < 2)
      stop("Need at least two classes to do classification.")
    
    n <- ncol(y)           # number of samples
    
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
    
      rfout <- .C("ptmJRF_regRF",
                  x,
                  y, ww,
                  as.integer(c(totsize, p,mptm)),
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
                  as.integer(numptm), as.integer(locptm))[c(16:28, 36:41)]
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


# --- MAIN function
"ptmJRF_permutation" <-
  function(X, ntree=NULL,mtry=NULL,genes.name,ptm.name,seed,to.store=NULL) {

    p<-length(genes.name); ptm.p<-length(ptm.name)
    nclasses<-length(X)
    sampsize<-rep(0,nclasses)
    for (j in 1:nclasses) {  X[[j]] <- t(apply(X[[j]], 1, function(x) { (x - mean(x)) / sd(x) } ))    
                             sampsize[j]<-dim(X[[j]])[2] }
    
    set.seed(seed); sample.perm<-list() # -- permuted sample for each class      
    for (j in 1:nclasses) {sample.perm[[j]]<-sample(sampsize[j])}
    
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
    
    imp<-array(0,c(p,length(genes.name),nclasses))
    
    if (is.null(to.store)) {imp.final<-matrix(0,p*(p-1)/2,nclasses);} else {imp.final<-matrix(0,to.store,nclasses);}
    index<-seq(1,p)
    imp<-array(0,c(p,ptm.p,nclasses))
    
    for (j in 1:ptm.p){
      set.seed((seed-1)*nclasses+1)
      covar<-matrix(0,ptm.p*nclasses,max(sampsize))             
      y<-matrix(0,nclasses,max(sampsize))             
      
      for (c in 1:nclasses)  {
        if (c==1) {
          y[c,seq(1,sampsize[c])]<-as.matrix(X[[c]][j,sample.perm[[c]]])
          covar[seq(1,ptm.p-numptm[genes.name==ptm.name[j]]),seq(1,sampsize[c])]<-X[[c]][-seq(locptm[genes.name==ptm.name[j]],locptm[genes.name==ptm.name[j]]+numptm[genes.name==ptm.name[j]]-1),]
          n.covar<-ptm.p-numptm[genes.name==ptm.name[j]]         } else {
            y[c,seq(1,sampsize[c])]<-as.matrix(X[[c]][genes.name==ptm.name[j],sample.perm[[c]]])
            covar[seq(n.covar+1,n.covar+p-1),seq(1,sampsize[c])]<-X[[c]][genes.name!=ptm.name[j],]
            n.covar<-n.covar+p-1
          }      
      }
      covar<-covar[seq(1,n.covar),]
      numptm.j<-numptm[genes.name!=ptm.name[j]]
      index<-seq(1,length(locptm))
      index<-index[genes.name==ptm.name[j]]
      locptm.j<-locptm; 
      if (index != p) locptm.j[seq(index+1,length(locptm))]<-locptm.j[seq(index+1,length(locptm))]-numptm[index]
      locptm.j<-locptm.j[-index]
      
      rfout<-ptmJRF_onetarget(x=covar,y=y,p=(p-1),mptm=ptm.p-numptm[genes.name==ptm.name[j]],
                              mtry=sqrt(p-1),importance=TRUE,sampsize=sampsize,nclasses=nclasses,
                              ntree=ntree,numptm=numptm.j,locptm=locptm.j)
      
      imp.rfout<-importance(rfout)
      for (s in 1:nclasses) imp[genes.name!=ptm.name[j],j,s]<-imp.rfout[seq((p-1)*(s-1)+1,(p-1)*(s-1)+p-1)] 
    }
    
    imp.new<-array(0,c(p,p,nclasses))
    for (j in 1:p){
      if (numptm[j]==1){  
        for (c in 1:nclasses) imp.new[,j,c]<-imp[,locptm[j],c]
      } else {  
        for (c in 1:nclasses) imp.new[,j,c]<-apply(imp[,seq(locptm[j],locptm[j]+numptm[j]-1),c], 1, function(x) { mean(x) } )
      }
    }
    
    # --- Derive importance score for each interaction 
    for (s in 1:nclasses){ 
      imp.s<-imp.new[,,s]; t.imp<-t(imp.s)
      if (is.null(to.store)) imp.final[,s]<-(imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2        
      if (is.null(to.store)==FALSE) imp.final[,s]<-sort((imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2,decreasing=TRUE)[seq(1,to.store)]        
    }
    
    return(imp.final)
    
  }
