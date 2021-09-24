# --- Functions Called by Main Function Run_permutation

"iJRF_permutation" <-
  function(X, W=NULL,ntree=NULL,mtry=NULL,genes.name=NULL,seed,to.store=NULL) {
    
    nclasses<-length(X)
    sampsize<-rep(0,nclasses)
    p<-dim(X[[1]])[1];
    
    if (is.null(mtry)) mtry=sqrt(p)
    if (is.null(ntree)) ntree=1000
    if (is.null(genes.name)) genes.name=paste("G",seq(1,p),sep="")
    
    imp<-array(0,c(p,length(genes.name),nclasses))
    if (is.null(to.store)) {imp.final<-matrix(0,p*(p-1)/2,nclasses)} else {imp.final<-matrix(0,to.store,nclasses);}
    index<-seq(1,p)
    
    for (j in 1:nclasses) {  X[[j]] <- t(apply(X[[j]], 1, function(x) { (x - mean(x)) / sd(x) } ))    
                             sampsize[j]<-dim(X[[j]])[2] }
    tot<-max(sampsize);
    
    if (is.null(W)) { # -- implement standard JRF
      
      for (j in 1:length(genes.name)){
        set.seed((seed-1)*nclasses+1)
        covar<-matrix(0,(p-1)*nclasses,tot)             
        y<-matrix(0,nclasses,tot)             
        
        for (c in 1:nclasses)  {
          y[c,seq(1,sampsize[c])]<-as.matrix(X[[c]][j,sample(sampsize[c])])
          covar[seq((c-1)*(p-1)+1,c*(p-1)),seq(1,sampsize[c])]<-X[[c]][-j,]
        }
        jrf.out<-JRF_onetarget(x=covar,y=y,mtry=mtry,importance=TRUE,sampsize=sampsize,nclasses=nclasses,ntree=ntree)
        for (s in 1:nclasses) imp[-j,j,s]<-importance(jrf.out,scale=FALSE)[seq((p-1)*(s-1)+1,(p-1)*(s-1)+p-1)]  #- save importance score for net1
        
      }
    } else { # -- implement iJRF (integrative JRF)
      
      for (j in 1:length(genes.name)){
        set.seed((seed-1)*nclasses+1)
        weights.rf<-as.matrix(W[,j]); 
        weights.rf[j]<-0
        weights.rf<-weights.rf/sum(weights.rf);
        
        w.sorted<-sort(weights.rf,decreasing = FALSE,index.return=T)
        index<-w.sorted$ix
        w.sorted<-w.sorted$x
        
        covar<-matrix(0,p*nclasses,tot)             
        y<-matrix(0,nclasses,tot)             
        for (c in 1:nclasses)  {
          y[c,seq(1,sampsize[c])]<-X[[c]][j,sample(sampsize[c])]
          covar[seq((c-1)*(p)+1,c*p),seq(1,sampsize[c])]<-X[[c]][index,]
        }
        
        jrf.out<-iJRF_onetarget(x=covar,y=y,mtry=mtry,importance=TRUE,sampsize=sampsize,
                                nclasses=nclasses,ntree=ntree,sw=as.double(w.sorted))
        
        for (s in 1:nclasses) imp[index,j,s]<-importance(jrf.out,scale=FALSE)[seq(p*(s-1)+1,p*s)]  #- save importance score for net1
        
      }
      
    }
    
    # --- Derive importance score for each interaction 
    for (s in 1:nclasses){ 
      imp.s<-imp[,,s]; t.imp<-t(imp.s)
      if (is.null(to.store)) imp.final[,s]<-(imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2        
      if (is.null(to.store)==FALSE) imp.final[,s]<-sort((imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2,decreasing=TRUE)[seq(1,to.store)]        
    }
    
    out<-(imp.final)
    return(out)
    
  }


"iRafNet_permutation" <-
  function(X,W,ntree=NULL,mtry=NULL,genes.name,seed,to.store=NULL) {
    
    X<-t(X[[1]])
    p<-dim(X)[2]
    imp<-matrix(0,p,p)
    index<-seq(1,p)
    
    if (is.null(to.store)) {imp.final<-matrix(0,p*(p-1)/2,1);} else {imp.final<-matrix(0,to.store,1);}
    
    X <- (apply(X, 2, function(x) { (x - mean(x)) / sd(x) } ))    
    set.seed(seed)
    perm.sample<-sample(dim(X)[1])
    for (j in 1:p){ 
      y<-X[perm.sample,j];
      
      weights.rf<-as.matrix(W[,j]); 
      weights.rf[j]<-0
      weights.rf<-weights.rf/sum(weights.rf);
      
      w.sorted<-sort(weights.rf,decreasing = FALSE,index.return=T)
      index<-w.sorted$ix
      x.sorted<-X[,index]
      w.sorted<-w.sorted$x
      
      rout<-irafnet_onetarget(x=x.sorted,y=as.double(y),importance=TRUE,mtry=round(sqrt(p-1)),ntree=1000,
                              sw=as.double(w.sorted))
      
      imp[index,j]<-c(importance(rout))
    }
    
    # --- Return importance score for each interaction
    imp.s<-imp; t.imp<-t(imp.s)
    if (is.null(to.store)) imp.final<-(imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2        
    if (is.null(to.store)==FALSE) imp.final<-sort((imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2,decreasing=TRUE)[seq(1,to.store)]        
    
    return(imp.final)
    
    
  }

"iJRF_permutation" <-
  function(X, W=NULL,ntree=NULL,mtry=NULL,genes.name=NULL,seed,to.store=NULL) {
    
    nclasses<-length(X)
    sampsize<-rep(0,nclasses)
    p<-dim(X[[1]])[1];
    
    if (is.null(mtry)) mtry=sqrt(p)
    if (is.null(ntree)) ntree=1000
    if (is.null(genes.name)) genes.name=paste("G",seq(1,p),sep="")
    
    imp<-array(0,c(p,length(genes.name),nclasses))
    if (is.null(to.store)) {imp.final<-matrix(0,p*(p-1)/2,nclasses);} else {imp.final<-matrix(0,to.store,nclasses);}
    index<-seq(1,p)
    
    for (j in 1:nclasses) {  X[[j]] <- t(apply(X[[j]], 1, function(x) { (x - mean(x)) / sd(x) } ))    
                             sampsize[j]<-dim(X[[j]])[2] }
    tot<-max(sampsize);
    set.seed(seed); sample.perm<-list() # -- permuted sample for each class      
    for (j in 1:nclasses) {sample.perm[[j]]<-sample(sampsize[j])}
    
    if (is.null(W)) { # -- implement standard JRF
      
      for (j in 1:length(genes.name)){
        set.seed((seed-1)*nclasses+1)
        covar<-matrix(0,(p-1)*nclasses,tot)             
        y<-matrix(0,nclasses,tot)             
        
        for (c in 1:nclasses)  {
          y[c,seq(1,sampsize[c])]<-as.matrix(X[[c]][j,sample.perm[[c]]])
          covar[seq((c-1)*(p-1)+1,c*(p-1)),seq(1,sampsize[c])]<-X[[c]][-j,]
        }
        jrf.out<-JRF_onetarget(x=covar,y=y,mtry=mtry,importance=TRUE,sampsize=sampsize,nclasses=nclasses,ntree=ntree)
        for (s in 1:nclasses) imp[-j,j,s]<-importance(jrf.out,scale=FALSE)[seq((p-1)*(s-1)+1,(p-1)*(s-1)+p-1)]  #- save importance score for net1
        
      }
    } else { # -- implement iJRF (integrative JRF)
      
      for (j in 1:length(genes.name)){
        set.seed((seed-1)*nclasses+1)
        weights.rf<-as.matrix(W[,j]); 
        weights.rf[j]<-0
        weights.rf<-weights.rf/sum(weights.rf);
        
        w.sorted<-sort(weights.rf,decreasing = FALSE,index.return=T)
        index<-w.sorted$ix
        w.sorted<-w.sorted$x
        
        covar<-matrix(0,p*nclasses,tot)             
        y<-matrix(0,nclasses,tot)             
        for (c in 1:nclasses)  {
          y[c,seq(1,sampsize[c])]<-X[[c]][j,sample.perm[[c]]]
          covar[seq((c-1)*(p)+1,c*p),seq(1,sampsize[c])]<-X[[c]][index,]
        }
        
        jrf.out<-iJRF_onetarget(x=covar,y=y,mtry=mtry,importance=TRUE,sampsize=sampsize,
                                nclasses=nclasses,ntree=ntree,sw=as.double(w.sorted))
        
        for (s in 1:nclasses) imp[index,j,s]<-importance(jrf.out,scale=FALSE)[seq(p*(s-1)+1,p*s)]  #- save importance score for net1
        
      }
      
    }
    
    # --- Derive importance score for each interaction 
    for (s in 1:nclasses){ 
      imp.s<-imp[,,s]; t.imp<-t(imp.s)
      if (is.null(to.store)) imp.final[,s]<-(imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2        
      if (is.null(to.store)==FALSE) imp.final[,s]<-sort((imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2,decreasing=TRUE)[seq(1,to.store)]        
    }
    
    out<-(imp.final)
    return(out)
    
  }
# ---- MAIN Code    ----------------------------------------------------------------- #
"iJRFNet_permutation" <-
  function(X, W=NULL,ntree=NULL,mtry=NULL,genes.name=NULL,M=NULL,model,ptm.name=NULL,to.store=NULL) {

    if (model=="iJRF" | model =="iRafNet") p<-dim(X[[1]])[1]  
    if (model=="ptmJRF") p<-dim(X[[2]])[1]  
    
    if (is.null(model)) { print("Error: Specify Model") } else { 
    if (is.null(M)) M=100
    if (is.null(mtry)) mtry=sqrt(p)
    if (is.null(ntree)) ntree=1000
    if (is.null(genes.name)) genes.name=paste("G",seq(1,p),sep="")
      
    nclasses<-length(X)
    perm<-array(0,c((p^2-p)/2,M,nclasses)) 
    for (j in 1:M)  {
     if (model=="iJRF") perm[,j,]<-as.matrix(iJRF_permutation(X=X,W=W, ntree=ntree,mtry=mtry,genes.name=genes.name,seed=j,to.store=to.store))
     if (model=="iRafNet") perm[,j]<-as.matrix(iRafNet_permutation(X=X, W=W, ntree=ntree,mtry=mtry,genes.name=genes.name,seed=j,to.store=to.store))   
     if (model=="ptmJRF") perm[,j,]<-as.matrix(ptmJRF_permutation(X=X, ntree=ntree,mtry=mtry,genes.name=genes.name,ptm.name=ptm.name,seed=j,to.store=to.store))       
    }
    return(perm)
    
    }
  }
    