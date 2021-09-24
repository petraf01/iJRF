"FinalScore_parallel_permutation" <-
  function(importance, model,genes.name,to.store=NULL) {
    
    p<-dim(importance)[1]
    if (is.null(to.store)) to.store<-(p-1)*p/2

    if (is.na(dim(importance)[3])==FALSE){ 
      nclasses<-dim(importance)[3]
      imp.final<-matrix(0,to.store,nclasses); } else imp.final<-matrix(0,p*(p-1)/2,1);
        
    if (model=="iJRF" | model=="ptmJRF") { 
      for (s in 1:nclasses){ 
        imp.s<-importance[,,s]; t.imp<-t(imp.s)
        imp.final[,s]<-sort((imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2,decreasing=TRUE)[seq(1,to.store)]        
      }
      
      out<-imp.final
    }
    if (model=="iRafNet") {
      imp.s<-importance; t.imp<-t(imp.s)
      imp.final<-sort((imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2,decreasing=TRUE)[seq(1,to.store)]        
      out<-imp.final
    }
    
    return(out)
    
  }