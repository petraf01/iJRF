"FinalScore_parallel" <-
  function(importance, model,genes.name) {
     
    p<-dim(importance)[1]
    if (is.na(dim(importance)[3])==FALSE){ 
     nclasses<-dim(importance)[3]
     imp.final<-matrix(0,p*(p-1)/2,nclasses); } else imp.final<-matrix(0,p*(p-1)/2,1);

     vec1<-matrix(rep(genes.name,p),p,p)
     vec2<-t(vec1)
     vec1<-vec1[lower.tri(vec1,diag=FALSE)]
     vec2<-vec2[lower.tri(vec2,diag=FALSE)]
     
   if (model=="iJRF" | model=="ptmJRF") { 
     for (s in 1:nclasses){ 
       imp.s<-importance[,,s]; t.imp<-t(imp.s)
       imp.final[,s]<-(imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2        
     }
     out<-cbind(as.character(vec1),as.character(vec2),as.data.frame(imp.final),stringsAsFactors=FALSE)
     colnames(out)<-c(paste0('gene',1:2),paste0('importance',1:nclasses))
   }
   if (model=="iRafNet") {
     imp.s<-importance; t.imp<-t(imp.s)
     imp.final<-matrix(0,length(vec1),1)
     imp.final<-(imp.s[lower.tri(imp.s,diag=FALSE)]+t.imp[lower.tri(t.imp,diag=FALSE)])/2        
     out<-cbind(as.character(vec1),as.character(vec2),imp.final,stringsAsFactors=FALSE)
     colnames(out)<-c(paste0('gene',1:2),'importance')
   }
   
   return(out)
   
}