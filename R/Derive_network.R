Derive_network <- function(out.iJRFNet,out.perm,TH) {
 
  nclasses<-dim(out.perm)[3]
  if (is.null(nclasses)) nclasses=1; # -- if iRafNet

  M<-dim(out.perm)[2]; out<-list()
  
  for (net in 1:nclasses) { 
    j.np<-sort(out.iJRFNet[,2+net],decreasing=TRUE)
    FDR<-matrix(0,dim(out.perm)[1],1); 
    for (s in 1:length(j.np)){ 
      FP<-sum(sum(out.perm[,,net]>=j.np[s]))/M
      FDR[s]<-FP/s;
      if (FDR[s]>TH) {th<-j.np[s];  break;}
    }
    out[[net]]<-out.iJRFNet[out.iJRFNet[,2+net]>th,seq(1,2)]
  }
  return(out)
}