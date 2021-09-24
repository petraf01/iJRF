\name{FinalScore_parallel_permutation}
\alias{FinalScore_parallel_permutation}
\title{Derive final importance scores of interactions for one permuted data set run in parallel.}
\usage{
FinalScore_parallel_permutation(importance, model, genes.name,
                                 to.store=NULL)
                                 }
\arguments{
\item{importance}{A matrix containing importance scores. When model \code{iRafNet} is implemented, \code{importance} is a two dimensional matrix (\code{p} x \code{p}) with \code{p} being the total number of genes/proteins. When either function \code{iJRF} or \code{ptmJRF} is implemented, \code{importance} is a three dimensional matrix of importance scores (\code{p} x \code{p} x \code{C}) with \code{p} being the total number of genes/proteins and \code{C} the number of classes.}

\item{model}{ Variable indicating which iJRFNet model needs to be imlemented. Takes values in \code{c("iJRF", "iRafNet","ptmJRF")} }

\item{genes.name}{Vector containing genes name. The order needs to match the rows/columns of \code{importance}.}

\item{to.store}{Optional Integer. Total number of importance scores to be stored. When omitted, all importance scores will be returned. Note that to compute FDR and derive the final network via function \code{Derive_network} we do not need all \code{(p-p) x p / 2} importance scores where \code{p} is the total number of proteins/genes. A sufficiently large number would work. This number is usually chosen based on the number of nodes and is the maximum number of interactions that you would expect. }
}
\value{
A matrix with \code{to.store} rows and \code{C} columns where \code{to.store} is the numbere of top importance scores and \code{C} is the number of classes.  When \code{to.store} is omitted, all importance scores will be returned.
}
\description{
This function returns final importance scores of interactions for one permuted data set run in parallel.
}
\examples{

 # --- Generate data sets
 nclasses=2               # number of data sets / classes
 n1<-n2<-20               # sample size for each data sets
 p<-5                   # number of variables (genes/proteins)
 genes.name<-paste("G",seq(1,p),sep="")   # genes/proteins name
 
 data1<-matrix(rnorm(p*n1),p,n1)       # generate data1
 data2<-matrix(rnorm(p*n2),p,n1)       # generate data2

  ##---------------------------------------------------------------------------##
  ## ---  Run iJRF 
  
  # --- Run multiple jobs in parallel and combine them
   out.new<-array(0,c(p,p,nclasses))
   n.var=0
   for (k in 1:3){ 
      out<-iJRFNet_parallel_permutation(X=list(data1,data2),
      genes.name=genes.name, model="iJRF",parallel=c(k,2),seed=1)
      
      n.target<-dim(out$importance)[2]
      for (c in 1:nclasses) {
      out.new[,seq(n.var+1,n.var+n.target),c]<-out$importance[,,c];}
      n.var=n.var+n.target
    }
     
   # --- Derive interactions 
   FinalScore_parallel_permutation(importance=out.new,model="iJRF",
                                   genes.name=genes.name)


}    
\references{
Petralia, F., Song, W.M., Tu, Z. and Wang, P. (2016). New method for joint network analysis reveals common and different coexpression patterns among genes and proteins in breast cancer. \emph{Journal of proteome research}, \bold{15}(3), pp.743-754.

}
