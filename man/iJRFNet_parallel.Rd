\name{iJRFNet_parallel}
\alias{iJRFNet_parallel}
\title{Derive importance scores for a subset of target genes for functions of class iJRFNet.}
\usage{
iJRFNet_parallel(X, W=NULL, ntree=NULL, mtry=NULL, model=NULL, 
                genes.name, ptm.name=NULL,parallel)}
\arguments{
\item{X}{List object containing expression data for each class, \code{X=list(x_1,x_2, ... )} where \code{x_j} is a \code{(p x n_j)} matrix with rows corresponding to genes and columns to samples. Rows need to be the same across objects, while samples can vary. Missing values are not allowed. If \code{model="ptmJRF"}, the first object of the list must contain the expression of post translational modification variables. Only in this case, the number of variables in the first object might differ from that of other objects. Rows of \code{X[[1]]} does not need to be ordered in a specific way.
}

\item{W}{\code{(p x p)} Optional symmetric matrix containing sampling scores. When omitted, the standard JRF algorithm without weighted sampling scheme will be implemented. Element \code{(i,j)} contains interaction score \code{(i - j)}. Scores must be non-negative. Larger value of sampling score corresponds to higher likelihood of gene \code{i} interacting with gene \code{j}. Columns and rows of \code{W} must be in the same order as the columns of \code{X}. Sampling scores \code{W} are computed considering one prior data such as protein-protein interactions.}

\item{ntree}{Numeric value: number of trees. If omitted, 1000 trees are considered.}

\item{mtry}{Numeric value: number of predictors to be sampled at each node. If omitted, \code{mtry} is set to the square root of the number of variables.}

\item{model}{Variable indicating which iJRFNet model needs to be run. Takes values in \code{c("iJRF", "iRafNet","ptmJRF")} }

\item{genes.name}{Vector containing genes name. The order needs to match the rows of \code{x_j}.}

\item{ptm.name}{List of post translational modification variables in protein domain. This list must be ordered as rows of \code{X[[1]]}. This is required only if function \code{ptmJRF} is implemented.}

\item{parallel}{Vector containing two elements \code{c(num.job,num.targets)}. The first element is the job number that is implemented, target genes will be divided in \code{J} jobs each containing a specific number of target genes. The second element contains the number of target genes considered in each job.}

}
\value{
List object containing:

\item{num.par}{Integer. Parallel batch implemented.}

\item{model}{Variable indicating which iJRFNet model needs to be run. Takes values in \code{c("iJRF", "iRafNet","ptmJRF")} }

\item{importance}{A matrix containing importance score. When function \code{iRafNet}, this is a two dimensional matrix (\code{p} x \code{num.targets}) with \code{num.targets} being the number of targets considered for this parallel batch and \code{p} the number of genes. When function \code{iJRF} or \code{ptmJRF} is implemented, this is a three dimensional matrix of importance scores (\code{p} x \code{num.targets} x \code{C}) with \code{num.targets} being the number of targets considered in each batch, \code{p} the total number of genes/proteins and \code{C} the number of classes.}
}
\description{
This function computes importance score in parallel for a subset of target genes. 
}
\examples{

 # --- Generate data sets
 nclasses=2               # number of data sets / classes
 n1<-n2<-20               # sample size for each data sets
 p<-5                   # number of variables (genes)
 genes.name<-paste("G",seq(1,p),sep="")   # genes name
 
 data1<-matrix(rnorm(p*n1),p,n1)       # generate data1
 data2<-matrix(rnorm(p*n2),p,n1)       # generate data2
 
 # --- Run moultiple jobs and combine them for each function
  
   # -- function iJRF
   out.new<-array(0,c(p,p,nclasses))
   n.var=0
   for (k in 1:3){ 
      out<-iJRFNet_parallel(X=list(data1,data2),genes.name=genes.name,
      model="iJRF",parallel=c(k,2))
      
      n.target<-dim(out$importance)[2]
      for (c in 1:nclasses) {
      out.new[,seq(n.var+1,n.var+n.target),c]<-out$importance[,,c];}
      n.var=n.var+n.target
    }

   # -- function iRafNet
   W<-abs(matrix(rnorm(p*p),p,p))    # generate weights for interactions
   for (k in 1:3){ 
      out<-iJRFNet_parallel(X=list(data1),W=W,genes.name=genes.name,
      model="iRafNet",parallel=c(k,2))
      print(dim(out$importance))
      if (k==1) out.new<-out$importance
      if (k >2) out.new<-cbind(out.new,out$importance)
   }

    # -- function ptmJRF
    genes.name<-paste("G",seq(1,p),sep="")   # genes name
    ptm.name<-c("G1","G2","G3","G3","G4","G5","G1")   # ptm name
    p.ptm<-length(ptm.name)
 
    data1<-matrix(rnorm(p.ptm*n2),p.ptm,n1)       # generate PTM data
    data2<-matrix(rnorm(p*n1),p,n1)       # generate global proteomics data
 
    out.new<-array(0,c(p,p,nclasses)) # -- p x p matrix of importance scores
    n.var=0
    for (k in 1:3){ 

      out<-iJRFNet_parallel(X=list(data1,data2),genes.name=genes.name,
                     ptm.name=ptm.name,model="ptmJRF",parallel=c(k,2))
      
      n.target<-dim(out$importance)[2]
      for (c in 1:nclasses) {
      out.new[,seq(n.var+1,n.var+n.target),c]<-out$importance[,,c];}
      n.var=n.var+n.target
    }

  }
\references{
Petralia, F., Song, W.M., Tu, Z. and Wang, P. (2016). New method for joint network analysis reveals common and different coexpression patterns among genes and proteins in breast cancer. \emph{Journal of proteome research}, \bold{15}(3), pp.743-754.

A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. \emph{R News} \bold{2}, 18--22.
}
