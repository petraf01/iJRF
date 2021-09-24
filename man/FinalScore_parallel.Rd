\name{FinalScore_parallel}
\alias{FinalScore_parallel}
\title{Derive final importance scores for object of class iJRFNet_parallel.}
\usage{
FinalScore_parallel(importance, model, genes.name)}
\arguments{
\item{importance}{A matrix containing importance scores. When model \code{iRafNet} is implemented, \code{importance} is a two dimensional matrix (\code{p} x \code{p}) with \code{p} being the total number of genes/proteins. When either function \code{iJRF} or \code{ptmJRF} is implemented, \code{importance} is a three dimensional matrix of importance scores (\code{p} x \code{p} x \code{C}) with \code{p} being the total number of genes/proteins and \code{C} the number of classes.}

\item{model}{Variable indicating which \code{iJRFNet} model will be implemented. Takes values in \code{c("iJRF", "iRafNet", "ptmJRF")} }

\item{genes.name}{Vector containing genes name. The order needs to match the rows/columns of \code{importance}.}

}
\value{
A matrix with \code{I} rows and \code{C + 2} columns where \code{I} is the total number of gene-gene (protein-protein) interactions and \code{C} is the total number of classes. The first two columns contain gene names for each interaction while the remaining columns contain importance scores for different classes. When model \code{iRafNet} is implemented, the number of classes is \code{1} and therefore only three columns will be returned.
}
\description{
This function returns importance score for each gene-gene (protein-protein) interaction.
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
  ## ---  Run iJRFNet 

  ## --- Run multiple jobs in parallel and combine them
  
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
     
   ## --- Derive interactions 
   FinalScore_parallel(importance=out.new,model="iJRF",genes.name=genes.name)



}    
\references{
Petralia, F., Song, W.M., Tu, Z. and Wang, P. (2016). New method for joint network analysis reveals common and different coexpression patterns among genes and proteins in breast cancer. \emph{Journal of proteome research}, \bold{15}(3), pp.743-754.

A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. \emph{R News} \bold{2}, 18--22.
}
