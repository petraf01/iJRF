\name{Derive_network}
\alias{Derive_network}
\title{Compute permutation-based FDR of importance scores and return estimated interactions.}
\usage{
Derive_network(out.iJRFNet,out.perm,TH)
}
\arguments{
\item{out.iJRFNet}{Output from object of class \code{iJRFNet}.}

\item{out.perm}{Output from object of class \code{Run_permutation}.}

\item{TH}{Threshold for FDR.}
}

\value{
List of estimated interactions. 
}
\description{
This function computes permutation-based FDR of importance scores and returns gene-gene interactions.  
}
\examples{


# --- Generate data sets
 nclasses=2               # number of data sets / classes
 n1<-n2<-20               # sample size for each data set
 p<-5                   # number of variables (genes/proteins)
 genes.name<-paste("G",seq(1,p),sep="")   # genes name

 data1<-matrix(rnorm(p*n1),p,n1)       # generate data1
 data2<-matrix(rnorm(p*n2),p,n1)       # generate data2

  # --- Run iJRFNet and obtain importance score of interactions
  out.iJRFNet<-iJRFNet(X=list(data1,data2),genes.name=genes.name,model="iJRF")

 # --- Obtain importance scores for M permuted data sets
  out.perm<-iJRFNet_permutation(X=list(data1,data2), ntree=1000,mtry=sqrt(5),
  genes.name=genes.name,M=5,model="iJRF")

  # --- Derive final networks
  final.net<-Derive_network(out.iJRFNet,out.perm,0.001)

}
\references{
Petralia, F., Song, W.M., Tu, Z. and Wang, P. (2016). New method for joint network analysis reveals common and different coexpression patterns among genes and proteins in breast cancer. \emph{Journal of proteome research}, \bold{15}(3), pp.743-754.

A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. \emph{R News} \bold{2}, 18--22.

Xie, Y., Pan, W. and Khodursky, A.B., 2005. A note on using permutation-based false discovery rate estimates to compare different analysis methods for microarray data. \emph{Bioinformatics}, \bold{21}(23), pp.4280-4288.
}
