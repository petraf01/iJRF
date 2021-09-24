\name{Plot_Modules}
\alias{Plot_Modules}
\title{Derive Network Modules.}
\usage{
Plot_Modules(net.final,genes.name)
}
\arguments{
\item{net.final}{Network to plot.}
\item{genes.name}{A vector containing gene names.}
}

\value{Return list of modules.}
\description{
This function returns the list of modules.  
}
\examples{

 # --- Generate data sets
 nclasses=2               # number of data sets / classes
 n1<-n2<-20               # sample size for each data sets
 p<-40                   # number of variables (genes/proteins)
 genes.name<-paste("G",seq(1,p),sep="")   # genes/proteins name
 
 data1<-matrix(rnorm(p*n1),p,n1)       # generate data1
 data2<-matrix(rnorm(p*n2),p,n1)       # generate data2

  # --- Run iJRF and obtain importance score of interactions
  out.iJRFNet<-iJRFNet(X=list(data1,data2),genes.name=genes.name,
                      model="iJRF")

  # --- Degree plot
  final.net<-Plot_Modules(out.iJRFNet[sample(dim(out.iJRFNet)[1],200),c(1,2)]
               ,genes.name=genes.name)

}
\references{
Petralia, F., Song, W.M., Tu, Z. and Wang, P. (2016). New method for joint network analysis reveals common and different coexpression patterns among genes and proteins in breast cancer. \emph{Journal of proteome research}, \bold{15}(3), pp.743-754.

A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. \emph{R News} \bold{2}, 18--22.

Xie, Y., Pan, W. and Khodursky, A.B., 2005. A note on using permutation-based false discovery rate estimates to compare different analysis methods for microarray data. \emph{Bioinformatics}, \bold{21}(23), pp.4280-4288.
}
