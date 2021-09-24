\name{Hubgenes_barplot}
\alias{Hubgenes_barplot}
\title{Degree plot of hubgenes.}
\usage{
Hubgenes_barplot(net.final,net1=NULL,net2=NULL,genes.name,
                name.net1=NULL,name.net2=NULL,num.hub=NULL)
}
\arguments{
\item{net.final}{List object from function \code{Derive_network}. The \code{c}th element of the list contains estimated interaction for the \code{c}th class.}

\item{net1}{Optional integer. The network whose hubgenes must be plot. This is basically the element of the list \code{net.final} that needs to be considered. If omitted, it is set to 1.}

\item{net2}{Optional integer. This is a network used for comparison and is the element of the list \code{net.final} that needs to be considered for comparison. When \code{net2} is not omitted, the plot shows the degree of hubgenes in \code{net1} indicating for each of them how many edges are shared with \code{net2} and how many edges are \code{net1}- and \code{net2}- specific.}

\item{genes.name}{A vector containing gene names.}

\item{name.net1}{Optional string containing the name of \code{net1}.}

\item{name.net2}{Optional string containing the name of \code{net2}.}

\item{num.hub}{Number of top hubgenes the function will plot.}
}

\value{
Degree plot of top \code{num.hub} hub-genes in \code{net1}. When \code{net2} is not omitted, for each hub-gene in \code{net1}, the plot will show the number of connecting edges shared with \code{net2} and the number of connecting edges that are \code{net1} and \code{net2} specific. 
}
\description{
This function returns a degree plot for the top hub-genes in a given network.  
}
\examples{

   # --- Generate data sets
 nclasses=2               # number of data sets / classes
 n1<-n2<-20               # sample size for each data sets
 p<-5                   # number of variables (genes/proteins)
 genes.name<-paste("G",seq(1,p),sep="")   # genes/proteins name
 
 data1<-matrix(rnorm(p*n1),p,n1)       # generate data1
 data2<-matrix(rnorm(p*n2),p,n1)       # generate data2

  # --- Run iJRF and obtain importance score of interactions
  out.iJRFNet<-iJRFNet(X=list(data1,data2),genes.name=genes.name,
                      model="iJRF")

  # --- Obtain importance scores for M permuted data sets
  out.perm<-iJRFNet_permutation(X=list(data1,data2), ntree=1000,
              mtry=sqrt(5),genes.name=genes.name,M=5,model="iJRF")

  # --- Derive final networks
  net.final<-Derive_network(out.iJRFNet,out.perm,0.001)

  # --- Degree plot
  final.net<-Hubgenes_barplot(net.final,genes.name=genes.name)

}
\references{
Petralia, F., Song, W.M., Tu, Z. and Wang, P. (2016). New method for joint network analysis reveals common and different coexpression patterns among genes and proteins in breast cancer. \emph{Journal of proteome research}, \bold{15}(3), pp.743-754.

A. Liaw and M. Wiener (2002). Classification and Regression by randomForest. \emph{R News} \bold{2}, 18--22.

Xie, Y., Pan, W. and Khodursky, A.B., 2005. A note on using permutation-based false discovery rate estimates to compare different analysis methods for microarray data. \emph{Bioinformatics}, \bold{21}(23), pp.4280-4288.
}
