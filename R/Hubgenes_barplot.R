Hubgenes_barplot <- function(net.final,net1=NULL,net2=NULL,genes.name,
                        name.net1=NULL,name.net2=NULL,num.hub=NULL) {

  p<-length(genes.name)
  if (is.null(name.net1)) name.net1<-"Network 1"
  if (is.null(name.net2)) name.net2<-"Network 2"
  if (length(net.final)==1 | is.null(net1)) net1=1
  if (is.null(num.hub)) num.hub=min(20,p)
  
  net1.edges<-net.final[[net1]]

  mat.net1<-mat.net2<-matrix(0,p,p)
  m.g1<-match(net1.edges[,1],genes.name);   m.g2<-match(net1.edges[,2],genes.name)
  cond<-(is.na(m.g1)==FALSE & is.na(m.g2)==FALSE)
  mat.net1[(m.g1[cond]-1)*p+m.g2[cond]]<-mat.net1[(m.g2[cond]-1)*p+m.g1[cond]]<-1
  net1.specific<-colSums(mat.net1==1)
  i=sort(net1.specific,decreasing=TRUE,index=TRUE)
 
 Degree<-Gene<-Edges<-NULL
 d<-data.frame(Degree=net1.specific[i$ix[seq(1,num.hub)]],Gene=genes.name[i$ix[seq(1,num.hub)]])
  
  if (is.null(net2)==FALSE) { # -- compare two networks
    net2.edges<-net.final[[net2]]
    m.g1<-match(net2.edges[,1],genes.name);   m.g2<-match(net2.edges[,2],genes.name)
    cond<-(is.na(m.g1)==FALSE & is.na(m.g2)==FALSE)
    mat.net2[(m.g1[cond]-1)*p+m.g2[cond]]<-mat.net2[(m.g2[cond]-1)*p+m.g1[cond]]<-1 
  
    net1.specific<-colSums(mat.net1==1 & mat.net2==0)
    shared<-colSums(mat.net2==1 & mat.net1==1)
    net2.specific<-colSums(mat.net1==0 & mat.net2==1)
    
    tot.non.zero<-sum(net1.specific+shared)
    if (num.hub>tot.non.zero) num.hub=tot.non.zero  
    d<-data.frame(Degree=c(shared[i$ix[seq(1,num.hub)]],net1.specific[i$ix[seq(1,num.hub)]],
     net2.specific[i$ix[seq(1,num.hub)]]),Edges=factor(c(rep("Shared",num.hub),
     rep(name.net1,num.hub),rep(name.net2,num.hub))),Gene=factor(rep(genes.name[i$ix[seq(1,num.hub)]],3)))

}
  
  d$Gene <- factor(d$Gene, levels=genes.name[i$ix])
 
 if (is.null(net2)==FALSE)  {ggplot(d, aes(x = Gene,Degree,fill=Edges)) + geom_bar(stat = "identity",col="black")+
    scale_fill_manual(values = c("red","blue","chartreuse3"))+theme_bw()+
    theme(legend.position="top",axis.text.x=element_text(angle=-90),
          text = element_text(size=15),
          legend.text= element_text(size=10),  
          axis.text = element_text(colour = "black"),legend.direction="horizontal")+labs(fill="")+
    coord_flip()            }  else {
      
      ggplot(d, aes(x = Gene,Degree)) + geom_bar(stat = "identity",fill="red",col="black")+theme_bw()+
        theme(legend.position="top",axis.text.x=element_text(angle=-90),
      text = element_text(size=15),
      legend.text= element_text(size=10),  
      axis.text = element_text(colour = "black"),legend.direction="horizontal")+labs(fill="")+
      coord_flip()  
      
      
    }
}