Plot_Modules<-function(net.final,genes.name) {
  
  p<-length(genes.name)
  truth<-matrix(0,p,p)
  m.g1<-match(net.final[,1],genes.name);   m.g2<-match(net.final[,2],genes.name)
  cond<-(is.na(m.g1)==FALSE & is.na(m.g2)==FALSE)
  truth[(m.g1[cond]-1)*p+m.g2[cond]]<-truth[(m.g2[cond]-1)*p+m.g1[cond]]<-1
  
  colnames(truth)<-rownames(truth)<-genes.name
  new.j<-graph.adjacency(truth, mode="undirected")
  
  l <-layout.fruchterman.reingold(new.j)
  
  # -- remove modules with less than 4 vertices
  C<-clusters(new.j, mode=c("weak"))
  bad.vs<-V(new.j)[C$csize[C$membership]<4]     # remove components with less than 4 vertices
  new.j<-delete.vertices(new.j, bad.vs) # exclude them from the graph
  l <-layout.fruchterman.reingold(new.j)
  
  # -- derive modules
  all_eb<-edge.betweenness.community(new.j,directed=F)
  cl<-membership(all_eb)
  u.cl<-unique(cl)
  colors = rainbow(max(cl))
  
  plot(new.j, vertex.color=colors[cl],edge.width=.5,vertex.size=3,vertex.frame.color="white",layout=l,
       vertex.label=V(new.j)$name,    #specifies the lables of the vertices. in this case the 'name' attribute is used
       vertex.label.cex=.1,vertex.label.color="black",vertex.label.font=2) # 
 
  out<-cbind(as.character(V(new.j)$name),as.data.frame(as.numeric(cl)),stringsAsFactors=FALSE)
  colnames(out)<-c("Vertex","Cluster Membership")
  return(out)
}