MaxCon <-
function(PruneNet.o){
library(igraph)
 PruneNet.o$pradj->pradj.1
 p.sign1<-sign(PruneNet.o$sign)
 graph.adjacency(pradj.1, mode="undirected",weight=NULL)->g1.n
 c.n <-clusters(g1.n)
which(c.n$membership==(which.max(c.n$csize)-1))->pradj.idx
#print(length(pradj.idx))
 pradj.max<-pradj.1[pradj.idx,pradj.idx]
 p.sign<-p.sign1[pradj.idx]
 pradj.max[lower.tri( pradj.max)]<-0
 return(list(pradj=pradj.max,sign=p.sign));
 } ### END of FUNCTION

