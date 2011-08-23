PruneNet <-
function(evalNet.o){

######################################
adj.m <- evalNet.o$adj;
netsignedge.m <- evalNet.o$netsign;
sign.v <- evalNet.o$s;
netcons.v <- evalNet.o$netcons
######################################

## first need to define inverse map
ng <- nrow(adj.m);
imap.m <- matrix(nrow=2,ncol=0.5*ng*(ng-1));
ie <- 1;
for(g1 in 1:(ng-1)){
 for(g2 in (g1+1):ng){
  imap.m[,ie] <- c(g1,g2);
  ie <- ie+1;
 }
}
pruneE.idx <- which(netsignedge.m[1,]!=netsignedge.m[2,]);
pradj.m <- adj.m;
for(e in pruneE.idx){
  rowe <- imap.m[1,e];
  cole <- imap.m[2,e];
  pradj.m[rowe,cole] <- 0;
  pradj.m[cole,rowe] <- 0;
}


score <- sum(pradj.m)/sum(adj.m);

library(igraph);
gr.o <- graph.adjacency(pradj.m,mode="undirected");
clust.o <-clusters(gr.o);
maxc.idx <- which(clust.o$membership==(which.max(clust.o$csize)-1));
pradjMC.m <- pradj.m[maxc.idx,maxc.idx];
signMC.v <- sign.v[maxc.idx];

 return(list(pradj=pradj.m,sign=sign.v,score=score,netcons=netcons.v,pradjMC=pradjMC.m,signMC=signMC.v));
}

