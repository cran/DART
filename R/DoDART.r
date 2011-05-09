DoDART <-
function(data.m,sign.v,fdr){

  ### find signature genes in data matrix
  match(names(sign.v),rownames(data.m)) -> map.idx;
  rep.idx <- which(is.na(map.idx)==FALSE);
  print(paste("Found ",100*round(length(rep.idx)/length(sign.v),2),"% of signature genes in data matrix",sep=""));
  tmp.m <- data.m[map.idx[rep.idx],];
  ### compute correlations between signature genes
  cor.m <- cor(t(tmp.m),method="pearson",use="pairwise.complete.obs");
  ### set trivial correlations to zero
  diag(cor.m) <- 0;
  z.m <- 0.5*log( (1+cor.m)/(1-cor.m) );
  std <- 1/sqrt(ncol(tmp.m)-3);
  pv.m <-  2*pnorm(abs(z.m),0,std,lower.tail=FALSE);
  pv.v <- as.vector(pv.m[upper.tri(pv.m)]);
  ### now correct for multiple testing
  qv.v <- p.adjust(pv.v,method="BH");
    qv.m <- pv.m;
  qv.m[cor.m<3] <- 0
 qv.m[upper.tri(qv.m)] <- qv.v;
 qv.m<-qv.m+t(qv.m)
  ### define relevance network
  adj.m <- cor.m;
adj.m[cor.m<3] <- 0; ### a trick to set all elements to zero
  adj.m[qv.m < fdr] <- 1;
  ### test symmetry
  if(identical(adj.m,t(adj.m))==FALSE){
    print("PROBLEM:adjaceny matrix is not symmetric");
  }
   diag(adj.m)<-0
  #return(list(adj=adj.m,s=sign.v[rep.idx],sd=tmp.m,c=cor.m,d=data.m, rep.idx= rep.idx));


sign.v <-sign.v[rep.idx]
adj.m <- adj.m

########################################

netconsist.v <- vector(length=5);
names(netconsist.v) <- c("nG","nE","fE","fconsE","Pval(consist)");

ng <- length(sign.v);
netconsist.v[1] <- ng;
signPRIOR.v <- sign(sign.v);
signedgePRIOR.v <- rep(-1,0.5*ng*(ng-1));
ie <- 1;
for(n1 in 1:(ng-1)){
 for(n2 in (n1+1):ng){
   if(signPRIOR.v[n1]==signPRIOR.v[n2]){
     signedgePRIOR.v[ie] <- 1;
   }
   ie <- ie+1;
 }
}

signedgeADJ.v <- rep(NA,0.5*ng*(ng-1));
ie <- 1;
for(n1 in 1:(ng-1)){
 for(n2 in (n1+1):ng){
   if(adj.m[n1,n2]==1){
    signedgeADJ.v[ie] <- sign(cor.m[n1,n2]);
   }
   ie <- ie+1;
 }
}
edges.idx <- which(is.na(signedgeADJ.v)==FALSE);
ne <- length(edges.idx);
netconsist.v[2] <- ne;
netconsist.v[3] <- ne/(0.5*ng*(ng-1));
netconsist.v[4] <- length(which(signedgeADJ.v[edges.idx]==signedgePRIOR.v[edges.idx]))/ne;  

netsignedge.m <- rbind(signedgePRIOR.v,signedgeADJ.v);
rownames(netsignedge.m) <- c("PRIOR","ADJ");
## estimate significance of fraction of consistent edges: derive p-value
## randomly assign directionality values depending on data matrix
tmp.idx <- sample(1:nrow(data.m),min(c(1000,nrow(data.m))),replace=FALSE);
tmpC.m <- cor(t(data.m[tmp.idx,]));
w.v <- summary(factor(sign(as.vector(tmpC.m[upper.tri(tmpC.m)]))));
if(length(w.v)==1){
  w <- 1;
}else {
  w <- w.v[2]/sum(w.v);
}

permedgesADJ.v <- rep(NA,0.5*ng*(ng-1));
np <- 1000;
count <- 0;
for(p in 1:np){

 permedgesADJ.v[edges.idx] <- c(-1,1)[rbinom(length(edges.idx),1,w)+1];

  netconsist <- length(which(permedgesADJ.v[edges.idx]==signedgePRIOR.v[edges.idx]))/ne;   

  if(netconsist > netconsist.v[4]){   

   count <- count+1;
  }

}
netconsist.v[5] <- count/np;



# PruneNet.R

### DESCRIPTION
### Evaluates consistency of inferred networks with prior in vitro information

### INPUT
### evalNet.o : output object from EvalConsNet.R

### OUTPUT
### pradj: pruned network
### sign: gene signature

#PruneNet <- function(evalNet.o){
#
#######################################
#adj.m <- evalNet.o$adj;
#netsignedge.m <- evalNet.o$netsign;
#sign.v <- evalNet.o$s;
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
diag(pradj.m)<-0


(sum(pradj.m))/sum(adj.m)->  consist.score



pradj.m->pradj.1
 p.sign1<-sign(sign.v)
 graph.adjacency(pradj.1, mode="undirected",weight=NULL)->g1.n
 c.n <-clusters(g1.n)
which(c.n$membership==(which.max(c.n$csize)-1))->pradj.idx
#print(length(pradj.idx))
 pradj.max<-pradj.1[pradj.idx,pradj.idx]
 p.sign<-p.sign1[pradj.idx]
 pradj.max[lower.tri( pradj.max)]<-0

#return(list(pradj=pradj.m,sign=sign.v,consist.score=consist.score));
return(list(pradj=pradj.max,netcons=netconsist.v,netsign=netsignedge.m,sign=p.sign,c=cor.m,consist.score=consist.score));
} ### END of FUNCTION

