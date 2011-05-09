EvalConsNet <-
function(buildRN.o){

#######################################
sign.v <- buildRN.o$s;
adj.m <- buildRN.o$adj;
cor.m <- buildRN.o$c;
data.m <- buildRN.o$d;
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
netconsist.v[4] <- length(which(signedgeADJ.v[edges.idx]==signedgePRIOR.v[edges.idx]))/ne;  #AET original one
#netconsist.v[4] <- length(which(signedgeADJ.v[edges.idx]==signedgePRIOR.v[edges.idx]))/(0.5*ng*(ng-1)) ;      #yan modified
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
#round(w*length(edges.idx))->w1     #yan
#length(edges.idx)-w1->w2        #yan
 permedgesADJ.v[edges.idx] <- c(-1,1)[rbinom(length(edges.idx),1,w)+1];#AET original one
  #permedgesADJ.v[edges.idx] <- sample(c(rep(1,w1),rep(-1,w2)));
  netconsist <- length(which(permedgesADJ.v[edges.idx]==signedgePRIOR.v[edges.idx]))/ne;   # AET original one
  #netconsist <- length(which(permedgesADJ.v[edges.idx]==signedgePRIOR.v[edges.idx]))/(0.5*ng*(ng-1)) ;    #yan modified
  if(netconsist > netconsist.v[4]){   #AET original one
   #if(netconsist > w){
   count <- count+1;
  }
# print(netconsist)
}
netconsist.v[5] <- count/np;

return(list(netcons=netconsist.v,netsign=netsignedge.m,adj=adj.m,s=sign.v,c=cor.m));
}

