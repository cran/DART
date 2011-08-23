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
netconsist.v[4] <- length(which(signedgeADJ.v[edges.idx]==signedgePRIOR.v[edges.idx]))/ne;  
netsignedge.m <- rbind(signedgePRIOR.v,signedgeADJ.v);
rownames(netsignedge.m) <- c("PRIOR","ADJ");

## estimate significance of fraction of consistent edges: derive p-value
## first estimate probability that a correlation is positive

tmp.idx <- sample(1:nrow(data.m),min(c(1000,nrow(data.m))),replace=FALSE);
tmpC.m <- cor(t(data.m[tmp.idx,]));
w.v <- summary(factor(sign(as.vector(tmpC.m[upper.tri(tmpC.m)]))));
w <- w.v[match(1,names(w.v))]/sum(w.v);


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

return(list(netcons=netconsist.v,netsign=netsignedge.m,adj=adj.m,s=sign.v,c=cor.m));
}

