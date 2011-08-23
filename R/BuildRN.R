BuildRN <-
function(data.m,sign.v,fdr=0.05){

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
  qv.m[cor.m<3] <- 0; ## set all elements to zero
  qv.m[upper.tri(qv.m)] <- qv.v;
  qv.m <- qv.m+t(qv.m)
  diag(qv.m) <- 1;
  ### define relevance network
  adj.m <- cor.m;
  adj.m[cor.m<3] <- 0; ### set all elements to zero
  adj.m[qv.m < fdr] <- 1;
  ### test symmetry
  if(identical(adj.m,t(adj.m))==FALSE){
    print("PROBLEM:adjaceny matrix is not symmetric");
  }
  diag(adj.m)<- 0;
  return(list(adj=adj.m,s=sign.v[rep.idx],sd=tmp.m,c=cor.m,d=data.m,rep.idx=rep.idx));

}

