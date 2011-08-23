PredActScore <- function(pradjMC.m,signMC.v,data.m){

  if(length(which(is.na(data.m))) > 0){
    print("missing data in data.m must be imputed before proceeding");
  }
  ### next standardize rows of data matrix
  StdRow <- function(tmp.v){
    out.v <- (tmp.v-mean(tmp.v))/sqrt(var(tmp.v));
    return(out.v);
  }
  stdata.m <- data.m;
  for(r in 1:nrow(data.m)){
    stdata.m[r,] <- StdRow(data.m[r,]);
  }
  
    
  match(rownames(pradjMC.m),rownames(stdata.m)) -> idx;
  which(is.na(idx)) -> idx.na
  naF <- round(length(idx.na)/length(idx),2)
  print(paste("Found ",100-100*naF,"% of maximally connected pruned network genes in the data",sep=""));
  
  if (naF==1){
   print("no gene overlap between new data set and pruned network, can't perform prediction")
  }
  else {

  setdiff(1:length(idx),idx.na) -> idx.sel;
  tmp.m <- stdata.m[idx[idx.sel],];
  tmpA.m <- pradjMC.m[idx.sel,idx.sel];
  k.v <- apply(tmpA.m,1,sum);  
  genescores.m <- k.v*sign(signMC.v[idx.sel])*tmp.m;
  score.v <- apply(genescores.m,2,sum)/sqrt(sum(k.v*k.v));
  

  }


  return(list(adj=tmpA.m,sign=signMC.v[idx.sel],score=score.v,degree=k.v));
}
