PredActScore<- function(pradj.m,sign.v,test.data){

  
  
  #pradj.m<-Dart.o$pradj
  #sign.v<-Dart.o$sign
  match(as.numeric(rownames( pradj.m)),as.numeric(rownames(test.data)))->idx
  which(is.na(idx)) ->idx.na
  gene.percent<-round(length(idx.na)/length(idx),2)
  #paste("Found ",100*round(length(rep.idx)/length(sign.v),2),"% of signature genes in data matrix",sep="")
  print(paste("Found ",100-100*gene.percent,"%of signature genes in test data set"));
  if (gene.percent==1)
  {
  print("no gene overlap between test data set and prunned network, can't perform prediction")
  }else{
  setdiff(c(1:length(idx)),idx.na)->idx.sel
  test.data.sel<-test.data[idx[idx.sel],]
  pradj.m[idx.sel,idx.sel]-> pradj.sel
 score<- (colSums((pradj.sel)%*% ((sign.v[idx.sel])*test.data.sel)))
 n.node<-dim(pradj.sel)[1]-length(which( colSums(pradj.sel)==0))
edge.n<-colSums(pradj.m)
 score<-(sqrt(n.node)*score)/(sqrt(sum(edge.n*edge.n)))
  }


return(list(pradj=pradj.m,s=sign.v[idx.sel],score=score));

}
