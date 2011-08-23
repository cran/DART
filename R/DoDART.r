DoDART <-
function(data.m,sign.v,fdr){

  rn.o <- BuildRN(data.m,sign.v,fdr);
  evalNet.o <- EvalConsNet(rn.o);
  prnet.o <- PruneNet(evalNet.o);

  pred.o <- PredActScore(prnet.o$pradjMC,prnet.o$signMC,data.m);

  return(list(netcons=prnet.o$netcons,adj=pred.o$adj,sign=pred.o$sign,score=pred.o$score,degree=pred.o$degree));
  
 
} ### END of FUNCTION

