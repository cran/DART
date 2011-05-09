\name{BuildRN}
\alias{BuildRN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Building relevance network
}
\description{This function builds a relevance network over a given data set and model pathway signature}
\usage{
BuildRN(data.m, sign.v, fdr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{data.m}{Data matrix: Rows label features, columns label samples. It is assumed that number of features is much larger than number of samples.}
  \item{sign.v}{Signature vector: elements are numeric, names of sign.v are gene names (must match rows of data.m)}
  \item{fdr}{Desired false discovery rate}
}

\value{
  A list with following entries:
  
  \item{adj}{Inferred relevance network adjacency matrix}
  \item{s}{Gene signature in data}
  \item{sd}{Gene signature data matrix}
  \item{c}{Correlations between signature genes}
  \item{d}{Data matrix}
  \item{rep.idx}{Index of the gene signature which could be found in data matrix}

}
\references{{Jiao Y, Lawler K, Patel GS, Purushotham A, Jones AF, Grigoriadis A, Ng T, Teschendorff AE. Denoising algorithm based on relevance network topology improves molecular pathway activity inference. Submitted.}

{Teschendorff AE, Gomez S, Arenas A, El-Ashry D, Schmidt M, et al. (2010) Improved prognostic classification of breast
cancer defined by antagonistic activation patterns of immune response pathway modules. BMC Cancer 10:604.}


}
\author{Andrew E Teschendorff, Yan Jiao}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{



## The function is currently defined as
function(data.m,sign.v,fdr=0.05){

  ### find signature genes in data matrix
  match(names(sign.v),rownames(data.m)) -> map.idx;
  rep.idx <- which(is.na(map.idx)==FALSE);
  print(paste("Found ",100*round(length(rep.idx)/length(sign.v),2),"\% of signature genes in data matrix",sep=""));
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
  return(list(adj=adj.m,s=sign.v[rep.idx],sd=tmp.m,c=cor.m,d=data.m, rep.idx= rep.idx));

  }### END OF FUNCTION


### Example
### load in example data
library(DARTData)
data(data.m)
### load in signature 
data(sign.v)
### do BuildRN
BuildRN.o <- BuildRN(data.m,sign.v,fdr=0.05)


  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
