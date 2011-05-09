\name{PredActScore}
\alias{PredActScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predicting activity score given test data and pruned network from running DART
}
\description{
Given the adjacency matrix and sign of the perturbation signature(both could be the output of function MaxCon) and test dataset, one could use PredActScore to obtain activity score for test dataset. 
}
\usage{
PredActScore(pradj.m, sign.v, test.data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pradj.m}{
Adjacency matrix for calculating score
}
  \item{sign.v}{
Gene signatures
}
  \item{test.data}{
Test data set
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{pradj}{Adjacency matrix for calculating score}
  \item{sign}{Gene signatures}
  \item{score}{Score for test data set}
%% ...
}
\references{
{Jiao Y, Lawler K, Patel GS, Purushotham A, Jones AF, Grigoriadis A, Ng T, Teschendorff AE. Denoising algorithm based on relevance network topology improves molecular pathway activity inference. Submitted.}

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
function(pradj.m,sign.v,test.data){

  match(as.numeric(rownames( pradj.m)),as.numeric(rownames(test.data)))->idx
  which(is.na(idx)) ->idx.na
  gene.percent<-round(length(idx.na)/length(idx),2)

  print(paste("Found ",100-100*gene.percent,"\%of signature genes in test data set"));
  if (gene.percent==1)
  {
  print("no gene overlap between test data set and prunned network, can't perform prediction")
  }else{
  setdiff(c(1:length(idx)),idx.na)->idx.sel
  test.data.sel<-test.data[idx[idx.sel],]
  pradj.m[idx.sel,idx.sel]-> pradj.sel
  score<- (colSums((pradj.sel)\%*\% ((sign.v[idx.sel])*test.data.sel)))
  n.node<-dim(pradj.sel)[1]-length(which( colSums(pradj.sel)==0))
  edge.n<-colSums(pradj.m)
  score<-(sqrt(n.node)*score)/(sqrt(sum(edge.n*edge.n)))
  }
  
return(list(pradj=pradj.m,sign=sign.v[idx.sel],score=score));


  }
  

### Example
### load in example data
library(DARTData)
data(data.m)
data(sign.v)

### do DART

dart.o <- DoDART(data.m,sign.v,fdr=0.05)

maxcon.o <- MaxCon(dart.o)
score<-PredActScore(maxcon.o$pradj, maxcon.o$sign, data.m)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
