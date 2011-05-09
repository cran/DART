\name{MaxCon}
\alias{MaxCon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructing maximally connected pruned network 
}
\description{
Construction of relevance correlation pruned network as the maximally connected component where all edges reflect correlations that are consistent with the prior information given by perturbation signature. Library igraph is require for running this function
}
\usage{
MaxCon(PruneNet.o)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PruneNet.o}{include a adjacency matrix and gene signature vector, they are from running DoDART function}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{pradj}{maximally connected pruned network}
 \item{sign}{gene signature}
%% ...
}
\references{
{Jiao Y, Lawler K, Patel GS, Purushotham A, Jones AF, Grigoriadis A, Ng T, Teschendorff AE. Denoising algorithm based on relevance network topology improves molecular pathway activity inference. Submitted.}

{Teschendorff AE, Gomez S, Arenas A, El-Ashry D, Schmidt M, et al. (2010) Improved prognostic classification of breast
cancer defined by antagonistic activation patterns of immune response pathway modules. BMC Cancer 10:604.}

}
\author{Yan Jiao,Andrew E Teschendorff}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## The function is currently defined as
function(PruneNet.o){
 PruneNet.o$pradj->pradj.1
 p.sign1<-sign(PruneNet.o$sign)
 graph.adjacency(pradj.1, mode="undirected",weight=NULL)->g1.n
 c.n <-clusters(g1.n)
which(c.n$membership==(which.max(c.n$csize)-1))->pradj.idx
#print(length(pradj.idx))
 pradj.max<-pradj.1[pradj.idx,pradj.idx]
 p.sign<-p.sign1[pradj.idx]
 pradj.max[lower.tri( pradj.max)]<-0
 return(list(pradj=pradj.max,sign=p.sign));
 } ### END of FUNCTION
 
 ###example 
  library(DARTData)
  data(data.m)
  data(sign.v)
### do DART
  PruneNet.o <- DoDART(data.m,sign.v,fdr=0.05) 
### do MaxCon
  MaxCon.o <- MaxCon(PruneNet.o)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
