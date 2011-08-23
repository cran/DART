\name{PruneNet}
\alias{PruneNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prunes relevance network to allow only edges that are consistent with the predictions of the model signature
}
\description{
Prunes relevance network to allow only edges that are consistent with the predictions of the model signature, and extracts the maximally connected component. This is the denoising step in DART.
}
\usage{
PruneNet(evalNet.o)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{evalNet.o}{Output object from EvalConsNet}
}

\value{
  A list with following entries:
   \item{pradj}{The adjacency matrix of the pruned i.e consistent network.}
   \item{sign}{The model signature vector of genes in pruned network.}
   \item{score}{The fraction of edges surviving the pruning/denoising.}
   \item{netconst}{Same output as for EvalConsNet.}
   \item{pradjMC}{The adjacency matrix of the maximally connected component of pruned network.}
   \item{signMC}{The model signature vector of the genes in the maximally connected component.}
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

###See Example for DoDART.
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
