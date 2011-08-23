\name{PredActScore}
\alias{PredActScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the DART activation score of the model signature in the samples of a data set.
}
\description{
This function computes the DART activation score of the model signature for the samples of a data matrix. This data matrix can be the same data matrix in which DART was applied or it could be a totally independent data set.
}
\usage{
PredActScore(pradjMC.m, signMC.v, data.m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pradjMC.m}{The adjacency matrix of the maximally connected pruned network inferred using DART.}
  \item{signMC.v}{The corresponding model signature vector.}
  \item{data.m}{A data matrix with rows labeling genes/probes and columns labeling samples in which the signature activity scores are desired.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A list with following elements:
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{adj}{Adjacency matrix of the DART network for the genes found in the data set given by \code{data.m}.}
  \item{sign}{The corresponding model signature vector for the genes found in the data set.}
  \item{score}{A vector of the signature/pathway activity scores for each sample in the data set.}
%% ...
}
\references{
{Jiao Y, Lawler K, Patel GS, Purushotham A, Jones AF, Grigoriadis A, Ng T, Teschendorff AE. Denoising algorithm based on relevance network topology improves molecular pathway activity inference. Submitted.}

{Teschendorff AE, Gomez S, Arenas A, El-Ashry D, Schmidt M, et al. (2010) Improved prognostic classification of breast cancer defined by antagonistic activation patterns of immune response pathway modules. BMC Cancer 10:604.}

}
\author{Andrew E Teschendorff, Yan Jiao}

\examples{

###See example for DoDART.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
