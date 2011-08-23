\name{BuildRN}
\alias{BuildRN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds the relevance correlation network
}
\description{This function builds the relevance correlation network for the genes in the model pathway signature in the given data set.}
\usage{
BuildRN(data.m, sign.v, fdr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{data.m}{Data matrix: Rows label features, columns label samples. It is assumed that number of features is much larger than number of samples. Rownames must be a valid gene or probe identifier.}
  \item{sign.v}{Model pathway signature vector: elements are numeric and corresponding to the regulatory weights, i.e the sign indicates if up or downregulated. Names of sign.v must be a gene name (probe) identifier which must match the gene (probe) identifier used for the rows of \code{data.m}.}
  \item{fdr}{Desired false discovery rate of significant edges in relevance correlation network.}
}

\value{
  A list with following entries:
  \item{adj}{Adjacency matrix of inferred relevance network}
  \item{s}{Model signature vector in data}
  \item{sd}{Gene signature data matrix}
  \item{c}{Correlations between signature genes}
  \item{d}{Data matrix}
  \item{rep.idx}{Indeces of the genes in signature which could be found in data matrix}

}
\references{{Jiao Y, Lawler K, Patel GS, Purushotham A, Jones AF, Grigoriadis A, Ng T, Teschendorff AE. Denoising algorithm based on relevance network topology improves molecular pathway activity inference. Submitted.}

{Teschendorff AE, Gomez S, Arenas A, El-Ashry D, Schmidt M, et al. (2010) Improved prognostic classification of breast cancer defined by antagonistic activation patterns of immune response pathway modules. BMC Cancer 10:604.}


}
\author{Andrew E Teschendorff, Yan Jiao}



\examples{

## See example for DoDART.

  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
