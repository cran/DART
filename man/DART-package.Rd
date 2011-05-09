\name{DART-package}
\alias{DART}
\docType{package}
\title{
Denoising Algorithm based on Relevance network Topology
}
\description{
(DART) Denoising Algorithm based on Relevance network Topology is an unsupervised algorithm which evaluates the consistency of a model pathway signature in molecular tumour profiles before estimating molecular pathway activity. 
This has been shown to improve estimates of pathway activity in clinical tumour specimens. 
Examples of model pathway signatures could be an in-vitro derived perturbation signature, or a pathway model of signal transduction, or a curated list of genes predicted to be up or down regulated in response to pathway activation/inhibition. 
Five internal functions implement the steps in DART and are provided as explicit functions to allow flexibility. 
DoDART is the main user function which will automatically and sequentially run through the internal functions. 
(1) BuildRN: this builds a relevance correlation network over the model pathway signature in the data set in which the pathway activity estimate is desired. 
(2) EvalConsNet: this evaluates the consistency of the inferred network with the prior information in the model pathway signature. Only if the consistency score is higher than expected by random chance (see Pval(consist) from output: netcons), can the model signature be used to infer pathway activity. 
(3) PruneNet: this obtains the pruned, i.e consistent, network, in which any edge represents a significant correlation in gene expression whose directionality agrees with that predicted by the prior information. 
(4) MaxCon: this function builds the maximally connected pruned network from the output of PruneNet, and this subnetwork is then used for pathway activity estimation using a metric that takes the degree distribution of the nodes (genes) in the subnetwork into account. 
Finally, given the adjacency matrix of this subnetwork and the corresponding model pathway signature, one can use the function PredActScore to obtain the pathway activity scores in an independent data set.
}
\details{
\tabular{ll}{
Package: \tab DART\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-05-09\cr
License:  \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}


\author{
Andrew E Teschendorff, Yan Jiao
Maintainer:  <y.jiao@ucl.ac.uk>
}

\references{
{Jiao Y, Lawler K, Patel GS, Purushotham A, Jones AF, Grigoriadis A, Ng T, Teschendorff AE. Denoising algorithm based on relevance network topology improves molecular pathway activity inference. Submitted.}

{Teschendorff AE, Gomez S, Arenas A, El-Ashry D, Schmidt M, et al. (2010) Improved prognostic classification of breast
cancer defined by antagonistic activation patterns of immune response pathway modules. BMC Cancer 10:604.}
}




\keyword{  network
pathway
consistency
}
\seealso{
}
\examples{

### Example

###In this example, the model pathway signature:sign.v is ERBB2, data.m is a subset of wang data set(ER- tumor samples),DoDART is applied,
###then learned network is used to predict ERBB2 pathway activity in the same samples, phenotypic information(Her2/Basal)is not used in the training process,
###boxplot shows that the network predicts pathway activity higher in ERBB2+ ER- samples which meets the expectation. 

### load in example data
library(DARTData)

data(data.m)
data(sign.v)

### doDART, which is the same as sequentially run BuildRN, then EvalConsNet, followed by PruneNet
###
dart.o <- DoDART(data.m,sign.v,fdr=0.05)

### predict activity score for test data set
score<-PredActScore(dart.o$pradj, dart.o$sign, data.m)  #using the same data set as test dataset
### Testing whether pathway activity higher in ERBB2+ ER- samples
idx.Basal<- which( colnames(data.m)=="Basal")
idx.Her2<- which( colnames(data.m)=="Her2")
boxplot(score$score[idx.Basal],score$score[idx.Her2],names=c("Basal","Her2"))

}