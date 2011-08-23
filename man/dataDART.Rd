\name{dataDART}
\alias{dataDART}
\docType{data}
\title{Example data for DART package}
\description{
  This data object consists of a gene expression data matrix over 67 oestrogen receptor negative breast cancers, of which 51 are basals and 16 are HER2+/ERBB2+. This classification is based on the intrinsic subtype classifier (see Hu Z,...et.al, BMC Genomics. 2006 Apr 27;7:96), and while not equivalent to IHC or the amplification status at the ERBB2 locus, it broadly matches these other classifications. The model pathway signature is an in-vitro derived perturbation signature reflecting ERBB2 activation (see Creighton CJ,...et.al, Cancer Res 2006 Apr 1;66(7):3903-11). Both data and signature are annotated with Entrez Gene IDs.
}
\usage{dataDART}
\format{
A list with 3 elements: (i) data: is the gene expression matrix described above, (ii) sign: is the ERBB2 perturbation signature, (iii) pheno: basal/ERBB2 status of samples.
}