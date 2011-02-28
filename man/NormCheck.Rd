\name{NormCheck}
\alias{NormCheck}
\title{
Normalization
}
\description{
Normalize the Beta Values and plot a PCA and a hierarchical Clustering of the samples using the noralized data
}

\usage{
NormCheck(Dir, platform, pval, ChrX, ClustMethod)
}

\arguments{
  \item{Dir}{A character string containing the location of the directory 
    in which the input file are}
  \item{platform}{The type of Illumina Infinium BeadChip methylation assay. This must be one of "Hum27" (Infinium HumanMethylation27 BeadChip) or "Hum450" (Infinium HumanMethylation450 BeadChip). WARNING: currently only "Hum27" is supported}
  \item{pval}{A p-value threshold to define which sample keep for the normalization and the following analyses}
  \item{ClustMethod}{The distance measure to be used for the clustering. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "correlation", "spearman" or "kendall"}
  \item{ChrX}{Logical value indicating whether the CpGs that belongs to chromosome X should be deleted from normalization and the following analyses. The default is FALSE}
}
\details{
Creates:
\itemize{
\item  a PCA of normalized Beta values
\item  a Cluster of normalized Beta values
\item  an Excel file with the normalized M-values (see lumi package for details).
}}
\author{ F.M. Mancuso <francesco.mancuso@crg.es> }
\examples{
\dontrun{
	Dir <- system.file("extdata/",package="HumMeth27QCReport")
	NormCheck(Dir,pva=0.05,ClustMethod="euclidean",ChrX=F)
}
}
\keyword{ hplot }
\keyword{ manip }
