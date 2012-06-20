\name{NormCheck}
\alias{NormCheck}
\title{
Normalization
}
\description{
Normalize the Beta Values and plot a PCA and a hierarchical Clustering of the samples using the noralized data
}

\usage{
NormCheck(ImportDataR, platform, pval, ChrX, ClustMethod, normMethod)
}

\arguments{
  \item{ImportDataR}{the object that comes from ImportData function}
  \item{platform}{The type of Illumina Infinium BeadChip methylation assay. This must be one of "Hum27" (Infinium HumanMethylation27 BeadChip) or "Hum450" (Infinium HumanMethylation450 BeadChip)}
  \item{pval}{A p-value threshold to define which sample keep for the normalization and the following analyses}
  \item{ClustMethod}{The distance measure to be used for the clustering. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "correlation", "spearman" or "kendall"}
  \item{ChrX}{Logical value indicating whether the CpGs that belongs to chromosome X should be deleted from normalization and the following analyses. The default is FALSE}
  \item{normMethod}{Commonly either "quantile" or "ssn". Determines the normalization method. See documentation for lumiMethyN() function in lumi package. Default is "quantile"}
}
\details{
Creates:
\itemize{
\item  a PCA of normalized Beta values
\item  a Cluster of normalized Beta values
}}
\value{
 data.frame with the normalized M-values and pdf file called ExplorativeAnalysis.pdf with all the plots. 
 }
\author{ F.M. Mancuso <francesco.mancuso@crg.es>, N. Oien <Niles.Oien@colorado.edu> }
\examples{
\dontrun{
	Dir <- system.file("extdata/",package="HumMeth27QCReport")
	ImportDataR <- ImportData(Dir)
	normMvalues <- NormCheck(ImportDataR,platform="Hum27",pval=0.05,ChrX=F,ClustMethod="euclidean",normMethod="quantile")
}
}
\keyword{ hplot }
\keyword{ manip }
