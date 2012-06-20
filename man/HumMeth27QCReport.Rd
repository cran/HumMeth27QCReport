\name{HumMeth27QCReport}
\alias{HumMeth27QCReport}
\docType{package}
\title{ Generate a QC Report for Illumina Infinium BeadChip methylation assay}
\description{
HumMeth27QCReport is a tool for quality control and preprocessing of Illuminas Infinium BeadChip methylation assay. 
The automated analysis pipeline comprises data import, normalization, quality diagnostics and data export. 
}
\usage{
HumMeth27QCReport(ImportDataR, platform, pval, ChrX, ClustMethod, quoteOutput, normMethod)
}
\arguments{
  \item{ImportDataR}{the object that comes from ImportData function}
  \item{platform}{The type of Illumina Infinium BeadChip methylation assay. This must be one of "Hum27" (Infinium HumanMethylation27 BeadChip) or "Hum450" (Infinium HumanMethylation450 BeadChip).}
  \item{pval}{A p-value threshold number to define which sample keep for the normalization and the following analyses}
  \item{ClustMethod}{The distance measure to be used for the clustering. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary", "pearson", "correlation", "spearman" or "kendall"}
  \item{ChrX}{Logical value indicating whether the CpGs that belongs to chromosome X should be deleted from normalization and the following analyses. The default is FALSE}
  \item{quoteOutput}{Logical value determining if non-numeric items in the output file of normalized data should be in quotes. Default is TRUE}
  \item{normMethod}{Commonly either "quantile" or "ssn". Determines the normalization method. See documentation for lumiMethyN() function in lumi package. Default is "quantile"}
}
\details{
 A wide variety of tools from different packages are used, to produce pdfs with various statistics that should help determine whether there are data quality problems with the data.
 This package only works for Illumina Infinium methylation data and it requires a working Perl installation
 Creates different plots to asses the quality of the samples:
 
\itemize{
\item a histogram for each internal control.
\item an Intensity Graph plot foreach sample recalling the "plotSampleIntensities" function of methylumi package.
\item a histogram with the percentage of non dectected CPG (that is the CPGs tha have a detection p-value bigger than 0.05 or 0.01.
\item a histogram with the average p-value for each sample.
\item a PCA of normalized Beta values.
\item a Cluster of normalized Beta values.  
  }
 As ulterior output, an Excel file is provided. It contains the normalized BetaValues, a summary of the Internal Controls and of the gene detection and different lists of non-detected CPGs.  
 
 More details are available in the vignette.  
}
\value{
 data.frame with the normalized M-values and three pdf files (InternalControl.pdf, QualityCheck.pdf, ExplorativeAnalysis.pdf) with all the plots described previously.
}
\author{ F.M. Mancuso <francesco.mancuso@crg.es> }
\examples{
\dontrun{
	Dir <- system.file("extdata/",package="HumMeth27QCReport")
    ImportDataR <- ImportData(Dir)
	normMvalues <- HumMeth27QCReport(ImportDataR, platform="Hum27", pval=0.05, ChrX=FALSE, ClustMethod="euclidean", quoteOutput=FALSE, normMethod="quantile")
}
}
\keyword{ hplot }
\keyword{ manip }
