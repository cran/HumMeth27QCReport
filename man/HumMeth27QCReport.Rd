\name{HumMeth27QCReport}
\alias{HumMeth27QCReport}

\docType{package}
\title{ Generate a QC Report for Illuminas Infinium HumanMethylation27 BeadChip methylation assay}
\description{
HumMeth27QCReport is a tool for quality control and preprocessing of Illuminas Infinium HumanMethylation27 BeadChip methylation assay. 
The automated analysis pipeline comprises data import, normalization, quality diagnostics and data export. 
}

\usage{
HumMeth27QCReport(Dir)
}

\arguments{
  \item{Dir}{A character string containing the location of the directory 
    in which the input file are}
}

\details{
 A wide variety of tools from different packages are used, to produce pdfs with various statistics that should help determine whether there are data quality problems with the data.
 This package only works for Illumina Infinium HumanMethylation27 data and it requires a working Perl installation
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

\author{ F.M. Mancuso <francesco.mancuso@crg.es> }

\keyword{ hplot }
\keyword{ manip }
