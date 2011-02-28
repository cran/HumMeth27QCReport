\name{getAssayControls}
\alias{getAssayControls}
\title{Illuminas Infinium HumanMethylation Internal Controls}
\description{
Creates all the histogram plots relative to the internal controls of the Illumina assay into pdf files.
}
\usage{
getAssayControls(ImportDataR, platform)
}
\arguments{
	\item{ImportDataR}{the object that comes from ImportData function}
    \item{platform}{The type of Illumina Infinium BeadChip methylation assay. 
		This must be one of "Hum27" (Infinium HumanMethylation27 BeadChip) or "Hum450" (Infinium HumanMethylation450 BeadChip)}
}


\details{
 This function creates all the histograms foreach internal control. The sample-independent controls let you evaluate the quality of specific
steps in the process flow; instead, the sample-dependent controls let you evaluate performance across samples. 

\itemize{
\item Figure 1: Barplot of DNP staining control 
This figure represents the ratio (%) between background and signal for Staining control in the red channel (DNP). Staining controls are used to examine the efficiency of the staining step in both the red and green channels. Staining controls have dinitrophenyl (DNP) or biotin attached to the beads. 
The ratios should result in low signal, indicating that the staining step was efficient.
\item Figure 2: Barplot of Biotin staining control
This figure represents the ratio (%) between background and signal for Staining control in the green channel (Biotin). These controls are independent of the hybridization and extension step. 
The ratios should result in low signal, indicating that the staining step was efficient.
\item Figure 3: Barplot of Hybridization control
This figure represents the ratio (%) between background and signal for Hybridization controls in the green channel for three concentrations. The hybridization controls test the overall performance of the entire assay using synthetic targets instead of amplified DNA. These synthetic targets complement the sequence on the array perfectly, allowing the probe to extend on the synthetic target as template. The synthetic targets are present in the hybridization buffer at three levels, monitoring the response from high-concentration (5 pM), medium-concentration (1 pM), and low-concentration (0.2 pM) targets. 
All bead type IDs should result in signal with various intensities, corresponding to the concentrations of the initial synthetic targets.
\item Figure 4: Barplot of Target Removal control
This figure represents the intensity value for Target removal controls in the green channel. Target removal controls test the efficiency of the stripping step after the extension reaction. The control oligos are extended using the probe sequence as template. This process generates labeled targets. The probe sequences are designed such that extension from the probe does not occur. 
All target removal controls should result in low signal, indicating that the targets were removed efficiently after extension.  Values < 3400 have been detected (108 samples). There is not a range specified from illumina, the value is based on previous experiments run in our facility. 
\item Figure 5: Barplot of Extension control: green channel
This figure represents the ratio (%) between background and signal for Extension control in the green channel (C,G). Extension controls test the extension efficiency of A, T, C, and G nucleotides from a hairpin probe, and are therefore sample-independent. 
The ratios should result in low signal, indicating that the extension was efficient.
\item Figure 6: Barplot of Extension control: red channel 
This figure represents the ratio (%) between background and signal for Extension control in the red channel (A,T). 
The ratios should result in low signal, indicating that the extension was efficient.
\item Figure 7: Barplot of Bisulfite control
This figure represents the ratio (%) between background and signal for Bisultife conversion control. The Bisulfite conversion Control asses the efficiency of bisulfite conversion of the genomic DNA. The Infinium Methylation probes query a [C/T] polymorphism created by bisulfite conversion of two different Hind III sites [AAGCTT] in the genome. If the bisulfite conversion reaction was successful, the "C" (Converted) probes will match the converted sequence and get extended. If the sample has unconverted DNA, the "U" (Unconverted) probes will get extended. There are no underlying C bases in the primer landing sites, except for the query site itself. Performance of bisulfite conversion controls should only be monitored in the Green channel.
The ratios should result in low signal, indicating that the Bisulfite conversion was efficient. 
\item Figure 8: Barplot of specificity control (mismatch 1) in red channel
This figure represents the ratio (%) between background (MM) and signal (PM) for Specificity controls in red channel. In the Infinium Methylation assay, the methylation status of a particular cytosine is carried out following bisulfite treatment of DNA by using query probes for unmethylated and methylated state of each CpG locus. In assay oligo design, the A/T match corresponds to the unmethylated status of the interrogated C, and G/C match corresponds to the methylated status of C. G/T mismatch controls check for non-specific detection of methylation signal over unmethylated background. Specificity controls are designed against non-polymorphic T sites. PM controls correspond to A/T perfect match and should give high signal. MM controls correspond to G/T mismatch and should give low signal.
The ratios should result in low signal, indicating that the performance of the assay was efficient.
\item Figure 9: Barplot of specificity control (mismatch 2) in green channel
This figure represents the ratio (%) between background (MM) and signal (PM) for Specificity controls in the green channel. PM controls correspond to A/T perfect match and should give high signal. MM controls correspond to G/T mismatch and should give low signal.
The ratios should result in low signal, indicating that the performance of the assay was efficient.
\item Figure 10: Barplot of Negative control
This figure represents the intensity value for the Negative control. Negative control probes are randomly permutated sequences that should not hybridize to the DNA template. Negative controls are particularly important for methylation studies because of a decrease in sequence complexity after bisulfite conversion. The mean signal of these probes defines the system background. This is a comprehensive measurement of background, including signal resulting from cross-hybridization, as well as non-specific extension and imaging system background. 
All target negative controls should result in low signal. Values < 2500 have been detected (108 samples). There is not a range specified from illumina, the value is based on previous experiments run in our facility. 
\item Figure 11: Barplot for green channel of Non-Polymorphic control
This figure represents the ratio (%) between background and signal for Non-Polymorphic control in the green channel. Non-polymorphic controls test the overall performance of the assay, from amplification to detection, by querying a particular base in a non-polymorphic region of the bisulfite genome. They let compare assay performance across different samples. One non-polymorphic control has been designed to query each of the four nucleotides (A, T, C and G). The target with the C base results from querying the opposite whole genome amplified strand generated from the converted strand.
The ratios should result in low signal, indicating that the performance of the assay was efficient.
\item Figure 12: Barplot for red channel of Non-Polymorphic control
This figure represents the ratio (%) between background and signal for Non-Polymorphic control in the red channel. 
The ratios should result in low signal, indicating that the performance of the assay was efficient.
}
}
\value{
data.frame with all the values of the internal controls and a pdf file called InternalControl.pdf with the barplots.  
}
\author{ F.M. Mancuso <francesco.mancuso@crg.es> }
\examples{
\dontrun{
	Dir <- system.file("extdata/",package="HumMeth27QCReport")
	platform <- "Hum27"
    ImportDataR <- ImportData(Dir)	
	ControlResults <- getAssayControls(ImportDataR, platform)
}
}
\keyword{ hplot }
\keyword{ manip }
