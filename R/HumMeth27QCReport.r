 #### HumMeth27QCReport v. 1.2.15



###############################################
##  Get the separator character for a file ####
getFileSepChar <- function(File) {

# Init the return value
  sepChar <- NULL

# Open the file, return if this fails.
  inFile <- file(File, "r")
  if (!(isOpen(inFile))){
      return(sepChar)
  }

# List the possible separators we might consider. The order is important, if two
# separator characters occur with equal frequency then the one listed first will be
# used.
  possibleSeps=c(",", "\t", " ", ";", ":")

# Initialize a vector of counts for each possible separator
  possibleSepCounts=c( 0, 0, 0, 0, 0)

# Look at the first 10 lines after the header line (arbitrary decision to do that)
  for (lineNum in 1:11){
    if(length(input <- readLines(inFile, n=1)) > 0){
      if (lineNum > 1){ # Skip the first line as it is a header
	ml=strsplit(input, character(0))
	lc <- ml[[1]]

# For each char in the line, go through our list of possible separator characters, increment
# the count for that possible separator if the char in the line matches the possible separator.
	maxLen=length(lc)
	if (maxLen > 8192) maxLen=8192
	for (i in 1:maxLen){
	  for (j in 1:length(possibleSeps)){
	    if (lc[i] == possibleSeps[j]){
	      possibleSepCounts[j] <- possibleSepCounts[j] + 1
	      break
	    }
	  }
	}
      }
    }
  }

# Close the input file.
  close(inFile)

# Find the separator char that has the maximum count.
  maxIndex <- -1 
  maxCount <- 0
  for (j in 1:length(possibleSeps)){
    if (possibleSepCounts[j] > maxCount){
      maxCount <- possibleSepCounts[j]
      maxIndex <- j
    }
  }

  if (maxIndex > 0){
    sepChar <- possibleSeps[maxIndex]
  }

  return(sepChar)

}

###############################################
###############   Import data   ###############
ImportData <- function(Dir) {

	setwd(Dir);

	WorkFiles <- list.files(path = ".",pattern = ".*.txt")
	AverageBeta <- WorkFiles[grep("AvgBeta", WorkFiles)]
	Ctrls <- WorkFiles[grep("Control", WorkFiles)]
	SamplesName <- WorkFiles[grep("Sample", WorkFiles)]
	try(Discarder <- WorkFiles[grep("Discard", WorkFiles)], silent=T)

	control <- read.table(Ctrls, header = TRUE, sep = "\t", comment.char = "") 
	AvBeta <- read.table(AverageBeta, header = TRUE, sep = "\t", flush=T, comment.char = "") 
	
	samps <- read.table(SamplesName, header = TRUE, sep = "\t", comment.char = "")
	DiscarderII <- NULL
	try(DiscarderII <- read.table(Discarder, header = F, sep = "\t"), silent=T)
	try(DiscarderII <- as.character(DiscarderII$V1), silent=T)

  Label <- samps$Sample.ID
  if (is.null(Label)==T) { Label <- samps$SampleID }
	samps$SampleLabel <-  Label
  
	colnames(samps)[1] <- "Index"
	colnames(samps)[2] <- "SampleID"

	nsample <- length(samps$Index)
	samps2 <- samps[with(samps, order(Index)), ]

	
  pdf(file="Sample.pdf",width=15,height=9, fonts="Times")  
  if((nsample %% 2 == 0) ==F) { # Odd number of samples
    SamplePDF <- data.frame(samps2$Index[1:(length(samps2$Index)/2+1)], as.character(samps2$SampleID[1:(length(samps2$SampleID)/2+1)]), rep(" ",(length(samps2$SampleID)/2)+1), c(samps2$Index[ceil((length(samps2$Index)/2)+1):(length(samps2$Index))],""), c(as.character(samps2$SampleID[((ceil(length(samps2$SampleID)/2+1))):length(samps2$SampleID)]),""))
    colnames(SamplePDF)<-c("Index","SampleID","","Index","SampleID")
  } else { # Even number of samples
    SamplePDF <- data.frame(samps2$Index[1:(length(samps2$Index)/2)], samps2$SampleID[1:(length(samps2$SampleID)/2)], rep(" ",(length(samps2$SampleID)/2)), samps2$Index[((length(samps2$Index)/2)+1):length(samps2$Index)], samps2$SampleID[((length(samps2$SampleID)/2)+1):length(samps2$SampleID)])
    colnames(SamplePDF)<-c("Index","SampleID","","Index","SampleID")
  } 
 
	textplot(SamplePDF, halign="center", valign="center", show.rownames = F)
	title("Sample List")
	dev.off()

	return(list(ctrl=control,samples=samps, AverageBeta=AverageBeta, Ctrls=Ctrls, nsample=nsample, DiscarderII=DiscarderII))
}


###################################################
###############  Internal Controls  ###############	


getAssayControls <- function(ImportDataR,platform) {
	if (platform == "Hum450")   {
		Ctrl <- read.table(system.file("extdata/InternalCtrls450.txt",package="HumMeth27QCReport"),header=TRUE,as.is=TRUE,sep="\t")
	} 
	if (platform == "Hum27")   { 
		Ctrl <- read.table(system.file("extdata/InternalCtrls27.txt",package="HumMeth27QCReport"),header=TRUE,as.is=TRUE,sep="\t")
	}
	
	control <- ImportDataR$ctrl
	samps <- ImportDataR$samples
	
	pdf(file="InternalControl.pdf",width=15,height=9)
	par(oma=c(1,1,1,1), mar=c(10.1,4.1,3.1,1.1))

	
##### QUALITY CHECK - Staining DNP
	
	DNP.med.Ctrl <- Ctrl[Ctrl$Purpose == "Staining DNP Sig", ]
	DNP.bgnd.Ctrl <- Ctrl[Ctrl$Purpose == "Staining DNP Bgnd", ]

	
	DNP.med <- t(control[control$ProbeID %in% DNP.med.Ctrl,seq(5,length(control),3)])
	DNP.bgnd <- t(control[control$ProbeID %in% DNP.bgnd.Ctrl,seq(5,length(control),3)])
	DNP.Ratio <- as.vector(DNP.bgnd / DNP.med)*100
	DNP <- data.frame("SAMPLE"=samps$SampleID,"DNPRatio"=DNP.Ratio)
	DNP <- DNP[with(DNP, order(SAMPLE)), ]
	
	if (max(DNP[,2])>50)   { 
		try(gap.barplot(DNP[,2], gap=c(50.0001,50.001), col=rep("red",length(DNP[,2])), xaxlab=DNP$SAMPLE,main="Staining DNP Control: Background (BGND) on Signal (MED)",las=2,ylim=c(4,61), ylab="%", ytics=c(seq(0,50,10)),xlab=""), silent=T)
		abline(h=seq(0,50,10),col="grey")	
	}   else   {
		barplot2(DNP[,2], col="red", names=DNP$SAMPLE, las=2, ylim=c(0,50), main="Staining DNP Control: Background (BGND) on Signal (MED)", ylab="%")
		abline(h=seq(0,50,5),col="grey")
	}
	
	
##### QUALITY CHECK - Staining Biotin
	Biotin.med.Ctrl <- Ctrl[Ctrl$Purpose == "Staining Biotin Sig", ]
	Biotin.bgnd.Ctrl <- Ctrl[Ctrl$Purpose == "Staining Biotin Bgnd", ]
	
	Biotin.med <- t(control[control$ProbeID %in% Biotin.med.Ctrl,seq(4,length(control),3)])
	Biotin.bgnd <- t(control[control$ProbeID %in% Biotin.bgnd.Ctrl,seq(4,length(control),3)])
	Biotin.Ratio <- as.vector(Biotin.bgnd / Biotin.med)*100
	Biotin <- data.frame("SAMPLE"=samps$SampleID,"BiotinRatio"=Biotin.Ratio)
	Biotin <- Biotin[with(Biotin, order(SAMPLE)), ]
	
	if (max(Biotin[,2])>50)   { 
		try(gap.barplot(Biotin[,2], gap=c(50.0001,50.001), col=rep("green",length(Biotin[,2])), xaxlab=Biotin$SAMPLE,main="Staining Biotin Control: Background (BGND) on Signal (MED)",las=2,ylim=c(4,61), ylab="%", ytics=c(seq(0,50,10)),xlab=""), silent=T)
		abline(h=seq(0,50,10),col="grey")	
	}   else   {
		barplot2(Biotin[,2], col="green", names=Biotin$SAMPLE, las=2, ylim=c(0,50), main="Staining Biotin Control: Background (BGND) on Signal (MED)", ylab="%")
		abline(h=seq(0,50,5),col="grey")
	}
	
	
##### QUALITY CHECK - Hybridization
	HybL.Ctrl <- Ctrl[Ctrl$Purpose == "Hybridization L", ]
	HybM.Ctrl <- Ctrl[Ctrl$Purpose == "Hybridization M", ]
	HybH.Ctrl <- Ctrl[Ctrl$Purpose == "Hybridization H", ]
	
	HybL.g <- t(control[control$ProbeID %in% HybL.Ctrl,seq(4,length(control),3)])
	HybL.r <- t(control[control$ProbeID %in% HybL.Ctrl,seq(5,length(control),3)])
	HybM.g <- t(control[control$ProbeID %in% HybM.Ctrl,seq(4,length(control),3)])
	HybM.r <- t(control[control$ProbeID %in% HybM.Ctrl,seq(5,length(control),3)])
	HybH.g <- t(control[control$ProbeID %in% HybH.Ctrl,seq(4,length(control),3)])
	HybH.r <- t(control[control$ProbeID %in% HybH.Ctrl,seq(5,length(control),3)])
	RatioHybL <- round((HybL.r / HybL.g)*100,1)
	RatioHybM <- round((HybM.r / HybM.g)*100,1)
	RatioHybH <- round((HybH.r / HybH.g)*100,1)
	
	Hyb <- data.frame("SAMPLE"=samps$SampleID,"L"=RatioHybL,"M"=RatioHybM,"H"=RatioHybH)
	Hyb <- Hyb[with(Hyb, order(SAMPLE)), ]
	HybMax <- apply(as.matrix(t(data.frame("L"=Hyb[,2],"M"=Hyb[,3],"H"=Hyb[,4]))),1,max)
	
 
	if (max(HybMax)>50)   { 
		barplot2(as.matrix(t(data.frame("L"=Hyb[,2],"M"=Hyb[,3],"H"=Hyb[,4]))), names=Hyb$SAMPLE, ylim=c(0,max(HybMax)+50),ylab="%", col=c(1,2,3), las=2, beside=T, legend = c("HybL","HybM","HybH"), main="Hibridization Control: Background (Red) on Signal (Green)", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
	}   else   {
		barplot2(as.matrix(t(data.frame("L"=Hyb[,2],"M"=Hyb[,3],"H"=Hyb[,4]))), names=Hyb$SAMPLE, ylim=c(0,50),ylab="%",
				 col=c(1,2,3), las=2, beside=T, legend = c("HybL","HybM","HybH"), main="Hibridization Control: Background (Red) on Signal (Green)", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
	}
		
	
##### QUALITY CHECK - Target Removal
	TMCtrl <- Ctrl[Ctrl$Purpose == "Target Removal", ]
	TM <- t(control[control$ProbeID %in% TMCtrl$ID,seq(4,length(control),3)])
	TMmean <- apply(as.matrix(TM), 1,mean)
	
	TMmeanC <- data.frame("SAMPLE"=samps$SampleID, "tm"=TMmean)
	TMmeanC <- TMmeanC[with(TMmeanC, order(SAMPLE)), ]
	
	barplot2(TMmeanC[,2], names=TMmeanC$SAMPLE, col=c("green"), las=2, ylim=c(0,max(TMmeanC[,2])+1000), main="Target Removal Control on Green Channel", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
	
	
##### QUALITY CHECK - Extension
	
	Extension.A.Ctrl <- Ctrl[Ctrl$Purpose == "Extension A", ]
	Extension.T.Ctrl <- Ctrl[Ctrl$Purpose == "Extension T", ]
	Extension.C.Ctrl <- Ctrl[Ctrl$Purpose == "Extension C", ]
	Extension.G.Ctrl <- Ctrl[Ctrl$Purpose == "Extension G", ]
	
	
	Extension.Ag <- t(control[control$ProbeID %in% Extension.A.Ctrl,seq(4,length(control),3)])
	Extension.Ar <- t(control[control$ProbeID %in% Extension.A.Ctrl,seq(5,length(control),3)])
	Extension.Tg <- t(control[control$ProbeID %in% Extension.T.Ctrl,seq(4,length(control),3)]) 
	Extension.Tr <- t(control[control$ProbeID %in% Extension.T.Ctrl,seq(5,length(control),3)]) 
	Extension.Gg <- t(control[control$ProbeID %in% Extension.C.Ctrl,seq(4,length(control),3)])
	Extension.Gr <- t(control[control$ProbeID %in% Extension.C.Ctrl,seq(5,length(control),3)])
	Extension.Cg <- t(control[control$ProbeID %in% Extension.G.Ctrl,seq(4,length(control),3)])
	Extension.Cr <- t(control[control$ProbeID %in% Extension.G.Ctrl,seq(5,length(control),3)])

	Extension.AT.g <- apply(as.matrix(data.frame(Extension.Ag,Extension.Tg)),1,sum)
	Extension.GC.g <- apply(as.matrix(data.frame(Extension.Gg,Extension.Cg)),1,sum)
	Extension.AT.r <- apply(as.matrix(data.frame(Extension.Ar,Extension.Tr)),1,sum)
	Extension.GC.r <- apply(as.matrix(data.frame(Extension.Gr,Extension.Cr)),1,sum)
	
	ExtRatio.g <- data.frame("SAMPLE"=samps$SampleID,"Ext.g"=(Extension.AT.g / Extension.GC.g)*100) 
	ExtRatio.g <- ExtRatio.g[with(ExtRatio.g, order(SAMPLE)), ]
	
	ExtRatio.r <- data.frame("SAMPLE"=samps$SampleID,"Ext.r"=(Extension.GC.r / Extension.AT.r)*100) 
	ExtRatio.r <- ExtRatio.r[with(ExtRatio.r, order(SAMPLE)), ]
	
	if (max(ExtRatio.g[,2])>50)   { 
		try(gap.barplot(ExtRatio.g[,2], gap=c(50.0001,50.001), col=rep("green",length(ExtRatio.g[,2])), xaxlab=ExtRatio.g$SAMPLE,main="Extension Control (green channel): Background (AT) on Signal (GC)",las=2, ylab="%", ylim=c(4,61), ytics=c(seq(0,50,10)),xlab=""), silent=T)
		abline(h=seq(0,50,10),col="grey") 
	}   else   {
		barplot2(ExtRatio.g[,2], col="green", names=ExtRatio.g$SAMPLE, las=2, ylim=c(0,50), main="Extension Control (green channel): Background (AT) on Signal (GC)", ylab="%")
		abline(h=seq(0,50,5),col="grey")
	}
	
	if (max(ExtRatio.r[,2])>50)   { 
		try(gap.barplot(ExtRatio.r[,2], gap=c(50.0001,50.001), col=rep("red",length(ExtRatio.r[,2])), xaxlab=ExtRatio.r$SAMPLE,main="Extension control (red channel): Background (GC) on Signal (AT)",las=2,ylim=c(4,61), ylab="%", ytics=c(seq(0,50,10)),xlab=""), silent=T)
		abline(h=seq(0,50,10),col="grey")
	}   else   {
		barplot2(ExtRatio.r[,2], col="red", names=ExtRatio.r$SAMPLE, las=2, ylim=c(0,50), main="Extension Control (red channel): Background (AT) on Signal (GC)", ylab="%")
		abline(h=seq(0,50,5),col="grey")
	}
	
	
##### QUALITY CHECK - Bisulfite Conversion
	BisulfiteCtrl_SigG <- Ctrl[Ctrl$Purpose == "Bisulfite conversion I C green", ]
	BisulfiteCtrl_BckgG <- Ctrl[Ctrl$Purpose == "Bisulfite conversion I U green", ]
	BisulfiteCtrl_SigR <- Ctrl[Ctrl$Purpose == "Bisulfite conversion I C red", ]
	BisulfiteCtrl_BckgR <- Ctrl[Ctrl$Purpose == "Bisulfite conversion I U red", ]
	BisulfiteCtrl_II <- Ctrl[Ctrl$Purpose == "Bisulfite conversion II", ]
	
# Bisulfite I on green cheannel
	Bisulfite_SigG <- t(control[control$ProbeID %in% BisulfiteCtrl_SigG$ID,seq(4,length(control),3)])
	Bisulfite_Sig_SumG <- apply(as.matrix(Bisulfite_SigG), 1,sum)
	Bisulfite_BckgG <- t(control[control$ProbeID %in% BisulfiteCtrl_BckgG$ID,seq(4,length(control),3)])
	Bisulfite_Bckg_SumG <- apply(as.matrix(Bisulfite_BckgG), 1,sum)

	BisulfiteRatioG <- (Bisulfite_Bckg_SumG / Bisulfite_Sig_SumG)*100
	BisulfiteG <- data.frame("SAMPLE"=samps$SampleID,"BisulfiteRatio"=BisulfiteRatioG)
	BisulfiteG <- BisulfiteG[with(BisulfiteG, order(SAMPLE)), ]
	
	if (max(BisulfiteG[,2])>50)   { 
		gap.barplot(BisulfiteG[,2], gap=c(50.0001,50.001), col=rep("green",length(BisulfiteG[,2])), xaxlab=BisulfiteG$SAMPLE,
					main="Bisulfite Control (green channel): Background (U) on Signal (C)",	las=2, ylab="%", ytics=c(seq(0,50,10)),xlab="", horiz=F, ylim=c(4,61)) 
		abline(h=seq(0,50,5),col="grey") 
	}   else   {
		barplot2(BisulfiteG[,2], col="green", names=BisulfiteG$SAMPLE, las=2, ylim=c(0,50), main="Bisulfite Control (green channel): Background (U) on Signal (C)", ylab="%")
		abline(h=seq(0,50,5),col="grey")
	}
	
	if (platform == "Hum450")   {
		
# Bisulfite I on red cheannel
	Bisulfite_SigR <- t(control[control$ProbeID %in% BisulfiteCtrl_SigR$ID,seq(5,length(control),3)])
	Bisulfite_Sig_SumR <- apply(as.matrix(Bisulfite_SigR), 1,sum)
	Bisulfite_BckgR <- t(control[control$ProbeID %in% BisulfiteCtrl_BckgR$ID,seq(5,length(control),3)])
	Bisulfite_Bckg_SumR <- apply(as.matrix(Bisulfite_BckgR), 1,sum)
	
	BisulfiteRatioR <- (Bisulfite_Bckg_SumR / Bisulfite_Sig_SumR)*100
	BisulfiteR <- data.frame("SAMPLE"=samps$SampleID,"BisulfiteRatio"=BisulfiteRatioR)
	BisulfiteR <- BisulfiteR[with(BisulfiteR, order(SAMPLE)), ]
	
	if (max(BisulfiteR[,2])>50)   { 
		gap.barplot(BisulfiteR[,2], gap=c(50.0001,50.001), col=rep("red",length(BisulfiteR[,2])), xaxlab=BisulfiteR$SAMPLE,	main="Bisulfite Control (red channel): Background (U) on Signal (C)",	las=2, ylab="%", ytics=c(seq(0,50,10)),xlab="", horiz=F, ylim=c(4,61)) 
		abline(h=seq(0,50,5),col="grey") 
	}   else   {
		barplot2(BisulfiteR[,2], col="red", names=BisulfiteG$SAMPLE, las=2, ylim=c(0,50), main="Bisulfite Control (red channel): Background (U) on Signal (C)", ylab="%")
		abline(h=seq(0,50,5),col="grey")
	}	
	
# Bisulfite II
	BisulfiteIIR <- t(control[control$ProbeID %in% BisulfiteCtrl_II$ID,seq(5,length(control),3)])
	BisIIR <- apply(as.matrix(BisulfiteIIR), 1,sum)
	BisulfiteIIG <- t(control[control$ProbeID %in% BisulfiteCtrl_II$ID,seq(4,length(control),3)])
	BisIIG <- apply(as.matrix(BisulfiteIIG), 1,sum)	
	
	BisulfiteIIRatio <- (BisIIG / BisIIR)*100
		
	BisIIC <- data.frame("SAMPLE"=samps$SampleID,"BisII"=BisulfiteIIRatio)
	BisIIC <- BisIIC[with(BisIIC, order(SAMPLE)), ]
	
  if (max(BisIIC[,2])>50)   { 
		gap.barplot(BisIIC[,2], gap=c(50.0001,50.001), col=rep("red",length(BisIIC[,2])), xaxlab=BisIIC$SAMPLE,
					main="Bisulfite II Control: Background on Signal",	las=2, ylab="%", ytics=c(seq(0,50,10)),xlab="", horiz=F, ylim=c(4,61)) 
		abline(h=seq(0,50,5),col="grey") 
	}   else   {
		barplot2(BisIIC[,2], col="red", names=BisIIC$SAMPLE, las=2, ylim=c(0,50), main="Bisulfite II Control: Background on Signal", ylab="%")
		abline(h=seq(0,50,5),col="grey")
	}
	
	}
	
	
##### QUALITY CHECK - Specificity
	
	SpecificityCtrl_MM1 <- Ctrl[Ctrl$Purpose == "Specificity I MM green", ]
	SpecificityCtrl_PM1 <- Ctrl[Ctrl$Purpose == "Specificity I PM green", ]
	SpecificityCtrl_MM2 <- Ctrl[Ctrl$Purpose == "Specificity I MM red", ]
	SpecificityCtrl_PM2 <- Ctrl[Ctrl$Purpose == "Specificity I PM red", ]
	SpecificityCtrl_II <- Ctrl[Ctrl$Purpose == "Specificity II", ]
	
# Specificity I - mismatch 1	
	Specificity.mis1.PM <- t(control[control$ProbeID %in% SpecificityCtrl_PM2$ID,seq(5,length(control),3)])
	SpecPM.mis1_Sum <- apply(as.matrix(Specificity.mis1.PM), 1,sum)
	Specificity.mis1.MM <- t(control[control$ProbeID %in% SpecificityCtrl_MM2$ID,seq(5,length(control),3)])
	SpecMM.mis1_Sum <- apply(as.matrix(Specificity.mis1.MM), 1,sum)
	
	Specificity.mis1.Ratio <- (SpecMM.mis1_Sum / SpecPM.mis1_Sum)*100
	Spec_mm1 <- data.frame("SAMPLE"=samps$SampleID,"Specificity.mis1.Ratio"=Specificity.mis1.Ratio)
	Spec_mm1 <- Spec_mm1[with(Spec_mm1, order(SAMPLE)), ]
	
	if (max(Spec_mm1[,2])>50)   { 
		try(gap.barplot(Spec_mm1[,2], gap=c(50.0001,50.001), col=rep("red",length(Spec_mm1[,2])), xaxlab=Spec_mm1$SAMPLE,main="Specificity Control mismatch 1: Background (MM) on Signal (PM)",
						las=2,ylim=c(4,61), ylab="%", ytics=c(seq(0,50,10)),xlab=""), silent=T)
		abline(h=seq(0,50,10),col="grey")
	}   else   {
		barplot2(Spec_mm1[,2], col="red", names=Spec_mm1$SAMPLE, las=2, ylim=c(0,50), main="Specificity Control mismatch 1: Background (MM) on Signal (PM)", ylab="%")
		abline(h=seq(0,50,5),col="grey")
	}
	
# Specificity I - mismatch 2
	Specificity.mis2.PM <- t(control[control$ProbeID %in% SpecificityCtrl_PM1$ID,seq(4,length(control),3)])
	SpecPM.mis2_Sum <- apply(as.matrix(Specificity.mis2.PM), 1,sum)
	Specificity.mis2.MM <- t(control[control$ProbeID %in% SpecificityCtrl_MM1$ID,seq(4,length(control),3)])
	SpecMM.mis2_Sum <- apply(as.matrix(Specificity.mis2.MM), 1,sum)
	
	Specificity.mis2.Ratio <- (SpecMM.mis2_Sum / SpecPM.mis2_Sum)*100
	Spec_mm2 <- data.frame("SAMPLE"=samps$SampleID,"Specificity.mis2.Ratio"=Specificity.mis2.Ratio)
	Spec_mm2 <- Spec_mm2[with(Spec_mm2, order(SAMPLE)), ]
	
	if (max(Spec_mm2[,2])>50)   { 
		try(gap.barplot(Spec_mm2[,2], gap=c(50.0001,50.001), col=rep("green",length(Spec_mm2[,2])), xaxlab=Spec_mm2$SAMPLE,main="Specificity Control mismatch 2: Background (MM) on Signal (PM)",
						las=2,ylim=c(4,61), ylab="%", ytics=c(seq(0,50,10)),xlab=""), silent=T)
		abline(h=seq(0,50,10),col="grey")
	}   else   {
		barplot2(Spec_mm2[,2], col="green", names=Spec_mm2$SAMPLE, las=2, ylim=c(0,50), main="Specificity Control mismatch 2: Background (MM) on Signal (PM)", ylab="%")
		abline(h=seq(0,50,5),col="grey")
	}
	
	if (platform == "Hum450")   {
		
# Specificity II	
	SpecificityIIR <- t(control[control$ProbeID %in% SpecificityCtrl_II$ID,seq(5,length(control),3)])
	SpecIIR <- apply(as.matrix(SpecificityIIR), 1,sum)
  SpecificityIIG  <- t(control[control$ProbeID %in% SpecificityCtrl_II$ID,seq(4,length(control),3)])
	SpecIIG <- apply(as.matrix(SpecificityIIG), 1,sum)
  
  SpecIIRatio <- (SpecIIG / SpecIIR)*100
  
	SpecIIC <- data.frame("SAMPLE"=samps$SampleID,"SpecII"=SpecIIRatio)
	SpecIIC <- SpecIIC[with(SpecIIC, order(SAMPLE)), ]
	
  if (max(SpecIIC[,2])>50)   { 
		try(gap.barplot(SpecIIC[,2], gap=c(50.0001,50.001), col=rep("red",length(SpecIIC[,2])), xaxlab=SpecIIC$SAMPLE,main="Specificity Control II: Background on Signal",
						las=2,ylim=c(4,61), ylab="%", ytics=c(seq(0,50,10)),xlab=""), silent=T)
		abline(h=seq(0,50,10),col="grey")
	}   else   {
		barplot2(SpecIIC[,2], col="red", names=SpecIIC$SAMPLE, las=2, ylim=c(0,50), main="Specificity Control II: Background on Signal", ylab="%")
		abline(h=seq(0,50,5),col="grey")
	  }
	}

	
##### QUALITY CHECK - Negative
	NEGATIVECtrl <- Ctrl[Ctrl$Purpose == "Negative", ]
	Negative <- t(control[control$ProbeID %in% NEGATIVECtrl$ID,seq(4,length(control),3)])
	Neg <- apply(as.matrix(Negative), 1,mean)
	
	NegC <- data.frame("SAMPLE"=samps$SampleID,"Neg"=Neg)
	NegC <- NegC[with(NegC, order(SAMPLE)), ]
	
	barplot2(NegC[,2], col="blue", names=NegC$SAMPLE, las=2, ylim=c(0,max(NegC[,2])+1000), main="Negative Control", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
	
	
##### QUALITY CHECK - Non-Polymorphic (A+T / G+C)
	NPctrl.A <- Ctrl[Ctrl$Purpose == "Non-Polymorphic A", ]
	NPctrl.T <- Ctrl[Ctrl$Purpose == "Non-Polymorphic T", ]
	NPctrl.C <- Ctrl[Ctrl$Purpose == "Non-Polymorphic C", ]
	NPctrl.G <- Ctrl[Ctrl$Purpose == "Non-Polymorphic G", ]
	
	NP.Ag <- t(control[control$ProbeID %in% NPctrl.A$ID,seq(4,length(control),3)])  #### nonpolymorphic	control to query A in green channel
	NP.Ar <- t(control[control$ProbeID %in% NPctrl.A$ID,seq(5,length(control),3)])  #### nonpolymorphic	control to query A in red channel
	NP.Tg <- t(control[control$ProbeID %in% NPctrl.T$ID,seq(4,length(control),3)])  #### nonpolymorphic	control to query T in green channel
	NP.Tr <- t(control[control$ProbeID %in% NPctrl.T$ID,seq(5,length(control),3)])  #### nonpolymorphic	control to query T in red channel
	NP.Gg <- t(control[control$ProbeID %in% NPctrl.G$ID,seq(4,length(control),3)])  #### nonpolymorphic	control to query G in green channel
	NP.Gr <- t(control[control$ProbeID %in% NPctrl.G$ID,seq(5,length(control),3)])  #### nonpolymorphic	control to query G in red channel
	NP.Cg <- t(control[control$ProbeID %in% NPctrl.C$ID,seq(4,length(control),3)])  #### nonpolymorphic	control to query C in green channel
	NP.Cr <- t(control[control$ProbeID %in% NPctrl.C$ID,seq(5,length(control),3)])  #### nonpolymorphic	control to query C in red channel
	
	AvNP.AT.g <- apply(as.matrix(data.frame(NP.Ag,NP.Tg)),1,sum)
	AvNP.GC.g <- apply(as.matrix(data.frame(NP.Cg,NP.Gg)),1,sum)	
	AvNP.AT.r <- apply(as.matrix(data.frame(NP.Ar,NP.Tr)),1,sum)
	AvNP.GC.r <- apply(as.matrix(data.frame(NP.Cr,NP.Gr)),1,sum)	
	
	NPRatio.g <- data.frame("SAMPLE"=samps$SampleID,"NP.g"=(AvNP.AT.g / AvNP.GC.g)*100) 
	NPRatio.g <- NPRatio.g[with(NPRatio.g, order(SAMPLE)), ]
	NPRatio.r <- data.frame("SAMPLE"=samps$SampleID,"NP.r"=(AvNP.GC.r / AvNP.AT.r)*100) 
	NPRatio.r <- NPRatio.r[with(NPRatio.r, order(SAMPLE)), ]
	
	if (max(NPRatio.g[,2]) > 50)   { 
		try(gap.barplot(NPRatio.g[,2], gap=c(50.0001,50.001), col=rep("green",length(NPRatio.g[,2])), xaxlab=NPRatio.g$SAMPLE,
						main="Non-Polymorphic Control: Background (AT) on Signal (GC) - Green Channel",las=2,ylim=c(4,61), ylab="%", ytics=c(seq(0,50,10)),xlab=""), silent=T)
		abline(h=seq(0,50,5),col="grey")
	}   else   { 
		barplot2(NPRatio.g[,2], col="green", names=NPRatio.g$SAMPLE, las=2, ylim=c(0,50), main="Non-Polymorphic Control: Background (AT) on Signal (GC) - Green Channel", ylab="%")
		abline(h=seq(0,50,5),col="grey")
	}
	
	if (max(NPRatio.r[,2]) > 50)   { 
		try(gap.barplot(NPRatio.r[,2], gap=c(50.0001,50.001), col=rep("red",length(NPRatio.r[,2])), xaxlab=NPRatio.r$SAMPLE, main="Non-Polymorphic Control: Background (AT) on Signal (GC) - Red Channel",las=2,ylim=c(4,61), ylab="%", ytics=c(seq(0,50,10)),xlab=""), silent=T)
		abline(h=seq(0,50,5),col="grey")
	} else   {
		barplot2(NPRatio.r[,2], col="red", names=NPRatio.r$SAMPLE, las=2, ylim=c(0,50), main="Non-Polymorphic Control: Background (AT) on Signal (GC) - Red Channel", ylab="%")
		abline(h=seq(0,50,5),col="grey")
	}
	
	dev.off()
	
	samps <- samps[with(samps, order(SampleID)), ]
	if (platform == "Hum27")   {
		infotab_intC <- data.frame("Index"=samps$Index, "Sample ID"= samps$SampleID, "Sample Group" = samps$Sample.Group, "Sentrix Barcode"=samps$Sentrix.Barcode, "Sample Section"=samps$Sample.Section, "Sample_Plate"=samps$Sample_Plate, "Sample_Well"=samps$Sample_Well, "DNP staining" = round(DNP[,2],1), "Biotin staining" = round(Biotin[,2],1),"Low Hybridization Control" = round(Hyb[,2],1),"Medium Hybridization Control" = round(Hyb[,3],1), "High Hybridization Control" = round(Hyb[,4],1), "Target Removal Control" = round(TMmeanC[,2],1),"Extension.Control.Green" = round(ExtRatio.g[,2],1),"Extension.Control.Red" = round(ExtRatio.r[,2],1), "Bisulfite Conversion Control" = round(BisulfiteG[,2],1), "Specificity mismatch 1" = round(Spec_mm1[,2],1), "Specificity mismatch 2" = round(Spec_mm2[,2],1), "Negative Control" = round(NegC[,2],1), "Non.Polymorphic.Green" = round(NPRatio.g[,2],1),"Non.Polymorphic.Red" = round(NPRatio.r[,2],1))
	}
	
	if (platform == "Hum450")   {
		infotab_intC <- data.frame("Index"=samps$Index, "Sample ID"= samps$SampleID, "Sample Group" = samps$Sample.Group, "Sentrix Barcode"=samps$Sentrix.Barcode, "Sample Section"=samps$Sample.Section, "Sample_Plate"=samps$Sample_Plate, "Sample_Well"=samps$Sample_Well,"DNP staining" = round(DNP[,2],1), "Biotin staining" = round(Biotin[,2],1),"Low Hybridization Control" = round(Hyb[,2],1),"Medium Hybridization Control" = round(Hyb[,3],1), "High Hybridization Control" = round(Hyb[,4],1), "Target Removal Control" = round(TMmeanC[,2],1),"Extension.Control.Green" = round(ExtRatio.g[,2],1),"Extension.Control.Red" = round(ExtRatio.r[,2],1), "Bisulfite Conversion Control Green" = round(BisulfiteG[,2],1), "Bisulfite Conversion Control Red" = round(BisulfiteR[,2],1), "Bisulfite Conversion II" = round(BisIIC[,2],1), "Specificity mismatch 1" = round(Spec_mm1[,2],1), "Specificity mismatch 2" = round(Spec_mm2[,2],1), "Specificity II" = round(SpecIIC[,2],1), "Negative Control" = round(NegC[,2],1), "Non.Polymorphic.Green" = round(NPRatio.g[,2],1), "Non.Polymorphic.Red" = round(NPRatio.r[,2],1))
	}
	
	infotab_intC <- infotab_intC[with(infotab_intC, order(Index)), ]
	return(infotab_intC)
}


#################################################
###############   Quality Check   ###############

QCCheck <- function(ImportDataR, pval)   {

	if (missing(pval) || is.null(pval)) {
		pval <- 0.05
	}

	samps <- ImportDataR$samples
	AverageBeta <- ImportDataR$AverageBeta
	Ctrls <- ImportDataR$Ctrls
	nsample <- ImportDataR$nsample

	mldat <- methylumiR(AverageBeta, qcfile = Ctrls, sampleDescriptions = samps, sep=getFileSepChar(AverageBeta))
	pvalSample_CpG <- as.data.frame(pvals(mldat))
	
	pdf(file="QualityCheck.pdf",width=15,height=9)
	
##### Intensity Graphs
  if((nsample %% 2 == 0) ==F) {
    intNumber <- ceiling(nsample / 12)
    }	
  if((nsample %% 2 == 0) ==T) {
    intNumber <- round(nsample / 12,0)
    }
	for (k in 0:(intNumber-1))   {
		fileN <- paste("int",k+1,".pdf", sep = "", collapse = NULL)

		par(mfrow=c(4,3),oma=c(0.1,0.1,0.1,0.1), mar=c(8.1,2.1,2.1,2.1))
		for(i in ((1+12*k):(12*(k+1))))   { 
			try(plotSampleIntensities(mldat, s=i), silent=T)
    }
	}

	
##### No Detected
	NoDetect005 <- ""
	NoDetect001 <- ""
	for (i in 1:length(as.data.frame(pvals(mldat))))   {
		pvalSample <- as.data.frame(pvals(mldat)[,i])
		colnames(pvalSample) <- "pv"
		lenSample <- length(pvalSample$pv)
		NoDetect005[i] <-round(length(pvalSample[pvalSample$pv > 0.05,]) / lenSample * 100,2) 
		NoDetect001[i] <-round(length(pvalSample[pvalSample$pv > 0.01,]) / lenSample * 100,2) 	
	}

	samps2 <- samps[with(samps, order(Index)), ]
	
	NoDetectedAll <- data.frame("SampleID"=samps2$SampleID,"pv001"=as.numeric(NoDetect001),"pv005"=as.numeric(NoDetect005))
	NoDetectedAll <- NoDetectedAll[with(NoDetectedAll, order(SampleID)), ]

	par(mfrow=c(1,1),oma=c(1,1,1,1), mar=c(10.1,4.1,3.1,1.1))
	barplot2(rbind(as.numeric(NoDetectedAll[,2]),as.numeric(NoDetectedAll[,3])),beside=T, las = 2, 
			 col= c("blue","red"), names=NoDetectedAll$SampleID, axes=T, main="Percentage of non detected genes", ylab="%", 
			 ylim=c(0, 100), legend=c("% Not detected Genes (0.01)","% Not detected Genes (0.05)"), plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")


#### Average P-value
	avgPval <- as.data.frame(colMeans(pvals(mldat),na.rm = T))
	avgPval <- data.frame(rownames(avgPval),avgPval)
	colnames(avgPval) <- c("SampleID","avgPval")
	avgPval <- avgPval[with(avgPval, order(SampleID)), ]
	
	barplot2(avgPval[,2], ylab = "Average P-Value", las = 2, ylim=c(0,max(avgPval[,2])+0.1),col=4, main="Average p-value", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey",names=avgPval$SampleID)
	abline(h=pval, col="red",lty=2)
	dev.off()
	
	tempQC1 <- merge(samps,NoDetectedAll)
	tempQC2 <- merge(tempQC1,avgPval)
	
	
	Detection <- data.frame("Index"=tempQC2$Index, "Sample ID"= tempQC2$SampleID, "Sample Group" = tempQC2$Sample.Group, "Sentrix Barcode"=tempQC2$Sentrix.Barcode,
										"Sample Section"=tempQC2$Sample.Section, "Sample_Plate"=tempQC2$Sample_Plate, "Sample_Well"=tempQC2$Sample_Well, "Percentage of non detected genes (pval>0.01)"=tempQC2$pv001,
										"Percentage of non detected genes (pval>0.05)"=tempQC2$pv005, "Average Detection p-value"=tempQC2$avgPval)
	
	BadCpG_001 <- apply(pvalSample_CpG,1,function(y)sum(y > 0.01))		
	BadCpG_001_percent <- BadCpG_001/length(colnames(pvalSample_CpG))*100
	BadCpG_001_more5 <- as.data.frame(BadCpG_001_percent[BadCpG_001_percent>5])		
	BadCpG_001_more5 <- data.frame("CpG"=rownames(BadCpG_001_more5),"Percentage Bad CPG - 0.01"=BadCpG_001_more5) 
	
	BadCpG_005_count <- apply(pvalSample_CpG,1,function(y)sum(y > 0.05))
	BadCpG_005_percent <- BadCpG_005_count/length(colnames(pvalSample_CpG))*100
	BadCpG_005_more5 <- as.data.frame(BadCpG_005_percent[BadCpG_005_percent>5])		
	BadCpG_005_more5 <- data.frame("CpG"=rownames(BadCpG_005_more5),"Percentage Bad CPG - 0.05"=BadCpG_005_more5) 
	
	return(list("Detection"=Detection, "BadCpG_001_more5"=BadCpG_001_more5, "BadCpG_005_more5"=BadCpG_005_more5))
}

#############################################
################ PCA ANALYSIS ###############

pca.samples.plot<-function(data,int.col=(1:ncol(data)),main.str="Principal Component Analysis")
{
    data.nona<-na.delete(data[,int.col])
    pca.res<-prcomp(t(data.nona),tol=0.1,na.action=na.omit,center=T,scale=T)
    plot(pca.res$x,xlim=c(-200,200),xlab="PC 1",ylab="PC 2",main=main.str,pch=".")
    text(pca.res$x,colnames(data[,int.col]),cex=0.7)
}

################################################
############### Cluster Analysis ###############

cluster.samples.plot<-function(data,int.col=(1:ncol(data)), main.str="", method = "")
{
	# Matrix reprensentation of hierarchical clustering of the samples
    
    data.nona<- na.delete(data[,int.col])
    num.data <- t(data.nona)
    
    dist.mat <- as.matrix(Dist(num.data, method=method, diag =TRUE, upper =FALSE))
    resclust <- hcluster(dist.mat)
    dist.mat <- dist.mat[resclust$order,resclust$order]

    par(mfrow=c(2,1))
    image(x=seq(1,nrow(num.data)),y=seq(1,nrow(num.data)),z=dist.mat,axes=FALSE,xlab="",ylab="",main=main.str)
    axis(1,at=seq(1,nrow(num.data)),labels=rownames(dist.mat),las=2,cex.axis=0.5)
    axis(2,at=seq(1,nrow(num.data)),labels=rownames(dist.mat),las=2,cex.axis=0.25)
    plot(resclust,main="",xlab="",ylab="",cex=.6,axes=FALSE)
    par(mfrow=c(1,1))
}


#################################################
###############   Normalization   ###############

NormCheck <- function(ImportDataR, platform, pval, ChrX, ClustMethod, normMethod)   {
	
	if (missing(pval) || is.null(pval)) {
		pval <- 0.05
	}
	if (missing(ClustMethod) || is.null(ClustMethod)) {
		ClustMethod <- "euclidean"
	}
	if (missing(ChrX) || is.null(ChrX)) {
		ChrX <- "FALSE"
	}
        if (missing(normMethod) || is.null(normMethod)) {
                ChrX <- "quantile"
        }


	AverageBeta <- ImportDataR$AverageBeta
	Ctrl <- ImportDataR$Ctrl
	try(DiscarderII <- data.frame("ToDiscard"=ImportDataR$DiscarderII), silent=F)

	if (platform == "Hum450")   {
		lumiMethy <- lumiMethyR(AverageBeta, lib="IlluminaHumanMethylation450k.db", sep=getFileSepChar(AverageBeta))
		x <- IlluminaHumanMethylation450kCHR
		lumiMethy.c.adj <- lumiMethyC(lumiMethy)
	} 
	
	if (platform == "Hum27")   { 
		lumiMethy <- lumiMethyR(AverageBeta, lib="IlluminaHumanMethylation27k.db", sep=getFileSepChar(AverageBeta))
		x <- IlluminaHumanMethylation27kCHR
		lumiMethy.c.adj <- lumiMethyC(lumiMethy)  #### color balance adjustment
	}

	avgPval <- colMeans(detection(lumiMethy.c.adj),na.rm = T)
	toKeep <- (avgPval < pval)
	try(toKeep[as.character(DiscarderII$ToDiscard)] <- FALSE, silent=T)

# Do plot of intesity prior to normalization
	pdf(file="intensityPreNorm.pdf",width=15,height=9)
        par(oma=c(1,1,1,1), mar=c(7.1,4.1,3.1,1.1))
	plotDensity(estimateIntensity(lumiMethy.c.adj[,toKeep]), main="Pre normalization intensities")
	dev.off()

	lumiMethy.norm <- lumiMethyN(lumiMethy.c.adj[,toKeep], method=normMethod)   #### Normalization based on color balance adjusted data

# Plot of intesity post normalization
        pdf(file="intensityPostNorm.pdf",width=15,height=9)
        par(oma=c(1,1,1,1), mar=c(7.1,4.1,3.1,1.1))
        plotDensity(estimateIntensity(lumiMethy.norm), main="Post normalization intensities")
        dev.off()

	if (ChrX == "TRUE" || ChrX == "T")   {
		fData(lumiMethy.norm)$chrom <- as.factor(unlist(mget(featureNames(lumiMethy), x, ifnotfound=NA)))
		lumiMethy.norm.f <- data.frame(exprs(lumiMethy.norm), "Chr"= fData(lumiMethy.norm)$chrom)
		lumiMethy.norm.f <- lumiMethy.norm.f[lumiMethy.norm.f$Chr != "X",]
		lumiMethy.norm.f <- lumiMethy.norm.f[!is.na(lumiMethy.norm.f$Chr),]
		lumiMethy.norm.f$Chr <- NULL
		NormBetasVal <- lumiMethy.norm.f
	}   
	if (ChrX == "FALSE" || ChrX == "F")   {
		NormBetasVal <- as.data.frame(exprs(lumiMethy.norm[toKeep,]))
	}
	
	pdf(file="ExplorativeAnalysis.pdf",width=15,height=9)
	par(oma=c(1,1,1,1), mar=c(7.1,4.1,3.1,1.1))
	pca.samples.plot(NormBetasVal)
	cluster.samples.plot(NormBetasVal,main.str="Hierarchical Clustering",method=ClustMethod)
	dev.off()
	
	lenNormBetasVal <- length(NormBetasVal)
	NormBetasVal$TargetID <- rownames(NormBetasVal)
	NormBetasVal <- NormBetasVal[,c((lenNormBetasVal+1),1:lenNormBetasVal)]
	
	return(NormBetasVal)
}



#######################################
###############   All   ###############

HumMeth27QCReport <- function(ImportDataR, platform, pval, ChrX, ClustMethod, quoteOutput, normMethod)   {

	options(warn=1)
	options(showNCalls=400)

	if (missing(pval) || is.null(pval)) {
		pval <- 0.05
	}
	if (missing(ClustMethod) || is.null(ClustMethod)) {
		ClustMethod <- "euclidean"
	}
	if (missing(ChrX) || is.null(ChrX)) {
		ChrX <- "FALSE"
	}
	if (missing(quoteOutput) || is.null(quoteOutput)){
		quoteOutput <- TRUE
	}
        if (missing(normMethod) || is.null(normMethod)){
                normMethod <- "quantile"
        }


##### Internal controls
	dataVal <- getAssayControls(ImportDataR,platform)    
	
##### Quality controls
	dataQC <- QCCheck(ImportDataR, pval)
	infotab_DetectedGenes <- dataQC$Detection
	BadCpG_001_more5 <- dataQC$BadCpG_001_more5
	colnames(BadCpG_001_more5) <- c("CpG","Percentage")
	BadCpG_005_more5 <- dataQC$BadCpG_005_more5
	colnames(BadCpG_005_more5) <- c("CpG","Percentage")
	
##### Explorative analysis
	normVal <- NormCheck(ImportDataR, platform, pval, ChrX, ClustMethod, normMethod)
	write.table(normVal,"NormalizedMvalues.txt",sep="\t",row.names=FALSE,quote=quoteOutput)
	
	WriteXLS(c("dataVal", "infotab_DetectedGenes", "BadCpG_001_more5", "BadCpG_005_more5"), 
			 ExcelFileName = "QC_Analysis.xls", SheetNames = c("Internal Control", "Detected genes", "BadCpG_001_more5%", "BadCpG_005_more5%"), 
			 verbose = FALSE, Encoding = c("UTF-8", "latin1"), perl = "perl",BoldHeaderRow=T)
	
	return(normVal)
}
