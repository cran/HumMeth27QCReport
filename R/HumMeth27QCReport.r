
###############################################
###############   Import data   ###############
ImportData <- function(Dir) {
setwd(Dir);

WorkFiles <- list.files(path = ".",pattern = ".*.txt")
AverageBeta <- WorkFiles[grep("AvgBeta", WorkFiles)]
Ctrl <- WorkFiles[grep("Control", WorkFiles)]
SamplesName <- WorkFiles[grep("Sample", WorkFiles)]
try(Discarder <- WorkFiles[grep("Discard", WorkFiles)], silent=T)

control <- read.table(Ctrl, header = TRUE, sep = "\t") 
AvBeta <- read.table(AverageBeta, header = TRUE, sep = "\t", flush=T) 
samps <- read.table(SamplesName, header = TRUE, sep = "\t")
DiscarderII <- ""
try(DiscarderII <- read.table(Discarder, header = F, sep = "\t"), silent=T)

samps$SampleLabel <-  samps$Sample.ID
colnames(samps)[1] <- "Index"
colnames(samps)[2] <- "SampleID"

write.table(samps,"sample.txt", sep="\t",row.names=F)
samps_mod <- read.table("sample.txt", header = TRUE, sep = "\t")

samps2 <- samps[with(samps, order(SampleID)), ]
nsample <- length(samps$Index)

pdf(file="Sample.pdf",paper="a4r", fonts="Times")
SamplePDF <- data.frame(samps_mod$Index[1:(length(samps_mod$Index)/2)], samps_mod$SampleID[1:(length(samps_mod$SampleID)/2)],rep(" ",(length(samps_mod$SampleID)/2)), 
	samps_mod$Index[((length(samps_mod$Index)/2)+1):length(samps_mod$Index)], samps_mod$SampleID[((length(samps_mod$SampleID)/2)+1):length(samps_mod$SampleID)])
colnames(SamplePDF)<-c("Index","SampleID","","Index","SampleID")
textplot(SamplePDF, halign="center", valign="center", show.rownames = F)
title("Sample List")
dev.off()

return(list(ctrl=control,samples=samps, AverageBeta=AverageBeta, Ctrl=Ctrl, samps_mod=samps_mod, nsample=nsample, samps2=samps2, DiscarderII))
}


####################################################
###############   Internal Control   ###############
QCRep <- function(Dir) {
	
dataFiles <- ImportData(Dir)
control <- dataFiles$ctrl
samps <- dataFiles$samples

##### QUALITY CHECK - Bisulfite Conversion
Bisulfite_Sig1 <- control[control$ProbeID == "5270706",seq(4,length(control),3)]
Bisulfite_Sig1 <- as.data.frame(t(Bisulfite_Sig1))
Bisulfite_Sig2 <- control[control$ProbeID == "4670278",seq(4,length(control),3)]
Bisulfite_Sig2 <- as.data.frame(t(Bisulfite_Sig2))
Bisulfite_Bckg1 <- control[control$ProbeID == "5290048",seq(4,length(control),3)]
Bisulfite_Bckg1 <- as.data.frame(t(Bisulfite_Bckg1))
Bisulfite_Bckg2 <- control[control$ProbeID == "4670484",seq(4,length(control),3)]
Bisulfite_Bckg2 <- as.data.frame(t(Bisulfite_Bckg2))
BisulfiteRatio <- ((Bisulfite_Bckg1+Bisulfite_Bckg2) / (Bisulfite_Sig1+Bisulfite_Sig2))*100
Bisulfite <- data.frame("SAMPLE"=samps$SampleID,"BisulfiteRatio"=BisulfiteRatio)
Bisulfite <- Bisulfite[with(Bisulfite, order(SAMPLE)), ]

pdf(file="BisulfiteControl.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,4.1,4.1,2.1))
try(gap.barplot(Bisulfite[,2], gap=c(101,110), col=rep("blue",length(Bisulfite[,2])), xaxlab=Bisulfite$SAMPLE,main="Bisulfite Control: Background (U) on Signal (C)",las=2,ylim=c(0,120), ylab="%", ytics=c(seq(0,100,10),round(max(Bisulfite[,2]),0)-100),xlab=""), silent=T)
abline(h=seq(0,100,10),col="grey")
dev.off()


##### QUALITY CHECK - Extension
Extension.Ag <- control[control$ProbeID == "360446",seq(4,length(control),3)]
Extension.Ag <- data.frame("Sample"= samps$SampleID,"A"=t(Extension.Ag))
Extension.Ar <- control[control$ProbeID == "360446",seq(5,length(control),3)]
Extension.Ar <- data.frame("Sample"= samps$SampleID,"A"=t(Extension.Ar))
Extension.Tg <- control[control$ProbeID == "520537",seq(4,length(control),3)] 
Extension.Tg <- data.frame("Sample"= samps$SampleID,"T"=t(Extension.Tg))
Extension.Tr <- control[control$ProbeID == "520537",seq(5,length(control),3)] 
Extension.Tr <- data.frame("Sample"= samps$SampleID,"T"=t(Extension.Tr))
Extension.Gg <- control[control$ProbeID == "1190050",seq(4,length(control),3)]
Extension.Gg <- data.frame("Sample"= samps$SampleID,"G"=t(Extension.Gg))
Extension.Gr <- control[control$ProbeID == "1190050",seq(5,length(control),3)]
Extension.Gr <- data.frame("Sample"= samps$SampleID,"G"=t(Extension.Gr))
Extension.Cg <- control[control$ProbeID == "2630184",seq(4,length(control),3)]
Extension.Cg <- data.frame("Sample"= samps$SampleID,"C"=t(Extension.Cg))
Extension.Cr <- control[control$ProbeID == "2630184",seq(5,length(control),3)]
Extension.Cr <- data.frame("Sample"= samps$SampleID,"C"=t(Extension.Cr))

Extension.AT.g <- merge(Extension.Ag,Extension.Tg)
Extension.AT.r <- merge(Extension.Ar,Extension.Tr)
Extension.GC.g <- merge(Extension.Gg,Extension.Cg)
Extension.GC.r <- merge(Extension.Gr,Extension.Cr)

AvExt.AT.g <- apply(as.matrix(Extension.AT.g[,2:3]),1,mean)
AvExt.GC.g <- apply(as.matrix(Extension.GC.g[,2:3]),1,mean)
ExtRatio.g <- (AvExt.AT.g / AvExt.GC.g)*100 

AvExt.AT.r <- apply(as.matrix(Extension.AT.r[,2:3]),1,mean)
AvExt.GC.r <- apply(as.matrix(Extension.GC.r[,2:3]),1,mean)
ExtRatio.r <- (AvExt.GC.r / AvExt.AT.r)*100 

pdf(file="Extension_Green.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,4.1,4.1,2.1))
barplot2(as.table(ExtRatio.g), names=Extension.AT.g$Sample,ylim=c(0,max(ExtRatio.g)+0.5), col=3,las=2, 
	main="Extension Control (green channel): Background (AT) on Signal (GC)", ylab="%", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
dev.off()

pdf(file="Extension_Red.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,4.1,4.1,2.1))
barplot2(as.table(ExtRatio.r), names=Extension.AT.g$Sample,ylim=c(0, max(ExtRatio.r)+0.5),col=2,las=2, 
	main="Extension control (red channel): Background (GC) on Signal (AT)", ylab="%", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
dev.off()


##### QUALITY CHECK - Hybridization
HybL.g <- t(control[control$ProbeID == "2450040",seq(4,length(control),3)])
HybL.r <- t(control[control$ProbeID == "2450040",seq(5,length(control),3)])
HybM.g <- t(control[control$ProbeID == "5690110",seq(4,length(control),3)])
HybM.r <- t(control[control$ProbeID == "5690110",seq(5,length(control),3)])
HybH.g <- t(control[control$ProbeID == "5690072",seq(4,length(control),3)])
HybH.r <- t(control[control$ProbeID == "5690072",seq(5,length(control),3)])
RatioHybL <- round((HybL.r / HybL.g)*100,1)
RatioHybM <- round((HybM.r / HybM.g)*100,1)
RatioHybH <- round((HybH.r / HybH.g)*100,1)

Hyb <- data.frame("SAMPLE"=samps$SampleID,"L"=RatioHybL,"M"=RatioHybM,"H"=RatioHybH)
Hyb <- Hyb[with(Hyb, order(SAMPLE)), ]

pdf(file="Hibridization.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,4.1,3.1,1.1))
barplot2(as.matrix(t(data.frame("L"=Hyb[,2],"M"=Hyb[,3],"H"=Hyb[,4]))), names=Hyb$SAMPLE, ylim=c(0,max(RatioHybL)+1),ylab="%", 
	col=c(1,2,3), las=2, beside=T, legend = c("HybL","HybM","HybH"), main="Hibridization Control: Background (Red) on Signal (Green)", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
dev.off()


##### QUALITY CHECK - Target Removal
TM <- t(control[control$ProbeID == "580035",seq(4,length(control),3)])
TM <- data.frame("Sample"=samps$SampleID, TM)
colnames(TM) = c("SAMPLE","tm")
TM <- TM[with(TM, order(SAMPLE)), ]

pdf(file="TargetRemoval.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,3.1,3.1,1.1))
barplot2(TM[,2], names=TM$SAMPLE, col=c("blue"), las=2, ylim=c(0,max(TM[,2])+100), main="Target Removal Control on Green Channel", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
dev.off()


##### QUALITY CHECK - Negative
N1 <- t(control[control$ProbeID == "50110",seq(4,length(control),3)])
N2 <- t(control[control$ProbeID == "360079",seq(4,length(control),3)])
N3 <- t(control[control$ProbeID == "430114",seq(4,length(control),3)])
N4 <- t(control[control$ProbeID == "460494",seq(4,length(control),3)])
N5 <- t(control[control$ProbeID == "540577",seq(4,length(control),3)])
N6 <- t(control[control$ProbeID == "610692",seq(4,length(control),3)])
N7 <- t(control[control$ProbeID == "610706",seq(4,length(control),3)])
N8 <- t(control[control$ProbeID == "670750",seq(4,length(control),3)])
N9 <- t(control[control$ProbeID == "1190458",seq(4,length(control),3)])
N10 <- t(control[control$ProbeID == "1500059",seq(4,length(control),3)])
N11 <- t(control[control$ProbeID == "1500167",seq(4,length(control),3)])
N12 <- t(control[control$ProbeID == "1500398",seq(4,length(control),3)])
N13 <- t(control[control$ProbeID == "1660097",seq(4,length(control),3)])
N14 <- t(control[control$ProbeID == "1770019",seq(4,length(control),3)])
N15 <- t(control[control$ProbeID == "1940364",seq(4,length(control),3)])
N16 <- t(control[control$ProbeID == "1990692",seq(4,length(control),3)])

Negative <- data.frame(N1,N2,N3,N4,N5,N6,N7,N8,N9,N10,N11,N12,N13,N14,N15,N16)
Neg <- apply(as.matrix(Negative), 1,mean)
NegC <- data.frame("SAMPLE"=samps$SampleID,"Neg"=Neg)
NegC <- NegC[with(NegC, order(SAMPLE)), ]

pdf(file="Negative.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,3.1,3.1,1.1))
barplot2(NegC[,2], col="blue", names=NegC$SAMPLE, las=2, ylim=c(0,max(Neg)+20), main="Negative Control on Green Channel", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
dev.off()	


##### QUALITY CHECK - Non-Polymorphic (A+T / G+C)
NonPolymorphic.Ag <- control[control$ProbeID == "1740025",seq(4,length(control),3)]
NonPolymorphic.Ag <- data.frame("Sample"= samps$SampleID,"A"=t(NonPolymorphic.Ag))
NonPolymorphic.Tg <- control[control$ProbeID == "2480348",seq(4,length(control),3)]
NonPolymorphic.Tg <- data.frame("Sample"= samps$SampleID,"T"=t(NonPolymorphic.Tg))
NonPolymorphic.Gg <- control[control$ProbeID == "110184",seq(4,length(control),3)]
NonPolymorphic.Gg <- data.frame("Sample"= samps$SampleID,"G"=t(NonPolymorphic.Gg))
NonPolymorphic.Cg <- control[control$ProbeID == "2810035",seq(4,length(control),3)] 
NonPolymorphic.Cg <- data.frame("Sample"= samps$SampleID,"C"=t(NonPolymorphic.Cg))

NonPolymorphic.AT.g <- merge(NonPolymorphic.Ag,NonPolymorphic.Tg)
AvNP.AT.g <- apply(as.matrix(NonPolymorphic.AT.g[,2:3]),1,mean)
NonPolymorphic.GC.g <- merge(NonPolymorphic.Gg,NonPolymorphic.Cg)
AvNP.GC.g <- apply(as.matrix(NonPolymorphic.GC.g[,2:3]),1,mean)
NPRatio.g <- (AvNP.AT.g / AvNP.GC.g)*100 

pdf(file="NonPolimorphic_Green.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,4.1,3.1,1.1))
barplot2(as.table(NPRatio.g), names=NonPolymorphic.GC.g$Sample, col=3, las=2, ylim=c(0,max(NPRatio.g)+2),
	main="Non-Polymorphic Control: Background (AT) on Signal (GC) - Green Channel", ylab="%", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
dev.off()	

NonPolymorphic.Ar <- control[control$ProbeID == "1740025",seq(5,length(control),3)]
NonPolymorphic.Ar <- data.frame("Sample"= samps$SampleID,"A"=t(NonPolymorphic.Ar))
NonPolymorphic.Tr <- control[control$ProbeID == "2480348",seq(5,length(control),3)]
NonPolymorphic.Tr <- data.frame("Sample"= samps$SampleID,"T"=t(NonPolymorphic.Tr))
NonPolymorphic.Gr <- control[control$ProbeID == "110184",seq(5,length(control),3)]
NonPolymorphic.Gr <- data.frame("Sample"= samps$SampleID,"G"=t(NonPolymorphic.Gr))
NonPolymorphic.Cr <- control[control$ProbeID == "2810035",seq(5,length(control),3)] 
NonPolymorphic.Cr <- data.frame("Sample"= samps$SampleID,"C"=t(NonPolymorphic.Cr))

NonPolymorphic.AT.r <- merge(NonPolymorphic.Ar,NonPolymorphic.Tr)
AvNP.AT.r <- apply(as.matrix(NonPolymorphic.AT.r[,2:3]),1,mean)
NonPolymorphic.GC.r <- merge(NonPolymorphic.Gr,NonPolymorphic.Cr)
AvNP.GC.r <- apply(as.matrix(NonPolymorphic.GC.r[,2:3]),1,mean)
NPRatio.r <- (AvNP.GC.r / AvNP.AT.r)*100 

pdf(file="NonPolimorphic_Red.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,4.1,3.1,1.1))
barplot2(as.table(NPRatio.r), names=NonPolymorphic.AT.r$Sample, col=2, las=2, ylim=c(0,max(NPRatio.r)+1), 
	main="Non-Polymorphic Control: Background (GC) on Signal (AT) - Red Channel", ylab="%", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
dev.off()	


##### QUALITY CHECK - Staining DNP
DNP.med <- t(control[control$ProbeID == "4200736",seq(5,length(control),3)])
DNP.bgnd <- t(control[control$ProbeID == "5340168",seq(5,length(control),3)])
DNP.Ratio <- as.vector(DNP.bgnd / DNP.med)*100
DNP <- data.frame("SAMPLE"=samps$SampleID,"DNPRatio"=DNP.Ratio)
DNP <- DNP[with(DNP, order(SAMPLE)), ]

pdf(file="DNP.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,4.1,3.1,1.1))
barplot2(DNP[,2], names=DNP$SAMPLE, col="red", las=2, ylim=c(0,max(DNP.Ratio)+0.5), 
	main="Staining DNP Control: Background (BGND) on Signal (MED)", ylab="%", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
dev.off()	


##### QUALITY CHECK - Staining Biotin
Biotin.med <- t(control[control$ProbeID == "4570020",seq(4,length(control),3)])
Biotin.bgnd <- t(control[control$ProbeID == "5050601",seq(4,length(control),3)])
Biotin.Ratio <- as.vector(Biotin.bgnd / Biotin.med)*100
Biotin <- data.frame("SAMPLE"=samps$SampleID,"BiotinRatio"=Biotin.Ratio)
Biotin <- Biotin[with(Biotin, order(SAMPLE)), ]

pdf(file="Biotin.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,4.1,3.1,1.1))
barplot2(Biotin[,2], names=Biotin$SAMPLE, col=3, las=2, ylim=c(0,max(Biotin.Ratio)+1), 
	main="Staining Biotin Control: Background (BGND) on Signal (MED)", ylab="%", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
dev.off()	


##### QUALITY CHECK - Specificity (mismatch 1)
Specificity.mis1.PM <- t(control[control$ProbeID == "4610725",seq(5,length(control),3)])
Specificity.mis1.MM <- t(control[control$ProbeID == "4610400",seq(5,length(control),3)])
Specificity.mis1.Ratio <- as.vector(Specificity.mis1.MM / Specificity.mis1.PM)*100
Spec_mm1 <- data.frame("SAMPLE"=samps$SampleID,"Specificity.mis1.Ratio"=Specificity.mis1.Ratio)
Spec_mm1 <- Spec_mm1[with(Spec_mm1, order(SAMPLE)), ]

pdf(file="Specificity_mm1.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,4.1,3.1,1.1))
barplot2(Spec_mm1[,2], names=Spec_mm1$SAMPLE, col=3, las=2, ylim=c(0,max(Specificity.mis1.Ratio)+1), 
	main="Specificity Control mismatch 1: Background (MM) on Signal (PM)", ylab="%", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
dev.off()


##### QUALITY CHECK - Specificity (mismatch 2)
Specificity.mis2.PM <- t(control[control$ProbeID == "3800154",seq(4,length(control),3)])
Specificity.mis2.MM <- t(control[control$ProbeID == "3800086",seq(4,length(control),3)])
Specificity.mis2.Ratio <- as.vector(Specificity.mis2.MM / Specificity.mis2.PM)*100
Spec_mm2 <- data.frame("SAMPLE"=samps$SampleID,"Specificity.mis2.Ratio"=Specificity.mis2.Ratio)
Spec_mm2 <- Spec_mm2[with(Spec_mm2, order(SAMPLE)), ]

pdf(file="Specificity_mm2.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,4.1,3.1,1.1))
barplot2(Spec_mm2[,2], names=Spec_mm2$SAMPLE, col=3, las=2, ylim=c(0,max(Specificity.mis1.Ratio)+1), 
	main="Specificity Control mismatch 2: Background (MM) on Signal (PM)", ylab="%", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
dev.off()

return(list(Bisulfite=Bisulfite, ExtRatio.g=ExtRatio.g, ExtRatio.r=ExtRatio.r, Extension.AT.g=Extension.AT.g, Hyb=Hyb, TM=TM, NegC=NegC, NPRatio.g=NPRatio.g, NPRatio.r=NPRatio.r,
	DNP=DNP, Biotin=Biotin, Spec_mm1=Spec_mm1, Spec_mm2=Spec_mm2))
}



#################################################
###############   Quality Check   ###############

QCCheck <- function(Dir)   {

dataFiles <- ImportData(Dir)
control <- dataFiles$ctrl
samps <- dataFiles$samples
AverageBeta <- dataFiles$AverageBeta
Ctrl <- dataFiles$Ctrl
samps_mod <- dataFiles$samps_mod
nsample <- dataFiles$nsample

mldat <- methylumiR(AverageBeta, qcfile = Ctrl, sampleDescriptions = samps_mod)

##### Intensity Graphs
intNumber <- nsample / 12;
for (k in 0:(intNumber-1)) 
	{
	fileN <- paste("int",k+1,".pdf", sep = "", collapse = NULL)

	pdf(file=fileN,width=15,height=9)
	par(mfrow=c(4,3),oma=c(0.1,0.1,0.1,0.1), mar=c(8.1,2.1,2.1,2.1))
	for(i in ((1+12*k):(12*(k+1)))){ 
		try(plotSampleIntensities(mldat,s=i), silent=T)      
		}
	dev.off()
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

NoDetectedAll <- data.frame("SAMPLE"=samps$SampleID,"pv001"=as.numeric(NoDetect001),"pv005"=as.numeric(NoDetect005))
NoDetectedAll <- NoDetectedAll[with(NoDetectedAll, order(SAMPLE)), ]

pdf(file="NoDetect.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,4.1,3.1,1.1))
barplot2(rbind(as.numeric(NoDetectedAll[,2]),as.numeric(NoDetectedAll[,3])),beside=T, las = 2, 
	col= c("blue","red"), names=NoDetectedAll$SAMPLE, axes=T, main="Percentage of non detected genes", ylab="%", 
	ylim=c(0, max(as.numeric(NoDetect001)+2)), legend=c("% Not detected Genes (0.05)","% Not detected Genes (0.01)"), plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
dev.off()


#### Average P-value
avgPval <- colMeans(pvals(mldat),na.rm = T)

pdf(file="Avg_Pvalue.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,4.1,3.1,1.1))
barplot2(avgPval, ylab = "Average P-Value", las = 2, ylim=c(0,max(avgPval)+0.1),col=4, main="Average p-value", plot.grid = TRUE, grid.lty = "solid",grid.col = "grey")
abline(h=0.05, col="red",lty=2)
dev.off()

return(list(betasV=mldat,avgPval=avgPval,NoDetect005=NoDetect005,NoDetect001=NoDetect001 ))
}


# PCA ANALYSIS
pca.samples.plot<-function(data,int.col=(1:ncol(data)),main.str="Principal Component Analysis")
{
    data.nona<-na.delete(data[,int.col])
    pca.res<-prcomp(t(data.nona),tol=0.1,na.action=na.omit)
    plot(pca.res$x,xlab="PC 1",ylab="PC 2",main=main.str,pch=".")
    text(pca.res$x,colnames(data[,int.col]),cex=0.7)
}


# Cluster Analysis on normalized data
cluster.samples.plot<-function(data,int.col=(1:ncol(data)),main.str="", method = "euclidian")
{
	# Matrix reprensentation of hierarchical clustering of the samples
    
    data.nona<-na.delete(data[,int.col])
    num.data <-t(data.nona)
    
    dist.mat <-as.matrix(Dist(num.data, method =method, diag =TRUE, upper =FALSE))
    resclust <-hcluster(dist.mat)
    dist.mat <-dist.mat[resclust$order,resclust$order]

    old.mar<-par()$mar
    par(mar=old.mar+c(2,1,-3,-2))    
    par(mfrow=c(2,1))
        image(x=seq(1,nrow(num.data)),y=seq(1,nrow(num.data)),z=dist.mat,axes=FALSE,xlab="",ylab="",main=main.str)
            axis(1,at=seq(1,nrow(num.data)),labels=rownames(dist.mat),las=2,cex.axis=0.5)
            axis(2,at=seq(1,nrow(num.data)),labels=rownames(dist.mat),las=2,cex.axis=0.25)
        plot(resclust,main="",xlab="",ylab="",cex=.6,axes=FALSE)
    par(mar=old.mar,mfrow=c(1,1))    
}


#################################################
###############   Normalization   ###############

# The function, normalizeMethyLumiSet does this normalization.
# Basically, it looks at the median intensitiesin the methylated and unmethylated channels at very low and very high beta values and sets these medians equal. 
# Using the transformed unmethylated and methylated values, new beta values are calculated using one of two "map" functions. 
# The ratio function is the default and is the same as used by Illumina in the BeadStudio software, but values using the atan selection should be similar.
####

NormCheck <- function(Dir)   {

dataFiles <- ImportData(Dir)
AverageBeta <- dataFiles$AverageBeta
Ctrl <- dataFiles$Ctrl
samps_mod <- dataFiles$samps_mod
try(DiscarderII <- dataFiles$DiscarderII, silent=T)

mldat <- methylumiR(AverageBeta, qcfile = Ctrl, sampleDescriptions = samps_mod)
avgPval <- colMeans(pvals(mldat),na.rm = T)
	
toKeep <- (avgPval < 0.01)
try(toKeep[as.character(DiscarderII$V1)]<-FALSE, silent=T)

mldat.norm <- normalizeMethyLumiSet(mldat[,toKeep])

pdf(file="PCA.pdf",width=15,height=9)
par(oma=c(1,1,1,1), mar=c(6.1,3.1,3.1,1.1))
pca.samples.plot(as.data.frame(betas(mldat.norm)))
dev.off()

pdf(file="Cluster.pdf",width=15,height=9)
cluster.samples.plot(as.data.frame(betas(mldat.norm)),main="Hierarchical Clustering - Euclidean")
dev.off()
return(list(mldat.norm=mldat.norm, toKeep=toKeep))

NormBetasVal <- as.data.frame(betas(mldat.norm[toKeep,]))
WriteXLS("NormBetasVal", ExcelFileName = "QC_Analysis.xls", SheetNames = "Normalized Beta Values", 
	verbose = FALSE, Encoding = c("UTF-8", "latin1"), perl = "perl",BoldHeaderRow=T)
}


#######################################
###############   All   ###############

HumMeth27QCReport <- function(Dir)   {

dataFiles <- ImportData(Dir)
AverageBeta <- dataFiles$AverageBeta
Ctrl <- dataFiles$Ctrl
samps <- dataFiles$samples
samps_mod <- dataFiles$samps_mod
nsample <- dataFiles$nsample
samps2 <- dataFiles$samps2
	
dataVal <- QCRep(Dir)
Bisulfite <- dataVal$Bisulfite
ExtRatio.g <- dataVal$ExtRatio.g
ExtRatio.r <- dataVal$ExtRatio.r
Extension.AT.g <- dataVal$Extension.AT.g
Hyb <- dataVal$Hyb
TM <- dataVal$TM
NegC <- dataVal$NegC
NPRatio.g <- dataVal$NPRatio.g
NPRatio.r <- dataVal$NPRatio.r
DNP <- dataVal$DNP
Biotin <- dataVal$Biotin
Spec_mm1 <- dataVal$Spec_mm1
Spec_mm2 <- dataVal$Spec_mm2

dataQC <- QCCheck(Dir)
mldat <- dataQC$betasV
avgPval <- dataQC$avgPval
NoDetect005 <- dataQC$NoDetect005
NoDetect001 <- dataQC$NoDetect001

normVal <- NormCheck(Dir)
mldat.norm <- normVal$mldat.norm
toKeep <- normVal$toKeep

#### Create the excel file with the summary 

NormBetasVal <- as.data.frame(betas(mldat.norm[toKeep,]))
lenNormBetasVal <- length(NormBetasVal)
NormBetasVal$TargetID <- rownames(NormBetasVal)
NormBetasVal <- NormBetasVal[,c((lenNormBetasVal+1),1:lenNormBetasVal)]

infotab_intC_part1 <- data.frame("Index"=samps2$Index, "Sample ID"= samps2$SampleID, "Sample Group" = samps2$Sample.Group, "Sentrix Barcode"=samps2$Sentrix.Barcode,
		"Sample Section"=samps2$Sample.Section, "Sample_Plate"=samps2$Sample_Plate, "Sample_Well"=samps2$Sample_Well, 
		"Bisulfite Conversion Control" = Bisulfite[,2], "Low Hybridization Control" = Hyb[,2], 
		"Medium Hybridization Control" = Hyb[,3], "High Hybridization Control" = Hyb[,4],
		"Target Removal Control" = TM[,2], "Negative Control" = NegC[,2],"DNP staining" = DNP[,2], "Biotin staining" = Biotin[,2], 
		"Specificity mismatch 1" = Spec_mm1[,2], "Specificity mismatch 2" = Spec_mm2[,2])

infotab_intC_part2 <- data.frame("Sample ID"=Extension.AT.g$Sample,"Extension Control - Green" =ExtRatio.g,"Extension Control - Red" = ExtRatio.r,
		"Non-Polymorphic - Green" = NPRatio.g,"Non-Polymorphic - Red" = NPRatio.r)
			
infotab_intC <- merge(infotab_intC_part1,infotab_intC_part2)
infotab_intC <- data.frame("Index"=infotab_intC$Index, "Sample ID"=infotab_intC$Sample.ID, "Sample Group"=infotab_intC$Sample.Group, "Sentrix Barcode"=infotab_intC$Sentrix.Barcode,
		"Sample Section"=infotab_intC$Sample.Section, "Sample_Plate"=infotab_intC$Sample_Plate, "Sample_Well"=infotab_intC$Sample_Well,
		"Bisulfite Conversion Control" = round(infotab_intC$Bisulfite.Conversion.Control,1),"Low Hybridization Control" = round(infotab_intC$Low.Hybridization.Control,1), 
		"Medium Hybridization Control" = round(infotab_intC$Medium.Hybridization.Control,1), "High Hybridization Control" = round(infotab_intC$High.Hybridization.Control,1),
		"Target Removal Control" = round(infotab_intC$Target.Removal.Control,1),"Negative Control" = round(infotab_intC$Negative.Control,1), "DNP staining" = round(infotab_intC$DNP.staining,1), 
		"Biotin staining" = round(infotab_intC$Biotin.staining,1), "Specificity mismatch 1" = round(infotab_intC$Specificity.mismatch.1,1), "Specificity mismatch 2" = round(infotab_intC$Specificity.mismatch.2,1),
		"Extension.Control.Green" = round(infotab_intC$Extension.Control...Green,1), "Extension.Control.Red" = round(infotab_intC$Extension.Control...Red,1),
		"Non.Polymorphic.Green" = round(infotab_intC$Non.Polymorphic...Green,1), "Non.Polymorphic.Red" = round(infotab_intC$Non.Polymorphic...Red,1))

infotab_intC <- infotab_intC[with(infotab_intC, order(Index)), ]

infotab_DetectedGenes <- data.frame("Index"=samps$Index, "Sample ID"= samps$SampleID, "Sample Group" = samps$Sample.Group, "Sentrix Barcode"=samps$Sentrix.Barcode,
		"Sample Section"=samps$Sample.Section, "Sample_Plate"=samps$Sample_Plate, "Sample_Well"=samps$Sample_Well, "Percentage of non detected genes (pval>0.01)"=NoDetect001,
		"Percentage of non detected genes (pval>0.05)"=NoDetect005, "Average Detection p-value"=avgPval)
		
pvalSample_CpG <- as.data.frame(pvals(mldat))		
BadCpG_001 <- apply(pvalSample_CpG,1,function(x)sum(x > 0.01))		
BadCpG_001_percent <- BadCpG_001/length(colnames(pvalSample_CpG))*100
BadCpG_001_5to10 <- as.data.frame(BadCpG_001_percent[BadCpG_001_percent<10 & BadCpG_001_percent>5])
BadCpG_001_5to10 <- data.frame("CpG"=rownames(BadCpG_001_5to10),BadCpG_001_5to10) 
BadCpG_001_more10 <- as.data.frame(BadCpG_001_percent[BadCpG_001_percent>10])		
BadCpG_001_more10 <- data.frame("CpG"=rownames(BadCpG_001_more10),BadCpG_001_more10) 

BadCpG_005_count <- apply(pvalSample_CpG,1,function(x)sum(x > 0.05))
BadCpG_005_percent <- BadCpG_005_count/length(colnames(pvalSample_CpG))*100
BadCpG_005_5to10 <- as.data.frame(BadCpG_005_percent[BadCpG_005_percent<10 & BadCpG_005_percent>5])
BadCpG_005_5to10 <- data.frame("CpG"=rownames(BadCpG_005_5to10),BadCpG_005_5to10) 
BadCpG_005_more10 <- as.data.frame(BadCpG_005_percent[BadCpG_005_percent>10])		
BadCpG_005_more10 <- data.frame("CpG"=rownames(BadCpG_005_more10),BadCpG_005_more10) 

WriteXLS(c("NormBetasVal","infotab_intC", "infotab_DetectedGenes", "BadCpG_001_5to10", "BadCpG_001_more10", "BadCpG_005_5to10", "BadCpG_005_more10"), 
	ExcelFileName = "QC_Analysis.xls", SheetNames = c("Normalized Beta Values","Internal Control", "Detected genes", "BadCpG_001_5to10", "BadCpG_001_more10", "BadCpG_005_5to10", "BadCpG_005_more10"), 
	verbose = FALSE, Encoding = c("UTF-8", "latin1"), perl = "perl",BoldHeaderRow=T)
}
