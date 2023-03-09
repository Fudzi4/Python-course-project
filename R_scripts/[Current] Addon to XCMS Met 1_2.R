#
#                             Modification to Metabolomics Script v1.2
#
#seting working directory to Bruker results####
#
if(winDialog(type = "yesno",
             "Set WD to 'Bruker results' (Yes) or select manualy (No)?") == "YES")    #Setting working directory
{WD <- "b"} else(WD <- "c")
if(WD=="b"){
  setwd("C:/Bruker results")}
if(WD=="c"){
  setwd(choose.dir(default = "", caption = "Select a folder to store results"))}
#
#Select Polarity####
if(winDialog(type = "yesno", "Positive mode (Yes) or negative mode (No) analysis?") == "YES") 
{Mode <- "positive"} else(Mode <- "negative")

if(winDialog(type = "yesno", "Do you have QC and IS (Yes) or (No)") == "YES") 
{QC.mode <- "yes"} else(QC.mode <- "no") 


###Loading dataset####
#Creating a data frame from XCMS output
output_Final<- read.table("output_Final.tsv", header = T, sep = " ", dec = ".")
#Creating a data frame from CAMERA output
CAMERA<-read.table("Result_CAMERA.csv", header = T, sep = ",", dec = ".")

#Add CAMERA-unique output to separate dataframe ####
CAMERA_unique<- data.frame(CAMERA$isotopes,CAMERA$adduct,CAMERA$pcgroup,CAMERA$mzmin)
colnames(CAMERA_unique)<-c("isotopes","adduct","pcgroup","mzmin")
#Merging dataframes####
output_Final<-output_Final[order(output_Final$mzmin),]
CAMERA_unique<-CAMERA_unique[order(CAMERA_unique$mzmin),]
XCMS_merged<-merge(output_Final,unique(CAMERA_unique),by="mzmin")

#Rearranging columns independently of number of samples####
XCMS_rearranged<-XCMS_merged[,c(2:7,1,8:10,(ncol(XCMS_merged)-2):
                                  ncol(XCMS_merged),11:(ncol(XCMS_merged)-3))]
#colnames(XCMS_rearranged)


if (QC.mode=="yes"){

#creating a subset of QCs
QCs<-XCMS_rearranged[,base::grep("QC_",colnames(XCMS_rearranged))]
#Creating a subset of IS for visualization
#subset of all samples
Samples<-XCMS_rearranged[,21:ncol(XCMS_rearranged)]
#Selecting Polarity mode for IS search
if(Mode=="positive")
{  # IS sampling updated for positive mode
  ############tyrosine search
  #Sampling all features fitting mz and RT criteria
  IS_tyrosine<-Samples[which(XCMS_rearranged$mzmed>192.1081&XCMS_rearranged$mzmed<192.1087&
                               XCMS_rearranged$rtmed>264&XCMS_rearranged$rtmed<312),]
  #Adding average column
  IS_tyrosine$AV<-rowMeans(IS_tyrosine, na.rm=T)
  #Overwrite data frame with only one feature with highest average signal
  IS_tyrosine<-IS_tyrosine[which(IS_tyrosine$AV==max(IS_tyrosine$AV)),]
  #remove average column
  IS_tyrosine<-IS_tyrosine[,-which(names(IS_tyrosine)=="AV")]
  #name row according to the IS
  rownames(IS_tyrosine)<-c("tyrosine")
  ############phenylalanine search
  IS_phenylalanine<-Samples[which(XCMS_rearranged$mzmed>176.1133 & XCMS_rearranged$mzmed<176.1137&
                                    XCMS_rearranged$rtmed>360&XCMS_rearranged$rtmed<420),]
  IS_phenylalanine$AV<-rowMeans(IS_phenylalanine, na.rm=T)
  IS_phenylalanine<-IS_phenylalanine[which(IS_phenylalanine$AV==max(IS_phenylalanine$AV)),]                        
  IS_phenylalanine<-IS_phenylalanine[,-which(names(IS_phenylalanine)=="AV")]
  rownames(IS_phenylalanine)<-c("phenylalanine")
  ############valine search
  IS_valine<-Samples[which(XCMS_rearranged$mzmed>123.1028 & XCMS_rearranged$mzmed<123.1032&
                             XCMS_rearranged$rtmed>110&XCMS_rearranged$rtmed<170),]
  IS_valine$AV<-rowMeans(IS_valine, na.rm=T)
  IS_valine<-IS_valine[which(IS_valine$AV==max(IS_valine$AV)),]                        
  IS_valine<-IS_valine[,-which(names(IS_valine)=="AV")]
  rownames(IS_valine)<-c("valine") 
}
if(Mode=="negative")
{  # IS sampling updated for negative mode
  ############tyrosine search
  #Sampling all features fitting mz and RT criteria
  IS_tyrosine<-Samples[which(XCMS_rearranged$mzmed>190.093 & XCMS_rearranged$mzmed<190.095&
                               XCMS_rearranged$rtmed>264&XCMS_rearranged$rtmed<312),]
  #Adding average column
  IS_tyrosine$AV<-rowMeans(IS_tyrosine, na.rm=T)
  #Overwrite data frame with only one feature with highest average signal
  IS_tyrosine<-IS_tyrosine[which(IS_tyrosine$AV==max(IS_tyrosine$AV)),]
  #remove average column
  IS_tyrosine<-IS_tyrosine[,-which(names(IS_tyrosine)=="AV")]
  #name row according to the IS
  rownames(IS_tyrosine)<-c("tyrosine")
  ############phenylalanine search
  IS_phenylalanine<-Samples[which(XCMS_rearranged$mzmed>174.098 & XCMS_rearranged$mzmed<174.100&
                                    XCMS_rearranged$rtmed>360&XCMS_rearranged$rtmed<420),]
  IS_phenylalanine$AV<-rowMeans(IS_phenylalanine, na.rm=T)
  IS_phenylalanine<-IS_phenylalanine[which(IS_phenylalanine$AV==max(IS_phenylalanine$AV)),]                        
  IS_phenylalanine<-IS_phenylalanine[,-which(names(IS_phenylalanine)=="AV")]
  rownames(IS_phenylalanine)<-c("phenylalanine")
  ############valine search
  IS_valine<-Samples[which(XCMS_rearranged$mzmed>121.088 & XCMS_rearranged$mzmed<121.090&
                             XCMS_rearranged$rtmed>110&XCMS_rearranged$rtmed<170),]
  IS_valine$AV<-rowMeans(IS_valine, na.rm=T)
  IS_valine<-IS_valine[which(IS_valine$AV==max(IS_valine$AV)),]                        
  IS_valine<-IS_valine[,-which(names(IS_valine)=="AV")]
  rownames(IS_valine)<-c("valine")
}

############Combining of IS into 1 dataframe
IS<-rbind(IS_tyrosine,IS_phenylalanine,IS_valine)
#transposing of IS-all samples dataset
IS_t<-data.frame(t(IS))
rownames(IS_t)<-colnames(IS)
colnames(IS_t)<-rownames(IS)
#Statistics of IS-all samples dataset
IS_stat_1<-sd(IS_t$tyrosine)/mean(IS_t$tyrosine)
IS_stat_2<-sd(IS_t$phenylalanine)/mean(IS_t$phenylalanine)
IS_stat_3<-sd(IS_t$valine)/mean(IS_t$valine)

IS_stat<-IS_t
IS_stat$tyrosine<-(IS_stat$tyrosine/mean(IS_stat$tyrosine))-1
IS_stat$phenylalanine<-(IS_stat$phenylalanine/mean(IS_stat$phenylalanine))-1
IS_stat$valine<-(IS_stat$valine/mean(IS_stat$valine))-1

#Boxplot of IS with statistics####
par(mfrow = c(1,1))
png("IS-all samples-boxplot.png", width = 1000, height = 700, res = 100)
boxplot(IS_stat,main="Internal standards in all samples",ylab="Residuals")
dev.off()
boxplot(IS_stat,main="Internal standards in all samples",ylab="Residuals")

#IS Scatterplot ####
#Adding elution order
IS_t$index<-rownames(IS_t)
IS_t$index=substr(IS_t$index,nchar(IS_t$index)-3,nchar(IS_t$index))
IS_t$index<-as.numeric(IS_t$index)
#Sorting acording to elution order
IS_t<-IS_t[order(IS_t$index),]
IS_t$index<-(1:nrow(IS_t))

library(ggplot2)
library(reshape2)
IS_plot<-melt(IS_t, id.vars = "index", variable.name = "series")
#Scatter plot of internal standards
plot(IS_plot$index,IS_plot$value,          #x variable, y variable
     col=IS_plot$series,                   #colour by IS compound name
     pch=16,                               #type of points to use
     cex=1.2,                              #size of points to use
     xlab="Elution order",                 #x axis label
     ylab="Intensity",                     #y axis label
     main="Internal standards stability")  #plot title
legend(x=0,y=(max(IS_plot$value)*1.03), legend = levels(IS_plot$series),col = c(1:3), pch = 16)
lines(predict(lm(IS_plot$value[IS_plot$series=="tyrosine"]~
                   log(IS_plot$index[IS_plot$series=="tyrosine"]))), col=1)
lines(predict(lm(IS_plot$value[IS_plot$series=="phenylalanine"]~
                   log(IS_plot$index[IS_plot$series=="phenylalanine"]))), col=2)
lines(predict(lm(IS_plot$value[IS_plot$series=="valine"]~
                   log(IS_plot$index[IS_plot$series=="valine"]))), col=3)
#Saving the plot  
png("IS-all samples-scatterplot.png", width = 1000, height = 700, res = 100)
plot(IS_plot$index,IS_plot$value,          #x variable, y variable
     col=IS_plot$series,                   #colour by IS compound name
     pch=16,                               #type of points to use
     cex=1.2,                              #size of points to use
     xlab="Elution order",                 #x axis label
     ylab="Intensity",                     #y axis label
     main="Internal standards stability")  #plot title
legend(x=0,y=(max(IS_plot$value)*1.03), legend = levels(IS_plot$series),col = c(1:3), pch = 16)
lines(predict(lm(IS_plot$value[IS_plot$series=="tyrosine"]~
                   log(IS_plot$index[IS_plot$series=="tyrosine"]))), col=1)
lines(predict(lm(IS_plot$value[IS_plot$series=="phenylalanine"]~
                   log(IS_plot$index[IS_plot$series=="phenylalanine"]))), col=2)
lines(predict(lm(IS_plot$value[IS_plot$series=="valine"]~
                   log(IS_plot$index[IS_plot$series=="valine"]))), col=3)
dev.off()

#Add QC_AV#####
XCMS_rearranged$QC_AV<-rowMeans(QCs)
#Add CV%####
library(matrixStats)
XCMS_rearranged$QC_CV<-rowSds(as.matrix(QCs))/XCMS_rearranged$QC_AV

#QC CV visualization
png("QC_CV.png", width = 1000, height = 700, res = 100)
par(mfrow = c(1,2))
boxplot(XCMS_rearranged$QC_CV,main="QC CV%")  
hist(XCMS_rearranged$QC_CV,main="QC CV%",xlab = "")
dev.off()
par(mfrow = c(1,2))
boxplot(XCMS_rearranged$QC_CV,main="QC CV%")  
hist(XCMS_rearranged$QC_CV,main="QC CV%",xlab = "")

#Rearranging columns to put QC_AV and QC_CV in proper position
XCMS_rearranged<-XCMS_rearranged[,c(1:18,(ncol(XCMS_rearranged)-1):
                                  ncol(XCMS_rearranged),19:(ncol(XCMS_rearranged)-2))]
#colnames(XCMS_rearranged)

#Writing csv for unfiltered data####
write.csv(XCMS_rearranged,"XCMS_full.csv", row.names = T)
#Filtering
#Creating a filtered dataframe
XCMS_filtered<-subset.data.frame(x=XCMS_rearranged,subset = rtmed>60&rtmed<1080&QC_AV>20000&QC_CV<0.3)
}

if (QC.mode=="no"){
  #Writing csv for unfiltered data####
  write.csv(XCMS_rearranged,"XCMS_full.csv", row.names = T)
  #Filtering
  #Creating a filtered dataframe
  XCMS_filtered<-subset.data.frame(x=XCMS_rearranged,subset = rtmed>60&rtmed<1080)
}

#Isotopic filtering 
####Selecting Mode

if(Mode=="positive")
{
  #Isotopic filtering for positive mode
  XCMS_isoM<-XCMS_filtered[base::grep("M]",XCMS_filtered$isotopes), ]
  XCMS_isoM2<-XCMS_isoM[base::grep("M]2",XCMS_isoM$isotopes),]
  XCMS_isoM3<-XCMS_isoM[base::grep("M]3",XCMS_isoM$isotopes),]
  XCMS_isoM4<-XCMS_isoM[base::grep("M]4",XCMS_isoM$isotopes),]
  
  if(nrow(XCMS_isoM2)>0)
  {XCMS_isoM<-XCMS_isoM[- base::grep("M]2",XCMS_isoM$isotopes),]}
  if(nrow(XCMS_isoM3)>0)
  {XCMS_isoM<-XCMS_isoM[- base::grep("M]3",XCMS_isoM$isotopes),]}
  if(nrow(XCMS_isoM4)>0)
  {XCMS_isoM<-XCMS_isoM[- base::grep("M]4",XCMS_isoM$isotopes),]}
  XCMS_isoN<-XCMS_filtered[XCMS_filtered$isotopes=="", ]
  XCMS_iso<-rbind(XCMS_isoM,XCMS_isoN)
}
if(Mode=="negative")
{
  #Isotopic filtering for negative mode
  XCMS_isoM<-XCMS_filtered[base::grep("M]",XCMS_filtered$isotopes), ]
  XCMS_isoM2<-XCMS_isoM[base::grep("M]2",XCMS_isoM$isotopes),]
  XCMS_isoM3<-XCMS_isoM[base::grep("M]3",XCMS_isoM$isotopes),]
  XCMS_isoM4<-XCMS_isoM[base::grep("M]4",XCMS_isoM$isotopes),]
  
  if(nrow(XCMS_isoM2)>0)
  {XCMS_isoM<-XCMS_isoM[- base::grep("M]2",XCMS_isoM$isotopes),]}
  if(nrow(XCMS_isoM3)>0)
  {XCMS_isoM<-XCMS_isoM[- base::grep("M]3",XCMS_isoM$isotopes),]}
  if(nrow(XCMS_isoM4)>0)
  {XCMS_isoM<-XCMS_isoM[- base::grep("M]4",XCMS_isoM$isotopes),]}
  XCMS_isoN<-XCMS_filtered[XCMS_filtered$isotopes=="", ]
  XCMS_iso<-rbind(XCMS_isoM,XCMS_isoN)
}

#Writing csv for filtered data####
write.csv(XCMS_iso,"XCMS_filtered.csv")

