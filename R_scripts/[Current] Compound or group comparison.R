#
#                             Script for Lauras compound detection and ploting
#
library(ggplot2)                                                                #Loading libraries
###seting working directory
if(winDialog(type = "yesno",
             "Set WD to 'Bruker results' (Yes) or select manualy (No)?") == "YES")    #Setting working directory
{WD <- "b"} else(WD <- "c")
if(WD=="b"){
  setwd("C:/Bruker results")}
if(WD=="c"){
  setwd(choose.dir(default = "", caption = "Select a folder to store results"))}

Main.table<-read.csv(choose.files(caption = "Select Standards.csv"),            #Manual loading of search output
                     sep = ",", header = T,dec = ".")
#Main.table<- read.table("Standards.csv", header = T, sep = ",", dec = ".")      #Loading compound search output
Main.table<-data.frame(t(Main.table))                                           #Tr of table for compound grouping

groups<-c("Clindamycin","Tipiracil","QC","Trifluridine")                        #ENTER DESIRED GROUPS HERE!!!

if(winDialog(type = "yesno",
             "Comparison between Groups (Yes) or Compounds (No)?") == "YES")    #Asking for Way of analysis
  {Mode <- "g"} else(Mode <- "c")

if(Mode=="c")
  {                                                                             #Compound comparison

  for (i in 1:length(groups))                                                   #Loop for grouping
    {
    name<-groups[i]
    temp<-data.frame(Main.table[grepl(name,row.names(Main.table)),])
    plot_name<-name                                                             #Boxplot
    png(paste("Cc - ",name,".png",sep = ""),width = 1500,height = 1000,res = 100)
    par(mfrow=c(1,length(temp)))
    #Loop for boxplots
    for (x in 1:length(temp)) 
      {
      temp2<-temp[,x]
      boxplot(temp2,xlab=names(temp[x]))
      }
    mtext(name,side = 3,line = -2,outer = TRUE)
    dev.off()
    }
  }

if(Mode=="g")                                                                    #Group comparison
  {
  for (i in 1:length(groups)) {                                                 #Addition column with groups
  Main.table$group[grepl(groups[i],row.names(Main.table))]<-groups[i]
  }
  
  for (i in 1:(length(Main.table)-1)){
                                                                                #Boxplot with individual intensities#####
   # plot_name<-name
   # png(paste(name,".png",sep = ""),width = 1500,height = 1000,res = 100)
   # par(mfrow=c(1,length(groups)))
    #Loop for boxplots
    #for (x in 1:length(groups)) {
    #  current_group<-groups[x]
    #boxplot(Main.table[grepl(current_group,row.names(Main.table)),x], 
    #xlab=current_group)
                                                                                #Boxplot all groups together
                                                                                #Standart boxplot####
    #plot_name<-colnames(Main.table[i])
    #png(paste(plot_name,".png",sep = ""),width = 1500,height = 1000,res = 100)
    #par(mfrow=c(1,1))
    #boxplot
    #boxplot(Main.table[,i]~Main.table$group,
    #        main=colnames(Main.table[i]),
    #        ylab = "intensity",xlab = "groups")
    
    #dev.off()
                                                                                #ggplot####
    plot_name<-colnames(Main.table[i])
    par(mfrow=c(1,1))
    
    p<-ggplot(Main.table, aes(x=group,y=Main.table[,i])) +                      #Plotting
      geom_boxplot() + 
      geom_jitter(shape=16,position = position_jitter(0.2),size=2)+
      labs(title = colnames(Main.table[i]),x="groups", y="intensity" )+
      theme_classic()+
      theme(plot.title=element_text(hjust=0.5,size=14, face = "bold"),
            axis.text = element_text(size = 10,color = "black"),
            axis.title = element_text(size = 10,face = "bold"))
    print(p)                                                                    #Visualizing plot
    
    png(paste("Gc - ",plot_name,".png",sep = ""),width = 500,height = 400,res = 100)    #saving png
    print(p)
    dev.off()
    }
  } 
winDialog(type = "ok", "R script execution is finished")                        #Message in the end
