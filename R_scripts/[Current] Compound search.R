#
#                             Control sample modification
###Intro#####
###seting working directory
if(winDialog(type = "yesno",
             "Set WD to 'Bruker results' (Yes) or select manualy (No)?") == "YES")    #Setting working directory
{WD <- "b"} else(WD <- "c")
if(WD=="b"){
  setwd("C:/Bruker results")}
if(WD=="c"){
  setwd(choose.dir(default = "", caption = "Select a folder to store results"))}

###Loading dataset
#Creating a data frame from XCMS output
Main.table<-read.csv(choose.files(caption = "Select XCMS output file"),         #Ugly window/but captions
                    sep = " ", header = T,dec = ".")
#Main.table<- read.table("output_Final.tsv", header = T, sep = " ", dec = ".")

###Subset of all samples
Samples<-Main.table[,17:ncol(Main.table)]

#Compound search range####
Mass.error<-0.05
RT.error <- 12

# Loading specs####
#winDialog(type = "ok", "Press OK and select specs file")                        #Message to select 
#Std.specs<-read.csv(file.choose(),                                              #Nice window/no captions
#                    sep = ";", header = FALSE,col.names = c("name","mz","RT"))
Std.specs<-read.csv(choose.files(caption = "Select specs file"),               #Ugly window/but captions
                    sep = ";", header = FALSE,col.names = c("name","mz","RT"))
#transforming RT from min to sec
Std.specs$RT <- Std.specs$RT*60

#Main loop####
temp <- data.frame()
temp2<-data.frame()
Std <- data.frame()
for (i in 1:nrow(Std.specs)) {
  
  name <- Std.specs[i,1]
  mz <- Std.specs[i,2] 
  RT = Std.specs[i,3]
  
  if(nrow(Samples[which(Main.table$mzmed>(mz-Mass.error) & Main.table$mzmed<(mz+Mass.error)&
                        Main.table$rtmed>(RT-RT.error)&Main.table$rtmed<(RT+RT.error)),])>0)
  {
    #Sampling####
    #Sampling all features fitting mz and RT criteria 
    temp <-  Samples[which(Main.table$mzmed>(mz-Mass.error) & Main.table$mzmed<(mz+Mass.error)&
                             Main.table$rtmed>(RT-RT.error)&Main.table$rtmed<(RT+RT.error)),]
    #subset for visualization
    temp_Int<-temp
    temp2<-data.frame(t(temp_Int))
    temp_mz<-Main.table$mzmed[which(Main.table$mzmed>(mz-Mass.error) & Main.table$mzmed<(mz+Mass.error)&
                                       Main.table$rtmed>(RT-RT.error)&Main.table$rtmed<(RT+RT.error))]
    temp_mz<-round(temp_mz,digits = 4)
    temp_rt<-Main.table$rtmed[which(Main.table$mzmed>(mz-Mass.error) & Main.table$mzmed<(mz+Mass.error)&
                                       Main.table$rtmed>(RT-RT.error)&Main.table$rtmed<(RT+RT.error))]/60
    temp_rt<-round(temp_rt,digits = 2)
    temp_name<-paste(temp_mz,"/",temp_rt,sep = "")
    colnames(temp2)<-temp_name
    #Adding average column
    temp$AV<-rowMeans(temp, na.rm=T)
    #Overwrite data frame with only one feature with highest average signal
    temp<-temp[which(temp$AV==max(temp$AV)),]
    #remove average column
    temp<-temp[,-which(names(temp)=="AV")] 
    #name row according to the CS compound
    rownames(temp)<-name
    Std <-rbind(Std, temp)  
    
    #Boxplot####
    plot_name<-paste(name,"/",c(RT/60),sep = "")  
    par(mfrow=c(1,1))
    
    ###Standard boxplot#####
    #png(paste(name,".png",sep = ""),width = 500,height = 400,res = 100)        
    #boxplot(temp2,main=plot_name,ylab="intensity",xlab="candidates")
    #dev.off()
    #boxplot(temp2,main=plot_name,ylab="intensity",xlab="candidates")
    
    
    ###ggplot####
    p<-ggplot(stack(temp2), aes(x=ind,y=values)) +
      geom_boxplot() + 
      geom_jitter(shape=16,position = position_jitter(0.2),size=2)+
      labs(title = plot_name,x="groups", y="intensity" )+
      theme_classic()+
      theme(plot.title=element_text(hjust=0.5,size=14, face = "bold"),
            axis.text = element_text(size = 10,color = "black"),
            axis.title = element_text(size = 10,face = "bold"))
    print(p)
    png(paste("S - ",name,".png",sep = ""),width = 500,height = 400,res = 100)  #Saving png
    print(p)
    dev.off()
    
   
  } else 
  {
    par(mar=c(0,0,0,0))
    png(paste("S - ",name,".png",sep = ""),width = 500,height = 400,res = 100)
    plot(c(0,1),c(0,1),ann=F,bty="n",type = "n",xaxt="n",yaxt="n")
    text(x=0.5,y=0.5,paste("For",name,"\n there was no matches", sep = " "),cex=1.5,col="black",font=2,adj=0.5)
    dev.off()
    plot(c(0,1),c(0,1),ann=F,bty="n",type = "n",xaxt="n",yaxt="n")
    text(x=0.5,y=0.5,paste("For",name,"\n there was no matches", sep = " "),cex=1.5,col="black",font=2,adj=0.5)
  }
  
}

write.table(Std,"Standards.csv", col.names = TRUE, row.names = TRUE, sep = ",")
winDialog(type = "ok", "R script execution is finished")                        #Message in the end
