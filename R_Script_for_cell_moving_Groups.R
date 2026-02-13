rm(list = ls())
library(sojourner)
# Specify folder with data
folder1=system.file("extdata","Scr_shRNA",package="sojourner")
folder2=system.file("extdata","shNDR1",package="sojourner")
folder3=system.file("extdata","shNDR2",package="sojourner")
folder4=system.file("extdata","shNDR12",package="sojourner")

# Create track list
trackll1<-createTrackll(folder=folder1, interact=FALSE, input=3, ab.track=FALSE, cores=1, frameRecord=TRUE)
trackll2<-createTrackll(folder=folder2, interact=FALSE, input=3, ab.track=FALSE, cores=1, frameRecord=TRUE)
trackll3<-createTrackll(folder=folder3, interact=FALSE, input=3, ab.track=FALSE, cores=1, frameRecord=TRUE)
trackll4<-createTrackll(folder=folder4, interact=FALSE, input=3, ab.track=FALSE, cores=1, frameRecord=TRUE)

# Take a look at the list by 
str(trackll1,1)

# Filter/choose tracks 10 frames or longer for all analysis
trackll.fi1<-filterTrack(trackll=trackll1, filter=c(min=10,max=Inf))
trackll.fi2<-filterTrack(trackll=trackll2, filter=c(min=10,max=Inf))
trackll.fi3<-filterTrack(trackll=trackll3, filter=c(min=10,max=Inf))
trackll.fi4<-filterTrack(trackll=trackll4, filter=c(min=10,max=Inf))

# Merge tracks from different image files in the folder
trackll.fi.me1=c(mergeTracks(folder1, trackll.fi1))
trackll.fi.me2=c(mergeTracks(folder2, trackll.fi2))
trackll.fi.me3=c(mergeTracks(folder3, trackll.fi3))
trackll.fi.me4=c(mergeTracks(folder4, trackll.fi4))
str(trackll.fi.me1,1)

# calculate MSD for all tracks longer than 10 frames
#resolution is the  ratio of pixel to uM, our is 1.3043um/pixel
#Folder-1
data_individual <- msd(trackll.fi.me1,dt=45,resolution=1.3043,summarize=F,cores=4,plot=TRUE,output=TRUE)
data_summarized <- msd(trackll.fi.me1,dt=45,resolution=1.3043,summarize=TRUE,cores=4,plot=TRUE,output=TRUE)
head(data_individual[[1]])
individual <- data_individual[[1]]
summarized <- data_summarized[[1]]
write.csv(individual,file = "Scr_shRNA_MSD_individual.csv")
write.csv(summarized,file = "Scr_shRNA_MSD_summarized.csv")

#Folder-2
data_individual2 <- msd(trackll.fi.me2,dt=45,resolution=1.3043,summarize=F,cores=4,plot=TRUE,output=TRUE)
data_summarized2 <- msd(trackll.fi.me2,dt=45,resolution=1.3043,summarize=TRUE,cores=4,plot=TRUE,output=TRUE)
head(data_individual2[[1]])
individual2 <- data_individual2[[1]]
summarized2 <- data_summarized2[[1]]
write.csv(individual2,file = "shNDR1_MSD_individual.csv")
write.csv(summarized2,file = "shNDR1_MSD_summarized.csv")

#Folder-3
data_individual3 <- msd(trackll.fi.me3,dt=45,resolution=1.3043,summarize=F,cores=4,plot=TRUE,output=TRUE)
data_summarized3 <- msd(trackll.fi.me3,dt=45,resolution=1.3043,summarize=TRUE,cores=4,plot=TRUE,output=TRUE)
head(data_individual3[[1]])
individual3 <- data_individual3[[1]]
summarized3 <- data_summarized3[[1]]
write.csv(individual3,file = "shNDR2_MSD_individual.csv")
write.csv(summarized3,file = "shNDR2_MSD_summarized.csv")

#Folder-4
data_individual4 <- msd(trackll.fi.me4,dt=45,resolution=1.3043,summarize=F,cores=4,plot=TRUE,output=TRUE)
data_summarized4 <- msd(trackll.fi.me4,dt=45,resolution=1.3043,summarize=TRUE,cores=4,plot=TRUE,output=TRUE)
head(data_individual4[[1]])
individual4 <- data_individual4[[1]]
summarized4 <- data_summarized4[[1]]
write.csv(individual4,file = "shNDR12_MSD_individua_new.csv")
write.csv(summarized4,file = "shNDR12_MSD_summarized_new.csv")


#	rsquare filter on Dcoef results. Default to be 0.8. Set value to 0 if rsquare filter is not desired.
#method=c('static','percentage','rolling.window')

trackll_compare=compareFolder(folders=c(folder1,folder2,folder3,folder4), input=3)

Dcoef(trackll=trackll_compare,dt=45,filter=c(min=7,max=Inf),rsquare=0.8,
      resolution=1.3043, binwidth=NULL,
      method=c('percentage'),
      plot=TRUE,output=T,t.interval=20,profile=NULL)

Dcoef(trackll=trackll_compare,dt=45,filter=c(min=7,max=Inf),rsquare=0.8,
      resolution=1.3043, binwidth=NULL,
      method=c('static'),
      plot=TRUE,output=T,t.interval=20,profile=NULL)


Dcoef(MSD=data_individual,dt=45,filter=c(min=7,max=Inf),rsquare=0.8,
      resolution=1.3043, binwidth=NULL,
      method=c('percentage'),
      plot=TRUE,output=FALSE,t.interval=20,profile=NULL)


#Dcoef(trackll=trackll.fi.me,dt=45, filter=c(min=6,max=Inf), method="static", plot=TRUE, output=TRUE)


## calculating displacement CDF,calculate cumulative distribution function of all displacement for individual trajectories.
#CDF-Cumulative distribution function 
# Calculate all dislacement frequency and displacement CDF for all tracks longer than X frames
displacementCDF(trackll= trackll.fi.me,dt=45,resolution=1.3043,plot=T,output=FALSE,
                bivar=T)























