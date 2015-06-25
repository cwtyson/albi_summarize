source('C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/Codes/albi_summarize/streamlined_functions.R')


WorDir<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/BBAL/individuals/"
outDir<-"C:/Temp"
filePattern="GPS"
Resamp<-180
timezone<-1
species<-"BBAL"
resmp<-1
summarize<-1


streamlined(WorDir,outDir,Resamp,timezone,species,resmp,summarize,filePattern)

