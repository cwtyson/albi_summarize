source('C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/Codes/albi_summarize/streamlined_functions.R')


#WorDir<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/WAAL/individuals/"
#outDir<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/BBAL/analysis/"

WorDir<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/BBAL/individuals/"
outDir<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/BBAL/analysis/"

filePattern="BBAL"
Res<-120
timezone<-1
species<-"BBAL"
resmp<-1
summarize<-1


streamlined(WorDir,outDir,Res,timezone,species,resmp,summarize,filePattern)



source('C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/Codes/albi_summarize/streamlined_functions.R')


WorDir<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/WAAL/individuals/"
outDir<-"C:/Users/Grant/Dropbox/GrantHumphriesBackup/Projects/Albatross/WAAL/analysis/"
filePattern="WAAL"
Resamp<-10
timezone<-2
species<-"WAAL"
resmp<-1
summarize<-1


streamlined(WorDir,outDir,Resamp,timezone,species,resmp,summarize,filePattern)
