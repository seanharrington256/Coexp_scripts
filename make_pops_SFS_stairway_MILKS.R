#### Script to make population assignment files to use in making SFS files for stairwayplot2
###     also creates a filtered vcf file that contains the vcf file pared down to just the
###     individuals in each indiviual population

#############################################################################
###  This script is specifically for the milk snake data from Frank, 
###      because this was handled differently than the rest 
#############################################################################




### load up relevant packages - copied over from _Pop_assignment.R - not all necessary
library(adegenet)
library(LEA)
library(vcfR)
library(plotrix)
library(mapdata)
# library(tess3r)
library(rworldmap)
# library(PipeMaster)


main_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/pop_assignment"
stairway2popsdir<-"~/Active_Research/Ecotone_genomics/GBS_Data/stairwayplot2/pop_files"
vcf_filtered_dir<-"~/Active_Research/Ecotone_genomics/GBS_Data/stairwayplot2/filtered_vcf"

setwd(main_dir)


# Set up assembly

species<-"milks_denovo-92"

# read in the file of coords and pop membership
popass<-read.table("~/Active_Research/Ecotone_genomics/GBS_Data/stairwayplot2/Data_D6_DAPC_TESS_assignments_4_taxa.txt", header=TRUE)
popass_red<-popass[which(popass$dapc2.assign %in% c(2, 3, 4)),] # reduce down to only individuals assigned to gentilis (3), triangulum (2), or elapsoides (4)


## Plot out these individuals and make sure things look good
col_plot<-cbind(popass_red$Sample, popass_red$dapc2.assign)  # make an object to assing colors
col_plot[which(col_plot[,2]==2),2]<-"white"
col_plot[which(col_plot[,2]==4),2]<-"black"
col_plot[which(col_plot[,2]==3),2]<-"blue"


map("worldHires", "Mexico", xlim=c(-125,-65), ylim=c(23,53),col="gray90", fill=TRUE)
map("state", xlim=c(-125,-65), ylim=c(23,53), add=TRUE,col="gray90", fill=TRUE)
points(x=popass_red[,"lon"], y=popass_red[,"lat"], pch=21, bg=col_plot[,2])



## Read in genetic data
path_vcf<-paste0(main_dir,"/", species,".vcf")
gendata_all<-read.vcfR(path_vcf) # read in all of the genetic data
gendata<-vcfR2genlight(gendata_all) # make the genetic data a biallelic matrix of alleles in genlight format
ind_names<-gendata@ind.names ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting


## Make sure that all of the samples we want to pull out are in the vcf file
which(!popass_red[,"Sample"] %in% ind_names) # looks good - nothing missing


## Give the populations meaningful names
fsc_pops<-col_plot
fsc_pops[which(fsc_pops[,2]=="white"),2]<-"tri"
fsc_pops[which(fsc_pops[,2]=="black"),2]<-"elap"
fsc_pops[which(fsc_pops[,2]=="blue"),2]<-"gent"



setwd(stairway2popsdir)

k<-3

## write the pops file - this has to be done in a convoluted way to prevent a newline character on the last line, which FSC reads as a final blank line and then will not run correctly--but will also not throw an informative error-I learned this the hard way
##   I think this is actually unnecessary for easySFS, but do it anyways just in case I wasnt to run FSC on these later
stairway_pops<-list()
for(i in 1:k){
  stairway_pops<-fsc_pops[which(fsc_pops[,2]==unique(fsc_pops[,2])[[i]]),]
  lines<-apply(stairway_pops, 1, function(x) paste(x[[1]], x[[2]], sep="\t"))
  lines[1:length(lines)-1]<-paste0(lines[1:length(lines)-1], "\n")
  cat(lines, file=paste0(species, "_pop", unique(fsc_pops[,2])[[i]], ".txt"), fill=FALSE, sep="")
}

setwd(vcf_filtered_dir)
## Write vcf files with only individuals in each pop
stairway_pops<-list()
for(i in 1:k){
  pop<-unique(fsc_pops[,2])[[i]] # name of pop
  inds<-fsc_pops[which(fsc_pops[,2]==pop),1] # indiividuals in that population
  vcf_filt<-gendata_all[,c("FORMAT", inds)] # filter down the vcf
  write.vcf(vcf_filt, file=paste0(species, "_pop", pop, ".vcf.gz")) # write the vcf out
}

