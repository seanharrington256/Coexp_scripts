#### Script to make population assignment files to use in making SFS files for stairwayplot2
###     also creates a filtered vcf file that contains the vcf file pared down to just the
###     individuals in each indiviual population

# Note that milk snakes, L. triangulum complex, are handled in a separate script because those data are published and handled differently
### !!!!! ALSO - note that if you run through this all from scratch, i.e., starting from your own sNMF runs,
###     then you'll need to make sure that the order of populations is correct - i.e., you'll have to go through one 
###     by one and make sure that designations like "east" actually correspond to the eastern population






### load up relevant packages - copied over from Pop_assignment.R - not all necessary
library(adegenet)
library(LEA)
library(vcfR)
library(plotrix)
library(mapdata)
library(rworldmap)

# set up directories
main_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad"
stairway2popsdir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad/stairwayplot2/pop_files"

vcf_filtered_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad/stairwayplot2/filtered_vcf"

setwd(main_dir)

# The assemblies to do this for - I had to do some manual stuff each time to make sure
#    that the population designations were correct (i.e., had to plot things out and correctly assign stuff list "east" to the eastern pop)
#    AS NOTED ABOVE, this will not necessarily work if you run this from your own sNMF output if population designations are shuffled - you will need to check that population names
#     go to the right populations

all_assemblies <- c(
  "Acontortrix_p123_v2_25miss",
  "Dpunctatus_p123_v3_25missEAST",
  "Lgetula_p123_v4_25miss",
  "Pguttatus_p123_v2_25miss",
  "Sdekayi_p123_v4_25miss",
  "erytro",
  "abacura_only",
  "Mflagellum_p123_v3_25missEast"
)


## read in coordinates of all samples
coords<-read.csv("all_coords_requested.csv", header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't



for(species in all_assemblies){
  if(species %in% c("Sdekayi_p123_v4_25miss", "erytro", "Mflagellum_p123_v3_25missEast")){
    k <- 1
  }
  if(species %in% c("Acontortrix_p123_v2_25miss", "Pguttatus_p123_v2_25miss", "abacura_only")){
    k <- 2
  }
  if(species %in% c("Lgetula_p123_v4_25miss", "Dpunctatus_p123_v3_25missEAST")){
    k <- 3
  }
  
  
  
  
  ## Read in genetic data just to get names
  path_vcf<-paste0(main_dir,"/", species,".vcf")
  gendata_all<-read.vcfR(path_vcf) # read in all of the genetic data
  gendata<-vcfR2genlight(gendata_all) # make the genetic data a biallelic matrix of alleles in genlight format
  ind_names<-gendata@ind.names ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting
  
  setwd(main_dir)
  
  if(k>1){
  
    obj.at<-load.snmfProject(paste0(species,".u.snmfProject")) # load up sNMF results
    
    #### For a selected value of k, get the best snmf run - k is specified up in if() statements showing file paths
    ce <- cross.entropy(obj.at, K = k)
    best.run <- which.min(ce) # find the run with the lowest cross validation error
    qmatrix <- Q(obj.at, K = k, run = best.run)   # Get the snmf Q matrix from the best run at the best k
    cluster <- apply(qmatrix, 1, which.max) # use the qmatrix to find which population cluster each sample has most membership in
    
    
    
    #### Combine individual names with the population membership - fsc is referenced because I originally wrote a version of this script to prep SFS for FastSimCoal
    fsc_pops<-cbind(ind_names, cluster)
    

    ## sort the coords into the right order
    fsc_pops[which(!fsc_pops[,1] %in% coords[,"number"]), 1]
    # match up the coordinates to the order of the individuals from snmf
    match_coords<-match(fsc_pops[,1], coords[,"number"])
    snmf_coords<-coords[match_coords,]
    
    
    col_plot <- fsc_pops # make this into an object to specify colors
    
    if(k==2){
      col_plot[which(col_plot[,2]==1),2]<-"white"
      col_plot[which(col_plot[,2]==2),2]<-"black"
    }
    if(k==3){
      col_plot[which(col_plot[,2]==1),2]<-"white"
      col_plot[which(col_plot[,2]==2),2]<-"black"
      col_plot[which(col_plot[,2]==3),2]<-"red"
    }
  
    
    map("worldHires", "Mexico", xlim=c(-125,-65), ylim=c(23,53),col="gray90", fill=TRUE)
    map("state", xlim=c(-125,-65), ylim=c(23,53), add=TRUE,col="gray90", fill=TRUE)
    points(x=snmf_coords[,"lon"], y=snmf_coords[,"lat"], pch=21, bg=col_plot[,"cluster"])
  
    
    if(species=="Acontortrix_p123_v2_25miss"){
      ## here, white = 1 = gut - so set this up:
      fsc_pops[which(fsc_pops[,2]==1),2]<-"east"
      fsc_pops[which(fsc_pops[,2]==2),2]<-"west"
    }
    
    if(species=="Dpunctatus_p123_v3_25missEAST"){
      ## set up the populatoins based on how it was colored:
      fsc_pops[which(fsc_pops[,2]==1),2]<-"south"
      fsc_pops[which(fsc_pops[,2]==2),2]<-"central"
      fsc_pops[which(fsc_pops[,2]==3),2]<-"north"
    }
    
    
    if(species=="Lgetula_p123_v4_25miss"){
      ## set up the populatoins based on how it was colored:
      fsc_pops[which(fsc_pops[,2]==1),2]<-"east"
      fsc_pops[which(fsc_pops[,2]==2),2]<-"cent"
      fsc_pops[which(fsc_pops[,2]==3),2]<-"holb"
    }
    
    if(species=="Pguttatus_p123_v2_25miss"){
      ## here, white = 1 = gut - so set this up:
      fsc_pops[which(fsc_pops[,2]==1),2]<-"gut"
      fsc_pops[which(fsc_pops[,2]==2),2]<-"emor"
    }
    
    if(species=="abacura_only"){
      ## here, white = 1 = west - so set this up:
      fsc_pops[which(fsc_pops[,2]==1),2]<-"west"
      fsc_pops[which(fsc_pops[,2]==2),2]<-"east"
    }
    
  }else{
    fsc_pops<-cbind(ind_names, strsplit(species, split="_")[[1]][[1]])
  }
  
  setwd(stairway2popsdir)
  
  ## write the pops file - this has to be done in a convoluted way to prevent a newline character on the last line, which FSC reads as a final blank line and then will not run correctly--but will also not throw an informative error-I learned this the hard way
  ##   I think this is actually unnecessary for easySFS, but do it anyways just in case I wasnt to run FSC on these later
  if(k>1){
    stairway_pops<-list()
    for(i in 1:k){
      stairway_pops<-fsc_pops[which(fsc_pops[,2]==unique(fsc_pops[,2])[[i]]),]
      lines<-apply(stairway_pops, 1, function(x) paste(x[[1]], x[[2]], sep="\t"))
      lines[1:length(lines)-1]<-paste0(lines[1:length(lines)-1], "\n")
      cat(lines, file=paste0(species, "_pop", unique(fsc_pops[,2])[[i]], ".txt"), fill=FALSE, sep="")
    }
  }else{
    lines<-apply(fsc_pops, 1, function(x) paste(x[[1]], x[[2]], sep="\t"))
    lines[1:length(lines)-1]<-paste0(lines[1:length(lines)-1], "\n")
    cat(lines, file=paste0(species, "_pop", strsplit(species, split="_")[[1]][[1]], ".txt"), fill=FALSE, sep="")
  }
  
  setwd(vcf_filtered_dir)
  ## Write vcf files with only individuals in each pop
  if(k>1){
    stairway_pops<-list()
    for(i in 1:k){
      pop<-unique(fsc_pops[,2])[[i]] # name of pop
      inds<-fsc_pops[which(fsc_pops[,2]==pop),1] # indiividuals in that population
      vcf_filt<-gendata_all[,c("FORMAT", inds)] # filter down the vcf
      write.vcf(vcf_filt, file=paste0(species, "_pop", pop, ".vcf.gz")) # write the vcf out
    }
  }else{  # if there is a single population, just write out the vcf as is
    write.vcf(gendata_all, file=paste0(species, "_pop", strsplit(species, split="_")[[1]][[1]], ".vcf.gz")) # write the vcf out
  }
}

