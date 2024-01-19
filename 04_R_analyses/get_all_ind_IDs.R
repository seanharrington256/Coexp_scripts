#### Script to do get all the indvidual IDs used in the paper
# and make a table of populations


# had to change YPM13949 in the coords to YP13949 to match the genetic data files
#    and YPM13969 to YP13969 for same reason

### load up relevant packages
library(vcfR)

## Set up an object to contain the path to the main directory with the data and then set that as the working directory
main_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad"
setwd(main_dir)




## Specify all of the assemblies, excluding milk snakes, already published
all_assemblies<-c(
  "Acontortrix_p123_v2_25miss",
  "Dpunctatus_p123_v3_25missEAST",
  "Lgetula_p123_v4_25miss",
  "Pguttatus_p123_v2_25miss",
  "Sdekayi_p123_v4_25miss",
  "erytro",
  "abacura_only",
  "Mflagellum_p123_v3_25missEast"
)





#### Some overall setup for mapping and plotting


######################################################################################################################
## Loop over the assemblies
######################################################################################################################

# empty list to put individual names into
ind_names <- list()

for(i in seq_along(all_assemblies)){ 
  ###########################################################
  ## Set up paths to input files
  ###########################################################
  species <- all_assemblies[[i]]
  
  # path to vcf file
  path_vcf<-paste0(main_dir,"/", species,".vcf")
  
  
  ## Read in genetic data 
  gendata_all <- read.vcfR(path_vcf) # read in all of the genetic data
  gendata <- vcfR2genlight(gendata_all) # make it a genlight for easy access to ind names
  ind_names[[i]] <- gendata@ind.names ## get the individual names
}

# write out the individuals
all_inds <- unlist(ind_names)
write.table(all_inds, file = "all_inds_coexp.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



##### After doing this, and then some stuff to grab the number of reads per individual
##   read that back in with other info to be combined into an appendix

all_inds <- read.table("all_inds_coexp.txt") # all the newly sequenced indiviudals
colnames(all_inds) <- "Sample_ID"
num_reads <- read.table("num_forw_reads_per_samp.txt")  # reads for individuals other than Farancia
colnames(num_reads) <- c("Sample_ID", "Reads")

# Merge these together:
inds_reads <- merge(all_inds, num_reads, all = TRUE)

# go to the directory with Stairwayplot population files and use these to fill in some extra info
setwd("~/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad/stairwayplot2/pop_files")

# Read in each file
popfiles <- list.files(pattern = ".txt")
popinfo <- list()

for(i in seq_along(popfiles)){
  popinfo[[i]] <- read.table(popfiles[[i]]) # get the pop info 
  popinfo[[i]]$Species <- strsplit(popfiles[[i]], split = "_")[[1]][[1]] # make a column of just the species name
}

popinfo_df <- do.call(rbind, popinfo)
colnames(popinfo_df) <- c("Sample_ID", "Population", "Species_or_complex")

# merge that with inds_reads
appendix_info <- merge(inds_reads, popinfo_df, all.x = TRUE)

# Do a little renaming:
appendix_info$Species_or_complex <- gsub("abacura", "Fabacura", appendix_info$Species_or_complex)
appendix_info$Species_or_complex <- gsub("erytro", "Ferytrogramma", appendix_info$Species_or_complex)

# write it out
write.table(appendix_info, file = "~/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad/Appendix_1.txt", row.names = FALSE, quote = FALSE)
