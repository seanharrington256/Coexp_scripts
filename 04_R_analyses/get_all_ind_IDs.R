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


# read in coords
coords<-read.csv("all_coords_requested.csv", header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't



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

## add in geographic coordinates
# Prune down and rename columns
coords <- coords[,c("number", "lon", "lat")]
colnames(coords) <- c("Sample_ID", "Longitude", "Latitude")
## make sure there aren't any individuals that don't have coordinates
appendix_info$Sample_ID[which(!appendix_info$Sample_ID %in% coords$number)]
# merge that:
appendix_info <- merge(appendix_info, coords, all.x = TRUE)

# Do a little renaming:
appendix_info$Species_or_complex <- gsub("abacura", "Fabacura", appendix_info$Species_or_complex)
appendix_info$Species_or_complex <- gsub("erytro", "Ferytrogramma", appendix_info$Species_or_complex)


# Add in the num reads for Farancia that I got from Ed
setwd(main_dir)
fara_reads <- read.csv("Farancia_raw_reads_info.csv")

# check an discrepancies - all good
fara_reads$specimen[!fara_reads$specimen %in% appendix_info$Sample_ID]

# add in the number of reads:
for(i in 1:nrow(fara_reads)){
  ID <- fara_reads$specimen[i] # get the ID of the sample from the datafram with number of reads
  row_in_app <- which(appendix_info$Sample_ID == ID) # get the row in the appendix corresponding to that sample
  appendix_info$Reads[row_in_app] <- fara_reads$reads_raw[i] # replace the num reads colum in that row in the appendix with the number of reads
}

# Get museum numbers for Farancia - everything else is field numbers
farancia_musnums<- read.csv("Farancia_museum_nums.csv")

farancia_musnums$Name_in_Seq_files[!farancia_musnums$Name_in_Seq_files %in% appendix_info$Sample_ID]

# add in museum numbers:
appendix_info$Farancia_mus_num <- NA
for(i in 1:nrow(farancia_musnums)){
  ID <- farancia_musnums$Name_in_Seq_files[i] # get the ID of the sample from the datafram with museum numbers
  row_in_app <- which(appendix_info$Sample_ID == ID) # get the row in the appendix corresponding to that sample
  appendix_info$Farancia_mus_num[row_in_app] <- farancia_musnums$Museum.Field.Series.Numbers[i] # replace the column in that row in the appendix with the museum number
}



# write it out
write.table(appendix_info, file = "~/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad/Appendix_1.txt", row.names = FALSE, quote = FALSE)
