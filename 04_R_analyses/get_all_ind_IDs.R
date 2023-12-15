#### Script to do get all the indvidual IDs used in the paper


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


all_inds <- unlist(ind_names)
write.table(all_inds, file = "all_inds_coexp.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
