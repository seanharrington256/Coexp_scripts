### Script to run generalised dissimilarity modelling (GDM) on GBS data
###    for manuscript looking at coexpanding eastern snakes

# load up relevant packages
library(adegenet)
library(gdm)
library(caret)

## Set up an object to contain the path to the main directory with the data and then set that as the working directory
main_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/pop_assignment"
setwd(main_dir)
coords<-read.csv("all_coords_requested.csv", header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't



## make a directory to put the output into
gdm_out_dir<-paste0(main_dir, "/Coexp_gdm_out")  # specify a full path to the directory
if(!dir.exists(gdm_out_dir)){ # check if the directory  exists and then only create it if it does not
  dir.create(gdm_out_dir)
}





### Set this up so that major code blocks only need to be written once, with subsitutions based on specific assemblies
######################################################################################################################
######################################################################################################################

assemblies_gdm<-c(
  "Acontortrix_p123_v2_25miss",
  "Dpunctatus_p123_v3_25missEAST",
  "Lgetula_p123_v4_25miss", # note that PTA & Stairway was made from Lgetula_p123_v2_25miss
  "Pguttatus_p123_v2_25miss",
  "Sdekayi_p123_v2_25miss",
  "erytro",
  "abacura_only",
  "Mflagellum_p123_v3_25missEast",
  "milks_denovo-92"
)


#####################################################################################################################
######################################################################################################################


for(species in assemblies_gdm){   ### if we want to loop over all species, this line and line starting "all_assemblies<-c" should be uncommented, as well as final "}"
  ###########################################################
  ## Set up paths to input files for each assembly
  ## these are not all required at present, but
  ## have left them in because they were alrady in the 
  ## _Pop_assignment.R script
  ###########################################################
  path_ugeno<-paste0(main_dir,"/", species,".ugeno")
  path_ustr<-paste0(main_dir,"/", species,".ustr")
  # path_usnps<-paste0(main_dir,"/", species,".usnps")
  path_vcf<-paste0(main_dir,"/", species,".vcf")


  ###########
  # ## read in the usnps file to get the number of individuals and snps for this assembly
  # nums_ind_snps<-as.numeric(unlist(strsplit(readLines(path_usnps)[[1]], split=" ")))
  
  # read in the geno file to get the number of individuals and snps for this assembly
  geno_txt<-readLines(path_ugeno)
  nums_snps<-length(geno_txt)
  num_ind<-length(strsplit(geno_txt[[1]], "")[[1]])
  


  ## quirk of read.structure function is that it requires the strucure file to have the file extension “.stru” - do some copying to make a new file with this extension
  path_stru<-gsub(".ustr", ".stru", path_ustr)  # Use a regular expression substitution to generate the new file name
  file.copy(path_ustr, path_stru) # make a copy of the file with the new name - if this returns FALSE, just means that this has already been done: should be fine unless we did it incorrectly previously
  # Now we can read in this file
  gendata<-read.structure(path_stru, n.ind=num_ind, n.loc=nums_snps, onerowperind = FALSE, col.lab=1, col.pop=0, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)
  ind_names<-rownames(gendata@tab) ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting
  
  
  ### For  down below, get the geographic coordinates sorted out
  ## make sure there aren't any individuals that don't have coordinates
  ind_names[which(!ind_names %in% coords[,"number"])]
  # match up the coordinates to the order of the individuals from snmf
  match_coords<-match(ind_names, coords[,"number"])
  sorted_coords<-coords[match_coords,]
  
  ## Calculate genetic distance among all individals
  Dgen<-dist(gendata$tab, diag=TRUE, upper=TRUE)
  gdm_dgen<-Dgen/max(Dgen) #normalize the genetic distances
  
  ## Read in environmental data (this was extracted from Bioclim variables in script "_IBD_IBE_2_plates.R")
  all_alt<-read.csv("Coexp_coords_altitude.csv", stringsAsFactors = FALSE, row.names=1)
  all_cur_bioclim_only<-read.csv("Coexp_coords_present_bioclim.csv", stringsAsFactors = FALSE, row.names=1)
  all_cur_envirem<-read.csv("Coexp_coords_present_envirem.csv", stringsAsFactors = FALSE, row.names=1)
  all_cur_bioclim<-cbind(all_cur_bioclim_only, all_cur_envirem) ## note that despite the name, this also includes Envirems
  
  
  
  # match up these data to the order of the individuals in genetic data (also prunes these data to just those for the assembly of interest)
  match_bioclim<-match(ind_names, rownames(all_cur_bioclim))
  bioclim<-all_cur_bioclim[match_bioclim,]
  match_alt<-match(ind_names, rownames(all_alt))
  alt<-all_alt[match_alt,]
  
  ### Remove any highly correlated environmental variables from bioclim
  # bioclimcor<-cor(bioclim) # get correlations of variables
  # drop_bio<-findCorrelation(bioclimcor, cutoff = 0.9,verbose = FALSE, names = FALSE, exact = TRUE) # identify variables to remove that have correlations with abs value > 0.9
  # bioclim_red<-bioclim[,-drop_bio]   # drop out the correlated columns to generate reduced set of bioclim
  
  # instead of the commented out automatic way above, I have manually done this in the _IBD_IBE.R script
  bio_to_drop<-c("wc2.1_30s_bio_10", "wc2.1_30s_bio_11", "wc2.1_30s_bio_6", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", "wc2.1_30s_bio_16", "wc2.1_30s_bio_17", "wc2.1_30s_bio_19", "wc2.1_30s_bio_3", "wc2.1_30s_bio_7")
  envirem_to_drop<-c("current_30arcsec_growingDegDays0", "current_30arcsec_growingDegDays5", "current_30arcsec_maxTempColdest", "current_30arcsec_monthCountByTemp10", "current_30arcsec_PETColdestQuarter", "current_30arcsec_embergerQ", "current_30arcsec_thermicityIndex", "current_30arcsec_continentality", "current_30arcsec_PETDriestQuarter", "current_30arcsec_PETWarmestQuarter", "current_30arcsec_PETWettestQuarter", "current_30arcsec_annualPET", "current_30arcsec_climaticMoistureIndex", "current_30arcsec_minTempWarmest")  
  vars_drop<-c(bio_to_drop, envirem_to_drop)
  
  bioclim_red<-bioclim[,!colnames(bioclim) %in% vars_drop] # reduced set of bioclim variables
  
  ## Slap some informative names onto the bioclims
  colnames(bioclim_red)<-gsub("wc2.1_30s_bio_1$", "AnnMeanT_B1", colnames(bioclim_red))
  colnames(bioclim_red)<-gsub("wc2.1_30s_bio_12", "AnnPrecip_B12", colnames(bioclim_red))
  colnames(bioclim_red)<-gsub("wc2.1_30s_bio_15", "PrecipSeas_B15", colnames(bioclim_red))
  colnames(bioclim_red)<-gsub("wc2.1_30s_bio_18", "PrecWarmQ_B18", colnames(bioclim_red))
  colnames(bioclim_red)<-gsub("wc2.1_30s_bio_2", "MDiurnRange_B2", colnames(bioclim_red))
  colnames(bioclim_red)<-gsub("wc2.1_30s_bio_4", "TSeas_B4", colnames(bioclim_red))
  colnames(bioclim_red)<-gsub("wc2.1_30s_bio_5", "MaxTWarmMon_B5", colnames(bioclim_red))
  colnames(bioclim_red)<-gsub("wc2.1_30s_bio_8", "MTWetQ_B8", colnames(bioclim_red))
  colnames(bioclim_red)<-gsub("wc2.1_30s_bio_9", "MTDryQ_B9", colnames(bioclim_red))
  colnames(bioclim_red)<-gsub("current_30arcsec_topoWet", "topoWet", colnames(bioclim_red))
  colnames(bioclim_red)<-gsub("current_30arcsec_tri", "TerrRough", colnames(bioclim_red))
  colnames(bioclim_red)<-gsub("current_30arcsec_aridityIndexThornthwaite", "aridityIndexThorn", colnames(bioclim_red))
  colnames(bioclim_red)<-gsub("current_30arcsec_PETseasonality", "PETseasonality", colnames(bioclim_red))
  
  
  # bind together the environmental data with individual names
  env_preds<-cbind(ind_names,  sorted_coords[c("lon", "lat")], bioclim_red, alt) # object of environmental predictors
  
  # bind together individual names with the genetic distance matrix  
  gdm_dgen_named<-cbind(ind_names, as.matrix(gdm_dgen))

  # Make the sitepair table necessary for gdm input (specific to gdm R package)
  gdm_sitepair<-formatsitepair(gdm_dgen_named, bioFormat=3, XColumn="lon", YColumn="lat", predData=env_preds, siteColumn="ind_names")

  # use this to run gdm
  gdm.1 <- gdm(gdm_sitepair, geo=T)
  
  # str(gdm.1)
  # summary(gdm.1)
  
  # plot
  setwd(gdm_out_dir)
  
  length(gdm.1$predictors) # get idea of number of panels
  pdf(file=paste0(species, "_gdm_splines.pdf"), width=13, height=20)
  plot(gdm.1, plot.layout=c(4,3))
  title(paste0(species, "\nGDM var explained", gdm.1$explained))
  dev.off()
  

  ## Get variable importance
  modTest_b <- gdm.varImp(gdm_sitepair, fullModelOnly = F, geo=T, nPerm=100, parallel=F)
  
  write.csv(modTest_b[[1]], paste0(species, "_GDM_p-value.csv"))
  write.csv(modTest_b[[2]], paste0(species, "_GDM_var_imp.csv"))
  write.csv(modTest_b[[3]], paste0(species, "_GDM_var_sig.csv"))

  pdf(file=paste0(species, "_GDM_var_imp.pdf"), width=6, height=6)
  par(mar=c(10,4,4,4))
  barplot(modTest_b[[2]][,1],main = paste0(species, "\nGDM full model var importance"), las=2)
  dev.off()

  setwd(main_dir)
  
}


