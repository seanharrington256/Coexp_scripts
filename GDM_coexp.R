### Script to run generalised dissimilarity modelling (GDM) on GBS data
###    for manuscript looking at coexpanding eastern snakes

# load up relevant packages
library(adegenet)
library(gdm)
library(caret)
library(vcfR)
library(reshape2)
library(ggplot2)

## Set up an object to contain the path to the main directory with the data and then set that as the working directory
main_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad"
setwd(main_dir)
coords<-read.csv("all_coords_requested.csv", header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't


## make a directory to put the output into
gdm_out_dir<-paste0(main_dir, "/Coexp_gdm_out")  # specify a full path to the directory
if(!dir.exists(gdm_out_dir)){ # check if the directory  exists and then only create it if it does not
  dir.create(gdm_out_dir)
}


## Read in environmental data (this was extracted from Bioclim variables in script "IBD_coexp.R")
all_alt<-read.csv("Coexp_coords_altitude.csv", stringsAsFactors = FALSE, row.names=1)
all_cur_bioclim_only<-read.csv("Coexp_coords_present_bioclim.csv", stringsAsFactors = FALSE, row.names=1)
all_cur_envirem<-read.csv("Coexp_coords_present_envirem.csv", stringsAsFactors = FALSE, row.names=1)
all_cur_bioclim<-cbind(all_cur_bioclim_only, all_cur_envirem) ## note that despite the name, this also includes Envirems







### Set this up so that major code blocks only need to be written once, with subsitutions based on specific assemblies
######################################################################################################################
######################################################################################################################

assemblies_gdm<-c(
  "Acontortrix_p123_v2_25miss",
  "Dpunctatus_p123_v3_25missEAST",
  "Lgetula_p123_v4_25miss",
  "Pguttatus_p123_v2_25miss",
  "Sdekayi_p123_v4_25miss",
  "erytro",
  "abacura_only",
  "Mflagellum_p123_v3_25missEast",
  "Milks_filtered_snps_taxa"
)


#####################################################################################################################
######################################################################################################################


for(species in assemblies_gdm){   ### if we want to loop over all species, this line and line starting "all_assemblies<-c" should be uncommented, as well as final "}"
  ###########################################################
  ## Set up paths to input file for each assembly
  ###########################################################
  path_vcf<-paste0(main_dir,"/", species,".vcf")

  
  ## Read in genetic data 
  gendata_all<-read.vcfR(path_vcf) # read in all of the genetic data
  gendata<-vcfR2genlight(gendata_all) # make the genetic data a biallelic matrix of alleles in genlight format
  if(species=="Milks_filtered_snps_taxa"){ # have to handle the milks slightly differenly because the names of the individuals in the genetic data have the species tacked onto the front
    ind_names<-sapply(gendata@ind.names, function(x) strsplit(x, "_")[[1]][[3]])
    gendata@ind.names <- ind_names
  }else{
    ind_names<-gendata@ind.names ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting
  }
  
  
  ### For  down below, get the geographic coordinates sorted out
  ## make sure there aren't any individuals that don't have coordinates
  ind_names[which(!ind_names %in% coords[,"number"])]
  # match up the coordinates to the order of the individuals from snmf
  match_coords<-match(ind_names, coords[,"number"])
  sorted_coords<-coords[match_coords,]
  
  ## Calculate genetic distance among all individals
  Dgen<-dist(gendata, diag=TRUE, upper=TRUE)
  gdm_dgen<-Dgen/max(Dgen) #normalize the genetic distances
  
  
  # match up environmental data to the order of the individuals in genetic data (also prunes these data to just those for the assembly of interest)
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
  gdm.1 <- gdm(gdm_sitepair, geo=TRUE)
  
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
  modTest_b <- gdm.varImp(gdm_sitepair, predSelect = FALSE, geo=T, nPerm=100, parallel=F)
  
  write.csv(modTest_b[[1]], paste0(species, "_GDM_p-value.csv"))
  write.csv(modTest_b[[2]], paste0(species, "_GDM_var_imp.csv"))
  write.csv(modTest_b[[3]], paste0(species, "_GDM_var_sig.csv"))

  pdf(file=paste0(species, "_GDM_var_imp.pdf"), width=6, height=6)
  par(mar=c(10,4,4,4))
  barplot(modTest_b[[2]][,1], main = paste0(species, "\nGDM full model var importance"), las=2, names.arg = rownames(modTest_b[[2]]))
  dev.off()

  setwd(main_dir)
  
}


####################################################################################################
#### Some summaries across species:
####################################################################################################

# Read in the csv files that tell if each model is significant
model_sig_files <- list.files(gdm_out_dir, pattern = "GDM_p-value.csv", full.names = TRUE)
model_sig <- lapply(model_sig_files, read.csv)
tax_names <- sapply(model_sig_files, function(x) strsplit(basename(x), "_GDM_")[[1]][1])
names(model_sig) <- tax_names

# Which are significant for the whole model?? - model p-value is row 3 in all, value is column 2
model_sig_values_only <- sapply(model_sig, function(x) x[3,2])
sig_models <- model_sig_values_only[model_sig_values_only < 0.05]
spec_sig_model <- names(sig_models) # species that are significant in the full model


# For each of these, find out which predictors are significant
var_sig_files <- list.files(gdm_out_dir, pattern = "GDM_var_sig.csv", full.names = TRUE)
var_sig <- lapply(var_sig_files, read.csv)
tax_names <- sapply(var_sig_files, function(x) strsplit(basename(x), "_GDM_")[[1]][1])
names(var_sig) <- tax_names


# if the full model isn't significant, I don't think there's any case where an individual variable should be called significant, but let's make sure
#     the opposite seems the be the case, it's possible for the full model to be significant even if no ineividual variable is
var_sig <- var_sig[names(var_sig) %in% spec_sig_model]


# Get out the names of predictors that are significant at < 0.05 for each species
sig_vars <- lapply(var_sig, function(x) x[x$All.predictors < 0.05, 1])


# Go in and pull out the variable importance for significant models, let's plot this as heatmaps
#     do it 2 ways - only include variables if they're significant for an assembly, or include all importances
#     for any species if the full model is significant

# get variable importance
var_imp_files <- list.files(gdm_out_dir, pattern = "GDM_var_imp.csv", full.names = TRUE)
var_imp <- lapply(var_imp_files, read.csv)
tax_names <- sapply(var_imp_files, function(x) strsplit(basename(x), "_GDM_")[[1]][1])
names(var_imp) <- tax_names

# If the full model is not significant, replace importances with NA
non_sig_models <- names(var_imp)[!names(var_imp) %in% spec_sig_model]
if(length(non_sig_models) > 0){
  for(i in 1:length(non_sig_models)){
    var_imp[[non_sig_models[i]]][,2] <- NA
  }
}



# fix up some names
names(var_imp) <- gsub("abacura_only", "F. abacura", names(var_imp))
names(var_imp) <- gsub("Acontortrix_p123_v2_25miss", "A. contortrix", names(var_imp))
names(var_imp) <- gsub("Dpunctatus_p123_v3_25missEAST", "D. punctatus", names(var_imp))
names(var_imp) <- gsub("erytro", "F. erytrogramma", names(var_imp))
names(var_imp) <- gsub("Lgetula_p123_v4_25miss", "L. getula", names(var_imp))
names(var_imp) <- gsub("Mflagellum_p123_v3_25missEast", "M. flagellum", names(var_imp))
names(var_imp) <- gsub("Milks_filtered_snps_taxa", "L. triangulum", names(var_imp))
names(var_imp) <- gsub("Pguttatus_p123_v2_25miss", "P. guttatus", names(var_imp))
names(var_imp) <- gsub("Sdekayi_p123_v2_25miss", "S. dekayi", names(var_imp))




# and for sig_vars
names(sig_vars) <- gsub("abacura_only", "F. abacura", names(sig_vars))
names(sig_vars) <- gsub("Acontortrix_p123_v2_25miss", "A. contortrix", names(sig_vars))
names(sig_vars) <- gsub("Dpunctatus_p123_v3_25missEAST", "D. punctatus", names(sig_vars))
names(sig_vars) <- gsub("erytro", "F. erytrogramma", names(sig_vars))
names(sig_vars) <- gsub("Lgetula_p123_v4_25miss", "L. getula", names(sig_vars))
names(sig_vars) <- gsub("Mflagellum_p123_v3_25missEast", "M. flagellum", names(sig_vars))
names(sig_vars) <- gsub("Milks_filtered_snps_taxa", "L. triangulum", names(sig_vars))
names(sig_vars) <- gsub("Pguttatus_p123_v2_25miss", "P. guttatus", names(sig_vars))
names(sig_vars) <- gsub("Sdekayi_p123_v2_25miss", "S. dekayi", names(sig_vars))





# If a variable is not significant for a given taxon, also make that importance NA
var_imp_sig_only <- var_imp
for(i in 1:length(sig_vars)){
  non_sig_vars <- !var_imp_sig_only[[names(sig_vars)[i]]][,1] %in% sig_vars[[i]]  # find the predictors that are not in the list of sig variables for each taxon
  var_imp_sig_only[[names(sig_vars)[i]]][non_sig_vars,2] <- NA # replace the importances for those with NA
}


### Make heatmaps of each of these - either only sig variables, or all sig variables for any model that is significant overall

# Formatting into a single dataframe each
# melt the list into a dataframe and ditch a column we don't need
var_imp_df <- melt(var_imp, id = "X")
var_imp_df <- var_imp_df[, !colnames(var_imp_df) %in% "variable"]
colnames(var_imp_df) <- c("variable", "importance", "species")

# I should've made this swap up higher somwehere, but I didn't, so here we are
var_imp_df$variable <- gsub("PrecipSeas_B15", "Prec. seas B15", var_imp_df$variable)
var_imp_df$variable <- gsub("AnnPrecip_B12", "Ann. prec. B12", var_imp_df$variable)
var_imp_df$variable <- gsub("PrecWarmQ_B18", "Prec. warmest Q B18", var_imp_df$variable)
var_imp_df$variable <- gsub("MTWetQ_B8", "M. T. wettest Q B8", var_imp_df$variable)
var_imp_df$variable <- gsub("MTDryQ_B9", "M. T. driest Q B9", var_imp_df$variable)
var_imp_df$variable <- gsub("topoWet", "Topo wetness", var_imp_df$variable)
var_imp_df$variable <- gsub("TerrRough", "Terr roughness", var_imp_df$variable)
var_imp_df$variable <- gsub("aridityIndexThorn", "Thorn's aridity", var_imp_df$variable)
var_imp_df$variable <- gsub("MDiurnRange_B2", "M diurn. range B2", var_imp_df$variable)
var_imp_df$variable <- gsub("TSeas_B4", "T. seas B4", var_imp_df$variable)
var_imp_df$variable <- gsub("MaxTWarmMon_B5", "Max T warmest mon. B5", var_imp_df$variable)
var_imp_df$variable <- gsub("AnnMeanT_B1", "Ann. mean T.
                            B1", var_imp_df$variable)
var_imp_df$variable <- gsub("PETseasonality", "PET seasonality", var_imp_df$variable)
var_imp_df$variable <- gsub("alt", "Altitude", var_imp_df$variable)



var_imp_sig_only_df <- melt(var_imp_sig_only, id = "X")
var_imp_sig_only_df <- var_imp_sig_only_df[, !colnames(var_imp_sig_only_df) %in% "variable"]
colnames(var_imp_sig_only_df) <- c("variable", "importance", "species")

# I should've made this swap up higher somwehere, but I didn't, so here we are
var_imp_sig_only_df$variable <- gsub("PrecipSeas_B15", "Prec. seas B15", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("AnnPrecip_B12", "Ann. prec. B12", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("PrecWarmQ_B18", "Prec. warmest Q B18", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("MTWetQ_B8", "M. T. wettest Q B8", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("MTDryQ_B9", "M. T. driest Q B9", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("topoWet", "Topo wetness", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("TerrRough", "Terr roughness", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("aridityIndexThorn", "Thorn's aridity", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("MDiurnRange_B2", "M diurn. range B2", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("TSeas_B4", "T. seas B4", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("MaxTWarmMon_B5", "Max T warmest mon. B5", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("AnnMeanT_B1", "Ann. mean T. B1", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("PETseasonality", "PET seasonality", var_imp_sig_only_df$variable)
var_imp_sig_only_df$variable <- gsub("alt", "Altitude", var_imp_sig_only_df$variable)




# make heatmap with ggplot
map_all_vars <- ggplot(var_imp_df, aes(species, variable, fill= importance, na.rm = TRUE)) + 
  geom_tile(na.rm = TRUE) +
  scale_fill_gradient2(mid = "white",  high = "red", na.value = NA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), axis.title.y = element_blank())

# make heatmap with ggplot for significant only
maps_sigVars <- ggplot(var_imp_sig_only_df, aes(species, variable, fill= importance, na.rm = TRUE)) + 
  geom_tile(na.rm = TRUE) +
  scale_fill_gradient2(mid = "white",  high = "red", na.value = NA) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.title.x = element_blank(), axis.title.y = element_blank())

map_all_vars
maps_sigVars


setwd(gdm_out_dir)

pdf(file="heatmap_all_gdm_vars.pdf", width = 6, height = 4)
print(map_all_vars)
dev.off()

pdf(file="heatmap_sig_gdm_vars.pdf", width = 6, height = 4)
print(maps_sigVars)
dev.off()



