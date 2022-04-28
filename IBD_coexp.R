#### Script to look at isolation by distance for GBS data

### load up relevant packages
library(adegenet)
library(MASS)
# library(vegan)
# library(ape)
library(fossil)
# library(adespatial)
# library(varhandle)
# library(logisticPCA)
library(raster)
# library(caret)
# library(LEA)
# library(psych)
library(vcfR)
library(dartR)
library(rSDM)


## Set up an object to contain the path to the main directory with the data and then set that as the working directory
main_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/pop_assignment"
worldclim_dir<-"/Users/harrington/worldclim_data"
bioclim_dir<-"/Users/harrington/worldclim_data/wc2.1_30s_bio"
envirem_dir<-"/Users/harrington/environ_layers_non_worldclim/Envirem_NAmerica"
setwd(main_dir)

### Set this up so that major code blocks only need to be written once, with subsitutions based on specific species
######################################################################################################################
######################################################################################################################


## Loop over all of the assemblies that we want to analyze
      

assemblies_mantel<-c(
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

# There is one loop set up right now, just to loop IBD plot and mantel tests over assemblies of interest
  # once the rest gets hammered out, I'll either loop it all together or make them seaprate loops

######################################################################################################################
######################################################################################################################


## make a directory to put the output plots into
ibd_out_dir<-paste0(main_dir, "/Coexp_ibd_out")  # specify a full path to the directory
if(!dir.exists(ibd_out_dir)){ # check if the directory  exists and then only create it if it does not
  dir.create(ibd_out_dir)
}



for(species in assemblies_mantel){   ### if we want to loop over all species, this line and line starting "all_assemblies<-c" should be uncommented, as well as final "}"
  ###########################################################
  ## Set up paths to input files
  ###########################################################
  path_ugeno<-paste0(main_dir,"/", species,".ugeno")
  path_ustr<-paste0(main_dir,"/", species,".ustr")
  # path_usnps<-paste0(main_dir,"/", species,".usnps")
  path_vcf<-paste0(main_dir,"/", species,".vcf")

  
  ## Read in coordinates
  setwd(main_dir)
  coords<-read.csv("all_coords_requested.csv", header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't
  ## read in the usnps file to get the number of individuals and snps for this assembly
  # nums_ind_snps<-as.numeric(unlist(strsplit(readLines(path_usnps)[[1]], split=" ")))
  
  
  ## Read in genetic data 
  gendata_all<-read.vcfR(path_vcf) # read in all of the genetic data
  gendata<-vcfR2genlight(gendata_all) # make the genetic data a biallelic matrix of alleles in genlight format
  ind_names<-gendata@ind.names ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting

  ### For  down below, get the geographic coordinates sorted out
  ## make sure there aren't any individuals that don't have coordinates
  ind_names[which(!ind_names %in% coords[,"number"])]
  # match up the coordinates to the order of the individuals from snmf
  match_coords<-match(ind_names, coords[,"number"])
  sorted_coords<-coords[match_coords,]
  
  ## Quick little Mantel test (not the best, but fast test)
  Dgen<-dist(gendata) # get the genetic distances
  Dgeo<-earth.dist(sorted_coords[,c("lon", "lat")]) # get the geographic distances
  ibd<-mantel.randtest(Dgen,Dgeo) # run the mantel test

  ## make a directory to put some plots into
  ## make a directory to put the output plots into
  ibd_out_dir<-paste0(main_dir, "/ibd_ibe_out")  # specify a full path to the directory
  if(!dir.exists(ibd_out_dir)){ # check if the directory  exists and then only create it if it does not
    dir.create(ibd_out_dir)
  }
  setwd(ibd_out_dir)
  
  ## pdf of plots
  pdf(file=paste0(species, "_Mantel_KD.pdf"), width=8, height=8)
  plot(ibd, main=paste0(species, "\nmantel p = ", ibd$pvalue)) # plot out the IBD significance
  ## make kernel density plot of genetic and geographic distances
  dens <- kde2d(as.numeric(Dgeo),as.numeric(Dgen), n=300)  # had to add as.numeric around the distances to make them numeric vectors--this wasn't necessary previously
  myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
  plot(Dgeo, Dgen, pch=20,cex=.5)
  image(dens, col=transp(myPal(300),.7), add=TRUE)
  abline(lm(as.numeric(Dgen)~as.numeric(Dgeo)))
  title(paste0(species, "\nIBD plot"))
  dev.off()
  
  setwd(main_dir)
}
  



#### Extract ecological data
###########################################################################

# Pull out altitude and climatic data for coordinates across all species. Do this only once and then write
#    a csv of the data to be used later
# bioclim data is downloaded from: https://www.worldclim.org/data/worldclim21.html  -- Bioclimatic variables at 30s -- specific DL link: https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_bio.zip
#    altitude also from https://www.worldclim.org/data/worldclim21.html -- specific download link: https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_elev.zip


## This is already done for all the coords, so don't do unless add in new coords - it takes a while

setwd(worldclim_dir)
alt_rast<-raster("wc2.1_30s_elev.tif")  # bring in altitude data
all_alt<-extract(alt_rast,coords[,c("lon", "lat")]) # get the altitiude data for each coordinate for all species
##  Downstream of this, found that there are NA values because of individual with offshore coords - map this to nearest non-NA cell
  na_index<-which(is.na(all_alt)) # get the row indices of NA values
  na_coords<-coords[na_index ,c("lon", "lat")] # extract out the NA coordinates from the coords object
  na_coords_form<-SpatialPoints(na_coords[,c("lon", "lat")], CRS(proj4string(alt_rast))) # format the NA coords into a SpatialPoints object
  new_points<-points2nearestcell(locs = na_coords_form, ras = alt_rast, move = TRUE,
                                 distance = NULL, showchanges = TRUE, showmap = TRUE, leaflet = FALSE)
  all_alt[na_index]<-extract(alt_rast, new_points)
rm(alt_rast) # remove alt_rast object now that data has been extracted from it
names(all_alt)<-coords[,"number"]
setwd(main_dir)
all_alt<-as.matrix(all_alt)
colnames(all_alt)<-"alt"
write.csv(all_alt, file="Coexp_coords_altitude.csv")
## now do the bioclim layers - current ecological data at 30s resolution
setwd(bioclim_dir)
bioclim_files <- list.files(pattern = "\\.tif$") # list out all files
cur_bioclim_stack<-stack(bioclim_files) # bring the data in as a raster stack
all_cur_bioclim<-extract(cur_bioclim_stack, coords[,c("lon", "lat")])
##  Downstream of this, found that there are NA values from samples offshore
colnames(all_cur_bioclim)[colSums(is.na(all_cur_bioclim)) > 0] ## looks like all of the envirem layers have some missing data
for(i in colnames(all_cur_bioclim)){  # loop over the columns to find the NA in each and then find the value from the nearest non-NA raster cell
  na_index<-which(is.na(all_cur_bioclim[,i])) # get the row indices of NA values
  na_coords<-coords[na_index ,c("lon", "lat")] # extract out the NA coordinates from the coords object
  na_coords_form<-SpatialPoints(na_coords[,c("lon", "lat")], CRS(proj4string(cur_bioclim_stack))) # format the NA coords into a SpatialPoints object
  new_points<-points2nearestcell(locs = na_coords_form, ras = cur_bioclim_stack, layer = i, move = TRUE,
                                 distance = NULL, showchanges = TRUE, showmap = TRUE, leaflet = FALSE)
  all_cur_bioclim[na_index,i]<-extract(cur_bioclim_stack[[i]], new_points)
}
rm(cur_bioclim_stack) # remove cur_bioclim_stack from environment now that data has been extracted from it
rownames(all_cur_bioclim)<-coords[,"number"]
setwd(main_dir)
write.csv(all_cur_bioclim, file="Coexp_coords_present_bioclim.csv")


## Do the same thing for Envirem layers, described here: http://envirem.github.io/
setwd(envirem_dir)
envirem_files <- list.files(pattern = "\\.tif$", recursive=TRUE) # list out all files
cur_envirem_stack<-stack(envirem_files) # bring the data in as a raster stack
all_cur_envirem<-extract(cur_envirem_stack, coords[,c("lon", "lat")])

##  Downstream of this, found that there are NA values in some of the Envirem vaiables
colnames(all_cur_envirem)[colSums(is.na(all_cur_envirem)) > 0] ## looks like all of the envirem layers have some missing data
for(i in colnames(all_cur_envirem)){  # loop over the columns to find the NA in each and then find the value from the nearest non-NA raster cell
  na_index<-which(is.na(all_cur_envirem[,i])) # get the row indices of NA values
  na_coords<-coords[na_index ,c("lon", "lat")] # extract out the NA coordinates from the coords object
  na_coords_form<-SpatialPoints(na_coords[,c("lon", "lat")], CRS(proj4string(cur_envirem_stack))) # format the NA coords into a SpatialPoints object
  new_points<-points2nearestcell(locs = na_coords_form, ras = cur_envirem_stack, layer = i, move = TRUE,
                                 distance = NULL, showchanges = TRUE, showmap = TRUE, leaflet = FALSE)
  all_cur_envirem[na_index,i]<-extract(cur_envirem_stack[[i]], new_points)
}

rm(cur_envirem_stack) # remove cur_envirem_stack from environment now that data has been extracted from it
rownames(all_cur_envirem)<-coords[,"number"]
setwd(main_dir)
write.csv(all_cur_envirem, file="Coexp_coords_present_envirem.csv")

  



# # pdf(file="all_samples_all_bioclim_corrs.pdf", width=50, height=50)
# # pairs.panels(all_cur_bioclim, scale=T) # For all bioclims across all points for all species, find correlations
# #   # do this once and then use just the same set of environmental variables for all - allows more direct comparison
# # dev.off()
# 
# 
# # to drop - drop any correlated above 0.8
# bio_to_drop<-c("wc2.1_30s_bio_10", "wc2.1_30s_bio_11", "wc2.1_30s_bio_6", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14", "wc2.1_30s_bio_16", "wc2.1_30s_bio_17", "wc2.1_30s_bio_19", "wc2.1_30s_bio_3", "wc2.1_30s_bio_7")
# all_spec_red_bioclim<-all_cur_bioclim[,!colnames(all_cur_bioclim) %in% bio_to_drop] # reduced set of bioclim variables, still for all species
# 
# # combine reduced bioclims with envirems
# all_specred_bioc_all_envirem<-cbind(all_spec_red_bioclim, all_cur_envirem)
# 
# 
# ## Now see to what degree envirem predictors are highly correlated with these or others  
# # pdf(file="all_samples_all_envirem_corrs.pdf", width=60, height=60)
# # pairs.panels(all_specred_bioc_all_envirem, scale=T) # For all bioclims across all points for all species, find correlations
# #   # do this once and then use just the same set of environmental variables for all - allows more direct comparison
# # dev.off()
# 
# # to drop - drop any correlated above 0.8
# envirem_to_drop<-c("current_30arcsec_growingDegDays0", "current_30arcsec_growingDegDays5", "current_30arcsec_maxTempColdest", "current_30arcsec_monthCountByTemp10", "current_30arcsec_PETColdestQuarter", "current_30arcsec_embergerQ", "current_30arcsec_thermicityIndex", "current_30arcsec_continentality", "current_30arcsec_PETDriestQuarter", "current_30arcsec_PETWarmestQuarter", "current_30arcsec_PETWettestQuarter", "current_30arcsec_annualPET", "current_30arcsec_climaticMoistureIndex", "current_30arcsec_minTempWarmest")  
# all_spec_bio_envir_red<-all_specred_bioc_all_envirem[,!colnames(all_specred_bioc_all_envirem) %in% envirem_to_drop]
# 
# 
# # double check that we don't have any high correlations here
# # pdf(file="all_samples_reduced_bio_evirem.pdf", width=60, height=60)
# # pairs.panels(all_spec_bio_envir_red, scale=T) # For all bioclims across all points for all species, find correlations
# #   # do this once and then use just the same set of environmental variables for all - allows more direct comparison
# # dev.off()
# ### Looks solid - nothing above 0.8









