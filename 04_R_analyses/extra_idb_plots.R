## A script to make IBD kernel density plots from the vcf files made for split up populations of
##    D. punctatus and A. contortrix


### load up relevant packages
library(adegenet)
library(MASS)
library(fossil)
library(raster)
library(vcfR)
library(dartR)
# library(rSDM)



main_dir <- "/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad/stairwayplot2/filtered_vcf"


all_assemblies <- c(
  "Dpunctatus_p123_v3_25missEAST_popsouth",
  "Dpunctatus_p123_v3_25missEAST_popnorth",
  "Dpunctatus_p123_v3_25missEAST_popcentral",
  "Acontortrix_p123_v2_25miss_popeast",
  "Acontortrix_p123_v2_25miss_popwest"
)

# make an object with the name of the species/complex for each assembly:
tax_names <- c(
  "D. punctatus South",
  "D. punctatus North",
  "D. punctatus Central",
  "A. contortrix East",
  "A. contortrix West"
)
ass_taxa <- data.frame(assembly = all_assemblies, taxon = tax_names)


coords<-read.csv("/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad/all_coords_requested.csv", header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't

## make a directory to put the output plots into
extra_ibd_out_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad/Coexp_ibd_out/extra_ibd_plots"
if(!dir.exists(extra_ibd_out_dir)){ # check if the directory  exists and then only create it if it does not
  dir.create(extra_ibd_out_dir)
}



# make a pdf with each page as the IBD plot for each species in loop just below
setwd(extra_ibd_out_dir)
pdf(file="Fig_S3_IBD_population_plots.pdf", width=8, height=8)



for(species in all_assemblies){
  # set the actual taxon name for pasting into figure titles
  taxon <- ass_taxa$taxon[ass_taxa$assembly == species]
  
  # read in file
  path_vcf<-paste0(main_dir,"/", species,".vcf")
  
  ## Read in genetic data 
  gendata_all<-read.vcfR(path_vcf) # read in all of the genetic data
  gendata<-vcfR2genlight(gendata_all) # make the genetic data a biallelic matrix of alleles in genlight format
  # if(species=="Milks_filtered_snps_taxa"){ # have to handle the milks slightly differenly because the names of the individuals in the genetic data have the species tacked onto the front
  #   ind_names<-sapply(gendata@ind.names, function(x) strsplit(x, "_")[[1]][[3]])
  #   gendata@ind.names <- ind_names
  # }else{
  ind_names<-gendata@ind.names ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting
  # }
  

  ### For  down below, get the geographic coordinates sorted out
  ## make sure there aren't any individuals that don't have coordinates
  ind_names[which(!ind_names %in% coords[,"number"])]
  # match up the coordinates to the order of the individuals from snmf
  match_coords<-match(ind_names, coords[,"number"])
  sorted_coords<-coords[match_coords,]
  
  
  ## Get genetic and geogrpahic distances, get kernel density, and make a color palette
  Dgen<-dist(gendata) # get the genetic distances
  Dgeo<-earth.dist(sorted_coords[,c("lon", "lat")]) # get the geographic distances
  dens <- kde2d(as.numeric(Dgeo),as.numeric(Dgen), n=300)  # had to add as.numeric around the distances to make them numeric vectors--this wasn't necessary previously
  myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
  
  
  
  ## pdf of plot
  plot(Dgeo, Dgen, pch=20,cex=.5)
  image(dens, col=transp(myPal(300),.7), add=TRUE)
  title(paste0("Fig. S2 IBD kernel density plot\n", taxon))
}
dev.off()




