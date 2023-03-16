#### Script to do population assignment on GBS data for manuscript
###     looking at coexpanding eastern snakes


# Prior to running this, I had to fix some names in the milk snake files
# ran these on command line in dir with the milk snake files:
# sed -i '.bak' 's/L_elapsoides_FTB2078/L_gentilis_FTB2078/g' *
# sed -i '.bak' 's/L_gentilis_FTB2109/L_elapsoides_FTB2109/g' *
# sed -i '.bak' 's/L_triangulum_FTB1538/L_triangulum_FTB1583/g' *
# had to also change YPM13949 in the coords to YP13949 to match the genetic data files
#    and YPM13969 to YP13969 for same reason

### load up relevant packages
library(adegenet)
library(LEA)
library(plotrix)
library(mapdata)
library(rworldmap)
library(ggplot2)
library(scatterpie)

## Set up an object to contain the path to the main directory with the data and then set that as the working directory
main_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad"
setwd(main_dir)


## make a directory to put the output plots into
sNMF_out_dir<-paste0(main_dir, "/Coexp_sNMF")  # specify a full path to the directory
if(!dir.exists(sNMF_out_dir)){ # check if the directory  exists and then only create it if it does not
  dir.create(sNMF_out_dir)
}


## Specify all of the assemblies that we want to run sNMF on 
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

## Get out the data for states and for Mexico and then combine them together
states <- map_data("state") # US states data
mex <- map_data("worldHires", "Mexico") # Mexico data
mex$group <- mex$group + length(states$group) # have to do this to get rid of weird lines that show up otherwise because of groups in Mexico already being group numbers in states
to_map <- rbind(states, mex) # combine these together
to_map <- dplyr::filter(to_map, lat > 23) # drop off souther coordinates that we don't need, only need northern Mexico

# make a list of colors:
colors_6<-c("V1" = "red", "V2" = "blue", "V3" = "white", "V4" = "purple", "V5" = "pink", "V6" = "yellow")

## Read in coordinates for plotting farther down
setwd(main_dir)
coords<-read.csv("all_coords_requested.csv", header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't



######################################################################################################################
## Loop to run sNMF over all assemblies
######################################################################################################################

for(species in all_assemblies){   ### if we want to loop over all assemblies, this line and line starting "all_assemblies<-c" should be uncommented, as well as final "}" if doing a single assembly, comment these lines instead
  ###########################################################
  ## Set up paths to input files
  ###########################################################
  setwd(main_dir)
  path_ugeno<-paste0(main_dir,"/", species,".ugeno")
  path_ustr<-paste0(main_dir,"/", species,".ustr")

  # snmf requires the geno file to have the extension .geno - the geno file of unlinked snps has ugeno
  #   as above, copy the geno and make one with the extension .u.geno
  path_geno<-gsub(".ugeno", ".u.geno", path_ugeno)  # Use a regular expression substitution to generate the new file name
  file.copy(path_ugeno, path_geno) # do the copying with the new name
  
  
  # Run sNMF using 1 to 10 ancestral populations and evaluate the fit of different k values to the data using cross entropy criterion
  # before running snmf, check if it's already been run
  if(dir.exists(gsub("geno", "snmf", basename(path_geno)))){
    obj.at<-load.snmfProject(gsub("geno", "snmfProject", basename(path_geno))) # if it has, just load up the results
  }else{ # otherwise, run sNMF
    obj.at <- snmf(input.file = path_geno,  # input file is the .geno format file. We set up the path to this above
                   K = 1:10, # we will test for k=1 through 10
                   ploidy = 2, 
                   entropy = T, # use the cross entropy criterion for assessing the best k value
                   repetitions = 10, # Run 10 independent replicate analyses
                   CPU = 1, 
                   project = "new", tolerance = 0.00001, iterations = 500)
  }
  
  setwd(sNMF_out_dir)
  
  # make pdf of cross-entropy plot
  plot(obj.at, col = "lightblue", cex = 1.2, pch = 19)
  
  # look at outstats
  outstats <- summary(obj.at)
  outstats # take a look
  
  # Plot k=2 through k=6 for all
  k_plot<-2:6

  
  ## This code block reads in the ustr file to get individual names in the order they show up
  ##    in data files, since geno files don't have ind names in them - this is not a great way
  ##    to do it, and is a holdover from when I used the ustr for other stuff that I've removed from this script,
  ##    but it works, so it stays - if I was building this ground-up again, I'd do this differently
  ##
  geno_txt<-readLines(path_ugeno)
  nums_snps<-length(geno_txt)
  num_ind<-length(strsplit(geno_txt[[1]], "")[[1]])
  ## quirk of read.structure function is that it requires the strucure file to have the file extension “.stru” - do some copying to make a new file with this extension
  path_stru<-gsub(".ustr", ".stru", path_ustr)  # Use a regular expression substitution to generate the new file name
  file.copy(path_ustr, path_stru) # make a copy of the file with the new name
  # Now we can read in this file
  ustr<-read.structure(path_stru, n.ind=num_ind, n.loc=nums_snps, onerowperind = FALSE, col.lab=1, col.pop=0, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)
  ind_names<-rownames(ustr@tab) ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting
  
  
  
  ### For  down below, get the geographic coordinates sorted out
  ## make sure there aren't any individuals that don't have coordinates
  ind_names[which(!ind_names %in% coords[,"number"])]
  # match up the coordinates to the order of the individuals from snmf
  match_coords<-match(ind_names, coords[,"number"])
  snmf_coords<-coords[match_coords,]
  
  
  pdf(file=paste0(species,"_SNMF_plots", ".pdf"), width=6, height=5)
  # put in the cross-entropy plot at the start
  plot(obj.at, col = "lightblue", cex = 1.2, pch = 19)
  #### use a loop to plot various different k values 
  for(i in k_plot){
    # confirm cross entropy values for K are consist. across runs
    ce <- cross.entropy(obj.at, K = i) 
    ce # pretty similar
    best.run <- which.min(ce) # find the run with the lowest cross validation error
    
    ## Get the snmf Q matrix from the best run at the best k
    qmatrix <- Q(obj.at, K = i, run = best.run)
    admix<-as.data.frame(qmatrix)
    
    # get the coordinate and admix data into a single dataframe
    for_pies <- cbind(snmf_coords, admix)
    
    # Get the right number of colors
    colors <- colors_6[1:ncol(admix)]
    
    
    ## plot it out
    snmf_plot <- ggplot(to_map, aes(long, lat, group = group)) + # map out the US & Mexico
      geom_polygon(data = to_map, fill = "grey90", color = "black", size = 0.2) + # make them polygons
      geom_scatterpie(data = for_pies, aes(x=lon, y=lat, group = number, r = 0.6), cols = grep("^V", colnames(for_pies), value = TRUE), size = 0.1) + # plot the pies - use grep to get the column names that start with V, these are the admix proportions
      scale_fill_manual(values = colors) +
      guides(fill="none") + # get rid of the legend for admixture
      theme_minimal() +
      labs(title=paste0(species,"_SNMF_K",i), x ="Longitude", y = "Latitude") +
      coord_map("moll") # Mollweide projection
    print(snmf_plot)
    
  }
  dev.off()
}



