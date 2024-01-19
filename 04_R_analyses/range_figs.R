## Script to plot out ranges and sNMF results for each of the 9 species/complexes

# maps for most from IUCN
#    milks from inaturalist: https://www.inaturalist.org/taxa/322431-Lampropeltis-triangulum/range.kml 
#    storeria from inaturalist: https://www.inaturalist.org/taxa/28562-Storeria-dekayi/range.kml 




### load up relevant packages
library(LEA)
library(mapdata)
library(rworldmap)
library(ggplot2)
library(scatterpie)
library(maps)
library(mapdata)
library(sf)
library(vcfR)
library(ggpubr)



## Set up an object to contain the path to the main directory with the data and then set that as the working directory
main_dir<-"/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad/"
setwd(main_dir)

# directory with shape files
poly_dir <- "/Users/harrington/Active_Research/Ecotone_genomics/GBS_Data/range_polygons/"


## make a directory to put the output plots into
plots_dir <- paste0(main_dir, "range_plots")  # specify a full path to the directory
if(!dir.exists(plots_dir)){ # check if the directory  exists and then only create it if it does not
  dir.create(plots_dir)
}

# Specify all of the assemblies that we want to plot out
all_assemblies<-c(
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


#### Some overall setup for mapping and plotting

# make a list of colors:
colors_6<-c("V1" = "red", "V2" = "blue", "V3" = "white", "V4" = "purple", "V5" = "pink", "V6" = "yellow")

## Read in coordinates for plotting farther down
setwd(main_dir)
coords<-read.csv("all_coords_requested.csv", header=TRUE, row.names=NULL) # coordinates of everything I sequenced and many I didn't


################################################# map data
sf_use_s2(FALSE) # turn off sphere geometry in sf - doesn't play nicely with the shape files I'm using

# Set up boundaries for map plotting:
xmin <- -115
xmax <- -65
ymin <- 25
ymax <- 50

# Define Albers Equal Area Conic projection
target_crs <- st_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# Set up basemap of eastern NA
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE)) # Get US states data
mexico <- st_as_sf(map("worldHires", "Mexico", plot = FALSE, fill = TRUE)) # Get Mexico data
canada <- st_as_sf(map("worldHires", "Canada", plot = FALSE, fill = TRUE)) # Get Canada data
to_map <- rbind(states, mexico, canada) # Combine the datasets

# crop the base map polygon - based on this: https://datascience.blog.wzb.eu/2019/04/30/zooming-in-on-maps-with-sf-and-ggplot2/
to_map_val <- st_make_valid(to_map)
cropped_map <- st_crop(to_map_val, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)

# make a ggplot object with the basemap
basemap <- ggplot() +
  geom_sf(data = cropped_map, fill = "gray95", color = "black") +
  theme_minimal()


for(species in all_assemblies){    
    ### If only one shape file, read it in, or combine multiple shape files if necessaruy:
    shapes_path <- paste0(poly_dir, species)
    all_shapes <- list.files(path = shapes_path, full.names = TRUE) # get list of each shape file
    
    if(length(all_shapes) > 1){ # if there are 2 of more polygons:
      all_get_poly <- lapply(all_shapes, st_read) # read in each shape file
      range_polygon <- st_union(all_get_poly[[1]], all_get_poly[[2]]) # merge together the first 2 polygons
      if(length(all_shapes) > 2){ # if there are more than 2
        for(i in 3:length(all_shapes)){ # for the 3rd through final shape:
          range_polygon <- st_union(range_polygon, all_get_poly[[i]]) # merge third or higher with previous
        }
      }
    }else{
      range_polygon <- st_read(all_shapes)
    }
    
    # crop the polygon to bounds we want
    range_poly_val <- st_make_valid(range_polygon)
    cropped_poly <- st_crop(range_poly_val, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
    
    
    #  set up value of k and size for plotting pies later:
    if(species == "Lgetula_p123_v4_25miss"){
      k <- 3
      pie_size <- 2.2
      sci_name <- "L. getula complex"
    }
    if(species == "Dpunctatus_p123_v3_25missEAST"){
      k <- 3
      pie_size <- 3.5
      sci_name <- "D. punctatus"
    }
    if(species == "Acontortrix_p123_v2_25miss"){
      k <- 2
      pie_size <- 2.75
      sci_name <- "A contortrix"
    }
    if(species == "Pguttatus_p123_v2_25miss"){
      k <- 2
      pie_size <- 1.5
      sci_name <- "P. guttatus complex"
    }
    if(species == "Sdekayi_p123_v4_25miss"){
      k <- 1
      pie_size <- 3
      sci_name <- "S. dekayi"
    }
    if(species == "erytro"){
      k <- 1
      pie_size <- 3
      sci_name <- "F. erytrogramma"
    }
    if(species == "abacura_only"){
      k <- 2
      pie_size <- 2.25
      sci_name <- "F. abacura"
    }
    if(species == "Mflagellum_p123_v3_25missEast"){
      k <- 1
      pie_size <- 3
      sci_name <- "M. flagellum"
    }
    if(species == "Milks_filtered_snps_taxa"){
      k <- 3
      pie_size <- 1.2
      sci_name <- "L. triangulum complex"
    }
    
    ## Read in vcf file to get individual names
    path_vcf<-paste0(main_dir, species,".vcf")
    
    ## Read in genetic data to get individual IDs for each dataset
    gendata_all<-read.vcfR(path_vcf) # read in all of the genetic data
    gendata<-vcfR2genlight(gendata_all) # make the genetic data a biallelic matrix of alleles in genlight format
    if(species=="Milks_filtered_snps_taxa"){ # have to handle the milks slightly differenly because the names of the individuals in the genetic data have the species tacked onto the front
      ind_names<-sapply(gendata@ind.names, function(x) strsplit(x, "_")[[1]][[3]])
      gendata@ind.names <- ind_names
    }else{
      ind_names<-gendata@ind.names ## get the individual names in the order that they show up in the various files - this is important farther down for getting coordinates into the right order for plotting
    }
    
    ### Get the geographic coordinates sorted out
    ## make sure there aren't any individuals that don't have coordinates
    ind_names[which(!ind_names %in% coords[,"number"])]
    # match up the coordinates to the order of the individuals from snmf
    match_coords <- match(ind_names, coords[,"number"])
    snmf_coords <- coords[match_coords,]
    
    # Milk snakes need to be handled differently from the rest
    if(species == "Milks_filtered_snps_taxa"){
      
      # read in the population designations from Burbrink et al.
      milk_pops<-read.table(paste0(main_dir, "/Data_D6_DAPC_TESS_assignments_4_taxa.txt"), header=TRUE)
      milk_pops_red<-milk_pops[which(milk_pops$dapc2.assign %in% c(2, 3, 4)),] # reduce down to only individuals assigned to gentilis (3), triangulum (2), or elapsoides (4)
      for_pies <- milk_pops_red[,c("Sample", "lon", "lat", "triangulum", "gentilis_dapc", "elapsoides_dapc")] # get just the columns for membership in triangulum, gentilis, elapsoides
      
      # name it for consistency with the other datasets
      colnames(for_pies)[4:6] <- c("V1", "V2", "V3")
      
    }else{ # for all other species/complexes
      # read in sNMF results
      path_snmf <- paste0(main_dir, species, ".u.snmfProject")
      obj.at<-load.snmfProject(path_snmf)
      
      # get cross entropy scores from sNMF:
      ce <- cross.entropy(obj.at, K = k) 
      best.run <- which.min(ce) # find the run with the lowest cross validation error
      
      ## Get the snmf Q matrix from the best run at this k
      qmatrix <- Q(obj.at, K = k, run = best.run)
      admix<-as.data.frame(qmatrix)
      
      # get the coordinate and admix data into a single dataframe
      for_pies <- cbind(snmf_coords, admix)
      
    }
    
    # convert the coordinates into the target projection
    for_pies_sf <- st_as_sf(for_pies, coords = c("lon", "lat"), crs = 4326) # Create an sf object with WGS84 coordinates
    for_pies_sf_transformed <- st_transform(for_pies_sf, target_crs) # Transform to the target CRS
    for_pies_transformed_df <- st_coordinates(for_pies_sf_transformed) %>% # Convert transformed sf object back to dataframe
      as.data.frame() %>%
      setNames(c("lon", "lat"))
    for_pies_transformed_df <- cbind(for_pies[, setdiff(names(for_pies), c("lon", "lat"))], for_pies_transformed_df) # Combine with other columns from the original dataframe
    
    
    # Get the right number of colors
    num_cols <- length(grep("^V", colnames(for_pies_transformed_df))) # how many V columns containing proportion of cluster memberships?
    colors <- colors_6[1:num_cols]
    
    
    if(k > 1){ # plot pie charts if K > 1
    range_snmfplot <- basemap +
      geom_sf(data = cropped_poly, fill = "gray50", color = NA, alpha = 0.85) + 
      geom_scatterpie(data = for_pies_transformed_df, aes(x=lon, y=lat), cols = grep("^V", colnames(for_pies_transformed_df), value = TRUE), size = 0.01, pie_scale = pie_size) + # plot the pies - use grep to get the column names that start with V, these are the admix proportions
      scale_fill_manual(values = colors) +
      guides(fill="none") + # get rid of the legend for admixture
      coord_sf(crs = target_crs) + # set target projection
      ggtitle(sci_name) + # set title to the species
      theme(plot.title = element_text(face="bold.italic", hjust = 0.5), # format the title
            axis.title.x = element_blank(),  # Remove x-axis label
            axis.title.y = element_blank())  # Remove y-axis label
    }else{ # at K 1, just plot points, not pies
      range_snmfplot <- basemap +
        geom_sf(data = cropped_poly, fill = "gray50", color = NA, alpha = 0.85) + 
        geom_point(data = for_pies_transformed_df, aes(x=lon, y=lat), shape = 21, fill = "red", size = pie_size) + # plot points if K = 1
        scale_fill_manual(values = colors) +
        guides(fill="none") + # get rid of the legend for admixture
        coord_sf(crs = target_crs) + # set target projection
        ggtitle(sci_name) + # set title to the species
        theme(plot.title = element_text(face="bold.italic", hjust = 0.5), # format the title
              axis.title.x = element_blank(),  # Remove x-axis label
              axis.title.y = element_blank())  # Remove y-axis label
    }
    assign(paste0("figplot_", species), range_snmfplot) # name the object with the species/complex name
}

# list out all of them
all_to_plot <- ls(pattern = "figplot")
all_plots <- mget(all_to_plot)

# arrange them in columns & rows
grid_plot <- ggarrange(plotlist = all_plots, ncol = 3, nrow = 3)


# plot it out to pdf & add in lat & long labels
pdf(file = paste0(plots_dir, "/ranges_pops.pdf"), width = 10, height = 10)
annotate_figure(grid_plot,
                bottom = text_grob("Longitude", face = "bold"),
                left = text_grob("Latitude", rot = 90, face = "bold"))
dev.off()


