library(dplyr)
library(purrr)
library(ggplot2)
library(Polychrome)


##go to directory with all of the stiarway plot table showing median Ne, Year etc.
setwd("~/Active_Research/Ecotone_genomics/GBS_Data/Coexp_dryad/stairwayplot2/StairwayOuts_2023_09_12")

# get the summary files
sum_files <- list.files(pattern = "summary")

# read them all in
results <- lapply(sum_files, read.table, header=TRUE)
names(results) <- gsub(".final.summary", "", sum_files) # give population name as name of each element of the list

# Do some kludgy renaming:
names(results) <- gsub("Acontortrix_east", "A. contortrix - east", names(results))
names(results) <- gsub("Acontortrix_west", "A. contortrix - west", names(results))
names(results) <- gsub("Dpunctatus_central", "D. punctatus - central", names(results))
names(results) <- gsub("Dpunctatus_north", "D. punctatus - north", names(results))
names(results) <- gsub("Dpunctatus_south", "D. punctatus - south", names(results))
names(results) <- gsub("F_abacura__west", "F. abacura - west", names(results))
names(results) <- gsub("F_abacura_east", "F. abacura - east", names(results))
names(results) <- gsub("F_erytrogramma", "F. erytrogramma", names(results))
names(results) <- gsub("Lelapsoides", "L. elapsoides", names(results))
names(results) <- gsub("Lgetula_central", "L. getula nigra", names(results))
names(results) <- gsub("Lgetula_east", "L. getula getula", names(results))
names(results) <- gsub("Lgetula_holb", "L. getula holbrooki", names(results))
names(results) <- gsub("Ltriangulum", "L. triangulum", names(results))
names(results) <- gsub("Mflagellum", "M. flagellum", names(results))
names(results) <- gsub("Pemoryi", "P. emoryi", names(results))
names(results) <- gsub("Pguttatus", "P. guttatus", names(results))
names(results) <- gsub("Sdekayi", "S. dekayi", names(results))



# For each element of the list, add a column that is just the population name
named_res <-  imap(results, ~mutate(.x, Lineage = .y))

# collapse it into a single dataframe
result_df <- bind_rows(named_res)
result_df$Lineage <- factor(result_df$Lineage) # make the Lineage a factor

# Do some transformation
result_df <- mutate(result_df, kya = year/1000)
result_df <- mutate(result_df, logNe = log(Ne_median))

# create your own color palette based on `seedcolors`
colX = createPalette(36,  c("#ff0000", "#00ff00", "#0000ff"))
names(colX) <- NULL # remove names



# plot  Ne - x axis is thousand years ago and y axis is median Ne
ggplot(data = result_df, mapping = aes(x = kya, y = Ne_median, color = Lineage)) +
  geom_line() +
  scale_color_manual(values = colX) +
  labs(x = "Thosand years ago", y = "median Ne") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"))


# plot logged Ne - x axis is thousand years ago and y axis is logged median Ne
ggplot(data = result_df, mapping = aes(x = kya, y = logNe, color = Lineage)) +
  geom_line() +
  scale_color_manual(values = colX) +
  labs(x = "Thosand years ago", y = "log(median Ne)") +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"))

### TRUNCATED VERSION
# plot logged Ne - x axis is thousand years ago and y axis is logged median Ne
yarn_log_trunc <- ggplot(data = result_df, mapping = aes(x = kya, y = logNe, color = Lineage)) +
  geom_line() +
  scale_color_manual(values = colX) +
  labs(x = "Thousand years ago", y = "log(median Ne)") +
  ylim(9.5, 16) +
  theme_classic() +
  theme(legend.text = element_text(face = "italic"))


pdf(file = "yarn_log_trunc.pdf", width = 10, height = 7)
print(yarn_log_trunc)
dev.off()

### make columns for log of the upper and lower confidence limits
result_df <- mutate(result_df, logNe975 = log(Ne_97.5.))
result_df <- mutate(result_df, logNe25 = log(Ne_2.5.))



### Make a plot for each species with the error around it
# open the pdf
pdf(file = "Fig_S4_all_stairway.pdf", width = 8, height = 6)
for(i in names(results)){ # loop over the results object by name
  # i <- "A. contortrix - east" #TESTING
  # i <- "F. erytrogramma"  
  loop_df <- filter(result_df, Lineage == i) # filter to a single lineage designated by i
  
  stair_w_error <- ggplot(data = loop_df) +
    geom_line(mapping = aes(x = kya, y = logNe), color = "darkred") +
    geom_ribbon(aes(x = kya, ymin = logNe25, ymax = logNe975), fill = "hotpink", alpha = 0.5) +
    labs(x = "Thousand years ago", y = "log(Ne)") +
    theme_classic() +
    labs(title = "Fig. S4. staiwayplot2 output", subtitle = bquote(italic(.(i)))) +
    theme(plot.title = element_text(hjust = 0.5, size = 15),
          plot.subtitle = element_text(hjust = 0.5, size = 15))
  print(stair_w_error)
}
dev.off()

