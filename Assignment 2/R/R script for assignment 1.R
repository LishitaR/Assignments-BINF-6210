#Assignment 1: Comparing Ixodes Tick BIN Diversity Between North America and Europe Using DNA Barcodes

##***** Name: Lishita Rowjee
##***** Class: BINF 6210
##***** Genus: Ixodes
##***** October 10, 2025

#Start-up and packages used
# Packages: tidyverse, ggplot2, vegan and dplyr

rm(list=ls())

# ==== EDIT START (Fangyi): Silent package setup using requireNamespace() ====
# EDIT: Replaced manual library() calls with automated, quiet package setup.
#       This version checks if each package is installed using requireNamespace(),
#       installs it if missing, an then loads it silently.
#       Ensures smooth execution on a clean system without startup messages.

# Old:
# library (tidyverse)
# library (ggplot2)
# library (vegan)
# library (dplyr)

packages_needed <- c("tidyverse", "ggplot2", "vegan", "dplyr", "tidyr")

for (pkg in packages_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, quite = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ==== EDIT END ====

# ==== EDIT (Fangyi): Explicit Namespaces Applied ====
# To improve clarity and reproducibility, and to avoid 'function not found' 
# error, all non-based functions in this script have been updated to use 
# explicit namespaces (pkg::function) - e.g.,
# dplyr::mutate(), readr::read_tsv(), ggplot2::ggplot(), vegan::vegdist()
# ===================================================

# Loading Ixodes data on this script

ixodes_data <- readr::read_tsv ("../data/Ixodes_BOLD_data.tsv")
View(ixodes_data)


#Finding column name where different countries in data set is located and looking at list of different countries

colnames(ixodes_data)
unique(ixodes_data$`country/ocean`)

# Defining countries per continent
north_america <- c("Canada", "United States", "Mexico", "Panama", "Costa Rica")
europe <- c("Germany", "Norway", "Sweden", "United Kingdom", "Poland", 
            "Russia", "Slovakia", "Luxembourg", "Bulgaria")

# Assigning countries to continent names using mutate function and removing NA values

country_to_continent <- ixodes_data %>%
  dplyr::mutate(region = case_when(
    `country/ocean` %in% north_america ~ "North America",
    `country/ocean` %in% europe ~ "Europe",
    TRUE ~ NA_character_
  )) %>%
  dplyr::filter(!is.na(region))

print(country_to_continent)

# Quality check – Filter original data by continent and confirm that correct countries were filtered out

ixodes_NorthAmerica <- ixodes_data %>%
  dplyr::filter(`country/ocean` %in% north_america)
unique(ixodes_NorthAmerica$`country/ocean`)

ixodes_europe <- ixodes_data %>%
  dplyr::filter(`country/ocean` %in% europe)
unique(ixodes_europe$`country/ocean`)


#Count unique BINs per continent (excluding NAs)

bins_NorthAmerica <- unique(na.omit(ixodes_NorthAmerica$bin_uri))
print(bins_NorthAmerica)
length(bins_NorthAmerica)

bins_europe <- unique(na.omit(ixodes_europe$bin_uri))
print(bins_europe)
length(bins_europe)

#Bin counts by continent after putting the bin_counts by region in one dataframe instead of 2 separate ones as above in order to make barplot- Figure 1 (remove NAs)

country_to_continent_clean <- country_to_continent %>%
  dplyr::filter(!is.na(bin_uri))

distinct_bins <- country_to_continent_clean %>%
  dplyr::distinct(region, bin_uri)

bin_counts_by_region <- distinct_bins %>%
  dplyr::count(region, name = "unique_bins")
print(bin_counts_by_region)

#same values (17 for Europe and 30 for North America obtained as above)

ggplot2::ggplot(bin_counts_by_region, ggplot2::aes(x = region, y = unique_bins, fill = region)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(
    title = "Unique BIN Counts by Continent",
    x = "Continent",
    y = "Number of Unique BINs") 

# ==== EDIT START (Fangyi): Correct Jaccard Calculation ====
# EDIT: Updated Jaccard computation to use a true presence-absence matrix (0/1)
#       and added binary = TRUE for correct similarity between continents.
#       This ensures that Jaccard reflects community composition rather than abundance.

#Building a presence-absence matrix (common format in ecology) to show BINS of Ixodes present or absent in the 2 continents

bin_continent_matrix <- table(country_to_continent$bin_uri, country_to_continent$region)

print(bin_continent_matrix)


#Converting presence-absence matrix to binary format where 1 indicates presence and 0 indicates absence

bin_continent_matrix_df <- as.data.frame.matrix(bin_continent_matrix)
print(bin_continent_matrix_df)

# EDIT (Fangyi): Renaming binary matrix for clarity
bin_continent_matrix_df_binary <- ifelse(bin_continent_matrix_df > 0, 1, 0)
print(bin_continent_matrix_df_binary)

#Calculating Jaccard distance of the 2 continents (statistic used in determining similarity and diversity of sample sets). Using the vegdist option from vegan package

distance_matrix <- vegan::vegdist(t(bin_continent_matrix_df), method = "jaccard")
print(distance_matrix)

# Value obtained is 0.9976303 which is close to 1 meaning that the BIN communities in Europe and North America are completely different

#Checking if same value would be obtained if transposing not the bin_continent_matrix_df (dataframe) but bin_continent_matrix (table format)-both should contain same values.

distance_matrix1 <- vegan::vegdist(t(bin_continent_matrix), method = "jaccard")
print(distance_matrix1)
rm(distance_matrix1)

#same jaccard distance obtained

# EDIT (Fangyi): 'bin_continent_matrix_df' and 'bin_continent_matrix' produce 
# the same Jaccard distance (0.9976303) because both represent quantitative data 
# (counts) and are effectively transposes of each other.
# However, to properly compare community composition, we need to compute Jaccard
# on a binary (presence-absence) matrix instead, where 1 indicates BIN presence
# and 0 indicates absence.

distance_matrix_binary <- vegan::vegdist(t(bin_continent_matrix_df_binary), 
                                         method = "jaccard", 
                                         binary = TRUE)
print(distance_matrix_binary)

# EDIT (Fangyi): The binary Jaccard dissimilarity (0.9782609) differs slightly 
# from the quantitative version (0.9976303) because it focuses solely on shared 
# BIN presence rather than abundance.
# This result still indicates that the BIN communities in Europe and North America
# are almost completely distinct, with very few shared taxa.

# ==== EDIT END ====

# Split the coordinates into latitude and longitude data. Use only latitude data for downstream analysis

unique(country_to_continent$`country/ocean`) #to confirm countries in my continents are still in country_to_continent data frame

country_to_continent_latitudeclean <- country_to_continent %>%
  dplyr::filter(!is.na(coord))


df_coords <- country_to_continent_latitudeclean %>% 
  tidyr::separate (coord, into = c("latitude", "longitude"), sep = ",", convert = TRUE)

latitudes <- df_coords$latitude
print(latitudes)

#remove "[" and non-numeric data except 0-9,decimal and negative sign from latitude data and make it numeric 

latitude_clean <- gsub("[^0-9.-]", "", latitudes)
print(latitude_clean)

latitudes_clean1 <- as.numeric(latitude_clean)
print(latitudes_clean1)


#Figure 2 to examine how bin richness vary by latitude at interval of 5 in North America and Europe 

df_latitude <- df_coords %>% 
  dplyr::mutate(lat_bin = cut(latitudes_clean1, breaks = seq(-90, 90, by = 5)))

unique(df_latitude$`country/ocean`) #confirming countries I selected are still there

# Count unique BINs per latitude bin at interval of 5 and region

bin_richness_latitude <- df_latitude %>%
  dplyr::group_by(region, lat_bin) %>%
  dplyr::summarise(bin_richness = dplyr::n_distinct(bin_uri))
print (bin_richness_latitude)

ggplot2::ggplot(bin_richness_latitude, 
                ggplot2::aes(x = lat_bin, y = bin_richness, fill = region)) +
  ggplot2::geom_col(position = "dodge") + 
  ggplot2::labs(title = "Ixodes BIN Diversity Across Latitude Bins",
                x = "Latitude Bin (5-degree intervals)",
                y = "Number of Unique BINs")

#Figure 3 to investigate whether Ixodes BINS occurring at higher latitude from these 2 continents have broader geographic ranges

df_latitude_bin <- df_coords %>% 
  dplyr::mutate(latitude_num = latitudes_clean1)

bin_latitude_summary <- df_latitude_bin %>%
  dplyr::group_by(region,bin_uri) %>%
  dplyr::summarise(
    mean_latitude = mean(latitude_num, na.rm = TRUE),
    latitude_range = max(latitude_num, na.rm = TRUE) - min(latitude_num, na.rm = TRUE))

head(bin_latitude_summary)

ggplot2::ggplot(bin_latitude_summary, 
                ggplot2::aes(x = mean_latitude, y = latitude_range)) +
  ggplot2::geom_point() + 
  ggplot2::geom_smooth(method = "lm", se = TRUE) + 
  ggplot2::labs(title = "Ixodes BIN Geographic (Latitudinal) Range vs. Mean Latitude",
                x = "Mean Latitude in degrees",
                y = "Latitude Range in degrees") 

#Testing the correlation between mean_latitude and latitude_range 

cor_test <- cor.test(bin_latitude_summary$mean_latitude, 
                     bin_latitude_summary$latitude_range, 
                     method = "spearman")
print(cor_test)






