# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OCE-TIMS Retailer Inspection Data
# Mar 3, 2022
# Sandeep Shetty
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# >>> Load packages <<<
# .......Packages needed
want <- c(
  "tidyverse", "lubridate",
  "geosphere", "tidymodels",
  "doParallel", "tibble",
  "skimr"
)
# ......Check, Install, Load packages
have <- want %in% rownames(installed.packages())
if (any(!have)) {
  install.packages(want[!have])
}
sapply(want, library, character.only = TRUE)
rm(have, want)

# >>> Path to data <<<
path_to_data <- "data/floridaOCE.RDS"

# >>> Read and Clean Data <<<
clean_data <- function(file_path) {
  # ......................................
  # Function to set up starting data
  # to ease resetting data upon errors
  # specific cleaning for the Florida TIMS
  # ......................................
  dats <- readRDS(file_path)
  dats <- dats[dats$MATCH == "1", ]
  names(dats) <- tolower(names(dats)) %>%
    stringr::str_replace_all(., "\\.", "_")

  # ... Keep only FL
  dats <- dats %>% filter(grepl("FL", rei))

  # ... Variable conversions
  dats$sale_to_minor <- as.integer(dats$sale_to_minor)
  dats$inspect_dt <- lubridate::as_date(dats$inspection_date,
    format = "%m/%d/%Y"
  )
  # ... Subset variables
  keep_vars <- c(
    "rei", "type", "outcome",
    "category", "sale_to_minor",
    "city", "longitude", "latitude",
    "inspect_dt", "zip"
  )
  dats <- dats %>% select(all_of(keep_vars))

  # ... Re-arrange, create some variables
  dats <- dats %>% mutate(inspect_year = lubridate::
  year(inspect_dt))
  dats <- dats %>%
    group_by(rei, inspect_year) %>%
    mutate(insp_in_year = n())
  dats <- dats %>%
    group_by(zip, inspect_year) %>%
    mutate(insp_in_zip_year = n())

  # Create days between inspections within a zip code
  insp_min_date <- dats %>%
    group_by(zip, inspect_year) %>%
    summarise(
      min_dat = min(inspect_dt),
      .groups = "drop"
    )

  # ... (average) duration between inspections in {zip-year}
  dats <- merge(dats, insp_min_date, by = c("zip", "inspect_year"))
  dats <- dats %>% mutate(insp_days = (inspect_dt - min_dat))
  dats <- dats %>%
    group_by(zip, inspect_year) %>%
    mutate(avg_days_insp_zip_yr = mean(insp_days))

  # ... convert lat and long into numeric
  cols_num <- c("latitude", "longitude")
  dats[cols_num] <- sapply(dats[cols_num], as.numeric)
  dats$rei <- as.factor(dats$rei)
  dats$zip <- as.factor(dats$zip)
  dats <- dats %>% filter(inspect_year != 2014)
  return(dats)
}

# >>> Run function to import data
fl <- clean_data(path_to_data)

# >>> Calculate Distance Between Retailers
# Unique REIs only
fl_uniq <- fl[!duplicated(fl$rei), ]

# >>> keep if more than 1 retailers in zip code
fl_uniq <- fl_uniq %>%
  group_by(zip) %>%
  mutate(totn = n()) %>%
  ungroup()

fl_uniq <- fl_uniq %>% filter(totn > 1)
fl_uniq$index_no <- seq(1, nrow(fl_uniq))

# >>> Calculating "Haversine" distance between retailers
distance_calc <- function(dat) {
  # ...................................................
  # To minimize processing, distance is calculated only between retailers
  # within the same zip code.
  # Input: data with retailer id, zip code, latitude and longitude
  # Output: distance (meters) matrix between retailers within same zip code
  # ....................................................
  dim1 <- nrow(dat)
  distance_matrix <- matrix(nrow = dim1, ncol = dim1)
  dat1 <- dat %>% select(index_no, rei, zip, longitude, latitude)
  for (zyp in levels(dat1$zip)) {
    geo_data <- dat1 %>% filter(zip == zyp)
    for (i in geo_data$index_no) {
      lati <- geo_data$latitude[which(geo_data$index_no == i)]
      loni <- geo_data$longitude[which(geo_data$index_no == i)]
      for (j in geo_data$index_no) {
        if (i > j) {
          latj <- geo_data$latitude[which(geo_data$index_no == j)]
          lonj <- geo_data$longitude[which(geo_data$index_no == j)]
          distance_matrix[i, j] <- geosphere::distm(
            c(loni, lati),
            c(lonj, latj)
          )
        }
      }
    }
  }
  return(distance_matrix)
}

# >>> Calculate distance between retailers
distance_matrix <- distance_calc(dat = fl_uniq)
distance_matrix <- distance_matrix / 1604 # in miles

# ... Row and column names for look-up using REI
rownames(distance_matrix) <- fl_uniq$rei
colnames(distance_matrix) <- fl_uniq$rei

# >>>  Indicator for retailers within a 'x' miles of a retailer
no_ret_miles <- function(miles, dist_matrix) {
  # .............................................
  # Input: distance matrix assumed
  # Input: enter miles to use from each retailer
  # output a matrix with 1/0 to indicate retailers is within x-mile
  # distance of a retailer
  # ..................................
  neighbor_in_x <- apply(dist_matrix, 2, function(x) (x <= miles) / 1)
  return(neighbor_in_x)
}

# >>>  Roll up stats for each retailer-year
# ....for violations in a x-mile neighborhood
# ....to output table --> {REI|Year|#violations|#count}
# ....merge this o/p to the main data set

# >>>  Neighbors within x=1 mile distance of a given REI
neighbor_in_x <- no_ret_miles(
  miles = 1,
  dist_matrix = distance_matrix
)

# ... Row and column names for look-up using REI
rownames(neighbor_in_x) <- fl_uniq$rei
colnames(neighbor_in_x) <- fl_uniq$rei

# >>> Calculates the number of inspection and violations
rei_ng <- function(reid) {
  # ......................................................
  # Calculates for the number of inspection and violations
  # with a given radius of a retailer for each year
  # Output: Dataframe
  # .....................................................
  select_reis <- names(which(neighbor_in_x[, reid] == 1))
  if (length(select_reis) > 0) {
    temp_dat <- fl %>% filter(rei %in% select_reis)
    new_temp_dat <- temp_dat %>%
      group_by(inspect_year) %>%
      # no of inspections
      mutate(
        neigh_insp_cnt = n(),
        # violations
        neigh_viol_cnt = sum(sale_to_minor)
      ) %>%
      distinct(
        inspect_year,
        neigh_insp_cnt,
        neigh_viol_cnt
      )
    new_temp_dat$rei <- reid # add REI id
    return(new_temp_dat)
  }
}

# >>> Create a data frame with violations in neighboring retailers
neigh_viols <- function() {
  # ...............................................
  # calculate number of inspections and violation within a defined
  # radius of each retailer by year
  # See function rei_ng for miles inputed (def = 1 miles)
  # output is a summary of violations and inspections
  # within a mile of a given retailer (grouped by year)
  # ................................................
  newdat <- data.frame()
  for (rei_id in fl_uniq$rei) {
    datInt <- rei_ng(rei_id)
    newdat <- rbind(newdat, datInt)
  }
  return(newdat)
}

# >>> save this data as it is not efficiently executed
if (file.exists("data/neighViolations.RDS")) {
  neigh_viols_data <- readRDS("data/neighViolations.RDS")
} else {
  neigh_viols_data <- neigh_viols()
  saveRDS(object = neigh_viols_data, file = "data/neighViolations.RDS")
}

# >>> combine neighborhood violations
fl_neigh_viol <- fl %>% left_join(neigh_viols_data,
  by = c("rei", "inspect_year")
)
# >>> replace "NA" in violation and inspections as O
fl_neigh_viol <- fl_neigh_viol %>%
  dplyr::mutate(
    neigh_insp_cnt = replace_na(neigh_insp_cnt, 0),
    neigh_viol_cnt = replace_na(neigh_viol_cnt, 0)
  )

# >>> average distance between retailers within zip
avg_distance_ret_zip <- colSums(distance_matrix,
  na.rm = TRUE
) / colSums(!is.na(distance_matrix))

avg_distance_ret_zip <- avg_distance_ret_zip %>%
  as_tibble() %>%
  mutate(rei = fl_uniq$rei)

# >>> merge with
fl_neigh_viol <- fl_neigh_viol %>%
  left_join(avg_distance_ret_zip, by = c("rei")) %>%
  rename(avg_dist_zip = value)

# >>> who was the first to be inspected
fl_neigh_viol <- fl_neigh_viol %>%
  mutate(
    first_inspected = case_when(
      insp_days == 0 ~ 1, TRUE ~ 0
    )
  )


# ... inspected in cluster (a.k.a. same day with many retailers)
fl_neigh_viol <- fl_neigh_viol %>%
  group_by(zip, inspect_year, inspect_dt) %>%
  mutate(clstinsptemp = n()) %>%
  mutate(clust_insp = case_when(
    clstinsptemp > 1 ~ 1, TRUE ~ 0
  ))

# ... one retailer or a group of retailers inspected first in the zip
first_insp <- fl_neigh_viol %>%
  filter(first_inspected == 1) %>%
  group_by(zip, inspect_year) %>%
  select(rei) %>%
  ungroup()

# ... distance from a retailer inspected nearest in days (ago)

fl_neigh_viol <- fl_neigh_viol %>% ungroup()

# focusing only non 1 timers
fl_not_first <- fl_neigh_viol %>%
  filter(first_inspected == 0)

# function that applies only to rei's that are not the first ones
# inspected in their {zip-year}

nneigh_insp_test2 <- function(rei_id, zip1, yr_insp, days_before) {
  collect_rei <- fl_neigh_viol[
    which(fl_neigh_viol$insp_days < days_before &
      fl_neigh_viol$inspect_year == yr_insp &
      fl_neigh_viol$zip == zip1 &
      fl_neigh_viol$rei != rei_id),
    c("rei")
  ]
  collect_rei <- unique((as.vector(collect_rei[[1]])))
  print(is.vector(collect_rei))
  if (length(collect_rei) > 0) {
    # lower tri-matrix - subsetting
    dist_neigh <- distance_matrix[
      collect_rei,
      rei_id
    ]
    print(dist_neigh)
    if (sum(!is.na(dist_neigh)) > 0) {
      return(min(dist_neigh, na.rm = TRUE))
    }
  }
  return(NA)
}
# nneigh_insp <- function(rei_id, zip1, yr_insp, days_before) {
#   #if (days_before > 0) {
#    collect_rei <- fl_neigh_viol %>%
#       filter(
#         insp_days < days_before,
#         inspect_year == yr_insp,
#         zip == zip1,
#         rei != rei_id
#       ) %>%
#       select(rei) %>%
#       unique() %>%
#       pull()
#      if (length(collect_rei) > 0) {
#       #try(
#         #{
#       dist_neigh <- distance_matrix[collect_rei, rei_id] #lower tri matrix
#         if (sum(!is.na(dist_neigh)) > 0){
#           return(min(dist_neigh, na.rm=TRUE))}
#         }
#         #)
#    return(NA)
#     }
# }
# return(NA)
# }

library(parallel)
nrst_ngh_insptd_bef <- mcmapply(
  nneigh_insp_test2, fl_neigh_viol$rei,
  fl_neigh_viol$zip,
  fl_neigh_viol$inspect_year,
  fl_neigh_viol$insp_days
)

fl_neigh_viol <- cbind(fl_neigh_viol, nrst_ngh_insptd_bef)

# ...
# fl_neigh_viol <-
#  fl_neigh_viol %>%
#  dplyr::rowwise() %>%
#  mutate(nrst_ngh_insptd = nneigh_insp(rei, zip, inspect_year, insp_days))

# split the column into two separate columms
# fl_neigh_viol <- fl_neigh_viol %>% tidyr::unnest_wider(newvar)

# >>> Pending Feature Engineering
# ... mean distance of those inspected before
# ... minimum distance of those insp+violated before
# ... mean distance of those insp+violated before

# >>> data type conversion - implement in analysis file
fl_neigh_viol$insp_days <- fl_neigh_viol$insp_days %>% as.numeric()
fl_neigh_viol$avg_insp_zip_yr <- fl_neigh_viol$avg_insp_zip_yr %>%
  as.numeric() %>%
  replace(is.na(.), 0)
fl_neigh_viol <- fl_neigh_viol %>% na.omit()
fl_neigh_viol <- fl_neigh_viol %>% filter(!is.na(avg_dist_zip))


# >>> save file
saveRDS(fl_neigh_viol, file = "data/floridaOCEClean.RDS")
write.csv(fl_neigh_viol,
  file = "data/florida_tims_clean.csv",
  row.names = FALSE
)
