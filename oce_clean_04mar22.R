# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# OCE-TIMS Data Set Up
# Spatial effects of retailer inspection (local neighborhood)
# Mar 3, 2022
# Sandeep Shetty
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~ Load packages
source("code/packages_needed.R")

# ~~~~~~ Read Data
gen_start_data <- function(file_path = "data/floridaOCE.RDS") {
  # ......................................
  # Function to set up starting data
  # to ease resetting data upon errors
  # specific cleaning for the Florida TIMS
  # ......................................
  dats <- readRDS(file_path)
  dats <- dats[dats$MATCH == "1", ]
  names(dats) <- tolower(names(dats))
  # ... Clean some variables
  dats$sale.to.minor <- as.integer(dats$sale.to.minor)
  dats$inspect.dt <- lubridate::as_date(dats$inspection.date,
    format = "%m/%d/%Y"
  )
  # ... Subset variables
  keep_vars <- c(
    "rei", "type", "outcome",
    "category",
    "sale.to.minor", "city",
    "longitude", "latitude",
    "inspect.dt", "zip"
  )
  dats <- dats %>% select(all_of(keep_vars))

  # ... Re-arrange, create some variables
  dats <- dats %>% mutate(inspect.year = lubridate::
  year(inspect.dt))
  dats <- dats %>%
    group_by(rei, inspect.year) %>%
    mutate(insp.in.year = n())
  dats <- dats %>%
    group_by(zip, inspect.year) %>%
    mutate(insp.in.zip.year = n())

# Create days between inspections within a zip code
  insp_min_date <- dats %>%
    group_by(zip, inspect.year) %>%
    summarise(
      min_dat = min(inspect.dt),
      .groups = "drop"
    )

  # ... (average) duration between inspections in {zip-year}  
  dats <- merge(dats, insp_min_date, by = c("zip", "inspect.year"))
  dats <- dats %>% mutate(insp_days = (inspect.dt - min_dat))
  dats <- dats %>%
    group_by(zip, inspect.year) %>%
    mutate(avg_insp_zip_yr = mean(insp_days))

  # ... convert lat and long into numeric
  cols.num <- c("latitude", "longitude")
  dats[cols.num] <- sapply(dats[cols.num], as.numeric)
  dats$rei <- as.factor(dats$rei)
  dats$zip <- as.factor(dats$zip)
  dats <- dats %>% filter(inspect.year != 2014)
  return(dats)
}

# ~~~~~~ Run function to import data
fl <- gen_start_data()

# ~~~~~~ Calculate Distance Between Retailers
# Unique REIs only
fl.uniq <- fl[!duplicated(fl$rei), ]

# .. keep if more than 1 retailers in zip code
fl.uniq <- fl.uniq %>%
  group_by(zip) %>%
  mutate(totn = n()) %>%
  ungroup()

fl.uniq <- fl.uniq %>% filter(totn > 1)
fl.uniq <- fl.uniq %>% filter(!fl.uniq$rei == "NJ707536")
fl.uniq$index_no <- seq(1, nrow(fl.uniq))

#~~~~~~ Calculating "Haversine" distance between retailers
distance_calc <- function(dat) {
  # ...................................................
  # To minimize processing, distance is calculated only between retailers
  # within the same zip code.
  # Input: data with retailer id, zip code, latitude and longitude
  # Output: distance (meters) matrix between retailers within same zip code
  # ....................................................
  dim1 <- nrow(dat)
  distanceMatrix <- matrix(nrow = dim1, ncol = dim1)
  dat1 <- dat %>% select(index_no, rei, zip, longitude, latitude)
  # dat1 <- dat1 %>% filter()
  for (zyp in levels(dat1$zip)) {
    LongLat <- dat1 %>% filter(zip == zyp)
    for (i in LongLat$index_no) {
      lati <- LongLat$latitude[which(LongLat$index_no == i)]
      loni <- LongLat$longitude[which(LongLat$index_no == i)]
      for (j in LongLat$index_no) {
        if (i > j) {
          latj <- LongLat$latitude[which(LongLat$index_no == j)]
          lonj <- LongLat$longitude[which(LongLat$index_no == j)]
          distanceMatrix[i, j] <- geosphere::distm(
            c(loni, lati),
            c(lonj, latj)
          )
        }
      }
    }
  }
  return(distanceMatrix)
}

# Calculate distance between retailers
distanceMatrix <- distance_calc(dat = fl.uniq)
distanceMatrix <- distanceMatrix / 1604 # in miles

# ~~~~~~~ Indicator for retailers within a 'x' miles of a retailer
no_ret_miles <- function(miles, distMatrix) {
  # .............................................
  # Input: distance matrix assumed
  # Input: enter miles to use from each retailer
  # output a matrix with 1/0 to indicate retailers is within x-mile
  # distance of a retailer
  # ..................................
  neighbor.in.x <- apply(distMatrix, 2, function(x) (x <= miles) / 1)
  return(neighbor.in.x)
}

# ~~~~~~~~ Roll up stats for each retailer-year
# ....for violations in a x-mile neighborhood
# ....to output table --> {REI|Year|#violations|#count}
# ....merge this o/p to the main data set

# Neighbors within x=1 mile distance of a given REI
neighbor.in.x <- no_ret_miles(
  miles = 1,
  distMatrix = distanceMatrix
)

# Row and column names for easy indexing/look-up
rownames(neighbor.in.x) <- fl.uniq$rei
colnames(neighbor.in.x) <- fl.uniq$rei

rei_ng <- function(reid) {
  #......................................................
  # Calculates for the number of inspection and violations
  # with a given radius of a retailer for each year
  # Output: Dataframe
  # .....................................................
  select.reis <- names(which(neighbor.in.x[, reid] == 1))
  if (length(select.reis) > 0) {
    temp.dat <- fl %>% filter(rei %in% select.reis)
    new.temp.dat <- (temp.dat %>%
      group_by(inspect.year) %>%
      # no of inspections
      mutate(
        neigh.insp.cnt = n(),
        # violations
        neigh.viol.cnt = sum(sale.to.minor)
      ))
    new.tmp.dat1 <- new.temp.dat %>% distinct(
      inspect.year,
      neigh.insp.cnt,
      neigh.viol.cnt
    )
    new.tmp.dat1$rei <- reid # add REI id
    return(new.tmp.dat1)
  }
}

#~~~~~~ Create a data frame with violations in neighboring retailers
neigh.viols <- function() {
  # ...............................................
  # calculate number of inspections and violation within a defined
  # radius of each retailer by year
  #................................................
  newdat <- data.frame()
  for (rei_id in fl.uniq$rei) {
    ## print(rei_id)
    datInt <- rei_ng(rei_id)
    newdat <- rbind(newdat, datInt)
  }
  return(newdat)
}

### save this data as it is not efficiently executed
## ....
if (file.exists("data/neighViolation.RDS")) {
  newdat <- readRDS(intdat, "data/neighViolations.RDS")
} else {
  newdat <- neigh.viols()
}


# average distance between retailers within zip



# who was the first to be inspected

# which group was inspected

# one retailer or a group of retailers inspected first



# Understand relationship of inspection and violation within zip-level
# How should be estimate the impact of neighborhood inspections on
# own violations.
# Need to consider the timing of these inspections


# need to figure out timing on inspections,
# follow-up information
# indicator for who is inspected before or after
# some more indicators
# then take to the model
# then think
