#Script pour utiliser les fonctions CESI fun

#library(tidyverse)
library(readxl)    # for importing excel file
library(tidyhydat) # for interacting with hydat database
library(dplyr)     # For data tidying and data tables
library(zoo)       # For working with time series
library(evir)      # For POT
library(purrr)     # function 'possibly' is easiest way to deal with errors
library(lfstat)    # package for low flow statistics
library(tibble)    # to work with tables
library(zyp)       # for function in Mann-Kendall test
library(MASS)      # for function in Negative Binomial & hurdle tests
library(countreg)  # for hurdle test
source('mainS.R')

#-------------------------------------------------------------------------------
# set up database
chk.hydat()

#-------------------------------------------------------------------------------
# import list of RHBN stations
url <- "https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/RHBN/RHBN_Metadata.xlsx"
destfile <- "./Output/RHBN_Metadata.xlsx"
download.file(url, destfile, method="curl")
stations <- read_xlsx(destfile, range="A3:P1286")
stations <- stations[-(1),]
stations <- filter(stations, Evaluation_Year==2020, DATA_TYPE=="Q")
yrs.of.int <- c(2001:2019)

#-------------------------------------------------------------------------------
# calculate annual mean flow
annual_mean_flow <- list()
for (l in 1:length(stations$STATION_NUMBER)){
ann_mean_flow <- Get_ann_mean_flow(stations$STATION_NUMBER[l])
annual_mean_flow[[l]] <- ann_mean_flow
}
annual_mean_flow <- bind_rows(annual_mean_flow)

# summarize mean annual flow with statistics on mean annual yield
snap <- list()
yield <- list()
for (k in 1:length(stations$STATION_NUMBER)){
  stn <- as.character(stations[k,1])
  data <- annual_mean_flow %>% filter(STATION_NUMBER==stn)
  data <- unique(data[order(data$Year),])
  print(stn)

  # set up table
  q <- data.frame(matrix(ncol=4, nrow=1,
                         dimnames= list(NULL, c("STATION_NUMBER", "q25", "q75",
                                                "median"))))
  q$STATION_NUMBER <- stn

  # Scale flow by watershed size to get yield. Units will now be in mm.
  hydat_info <- hy_stations(stn)
  stn_area <- hydat_info$DRAINAGE_AREA_GROSS

    time <-  365*24*3600
  data$ann_mean_yield <- data$ann_mean_flow*time/(stn_area*10^3)
  yield[[k]] <-  data

  # calculate quantiles and median variable based on 1981-2010 period, provided
  # there are at least 20 years of data within that period.
  data.ref <- data %>% filter(Year>=1981, Year<=2010)
  if (sum(!is.na(data.ref[["ann_mean_yield"]]))>=20){
    q.var <- round(quantile(data.ref[["ann_mean_yield"]], c(0.25, 0.75), na.rm=TRUE), 2)
    median <- round(median(data[["ann_mean_yield"]], na.rm=TRUE), 2)
  } else {
    q.var <- c(NA, NA)
    median <- NA
  }
  q[, c("q25", "q75")] <- q.var
  q[, "median"] <median

  # calculate percentile and ranks quantiles and ranks for years of interest
  for (refyear in yrs.of.int){
    if (refyear %in% data$Year){
      latest <- round(data[data$Year == refyear, "ann_mean_yield"], 2)
      print(refyear)
    } else {
      latest <- NA
    }
    if (!is.na(latest) & !is.na(median)){
      q[, paste0("rk.", refyear)] <-
        as.integer(percentile.ranked(data.ref[["ann_mean_yield"]], latest))
    } else {
      q[, paste0("rk.", refyear)] <- NA
    }
    q <- add_column(q, latest, .before = paste0("rk.",yrs.of.int[1]))
    names(q)[names(q) == 'latest'] <- paste0("res.", refyear)
  }
# combine all the stations and write output in a file
  snap[[k]] <- q
}
snap.all <- bind_rows(snap)
annual_mean_flow<-bind_rows(yield)
annual_mean_flow$ann_mean_yield <- round(annual_mean_flow$ann_mean_yield,2)
write.csv(snap.all, paste0("./Output/Summary.AnnMeanYield.csv"), row.names = FALSE)


#-------------------------------------------------------------------------------
# calculate flood metrics, specifically, for each year of interest:
#     - maximum flow,
#     - a threshold at 95th percentile of flow in reference period,
#     - number of days over threshold flow,
#     - number of events, and
#     - maximum duration of events

flood_metrics <- list()
for (i in 977:length(stations$STATION_NUMBER)) {
f_m <- Get_flood_metrics(stations$STATION_NUMBER[i], year="all")
flood_metrics[[i]] <- f_m
}
flood_metrics <- bind_rows(flood_metrics)

#-------------------------------------------------------------------------------
# calculate drought metrics, specifically, for each year of interest:
#     - minimum 7-day rolling average within the bdays
#     - a threshold at 5th percentile of flow in reference period,
#     - number of days under threshold flow,
#     - number of events, and
#     - maximum duration of events

drought_metrics <- list()
for (j in 1:length(stations$STATION_NUMBER)) {
    d_m <- Get_drought_metrics(stations$STATION_NUMBER[j], year="all")
  drought_metrics[[j]] <- d_m
}
drought_metrics <- bind_rows(drought_metrics)


#-------------------------------------------------------------------------------
# output flood and drought metrics in file
metrics <- merge(annual_mean_flow,flood_metrics, all=TRUE)
metrics <- merge(metrics, drought_metrics, all=TRUE)
write.csv(metrics,"./Output/Summary.FandDmetrics.csv", row.names = FALSE)


#-------------------------------------------------------------------------------
# calculate trends

# define metrics to calculate trends
var_list <- c("ann_mean_yield", "pot_days")
result_list <- paste0(var_list, "_trend")

# set up output files

for (j in 1:length(var_list)){
  var.t = var_list[j]
  print(var.t)
  output_name = paste0("./Output/Summary_", var.t, "_trends.csv")
  snap <- list()
  for (i in 1:length(stations$STATION_NUMBER)){
    stn.id <- stations$STATION_NUMBER[i]
    print(stn.id)

    # defaults to NA if data requirements aren't met
    trend <- data.frame(NA, NA, NA, NA, NA)
    colnames(trend) <- c("slope", "intercept",  "CATTrend", "years.for.trend",
                         "test")
    mapslope <- NA
    hurdlechk <- NA
    data <- metrics[(metrics$STATION_NUMBER==stn.id),]
    #read.csv(output1, header = TRUE)
    if (sum(!is.na(data[[var.t]]))>=30){
      data <- data[!is.na(data[[var.t]]),]
      data.p <- data[data$Year>=1970,]
      data.p <- data.p[data.p$Year<=2019,] # cap data range for 2022 CESI release

      # Data requirements: some data 1970-1975, >=30 points, no gap over 10 years
      goodyears <- data.p$Year[!is.na(data.p[[var.t]])]
      gap.check <- na.omit(goodyears - lag(goodyears))

      if (all(any(goodyears %in% c(1970:1975)), (length(goodyears) >= 30),
              (max(gap.check) <= 11))){

         # Are there any zero values?
        if((sum(data.p[[var.t]]==0)<=1)&(var.t %in% c("X1_day_max", "ann_mean_yield",
                                                      "X7_day_min"))){
          print("Mann-Kendall test")
          trend <-mk_test(data.p[,var.t],data.p[,"Year"],keep_z=TRUE)
          trend$test <- "Mann-Kendall"
        }else{
          # Is hurdle necessary?
          hurdle <- FALSE #default to False
          if (sum(data.p[[var.t]]==0) >= 3){
            model2 <- tryCatch(hurdle(data.p[[var.t]]~data.p$Year, data.p,
                                dist="negbin", zero.dist = "negbin"),
                                error=function(e){return("A")},
                                warning=function(w){return("B")})
            if(!is.character(model2)){
              hurdlechk <- hurdletest(model2)[2,4]
              if(!is.nan(hurdlechk)){
                hurdle <- ifelse(hurdlechk>0.1, TRUE, FALSE)
              }
            }
          }

          if(hurdle){
            #Apply the hurdle model
            print("Hurdle test")
            trend <- hurdle_test(data.p[,var.t],data.p[,"Year"])
            trend$test <- "Hurdle"
          } else {
            #Apply the negative binomial model
            print("Negative Binomial test")
            trend <-negbin(data.p[[var.t]],data.p[["Year"]])
            trend$test <- "Negative Binomial"
          }
        }
      }
    }

    # mapslope field only has values for likely/confident trends for mapping
    mapslope <- case_when( grepl("Likely", trend$CATTrend)    ~ trend$slope,
                           grepl("Confident", trend$CATTrend) ~ trend$slope)
    # Load data and subset
    snap[[i]] <- data.frame(STATION_NUMBER=stn.id, slope=trend$slope,
                            intercept=trend$intercept,
                            years.for.trend=trend$years.for.trend,
                            CATTrend=trend$CATTrend, test=trend$test,
                            hurdlechk = hurdlechk, mapslope=mapslope)
  }

  snap.all <- bind_rows(snap)

  if (file.exists(output_name)){
    file.remove(output_name)
  }
  write.csv(snap.all, output_name, row.names = FALSE)
}
