#' Function currently under test.
#' @param rv rv void
#' @export
#'

rv_info <- function(){
  rvi <- readRDS("/Users/DL/Documents/R project/My Package/CNV/rv_info.rds")
  rvi
}


#'Building a Model with Top Ten Features
#'
#' This funciton develops a prediction algorithm based on the top10 features in 'x' that are most predictive of 'y'.
#'
#' @param data a dataframe of the data that is pased
#' @param keyMap_std a dataframe containing the location information of standards used.
#' @param calib a dataframe containing calibrator information
#' @return a list of results for data summary, std summary, differnce between estimated and rounded copy numbers and sample location.
#'
#' @author Dong Liang
#' @details
#' This function runs a multivaritate regression of y on each predictor in x and calculates a p-value indicating the significance of the association. The final set of 10 predictors is taken from the 10 smallest p-values.
#' @seealso \code lm
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr tbl_df
#' @importFrom dplyr arrange
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom magrittr "%>%"
#' @importFrom stats lm
#' @import ggplot2
#' @importFrom xtable xtable
#' @importFrom knitr kable
#' @import reshape2
#' @export

makeBatch <- function(data = data, keyMap_std = "default", calib, locMap_plate_96 = "rivur" ){

  ###############################################
  ## set default keyMap_std
  ###############################################
  keyMap_std_default <- data.frame(
    pos = c("B13", "B14", "B15",
            "D13", "D14", "D15",
            "F13", "F14", "F15",
            "H13", "H14", "H15",
            "J13", "J14", "J15",
            "L13", "L14", "L15"),
    sample = c(rep("K106957-", 3),
               rep("K44781-", 3),
               rep("K96097-", 3),
               rep("K94705-", 3),
               rep("K84083-", 3),
               rep("K75881-", 3)),
    copy_Nbs = rep(c(1.2, 1.4, 1.2, 1.88, 1.62, 1.22), each = 3)
  )

  if (keyMap_std == "default"){ keyMap_std <- keyMap_std_default}

  ###############################################
  ## make location map for plate of interest - locMap_plate
  ###############################################
  rivurPlate = tbl_df(rivurPlate)
  rivurPlate$POSITION = strsplit(as.character(rivurPlate$POSITION), split = " ")%>%unlist # clean whitespace
  plate_code = unique(data$plate_ID)

  if (locMap_plate_96 == "rivur") {
    locMap_plate = filter(rivurPlate, Plate_Box_Inv_Code %in% plate_code)
  } else if (locMap_plate_96 == "columbia") {
    locMap_plate = filter(columbiaPlate, Plate_ID %in% plate_code)
    }


  ###############################################
  ## calibration setting
  ###############################################
  data = tbl_df(data)
  calib_df = filter(data, Pos %in% calib) %>% group_by(Pos, Sample.Name, plate_ID) %>% summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt)
  calib_dCt = mean(calib_df$dCt)
  calib_dCt = 0

  ###############################################
  ## data pre-proceesing
  ###############################################
  # data sumamary
  data_su = group_by(data, Pos, Sample.Name, plate_ID) %>%
    summarize(ddCt = Cp[1]-Cp[2]-calib_dCt, RelQ = 2^-ddCt)

  # merge key map and plate data
  data_384 <- merge(map_p96_p384, data_su, by.x = "p384", by.y = "Pos", all = T)%>% tbl_df

  if (data_su$plate_ID[1]  == "Columbia_P2"){
    data_384 <- merge(map_p96_p384_Columbia_P2, data_su, by.x = "p384", by.y = "Pos", all = T)%>% tbl_df
  }

  data_merge <- merge(data_384, locMap_plate, by.x = "p96", by.y = "POSITION", all.x = T, sort = T ) %>% tbl_df


  ###############################################
  ## data summary
  ###############################################
  data_m_su <- group_by(data_merge, SAMPLE, plate_ID) %>% # CHANGE IT TO RUID
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)))  %>%
    filter(SAMPLE != "NA") # CHANGE TO RUID


  ###############################################
  ## standard info
  ###############################################
  std_merge <- merge(keyMap_std, data_su, by = "Pos", all.x = T)

  std_m_su <- group_by(std_merge, plate_ID, Sample, CopyNum) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)))

  std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)

  ###############################################
  ## copy number calling
  ###############################################
  data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
  #         data_m_su = cbind(data_m_su, predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean), interval = "prediction")[, c(2,3)])
  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)
  diff <- abs(data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round)
  data_m_su$diff <- diff
#   file <- paste(c("DEFA1A3_Columbia-", plate_code, ".csv"), collapse = "")
#   write.csv(file = file, data_m_su)


  ###############################################
  ## return results
  ###############################################
  result = list()
  result$data_su = tbl_df(data_m_su)
  result$std_su = tbl_df(std_m_su)
  result$diff = diff
  result$sampleLoc = select(data_merge, c(1,2,4,8))
  re = result

  class(re) = "copyNbsCall"
  re
}



#'Building a Model with Top Ten Features
#'
#' This funciton develops a prediction algorithm based on the top10 features in 'x' that are most predictive of 'y'.
#'
#' @param data a dataframe of the data that is pased
#' @param keyMap_std a dataframe containing the location information of standards used.
#' @param calib a dataframe containing calibrator information
#' @return a list of results for data summary, std summary, differnce between estimated and rounded copy numbers and sample location.
#'
#' @author Dong Liang
#' @details
#' This function runs a multivaritate regression of y on each predictor in x and calculates a p-value indicating the significance of the association. The final set of 10 predictors is taken from the 10 smallest p-values.
#' @seealso \code lm
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom dplyr tbl_df
#' @importFrom dplyr arrange
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize
#' @importFrom magrittr "%>%"
#' @importFrom stats lm
#' @import ggplot2
#' @importFrom xtable xtable
#' @importFrom knitr kable
#' @import reshape2
#' @export
makeBatch_ctl <- function(data = data, keyMap_std = "default", calib ){

  # library -----------------------------------------------------------------
  library(knitr)
  library(xtable)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  # library -----------------------------------------------------------------


  # key-map -----------------------------------------------------------------
  map <- read.table("/Users/DL/Dropbox-work/My work-not sync/Protocol/Protocol Memphis/RIVUR - IGLL5/Map/RIVUR_DILUTION_PLATES.csv",sep = ",", header = T)
  p384 <- sapply(LETTERS[1:16], function(x) paste0(x, 1:24))%>%t
  p96 <- sapply(LETTERS[1:8], function(x) paste0(x, 1:12)) %>%t

  v96 = vector("character")
  v96_l <- lapply(seq_along(p96[,1]), function(i) {
    a <- paste(paste(p96[i,],p96[i,]), p96[i,])
    b <- strsplit(a, split = " ") %>% unlist
    b[37:48] = NA
    v96 = c(v96, b)
  })
  # key-map -----------------------------------------------------------------

  # make a 384-96 map -------------------------------------------------------
  keyMap <- data.frame(p96 = unlist(v96_l), p384 = as.vector(t(p384)))
  # make a 384-96 map -------------------------------------------------------

  # set keyMap_std ----------------------------------------------------------
  keyMap_std_default <- data.frame(
    pos = c("B13", "B14", "B15",
            "D13", "D14", "D15",
            "F13", "F14", "F15",
            "H13", "H14", "H15",
            "J13", "J14", "J15",
            "L13", "L14", "L15"),
    sample = c(rep("K106957-", 3),
               rep("K44781-", 3),
               rep("K96097-", 3),
               rep("K94705-", 3),
               rep("K84083-", 3),
               rep("K75881-", 3)),
    copy_Nbs = rep(c(NA, 1.5, NA, 3.5, 2, 1.2), each = 3)
  )

  if (keyMap_std == "default"){ keyMap_std <- keyMap_std_default}
  # set keyMap_std ----------------------------------------------------------

  # retrieve info from key map ----------------------------------------------
  map = tbl_df(map)
  map$POSITION = strsplit(as.character(map$POSITION), split = " ")%>%unlist # clean whitespace
  plate_code = unique(data$plate_ID)
  map_key = filter(map, Plate_Box_Inv_Code %in% plate_code)

  # retrieve info from key map -----------------------------------


  # calibration setting -----------------------------------------------------
  data = tbl_df(data)
  calib_df = filter(data, Pos %in% calib) %>% group_by(Pos, Sample.Name, plate_ID) %>% summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt)
  calib_dCt = mean(calib_df$dCt)
  calib_dCt = 0

  # calibration setting -----------------------------------------------------

  # data pre-proceesing ------------------------------------------------------------

  # data sumamary
  data_su = group_by(data, Pos, Sample.Name, plate_ID) %>%
    summarize(ddCt = Cp[1]-Cp[2]-calib_dCt, RelQ = 2^-ddCt)

  # merge key map and plate data
  data_384 <- merge(keyMap, data_su, by.x = "p384", by.y = "Pos", all = T)%>% tbl_df

  data_merge <- merge(data_384, map_key, by.x = "p96", by.y = "POSITION", all.x = T, sort = T ) %>% tbl_df

  # data pre-proceesin ------------------------------------------------------------

  # data summary ------------------------------------------------------------
  data_m_su <- group_by(data_merge, p96, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              conc = mean(Concentration_ng._L))  %>%
    filter(p96 != "NA") %>%
    select(-conc)
  # data summary ------------------------------------------------------------

  # standard info -----------------------------------------------------------
  std_merge <- merge(keyMap_std, data_su, by.x = "pos", by.y = "Pos", all.x = T)

  std_m_su <- group_by(std_merge, plate_ID, sample, copy_Nbs) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)))

  std.fit <- lm(copy_Nbs ~ RelQ.mean, data = std_m_su)
  # standard info -----------------------------------------------------------

  # copy number calling -----------------------------------------------------
  data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)
  diff <- abs(data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round)
  data_m_su$diff <- diff
  file <- paste(c("NRXN3-", plate_code, ".csv"), collapse = "")
  write.csv(file = file, data_m_su)

  # copy number calling -----------------------------------------------------

  # results -----------------------------------------------------------------
  result = list()
  result$data_su = tbl_df(data_m_su)
  result$std_su = tbl_df(std_m_su)
  result$diff = diff
  result$sampleLoc = select(data_merge, c(1,2,4,8))
  re = result
  # results -----------------------------------------------------------------


}





#' This funciton takes the file location and convert raw data into CNV class for downstream analysis.
#'
#' @param file_list
#' @param stdKeyMap
#' @param calib
#' @importFrom dplyr data_frame
#' @return a data bundle containing all the raw data information
#' @export

makeDataBatch <- function(file_list, stdKeyMap, calib, locMap_plate_96 = "rivur", BadData_remove = TRUE){

    rawDataBundle <- function(file_list, stdKeyMap, calib){
    para = dplyr::data_frame(data = file_list,
                      stdKeyMap = stdKeyMap,
                      calib = calib)
    para
  }

  removeBadData <- function(df, Cp_value = ifelse(BadData_remove, 35, 50)){
    Sample.Name_good <- df %>% filter(between(Cp, 15, Cp_value)) %>% .$Sample.Name %>% unique()
    df <- df %>% filter(Sample.Name %in% Sample.Name_good)
    df
  }
# file_list = rawData$file_list

  dataBundle <- rawDataBundle(file_list, stdKeyMap, calib)

  dataBatch <- lapply(1: nrow(dataBundle),
                  function(i) makeBatch(data = dataBundle[i,][[1]] %>% data.frame %>% removeBadData,
                                        keyMap_std = dataBundle[i,][[2]],
                                        calib = dataBundle[i,][[3]],
                                        locMap_plate_96 = locMap_plate_96))
# Wrap up data
  data_batch <- lapply(dataBatch, function(x) x$data_su)
  data <- do.call(rbind, data_batch)

  std_batch <- lapply(dataBatch, function(x) x$std_su)
  std <- do.call(rbind, std_batch)

#   sampleLoc_batch <- lapply(dataBatch, function(x) x$sampleLoc)
#   sampleLoc <- do.call(rbind, sampleLoc_batch)

  diff_batch <- lapply(dataBatch, function(x) x$diff)
  diff <- do.call(rbind, diff_batch)

  diff_batch <- lapply(dataBatch, function(x) x$diff)
  diff <- do.call(rbind, diff_batch)

  output = list()
  output$data = data
  output$std = std
  # output$sampleLoc = sampleLoc
  output$diff = diff
  output

}












