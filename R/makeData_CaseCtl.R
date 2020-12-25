#' Function currently under test for case control plate data
#' @param data_raw
#' @param keyMap_std
#' @param zeroCopy
#' @export makeData_CaseCtl
#' @export makeData_RIVUR
#' @export makeData_cutie
#' @export makeData
#' @export makeData_bp
#' @export makeData_dCt
#' @export makeData_Coriell
#'

makeData_CaseCtl <- function(data_raw, keyMap_std = "default", zeroCopy = NULL){
  ###############################################
  ## Data processing
  ###############################################
  # data_raw = rawData$data_raw

  # if(hasZero) {
  #     data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
  #     summarize(ddCt = Cp[1]-Cp[2], RelQ = 2^-ddCt, zeroCopy = zero[1])
  #   } else{
  #     data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
  #     summarize(ddCt = Cp[1]-Cp[2], RelQ = 2^-ddCt)
  #   }

  data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
    summarize(ddCt = Cp[1]-Cp[2], RelQ = 2^-ddCt)

  Ctl <- sapply(LETTERS[1:8], function(x) paste0(x, 1:6))
  Rivur <- sapply(LETTERS[1:8], function(x) paste0(x, 7:12))

  # creat old to new map for rivur10 and ctl/case plates
  C <- as.vector(Ctl)
  R <- as.vector(Rivur)
  map_LtoR = data.frame(old = R, rivur10_loc = C)

  rivur09 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013509")
  rivur10 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013510")
  rivur11 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013511")
  rivur12 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013512")
  rivur13 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013513")

  data_su = left_join(data_su, rivur10, by = c("Pos" = "p384")) %>% select(-RUID)

  data_su$Type = ifelse(data_su$p96 %in% Rivur, "RIVUR", "CTL")

  if(data_su$plate_ID.x == "50_50_A"){
    data_su = dplyr::left_join(data_su, map_LtoR, by = c("p96" = "old"))
  }else if(data_su$plate_ID.x == "50_50_B")
  {data_su$rivur10_loc = ifelse(data_su$Type == "RIVUR", as.character(data_su$p96), NA)
  }else {data_su$rivur10_loc = ifelse(data_su$Type == "RIVUR", as.character(data_su$p96), NA)}


  data_su = dplyr::left_join(data_su, rivur10[,c("p96", "RUID")]%>% unique, by = c("rivur10_loc" = "p96"))  # get RUID number
  # rivur10 %>% unique()

  if(!is.null(zeroCopy)){
    data_su = left_join(data_su, zeroCopy[,c("plate_ID.x","Pos.x", "zero")], by = c("plate_ID.x" = "plate_ID.x","Pos" = "Pos.x"))
  }else(data_su$zero = NA)

  data_m_su <- group_by(data_su, p96, plate_ID.x) %>%
      summarize(ddCt.mean = mean(ddCt),
                ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
                RelQ.mean = mean(RelQ),
                RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
                RUID = RUID[1],
                Type = Type[1],
                p384 = Pos[1],
                zeroCopy = zero[1])

  # data_m_su$Type = ifelse(data_m_su$p96 %in% Rivur, "RIVUR", "CTL")


  ###############################################
  ## Make standards
  ###############################################
  #

  std_su = left_join(keyMap_std, data_su, by = c("plate_ID" = "plate_ID.x", "Pos" = "Pos"))
  std_m_su <- group_by(std_su, Sample, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              p96 = unique(p96),
              CopyNum = mean(CopyNum)) %>%
    select(-p96) %>%
    separate(col = Sample, into = c("Sample", "Source"))


  ###############################################
  ## Copy Number Call
  ###############################################
  std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)

  data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
  #         data_m_su = cbind(data_m_su, predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean), interval = "prediction")[, c(2,3)])
  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)
  diff <- abs(data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round)
  data_m_su$diff <- diff



  output = list()
  output$data_su = data_su
  output$data_m_su = data_m_su

  output$std_su = std_su
  output$std_m_su = std_m_su
  output
}




# makeData_CaseCtl(rawData$data_raw)

##  testing ----------------------------

makeData_CaseCtl_old <- function(data_raw, keyMap_std = "default", zeroCopy = NULL){
  ###############################################
  ## Data processing
  ###############################################
  # data_raw = rawData$data_raw

  # if(hasZero) {
  #     data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
  #     summarize(ddCt = Cp[1]-Cp[2], RelQ = 2^-ddCt, zeroCopy = zero[1])
  #   } else{
  #     data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
  #     summarize(ddCt = Cp[1]-Cp[2], RelQ = 2^-ddCt)
  #   }

  data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
    summarize(ddCt = Cp[1]-Cp[2], RelQ = 2^-ddCt)

  Ctl <- sapply(LETTERS[1:8], function(x) paste0(x, 1:6))
  Rivur <- sapply(LETTERS[1:8], function(x) paste0(x, 7:12))

  # data_su$Type = ifelse(data_su$Pos %in% Rivur, "RIVUR", "CTL")


  rivur10 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013510")
  data_su = left_join(data_su, rivur10, by = c("Pos" = "p384"))




  # rivur10 %>% unique()

  if(!is.null(zeroCopy)){
    data_su = left_join(data_su, zeroCopy[,c("plate_ID.x","Pos.x", "zero")], by = c("plate_ID.x" = "plate_ID.x","Pos" = "Pos.x"))
  }else(data_su$zero = NA)

  data_m_su <- group_by(data_su, RUID, plate_ID.x) %>%
      summarize(ddCt.mean = mean(ddCt),
                ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
                RelQ.mean = mean(RelQ),
                RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
                p96 = unique(p96),
                p384 = Pos[1],
                zeroCopy = zero[1])  %>%
      filter(RUID != "NA")
  data_m_su$Type = ifelse(data_m_su$p96 %in% Rivur, "RIVUR", "CTL")


  ###############################################
  ## Make standards
  ###############################################
  #

  std_su = left_join(keyMap_std, data_su, by = c("plate_ID" = "plate_ID.x", "Pos" = "Pos"))
  std_m_su <- group_by(std_su, Sample, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              p96 = unique(p96),
              CopyNum = mean(CopyNum)) %>%
    select(-p96) %>%
    separate(col = Sample, into = c("Sample", "Source"))


  ###############################################
  ## Copy Number Call
  ###############################################
  std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)

  data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
  #         data_m_su = cbind(data_m_su, predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean), interval = "prediction")[, c(2,3)])
  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)
  diff <- abs(data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round)
  data_m_su$diff <- diff



  output = list()
  output$data_su = data_su
  output$data_m_su = data_m_su

  output$std_su = std_su
  output$std_m_su = std_m_su
  output
}




makeData_RIVUR_bp <- function(data_raw, keyMap_std = "default", zeroCopy = NULL, calib = NULL){
  ###############################################
  ## Data processing
  ###############################################
  # data_raw = rawData$data_raw

  # if(hasZero) {
  #     data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
  #     summarize(ddCt = Cp[1]-Cp[2], RelQ = 2^-ddCt, zeroCopy = zero[1])
  #   } else{
  #     data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
  #     summarize(ddCt = Cp[1]-Cp[2], RelQ = 2^-ddCt)
  #   }

  if(!is.null(calib)){
    data = tbl_df(data_raw)
    calib_df = filter(data, Pos %in% calib) %>%
      group_by(Pos, Sample.Name, plate_ID) %>%
      summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %>%
      group_by(plate_ID) %>%
      summarize(dCt.mean = mean(dCt), RelQ.mean = mean(RelQ))

    # calib_dCt = mean(calib_df$dCt)

    calibrator_dCt = as.character()

    calibrator_dCt["8000013509"] = calib_df$dCt.mean[1]
    calibrator_dCt["8000013510"] = calib_df$dCt.mean[2]
    calibrator_dCt["8000013511"] = calib_df$dCt.mean[3]
    calibrator_dCt["8000013512"] = calib_df$dCt.mean[4]
    calibrator_dCt["8000013513"] = calib_df$dCt.mean[5]

  } else {
    calibrator_dCt = as.character()
    calibrator_dCt["8000013509"] = 0
    calibrator_dCt["8000013510"] = 0
    calibrator_dCt["8000013511"] = 0
    calibrator_dCt["8000013512"] = 0
    calibrator_dCt["8000013513"] = 0
    calibrator_dCt["8000013513"] = 0
    calibrator_dCt["CORIELL1"] = 0
    calibrator_dCt["CORIELL2"] = 0
    }


  data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
    summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt)

  data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  data_su$RelQ = 2^-data_su$ddCt


#   data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
#     summarize(ddCt = Cp[1]-Cp[2]-calib_dCt, RelQ = 2^-ddCt)

#   Ctl <- sapply(LETTERS[1:8], function(x) paste0(x, 1:6))
#   Rivur <- sapply(LETTERS[1:8], function(x) paste0(x, 7:12))

  # creat old to new map for rivur10 and ctl/case plates
#   C <- as.vector(Ctl)
#   R <- as.vector(Rivur)
#   map_LtoR = data.frame(old = R, rivur10_loc = C)

  rivur09 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013509")
  rivur10 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013510")
  rivur11 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013511")
  rivur12 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013512")
  rivur13 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013513")
  rivur <- RIVURplate::locMap_rivur %>% filter(plate_ID != "8000013514")
  rivur$plate_ID = as.character(rivur$plate_ID)  # convert numeric to character to keep compatible types for left_join function(with data_su$plate_ID)

  data_su = left_join(data_su, rivur, by = c("Pos" = "p384", "plate_ID" = "plate_ID"))


  # data_su$Type = ifelse(data_su$p96 %in% Rivur, "RIVUR", "CTL")
  data_su$Type =  "RIVUR"

#   if(data_su$plate_ID.x == "50_50_A"){
#     data_su = dplyr::left_join(data_su, map_LtoR, by = c("p96" = "old"))
#   }else if(data_su$plate_ID.x == "50_50_B")
#   {data_su$rivur10_loc = ifelse(data_su$Type == "RIVUR", as.character(data_su$p96), NA)
#   }else {data_su$rivur10_loc = ifelse(data_su$Type == "RIVUR", as.character(data_su$p96), NA)}


#   data_su = dplyr::left_join(data_su, rivur10[,c("p96", "RUID")]%>% unique, by = c("rivur10_loc" = "p96"))  # get RUID number
#   # rivur10 %>% unique()

  if(!is.null(zeroCopy)){
    data_su = left_join(data_su, zeroCopy[,c("plate_ID.x","Pos.x", "zero")], by = c("plate_ID" = "plate_ID.x", "Pos" = "Pos.x"))
  }else(data_su$zero = NA)

  # ds <- data_su[!is.na(data_su$RUID),] %>% filter(plate_ID == "8000013510")



  data_m_su <- group_by(data_su, p96, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              RUID = RUID[1],
              Type = Type[1],
              p384 = Pos[1],
              zeroCopy = zero[1])

  # data_m_su$Type = ifelse(data_m_su$p96 %in% Rivur, "RIVUR", "CTL")


  ###############################################
  ## Make standards
  ###############################################
  #

  std_su = left_join(keyMap_std, data_su, by = c("plate_ID" = "plate_ID", "Pos" = "Pos"))
  std_m_su <- group_by(std_su, Sample, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              p96 = unique(p96),
              CopyNum = mean(CopyNum)) %>%
    select(-p96) %>%
    separate(col = Sample, into = c("Sample", "Source"))


  ###############################################
  ## Copy Number Call
  ###############################################
  std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)

  data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
  #         data_m_su = cbind(data_m_su, predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean), interval = "prediction")[, c(2,3)])
  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)
  diff <- abs(data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round)
  data_m_su$diff <- diff



  output = list()
  output$data_su = data_su
  output$data_m_su = data_m_su

  output$std_su = std_su
  output$std_m_su = std_m_su
  output
}
makeData_RIVUR <- function(data_raw, keyMap_std = "default", zeroCopy = NULL, calib = NULL, comment = NULL, mode = NULL){
  ###############################################
  ## Data processing
  ###############################################
#   rawData <- load_data(fileloc)
#   diag <- load_data(fileloc_diag)
#   diag <- do.call(rbind, diag$file_list)
#   correct <- function(target, ref){
#     plate_ID_targ = unique(target$plate_ID)
#     target_m <- merge(target[target$Type == "Target Unknown",], ref[toupper(ref$plate_ID) == toupper(plate_ID_targ),], by.x = "Sample.Name", by.y = "Name")
#     isZero <- is.na(target_m$Cp.y) & target_m$Cp.x != 0
#     transform(target_m, zero = ifelse(isZero, TRUE, FALSE))
#   }
# # #
#   zeroCopy <- lapply(rawData$file_list, function(x) correct(x, diag))
#   zeroCopy <- do.call(rbind, zeroCopy)
#
#   data_raw = rawData$data_raw

#   if(hasZero) {
#       data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
#       summarize(ddCt = Cp[1]-Cp[2], RelQ = 2^-ddCt, zeroCopy = zero[1])
#     } else{
#       data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
#       summarize(ddCt = Cp[1]-Cp[2], RelQ = 2^-ddCt)
#     }

  if(!is.null(calib)){
    data = tbl_df(data_raw)
    calib_df = filter(data, Pos %in% calib) %>%
      group_by(Pos, Sample.Name, plate_ID) %>%
      summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %>%
      group_by(plate_ID) %>%
      summarize(dCt.mean = mean(dCt), RelQ.mean = mean(RelQ))

    # calib_dCt = mean(calib_df$dCt)

    calibrator_dCt = as.character()

    calibrator_dCt["8000013509"] = calib_df$dCt.mean[calib_df$plate_ID == "8000013509"]
    calibrator_dCt["8000013510"] = calib_df$dCt.mean[calib_df$plate_ID == "8000013510"]
    calibrator_dCt["8000013511"] = calib_df$dCt.mean[calib_df$plate_ID == "8000013511"]
    calibrator_dCt["8000013512"] = calib_df$dCt.mean[calib_df$plate_ID == "8000013512"]
    calibrator_dCt["8000013513"] = calib_df$dCt.mean[calib_df$plate_ID == "8000013513"]

  } else {
    calibrator_dCt = as.character()
    calibrator_dCt["8000013509"] = 0
    calibrator_dCt["8000013510"] = 0
    calibrator_dCt["8000013511"] = 0
    calibrator_dCt["8000013512"] = 0
    calibrator_dCt["8000013513"] = 0
    calibrator_dCt["CORIELL1"] = 0
    calibrator_dCt["CORIELL2"] = 0
    }


  data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
    summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt)

  data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  data_su$RelQ = 2^-data_su$ddCt

  # group_by(data_raw, Pos, Sample.Name, plate_ID) %>% mutate(n = n())
  #


#   data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
#     summarize(ddCt = Cp[1]-Cp[2]-calib_dCt, RelQ = 2^-ddCt)

#   Ctl <- sapply(LETTERS[1:8], function(x) paste0(x, 1:6))
#   Rivur <- sapply(LETTERS[1:8], function(x) paste0(x, 7:12))

  # creat old to new map for rivur10 and ctl/case plates
#   C <- as.vector(Ctl)
#   R <- as.vector(Rivur)
#   map_LtoR = data.frame(old = R, rivur10_loc = C)

  rivur09 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013509")
  rivur10 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013510")
  rivur11 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013511")
  rivur12 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013512")
  rivur13 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013513")

  rivur <- RIVURplate::locMap_rivur %>% filter(plate_ID != "8000013514")
  rivur$plate_ID = as.character(rivur$plate_ID)  # convert numeric to character to keep compatible types for left_join function(with data_su$plate_ID)

  data_su = left_join(data_su, rivur, by = c("Pos" = "p384", "plate_ID" = "plate_ID"))


  # data_su$Type = ifelse(data_su$p96 %in% Rivur, "RIVUR", "CTL")
  data_su$Type =  "RIVUR"

#   if(data_su$plate_ID.x == "50_50_A"){
#     data_su = dplyr::left_join(data_su, map_LtoR, by = c("p96" = "old"))
#   }else if(data_su$plate_ID.x == "50_50_B")
#   {data_su$rivur10_loc = ifelse(data_su$Type == "RIVUR", as.character(data_su$p96), NA)
#   }else {data_su$rivur10_loc = ifelse(data_su$Type == "RIVUR", as.character(data_su$p96), NA)}


#   data_su = dplyr::left_join(data_su, rivur10[,c("p96", "RUID")]%>% unique, by = c("rivur10_loc" = "p96"))  # get RUID number
#   # rivur10 %>% unique()

  if(!is.null(zeroCopy)){
    data_su = left_join(data_su, zeroCopy[,c("plate_ID.x","Pos.x", "zero")], by = c("plate_ID" = "plate_ID.x", "Pos" = "Pos.x"))
  }else(data_su$zero = NA)

  # ds <- data_su[!is.na(data_su$RUID),] %>% filter(plate_ID == "8000013510")



  data_m_su <- group_by(data_su, p96, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              RUID = RUID[1],
              Type = Type[1],
              p384 = Pos[1],
              zeroCopy = zero[1])

  # data_m_su$Type = ifelse(data_m_su$p96 %in% Rivur, "RIVUR", "CTL")


  ###############################################
  ## Make standards
  ###############################################
  #

  std_su = left_join(keyMap_std, data_su, by = c("plate_ID" = "plate_ID", "Pos" = "Pos"))
  std_m_su <- group_by(std_su, Sample, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              p96 = unique(p96),
              CopyNum = mean(CopyNum)) %>%
    select(-p96) %>%
    separate(col = Sample, into = c("Sample", "Source"))
  std_m_su = std_m_su[!is.na(std_m_su$ddCt.sem),]


  ###############################################
  ## Copy Number Call
  ###############################################
#   std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)
#
#   data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
#   #         data_m_su = cbind(data_m_su, predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean), interval = "prediction")[, c(2,3)])

  plateID = std_m_su$plate_ID %>% unique

  fit = lapply(seq_along(plateID), function(p) lm(CopyNum ~ RelQ.mean, data = std_m_su[std_m_su$plate_ID == plateID[p],]))

  # std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)

  d = list()
  for(i in 1:length(fit)){
    d[[i]] = data_m_su[data_m_su$plate_ID == plateID[i],]
    d[[i]]$Copy.Nbs.esti = predict(fit[[i]], newdata = data.frame(RelQ.mean = d[[i]]$RelQ.mean))
  }

  data_m_su = rbind_all(d)

  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)

  #  codes specific for nrxn3  -----------------
  # data_m_su$Copy.Nbs.round = NA
  if(!is.null(comment)){
    data_m_su$Copy.Nbs.round[data_m_su$plate_ID == "8000013510"] = ifelse(data_m_su$Copy.Nbs.esti <= 1.15, 1, 2)
    data_m_su$Copy.Nbs.round[data_m_su$plate_ID == "8000013511"] = ifelse(data_m_su$Copy.Nbs.esti <= 1.45, 1, 2)
    data_m_su$Copy.Nbs.round[data_m_su$plate_ID == "8000013512"] = ifelse(data_m_su$Copy.Nbs.esti <= 1.4, 1, 2)
  }
  #  codes specific for nrxn3  -----------------


  #  codes specific for other genes  -----------------
  # data_m_su$Copy.Nbs.round = NA
#   if(comment == "preassumed"){
# #     data_m_su$Copy.Nbs.round[data_m_su$plate_ID == "8000013510"] = ifelse(data_m_su$Copy.Nbs.esti <= 1.15, 1, 2)
# #     data_m_su$Copy.Nbs.round[data_m_su$plate_ID == "8000013511"] = ifelse(data_m_su$Copy.Nbs.esti <= 1.45, 1, 2)
#     data_m_su$Copy.Nbs.round[data_m_su$plate_ID == "8000013512"] = ifelse(data_m_su$Copy.Nbs.esti <= 0.6, 1, 2)
#   }
  #  codes specific for other genes  -----------------



  if(mode == "simple"){
    data_m_su$Copy.Nbs.round = ifelse(data_m_su$RelQ.mean <= 0.49, 1, 2)
  }


  if(!is.na(data_m_su$zeroCopy[1])){
    data_m_su$Copy.Nbs.round = ifelse(data_m_su$zeroCopy, 0, data_m_su$Copy.Nbs.round)
    data_m_su$Copy.Nbs.esti = ifelse(data_m_su$zeroCopy, NA, data_m_su$Copy.Nbs.esti)
  }


  diff <- data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round
  data_m_su$diff <- diff



  output = list()
  output$data_su = data_su
  output$data_m_su = data_m_su

  output$std_su = std_su
  output$std_m_su = std_m_su
  output
}



makeData_dCt <- function(data_raw, keyMap_std = "default", zeroCopy = NULL, calib = NULL){
  ###############################################
  ## Data processing
  ###############################################
  # data_raw = rawData$data_raw

  if(!is.null(calib)){
    data = tbl_df(data_raw)
    calib_df = filter(data, Pos %in% calib) %>%
      group_by(Pos, Sample.Name, plate_ID) %>%
      summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %>%
      group_by(plate_ID) %>%
      summarize(dCt.mean = mean(dCt), RelQ.mean = mean(RelQ), sd = sd(RelQ))

    # calib_dCt = mean(calib_df$dCt)

    calibrator_dCt = as.character()

    calibrator_dCt["CORIELL1"] = calib_df[calib_df$plate_ID == "CORIELL1",][,2]
    calibrator_dCt["CORIELL2"] = calib_df[calib_df$plate_ID == "CORIELL2",][,2]
    calibrator_dCt["CORIELL3"] = calib_df[calib_df$plate_ID == "CORIELL3",][,2]
    calibrator_dCt["CORIELL4"] = calib_df[calib_df$plate_ID == "CORIELL4",][,2]
    calibrator_dCt["CORIELL5"] = calib_df[calib_df$plate_ID == "CORIELL5",][,2]

    calibrator_dCt["8000013509"] = calib_df$dCt.mean[calib_df$plate_ID == "8000013509"]
    calibrator_dCt["8000013510"] = ifelse(length(calib_df$dCt.mean[calib_df$plate_ID == "8000013510"]) == 0, NA, calib_df$dCt.mean[calib_df$plate_ID == "8000013510"])
    calibrator_dCt["8000013511"] = calib_df$dCt.mean[calib_df$plate_ID == "8000013511"]
    calibrator_dCt["8000013512"] = calib_df$dCt.mean[calib_df$plate_ID == "8000013512"]
    calibrator_dCt["8000013513"] = calib_df$dCt.mean[calib_df$plate_ID == "8000013513"]

  } else {
    calibrator_dCt = as.character()
    calibrator_dCt["CORIELL1"] = 0
    calibrator_dCt["CORIELL2"] = 0
    calibrator_dCt["CORIELL3"] = 0
    calibrator_dCt["CORIELL4"] = 0
    calibrator_dCt["CORIELL5"] = 0

    calibrator_dCt["8000013509"] = 0
    calibrator_dCt["8000013510"] = 0
    calibrator_dCt["8000013511"] = 0
    calibrator_dCt["8000013512"] = 0
    calibrator_dCt["8000013513"] = 0
  }


  data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
    summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt)

  # data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  # data_su$ddCt = data_su$dCt - 0
  data_su$RelQ = 2^-data_su$ddCt


  data_su$Type =  "RIVUR"

  if(!is.null(zeroCopy)){
    data_su = left_join(data_su, zeroCopy[,c("plate_ID.x","Pos.x", "zero")], by = c("plate_ID" = "plate_ID.x", "Pos" = "Pos.x"))
  }else(data_su$zero = NA)


  # rivur <- RIVURplate::locMap_rivur %>% filter(plate_ID != "8000013514")

#   rivur09 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013509")
#   rivur10 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013510")
#   rivur11 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013511")
#   rivur12 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013512")
#   rivur13 <- RIVURplate::locMap_rivur %>% filter(plate_ID == "8000013513")
  rivur <- RIVURplate::locMap_rivur %>% filter(plate_ID != "8000013514")
  rivur$plate_ID = as.character(rivur$plate_ID)  # convert numeric to character to keep compatible types for left_join function(with data_su$plate_ID)

  # data_su = left_join(data_su, map_p96_p384, by = c("Pos" = "p384"))

  data_su = left_join(data_su, rivur, by = c("Pos" = "p384", "plate_ID" = "plate_ID"))

#   data_su %>% group_by(plate_ID, is.na(RUID)) %>% summarise(n())
#   data_m_su %>% group_by(plate_ID, is.na(RUID)) %>% summarise(n())

  # data_su %>% filter(Pos %in% c("B1", "B2", "B3"))


  data_m_su <- group_by(data_su, p96, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              Type = Type[1],
              p384 = Pos[1],
              RUID = RUID[1],
              zeroCopy = zero[1])

  # data_m_su$Type = ifelse(data_m_su$p96 %in% Rivur, "RIVUR", "CTL")


  ###############################################
  ## Make standards
  ###############################################
  #

  std_su = left_join(keyMap_std, data_su, by = c("plate_ID" = "plate_ID", "Pos" = "Pos"))
  std_m_su <- group_by(std_su, Sample, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              p96 = unique(p96),
              CopyNum = mean(CopyNum)) %>%
    select(-p96) %>%
    separate(col = Sample, into = c("Sample", "Source"))


  ###############################################
  ## Copy Number Call
  ###############################################

  std_copyNbs = keyMap_std[keyMap_std$Pos %in% calib, ]$CopyNum[1]
  data_m_su$Copy.Nbs.esti = std_copyNbs * data_m_su$RelQ.mean
  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)

  data_m_su$Copy.Nbs.round = ifelse(data_m_su$zeroCopy, 0, data_m_su$Copy.Nbs.round)
  data_m_su$Copy.Nbs.esti = ifelse(data_m_su$zeroCopy, NA, data_m_su$Copy.Nbs.esti)

  diff <- abs(data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round)
  data_m_su$diff <- diff



  output = list()
  output$data_su = data_su
  output$data_m_su = data_m_su

  output$std_su = std_su
  output$std_m_su = std_m_su
  output
}


makeData_bp <- function(data_raw, keyMap_std = "default", zeroCopy = NULL, calib = NULL){
  ###############################################
  ## Data processing
  ###############################################
  # data_raw = rawData$data_raw

  if(!is.null(calib)){
    data = tbl_df(data_raw)
    calib_df = filter(data, Pos %in% calib) %>%
      group_by(Pos, Sample.Name, plate_ID) %>%
      summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %>%
      group_by(plate_ID) %>%
      summarize(dCt.mean = mean(dCt), RelQ.mean = mean(RelQ), sd = sd(RelQ))

    # calib_dCt = mean(calib_df$dCt)

    calibrator_dCt = as.character()

    calibrator_dCt["CORIELL1"] = calib_df[calib_df$plate_ID == "CORIELL1",][,2]
    calibrator_dCt["CORIELL2"] = calib_df[calib_df$plate_ID == "CORIELL2",][,2]
    calibrator_dCt["CORIELL3"] = calib_df[calib_df$plate_ID == "CORIELL3",][,2]
    calibrator_dCt["CORIELL4"] = calib_df[calib_df$plate_ID == "CORIELL4",][,2]
    calibrator_dCt["CORIELL5"] = calib_df[calib_df$plate_ID == "CORIELL5",][,2]
    calibrator_dCt["CORIELL6"] = calib_df[calib_df$plate_ID == "CORIELL6",][,2]

    calibrator_dCt["8000182997"] = calib_df[calib_df$plate_ID == "8000182997",][,2]
    calibrator_dCt["8000182998"] = calib_df[calib_df$plate_ID == "8000182998",][,2]

  } else {
    calibrator_dCt = as.character()
    calibrator_dCt["CORIELL1"] = 0
    calibrator_dCt["CORIELL2"] = 0
    calibrator_dCt["CORIELL3"] = 0
    calibrator_dCt["CORIELL4"] = 0
    calibrator_dCt["CORIELL5"] = 0
    calibrator_dCt["CORIELL6"] = 0
    calibrator_dCt["8000182997"] = 0
    calibrator_dCt["8000182998"] = 0
  }


  data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
    summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt)

  # data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  # data_su$ddCt = data_su$dCt - 0
  data_su$RelQ = 2^-data_su$ddCt

  data_su$plate_ID = ifelse(data_su$plate_ID == "CUTIE97", "8000182997",
                            ifelse(data_su$plate_ID == "CUTIE98", "8000182998", data_su$plate_ID))

  cutie = locMap_cutie

  data_su = left_join(data_su, cutie, by = c("Pos" = "p384", "plate_ID" = "Plate_Box_Inv_Code"))

  data_su$Type = "Plate"

  if(!is.null(zeroCopy)){
    data_su = left_join(data_su, zeroCopy[,c("plate_ID.x","Pos.x", "zero")], by = c("plate_ID" = "plate_ID.x", "Pos" = "Pos.x"))
  }else(data_su$zero = NA)


  # rivur <- RIVURplate::locMap_rivur %>% filter(plate_ID != "8000013514")

  data_su = left_join(data_su, map_p96_p384, by = c("Pos" = "p384"))

  # data_su %>% filter(Pos %in% c("B1", "B2", "B3"))


  data_m_su <- group_by(data_su, p96, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              Type = Type[1],
              p384 = Pos[1],
              zeroCopy = zero[1])

  # data_m_su$Type = ifelse(data_m_su$p96 %in% Rivur, "RIVUR", "CTL")


  ###############################################
  ## Make standards
  ###############################################
  #

  std_su = left_join(keyMap_std, data_su, by = c("plate_ID" = "plate_ID", "Pos" = "Pos"))
  std_m_su <- group_by(std_su, Sample, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              p96 = unique(p96),
              CopyNum = mean(CopyNum)) %>%
    select(-p96) %>%
    separate(col = Sample, into = c("Sample", "Source"))


  ###############################################
  ## Copy Number Call
  ###############################################
  std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)

  data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
  #         data_m_su = cbind(data_m_su, predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean), interval = "prediction")[, c(2,3)])
  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)
  diff <- abs(data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round)
  data_m_su$diff <- diff








  output = list()
  output$data_su = data_su
  output$data_m_su = data_m_su

  output$std_su = std_su
  output$std_m_su = std_m_su
  output
}


makeData_bp2 <- function(data_raw, keyMap_std = "default", zeroCopy = NULL, calib = NULL){
  ###############################################
  ## Data processing
  ###############################################
  # data_raw = rawData$data_raw

  if(!is.null(calib)){
    data = tbl_df(data_raw)
    calib_df = filter(data, Pos %in% calib) %>%
      group_by(Pos, Sample.Name, plate_ID) %>%
      summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %>%
      group_by(plate_ID) %>%
      summarize(dCt.mean = mean(dCt), RelQ.mean = mean(RelQ), sd = sd(RelQ))

    # calib_dCt = mean(calib_df$dCt)

    calibrator_dCt = as.character()

    calibrator_dCt["CORIELL1"] = calib_df[calib_df$plate_ID == "CORIELL1",][,2]
    calibrator_dCt["CORIELL2"] = calib_df[calib_df$plate_ID == "CORIELL2",][,2]
    calibrator_dCt["CORIELL3"] = calib_df[calib_df$plate_ID == "CORIELL3",][,2]
    calibrator_dCt["CORIELL4"] = calib_df[calib_df$plate_ID == "CORIELL4",][,2]
    calibrator_dCt["CORIELL5"] = calib_df[calib_df$plate_ID == "CORIELL5",][,2]
    calibrator_dCt["CORIELL6"] = calib_df[calib_df$plate_ID == "CORIELL6",][,2]

    calibrator_dCt["8000182997"] = calib_df[calib_df$plate_ID == "8000182997",][,2]
    calibrator_dCt["8000182998"] = calib_df[calib_df$plate_ID == "8000182998",][,2]

  } else {
    calibrator_dCt = as.character()
    calibrator_dCt["CORIELL1"] = 0
    calibrator_dCt["CORIELL2"] = 0
    calibrator_dCt["CORIELL3"] = 0
    calibrator_dCt["CORIELL4"] = 0
    calibrator_dCt["CORIELL5"] = 0
    calibrator_dCt["CORIELL6"] = 0
    calibrator_dCt["8000182997"] = 0
    calibrator_dCt["8000182998"] = 0
  }


  data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
    summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt)

  # data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  # data_su$ddCt = data_su$dCt - 0
  data_su$RelQ = 2^-data_su$ddCt

  data_su$plate_ID = ifelse(data_su$plate_ID == "CUTIE97", "8000182997",
                            ifelse(data_su$plate_ID == "CUTIE98", "8000182998", data_su$plate_ID))

  # cutie = locMap_cutie

  data_su = left_join(data_su, locMap_cutie, by = c("Pos" = "p384", "plate_ID" = "Plate_Box_Inv_Code"))

  data_su$Type = "Plate"

  if(!is.null(zeroCopy)){
    data_su = left_join(data_su, zeroCopy[,c("plate_ID.x","Pos.x", "zero")], by = c("plate_ID" = "plate_ID.x", "Pos" = "Pos.x"))
  }else(data_su$zero = NA)


  # rivur <- RIVURplate::locMap_rivur %>% filter(plate_ID != "8000013514")

  # data_su = left_join(data_su, map_p96_p384, by = c("Pos" = "p384"))

  # data_su %>% filter(Pos %in% c("B1", "B2", "B3"))


  data_m_su <- group_by(data_su, p96, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              Type = Type[1],
              p384 = Pos[1],
              RUID = RUID[1],
              zeroCopy = zero[1])

  # data_m_su$Type = ifelse(data_m_su$p96 %in% Rivur, "RIVUR", "CTL")


  ###############################################
  ## Make standards
  ###############################################
  #

  std_su = left_join(keyMap_std, data_su, by = c("plate_ID" = "plate_ID", "Pos" = "Pos"))
  std_m_su <- group_by(std_su, Sample, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              p96 = unique(p96),
              CopyNum = mean(CopyNum)) %>%
    select(-p96) %>%
    separate(col = Sample, into = c("Sample", "Source"))


  ###############################################
  ## Copy Number Call
  ###############################################
  std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)

  data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
  #         data_m_su = cbind(data_m_su, predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean), interval = "prediction")[, c(2,3)])
  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)


#   std_copyNbs = keyMap_std[keyMap_std$Pos %in% calib, ]$CopyNum[1]
#   data_m_su$Copy.Nbs.esti = std_copyNbs * data_m_su$RelQ.mean
#   data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)

  diff <- abs(data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round)
  data_m_su$diff <- diff








  output = list()
  output$data_su = data_su
  output$data_m_su = data_m_su

  output$std_su = std_su
  output$std_m_su = std_m_su
  output
}


makeData_bp3 <- function(data_raw, keyMap_std = "default", zeroCopy = NULL, calib = NULL){
  ###############################################
  ## Data processing
  ###############################################
  # data_raw = rawData$data_raw

  if(!is.null(calib)){
    data = tbl_df(data_raw)
    calib_df = filter(data, Pos %in% calib) %>%
      group_by(Pos, Sample.Name, plate_ID) %>%
      summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %>%
      group_by(plate_ID) %>%
      summarize(dCt.mean = mean(dCt), RelQ.mean = mean(RelQ), sd = sd(RelQ))

    # calib_dCt = mean(calib_df$dCt)

    calibrator_dCt = as.character()

    calibrator_dCt["CORIELL1"] = calib_df[calib_df$plate_ID == "CORIELL1",][,2]
    calibrator_dCt["CORIELL2"] = calib_df[calib_df$plate_ID == "CORIELL2",][,2]
    calibrator_dCt["CORIELL3"] = calib_df[calib_df$plate_ID == "CORIELL3",][,2]
    calibrator_dCt["CORIELL4"] = calib_df[calib_df$plate_ID == "CORIELL4",][,2]
    calibrator_dCt["CORIELL5"] = calib_df[calib_df$plate_ID == "CORIELL5",][,2]
    calibrator_dCt["CORIELL6"] = calib_df[calib_df$plate_ID == "CORIELL6",][,2]

    calibrator_dCt["8000182997"] = calib_df[calib_df$plate_ID == "8000182997",][,2]
    calibrator_dCt["8000182998"] = calib_df[calib_df$plate_ID == "8000182998",][,2]

  } else {
    calibrator_dCt = as.character()
    calibrator_dCt["CORIELL1"] = 0
    calibrator_dCt["CORIELL2"] = 0
    calibrator_dCt["CORIELL3"] = 0
    calibrator_dCt["CORIELL4"] = 0
    calibrator_dCt["CORIELL5"] = 0
    calibrator_dCt["CORIELL6"] = 0
    calibrator_dCt["8000182997"] = 0
    calibrator_dCt["8000182998"] = 0
  }


  data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
    summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt)

  # data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  # data_su$ddCt = data_su$dCt - 0
  data_su$RelQ = 2^-data_su$ddCt

  data_su$plate_ID = ifelse(data_su$plate_ID == "CUTIE97", "8000182997",
                            ifelse(data_su$plate_ID == "CUTIE98", "8000182998", data_su$plate_ID))

  # cutie = locMap_cutie

  data_su = left_join(data_su, locMap_cutie, by = c("Pos" = "p384", "plate_ID" = "Plate_Box_Inv_Code"))

  data_su$Type = "Plate"

  if(!is.null(zeroCopy)){
    data_su = left_join(data_su, zeroCopy[,c("plate_ID.x","Pos.x", "zero")], by = c("plate_ID" = "plate_ID.x", "Pos" = "Pos.x"))
  }else(data_su$zero = NA)


  # rivur <- RIVURplate::locMap_rivur %>% filter(plate_ID != "8000013514")

  # data_su = left_join(data_su, map_p96_p384, by = c("Pos" = "p384"))

  # data_su %>% filter(Pos %in% c("B1", "B2", "B3"))


  data_m_su <- group_by(data_su, p96, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              Type = Type[1],
              p384 = Pos[1],
              RUID = RUID[1],
              zeroCopy = zero[1])

  # data_m_su$Type = ifelse(data_m_su$p96 %in% Rivur, "RIVUR", "CTL")


  ###############################################
  ## Make standards
  ###############################################
  #

  std_su = left_join(keyMap_std, data_su, by = c("plate_ID" = "plate_ID", "Pos" = "Pos"))
  std_m_su <- group_by(std_su, Sample, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              p96 = unique(p96),
              CopyNum = mean(CopyNum)) %>%
    select(-p96) %>%
    separate(col = Sample, into = c("Sample", "Source"))


  ###############################################
  ## Copy Number Call
  ###############################################
  std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)

  data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
  #         data_m_su = cbind(data_m_su, predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean), interval = "prediction")[, c(2,3)])
  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)


  #   std_copyNbs = keyMap_std[keyMap_std$Pos %in% calib, ]$CopyNum[1]
  #   data_m_su$Copy.Nbs.esti = std_copyNbs * data_m_su$RelQ.mean
  #   data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)

  diff <- abs(data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round)
  data_m_su$diff <- diff


  output = list()
  output$data_su = data_su
  output$data_m_su = data_m_su

  output$std_su = std_su
  output$std_m_su = std_m_su
  output
}


makeData_cutie <- function(data_raw, keyMap_std = "default", zeroCopy = NULL, calib = NULL){
  ###############################################
  ## Data processing
  ###############################################
  # data_raw = rawData$data_raw

  if(!is.null(calib)){
    data = tbl_df(data_raw)
    calib_df = filter(data, Pos %in% calib) %>%
      group_by(Pos, Sample.Name, plate_ID) %>%
      summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %>%
      group_by(plate_ID) %>%
      summarize(dCt.mean = mean(dCt), RelQ.mean = mean(RelQ), sd = sd(RelQ))

    # calib_dCt = mean(calib_df$dCt)

    calibrator_dCt = as.character()

    calibrator_dCt["CORIELL1"] = calib_df[calib_df$plate_ID == "CORIELL1",][,2]
    calibrator_dCt["CORIELL2"] = calib_df[calib_df$plate_ID == "CORIELL2",][,2]
    calibrator_dCt["CORIELL3"] = calib_df[calib_df$plate_ID == "CORIELL3",][,2]
    calibrator_dCt["CORIELL4"] = calib_df[calib_df$plate_ID == "CORIELL4",][,2]
    calibrator_dCt["CORIELL5"] = calib_df[calib_df$plate_ID == "CORIELL5",][,2]
    calibrator_dCt["CORIELL6"] = calib_df[calib_df$plate_ID == "CORIELL6",][,2]

    calibrator_dCt["8000182997"] = calib_df[calib_df$plate_ID == "8000182997",][,2]
    calibrator_dCt["8000182998"] = calib_df[calib_df$plate_ID == "8000182998",][,2]

  } else {
    calibrator_dCt = as.character()
    calibrator_dCt["CORIELL1"] = 0
    calibrator_dCt["CORIELL2"] = 0
    calibrator_dCt["CORIELL3"] = 0
    calibrator_dCt["CORIELL4"] = 0
    calibrator_dCt["CORIELL5"] = 0
    calibrator_dCt["CORIELL6"] = 0
    calibrator_dCt["8000182997"] = 0
    calibrator_dCt["8000182998"] = 0
  }


  data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
    summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt)

  # data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  # data_su$ddCt = data_su$dCt - 0
  data_su$RelQ = 2^-data_su$ddCt

  data_su$plate_ID = ifelse(data_su$plate_ID == "CUTIE97", "8000182997",
                            ifelse(data_su$plate_ID == "CUTIE98", "8000182998", data_su$plate_ID))

  # cutie = locMap_cutie

  data_su = left_join(data_su, locMap_cutie, by = c("Pos" = "p384", "plate_ID" = "Plate_Box_Inv_Code"))

  data_su$Type = "CUTIE"

  if(!is.null(zeroCopy)){
    data_su = left_join(data_su, zeroCopy[,c("plate_ID.x","Pos.x", "zero")], by = c("plate_ID" = "plate_ID.x", "Pos" = "Pos.x"))
  }else(data_su$zero = NA)


  # rivur <- RIVURplate::locMap_rivur %>% filter(plate_ID != "8000013514")

  # data_su = left_join(data_su, map_p96_p384, by = c("Pos" = "p384"))

  # data_su %>% filter(Pos %in% c("B1", "B2", "B3"))


  data_m_su <- group_by(data_su, p96, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              Type = Type[1],
              p384 = Pos[1],
              RUID = RUID[1],
              zeroCopy = zero[1])

  # data_m_su$Type = ifelse(data_m_su$p96 %in% Rivur, "RIVUR", "CTL")


  ###############################################
  ## Make standards
  ###############################################
  #

  std_su = left_join(keyMap_std, data_su, by = c("plate_ID" = "plate_ID", "Pos" = "Pos"))
  std_m_su <- group_by(std_su, Sample, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              p96 = unique(p96),
              CopyNum = mean(CopyNum)) %>%
    select(-p96) %>%
    separate(col = Sample, into = c("Sample", "Source"))


  ###############################################
  ## Copy Number Call
  ###############################################
#   std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)
#
#   data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
  #         data_m_su = cbind(data_m_su, predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean), interval = "prediction")[, c(2,3)])


  plateID = std_m_su$plate_ID %>% unique

  fit = lapply(seq_along(plateID), function(p) lm(CopyNum ~ RelQ.mean, data = std_m_su[std_m_su$plate_ID == plateID[p],]))

  # std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)

  d = list()
  for(i in 1:length(fit)){
    d[[i]] = data_m_su[data_m_su$plate_ID == plateID[i],]
    d[[i]]$Copy.Nbs.esti = predict(fit[[i]], newdata = data.frame(RelQ.mean = d[[i]]$RelQ.mean))
  }

  data_m_su = rbind_all(d)



  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)
  data_m_su$Copy.Nbs.round = ifelse(data_m_su$zeroCopy, 0, data_m_su$Copy.Nbs.round)
  data_m_su$Copy.Nbs.esti = ifelse(data_m_su$zeroCopy, NA, data_m_su$Copy.Nbs.esti)

  #   std_copyNbs = keyMap_std[keyMap_std$Pos %in% calib, ]$CopyNum[1]
  #   data_m_su$Copy.Nbs.esti = std_copyNbs * data_m_su$RelQ.mean
  #   data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)

  diff <- data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round
  data_m_su$diff <- diff


  output = list()
  output$data_su = data_su
  output$data_m_su = data_m_su

  output$std_su = std_su
  output$std_m_su = std_m_su
  output
}


makeData_Coriell <- function(data_raw, keyMap_std = "default", zeroCopy = NULL, calib = NULL, mode = NULL){
  ###############################################
  ## Data processing
  ###############################################
  # data_raw = rawData$data_raw
#
#     rawData <- load_data(fileloc)
#     data_raw = rawData$data_raw
#     diag <- load_data(fileloc_diag)
#     diag <- do.call(rbind, diag$file_list)
#     correct <- function(target, ref){
#       plate_ID_targ = unique(target$plate_ID)
#       target_m <- merge(target[target$Type == "Target Unknown",], ref[toupper(ref$plate_ID) == toupper(plate_ID_targ),], by.x = "Sample.Name", by.y = "Name")
#       isZero <- is.na(target_m$Cp.y) & target_m$Cp.x != 0
#       transform(target_m, zero = ifelse(isZero, TRUE, FALSE))
#     }
#   #
#     zeroCopy <- lapply(rawData$file_list, function(x) correct(x, diag))
#     zeroCopy <- do.call(rbind, zeroCopy)
#     keyMap_std = stdInfo





  if(!is.null(calib)){
    data = tbl_df(data_raw)
    calib_df = filter(data, Pos %in% calib) %>%
      group_by(Pos, Sample.Name, plate_ID) %>%
      summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %>%
      group_by(plate_ID) %>%
      summarize(dCt.mean = mean(dCt), RelQ.mean = mean(RelQ), sd = sd(RelQ))

    # calib_dCt = mean(calib_df$dCt)
    calib_df$plate_ID = toupper(calib_df$plate_ID)
    calibrator_dCt = as.character()

    calibrator_dCt["CORIELL1"] = calib_df[calib_df$plate_ID == "CORIELL1",][,2]
    calibrator_dCt["CORIELL2"] = calib_df[calib_df$plate_ID == "CORIELL2",][,2]
    calibrator_dCt["CORIELL3"] = calib_df[calib_df$plate_ID == "CORIELL3",][,2]
    calibrator_dCt["CORIELL4"] = calib_df[calib_df$plate_ID == "CORIELL4",][,2]
    calibrator_dCt["CORIELL5"] = calib_df[calib_df$plate_ID == "CORIELL5",][,2]
    calibrator_dCt["CORIELL6"] = calib_df[calib_df$plate_ID == "CORIELL6",][,2]

    calibrator_dCt["8000182997"] = calib_df[calib_df$plate_ID == "8000182997",][,2]
    calibrator_dCt["8000182998"] = calib_df[calib_df$plate_ID == "8000182998",][,2]

  } else {
    calibrator_dCt = as.character()
    calibrator_dCt["CORIELL1"] = 0
    calibrator_dCt["CORIELL2"] = 0
    calibrator_dCt["CORIELL3"] = 0
    calibrator_dCt["CORIELL4"] = 0
    calibrator_dCt["CORIELL5"] = 0
    calibrator_dCt["CORIELL6"] = 0
    calibrator_dCt["8000182997"] = 0
    calibrator_dCt["8000182998"] = 0
  }


  data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
    summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt)

  # data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID%>%toupper])
  # data_su$ddCt = data_su$dCt - 0
  data_su$RelQ = 2^-data_su$ddCt

#   data_su$plate_ID = ifelse(data_su$plate_ID == "CUTIE97", "8000182997",
#                             ifelse(data_su$plate_ID == "CUTIE98", "8000182998", data_su$plate_ID))

  # cutie = locMap_cutie

  # data_su = left_join(data_su, locMap_cutie, by = c("Pos" = "p384", "plate_ID" = "Plate_Box_Inv_Code"))

  data_su$Type = "CORIELL"

  if(!is.null(zeroCopy)){
    data_su = left_join(data_su, zeroCopy[,c("plate_ID.x","Pos.x", "zero")], by = c("plate_ID" = "plate_ID.x", "Pos" = "Pos.x"))
  }else(data_su$zero = NA)


  # rivur <- RIVURplate::locMap_rivur %>% filter(plate_ID != "8000013514")

  data_su = left_join(data_su, map_p96_p384, by = c("Pos" = "p384"))

  # data_su %>% filter(Pos %in% c("B1", "B2", "B3"))


  data_m_su <- group_by(data_su, p96, plate_ID) %>%
  # data_m_su <- group_by(data_su, Sample.Name, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              Type = Type[1],
              p384 = Pos[1],
              # RUID = RUID[1],
              zeroCopy = zero[1])

  # data_m_su$Type = ifelse(data_m_su$p96 %in% Rivur, "RIVUR", "CTL")


  ###############################################
  ## Make standards
  ###############################################
  #
  data_su$plate_ID = toupper((data_su$plate_ID))
  std_su = left_join(keyMap_std, data_su, by = c("plate_ID" = "plate_ID", "Pos" = "Pos"))
  std_m_su <- group_by(std_su, Sample, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              p96 = unique(p96),
              CopyNum = mean(CopyNum)) %>%
    select(-p96) %>%
    separate(col = Sample, into = c("Sample", "Source"))


  ###############################################
  ## Copy Number Call
  ###############################################

  std_m_su = std_m_su[!is.na(std_m_su$ddCt.mean),]
  plateID = std_m_su$plate_ID %>% unique

  fit = lapply(seq_along(plateID), function(p) lm(CopyNum ~ RelQ.mean, data = std_m_su[std_m_su$plate_ID == plateID[p],]))

  # std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)

  d = list()
  for(i in 1:length(fit)){
    d[[i]] = data_m_su[toupper(data_m_su$plate_ID) == plateID[i],]
    d[[i]]$Copy.Nbs.esti = predict(fit[[i]], newdata = data.frame(RelQ.mean = d[[i]]$RelQ.mean))
  }

  data_m_su = rbind_all(d)

  # data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
  #         data_m_su = cbind(data_m_su, (std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean), interval = "prediction")[, c(2,3)])
  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)
  data_m_su$Copy.Nbs.round = ifelse(data_m_su$zeroCopy, 0, data_m_su$Copy.Nbs.round)
  data_m_su$Copy.Nbs.esti = ifelse(data_m_su$zeroCopy, NA, data_m_su$Copy.Nbs.esti)


  #   std_copyNbs = keyMap_std[keyMap_std$Pos %in% calib, ]$CopyNum[1]
  #   data_m_su$Copy.Nbs.esti = std_copyNbs * data_m_su$RelQ.mean
  #   data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)

  diff <- abs(data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round)
  data_m_su$diff <- diff


  output = list()
  output$data_su = data_su
  output$data_m_su = data_m_su

  output$std_su = std_su
  output$std_m_su = std_m_su
  output
}



makeData_Coriell_test <- function(){
  add_isZeroCopy = function(target, diag){
    plate_ID_targ = unique(target$plate_ID)
    # target = target[Type == "Target Unknown"]
    diag_targ = diag[plate_ID == plate_ID_targ]
    setkey(diag_targ, Pos)
    target[Type == "Target Unknown", Cp_diag := diag_targ[Pos , Cp]]
    target[Type == "Target Unknown", is_ZeroCopy := is.na(Cp_diag) & Cp !=0]
    target

  }

  data_transfromation <- function(ds, diag = NULL){
    library(stringr)
    names(ds) = make.names(names(ds))

    # Label with sample type info
    std_loc = map(.x = list("B", "D", "F", "H", "J", "L", "N"), .f = ~paste0(.x, 13:15)) %>% unlist
    blank_loc = c(map(list("F", "H", "J", "L", "N", "P"), .f = ~paste0(.x, 13:24)) %>% unlist,
                  map(list("B", "D"), .f = ~paste0(.x, 16:24)) %>% unlist)
    sample_loc = c(map(list("A", "C", "E", "G", "I", "K", "M", "O"), .f = ~paste0(.x, 1:24)) %>% unlist,
                   map(list("B", "D", "F", "H", "J", "L", "N", "P"), .f = ~paste0(.x, 1:12)) %>% unlist)
    setkey(ds, Pos)
    ds[, Category := "Unknown"]
    ds[blank_loc, Category := 'blank']
    ds[sample_loc, Category := 'sample']
    ds[std_loc, Category := 'std']

    if(!is.null(diag)) {ds = add_isZeroCopy(ds, diag)} else ds$is_ZeroCopy = FALSE
    # add_isZeroCopy(dataset[[1]], diag)[] -> ds

    # data transform - calculate dCT, RelQ etc.
    ds[, .(dCt = .SD[Type == "Target Unknown", Cp] - .SD[Type == "Reference Unknown", Cp],
           RelQ = 2^-(.SD[Type == "Target Unknown", Cp] - .SD[Type == "Reference Unknown", Cp]),
           Target.Name = .SD$Target.Name[1],
           is_ZeroCopy = is_ZeroCopy[1],
           plate_ID = plate_ID[1],
           Category = Category[1]),
       by = .(Pos, Sample.Name)] %>%
      .[, Label := ((str_extract(Sample.Name, "[0-9]*$") %>% as.numeric) / 3) %>% ceiling()] %>%
      .[, .(dCt.mean = mean(dCt),
            RelQ.mean = mean(RelQ),
            sd.mean = sd(dCt),
            Pos = Pos[1],
            plate_ID = plate_ID[1],
            is_ZeroCopy = sum(is_ZeroCopy) >=2 ,  # Determine the # of NAs being equal to or greater than 2 as zero copy.
            Type = Category[1]),
        by = .(Label, Target.Name)] -> ds1

    # Add posistion info for 96 well plate
    Pos_96 = map(LETTERS[1:8], ~ paste0(.x, 1:12)) %>% unlist
    setkey(ds1, Pos)
    ds2 = ds1[mixedsort(Pos)][Type == 'sample', Pos_96 := Pos_96]
    ds2
  }

  add_plateID <- function(x){
    get_plateID <- function(x) {
      l = names(x)
      # plate_loc_start = gregexpr("CORIELLE[0-9]+-", l)[[1]][1] + 5
      # plate_loc_end = attr(gregexpr("CORIELLE[0-9]+-", l)[[1]], 'match.length') + plate_loc_start - 7
      # plate_ID = substr(l, plate_loc_start-5, plate_loc_end)
      # plate_ID
      str_split(l, "-") %>% map_chr(., ~ .[3])
    }

    plate_ID = get_plateID(x)

    map2(x, plate_ID, ~ .x[, plate_ID := .y])
  }

}




makeData <- function(data_raw, keyMap_std = "default", zeroCopy = NULL, calib = NULL){
  ###############################################
  ## Data processing
  ###############################################
  # data_raw = rawData$data_raw

  if(!is.null(calib)){
    data = tbl_df(data_raw)
    calib_df = filter(data, Pos %in% calib) %>%
      group_by(Pos, Sample.Name, plate_ID) %>%
      summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt) %>%
      group_by(plate_ID) %>%
      summarize(dCt.mean = mean(dCt), RelQ.mean = mean(RelQ), sd = sd(RelQ))

    # calib_dCt = mean(calib_df$dCt)

    calibrator_dCt = as.character()

    calibrator_dCt["CORIELL1"] = calib_df[calib_df$plate_ID == "CORIELL1",][,2]
    calibrator_dCt["CORIELL2"] = calib_df[calib_df$plate_ID == "CORIELL2",][,2]
    calibrator_dCt["CORIELL3"] = calib_df[calib_df$plate_ID == "CORIELL3",][,2]
    calibrator_dCt["CORIELL4"] = calib_df[calib_df$plate_ID == "CORIELL4",][,2]
    calibrator_dCt["CORIELL5"] = calib_df[calib_df$plate_ID == "CORIELL5",][,2]
    calibrator_dCt["CORIELL6"] = calib_df[calib_df$plate_ID == "CORIELL6",][,2]

    calibrator_dCt["8000182997"] = calib_df[calib_df$plate_ID == "8000182997",][,2]
    calibrator_dCt["8000182998"] = calib_df[calib_df$plate_ID == "8000182998",][,2]

  } else {
    calibrator_dCt = as.character()
    calibrator_dCt["CORIELL1"] = 0
    calibrator_dCt["CORIELL2"] = 0
    calibrator_dCt["CORIELL3"] = 0
    calibrator_dCt["CORIELL4"] = 0
    calibrator_dCt["CORIELL5"] = 0
    calibrator_dCt["CORIELL6"] = 0
    calibrator_dCt["8000182997"] = 0
    calibrator_dCt["8000182998"] = 0
  }


  data_su = group_by(data_raw, Pos, Sample.Name, plate_ID) %>%
    summarize(dCt = Cp[1]-Cp[2], RelQ = 2^-dCt)

  # data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  data_su$ddCt = data_su$dCt - as.numeric(calibrator_dCt[data_su$plate_ID])
  # data_su$ddCt = data_su$dCt - 0
  data_su$RelQ = 2^-data_su$ddCt

  data_su$plate_ID = ifelse(data_su$plate_ID == "CUTIE97", "8000182997",
                            ifelse(data_su$plate_ID == "CUTIE98", "8000182998", data_su$plate_ID))

  # cutie = locMap_cutie

  data_su = left_join(data_su, locMap_cutie, by = c("Pos" = "p384", "plate_ID" = "Plate_Box_Inv_Code"))

  data_su$Type = "Plate"

  if(!is.null(zeroCopy)){
    data_su = left_join(data_su, zeroCopy[,c("plate_ID.x","Pos.x", "zero")], by = c("plate_ID" = "plate_ID.x", "Pos" = "Pos.x"))
  }else(data_su$zero = NA)


  # rivur <- RIVURplate::locMap_rivur %>% filter(plate_ID != "8000013514")

  # data_su = left_join(data_su, map_p96_p384, by = c("Pos" = "p384"))

  # data_su %>% filter(Pos %in% c("B1", "B2", "B3"))


  data_m_su <- group_by(data_su, p96, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              Type = Type[1],
              p384 = Pos[1],
              RUID = RUID[1],
              zeroCopy = zero[1])

  # data_m_su$Type = ifelse(data_m_su$p96 %in% Rivur, "RIVUR", "CTL")


  ###############################################
  ## Make standards
  ###############################################
  #

  std_su = left_join(keyMap_std, data_su, by = c("plate_ID" = "plate_ID", "Pos" = "Pos"))
  std_m_su <- group_by(std_su, Sample, plate_ID) %>%
    summarize(ddCt.mean = mean(ddCt),
              ddCt.sem = sd(ddCt)/sqrt(length(ddCt)),
              RelQ.mean = mean(RelQ),
              RelQ.sem = sd(RelQ)/sqrt(length(RelQ)),
              p96 = unique(p96),
              CopyNum = mean(CopyNum)) %>%
    select(-p96) %>%
    separate(col = Sample, into = c("Sample", "Source"))


  ###############################################
  ## Copy Number Call
  ###############################################
  #   std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)
  #
  #   data_m_su$Copy.Nbs.esti = predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean))
  #         data_m_su = cbind(data_m_su, predict(std.fit, newdata = data.frame(RelQ.mean = data_m_su$RelQ.mean), interval = "prediction")[, c(2,3)])


  plateID = std_m_su$plate_ID %>% unique

  fit = lapply(seq_along(plateID), function(p) lm(CopyNum ~ RelQ.mean, data = std_m_su[std_m_su$plate_ID == plateID[p],]))

  # std.fit <- lm(CopyNum ~ RelQ.mean, data = std_m_su)

  d = list()
  for(i in 1:length(fit)){
    d[[i]] = data_m_su[data_m_su$plate_ID == plateID[i],]
    d[[i]]$Copy.Nbs.esti = predict(fit[[i]], newdata = data.frame(RelQ.mean = d[[i]]$RelQ.mean))
  }

  data_m_su = rbind_all(d)



  data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)
  data_m_su$Copy.Nbs.round = ifelse(data_m_su$zeroCopy, 0, data_m_su$Copy.Nbs.round)
  data_m_su$Copy.Nbs.esti = ifelse(data_m_su$zeroCopy, NA, data_m_su$Copy.Nbs.esti)


  # if(mode == "simple"){
  #   data_m_su$Copy.Nbs.round = ifelse(data_m_su$RelQ.mean <= 0.7, 1, 2)
  # }
  #   std_copyNbs = keyMap_std[keyMap_std$Pos %in% calib, ]$CopyNum[1]
  #   data_m_su$Copy.Nbs.esti = std_copyNbs * data_m_su$RelQ.mean
  #   data_m_su$Copy.Nbs.round = round(data_m_su$Copy.Nbs.esti)

  diff <- data_m_su$Copy.Nbs.esti-data_m_su$Copy.Nbs.round
  data_m_su$diff <- diff


  output = list()
  output$data_su = data_su
  output$data_m_su = data_m_su

  output$std_su = std_su
  output$std_m_su = std_m_su
  output

}





