#' CN object
#'
#'
#'
#' @param data
#' @importFrom  magrittr %>%
#' @importFrom  R6 R6Class
#' @import dplyr
#' @export CN
#' @export CN.cutie
#' @export CN.rivur
#' @export CN.coriell


CN = R6::R6Class(
  "CN",
  public = list(
    rawData = NA,
    diag_list=NA,
    fileName = NA,
    diag = NA,
    zeroCopy = NA,
    keyMap_std = NULL,
    calib = NULL,
    comment = NULL,
    mode = NULL,

    initialize = function(fileloc, fileloc_diag, keyMap_std = NULL, calib = NULL, comment = NULL, mode = NA){   # dataFile, dataFile_diag
      self$rawData = RIVURplate::load_data(fileloc)
      self$fileName = (list.files(fileloc, pattern = "*.csv"))

      self$zeroCopy = NULL

      if(!is.null(fileloc_diag)){
        self$diag_list = RIVURplate::load_data(fileloc_diag)
        self$diag = do.call(rbind, self$diag_list$file_list)

        zc = lapply(self$rawData$file_list, function(x) private$correct_isZero(x, self$diag))
        self$zeroCopy = do.call(rbind, zc)
      }


      self$keyMap_std = keyMap_std
      self$calib = calib
      self$comment = comment
      self$mode = mode

    },

    dataBatch = function(data_raw = self$rawData$data_raw, keyMap_std = self$keyMap_std, zeroCopy = self$zeroCopy, calib = self$calib){
      RIVURplate::makeData(data_raw = data_raw, keyMap_std = keyMap_std, zeroCopy = zeroCopy, calib = calib)
    },

    setKeyMap_std = function(keyMap_std){
      self$keyMap_std = keyMap_std
    },

    setCalib = function(calib){
      self$calib = calib
    },

    getData = function(plateID = ""){
      data = self$dataBatch()$data_m_su
      if(plateID != ""){
        data = data[data$plate_ID == plateID, ]
      }

      data

    },

    getSTD = function(plateID = ""){
      std = self$dataBatch()$std_m_su
      if(plateID != ""){
        std = std[std$plate_ID == plateID, ]
      }

      std

    },

    export = function(){
      # write.csv()
      print("OK")
    },

    plot = function(type = "std"){
      switch(type,
             std = ggplot(self$getSTD(), aes(x = CopyNum, y = RelQ.mean, label = Sample)) +
               geom_point(size = 4) +
               ylim(c(0,2)) +
               geom_smooth(method = "lm") +
               facet_grid(plate_ID ~ .) +
               theme_bw() +
               iMap::theme_science(),

             stat = {
               data_m_su = self$getData()
               data_m_su_1 <- data_m_su %>%
                 subset(RelQ.sem < 0.1) %>%
                 subset(zeroCopy == F | is.na(zeroCopy)) %>%
                 subset(RelQ.mean >0.1)

               ggplot(data_m_su_1, aes(x = RelQ.mean)) +   # , fill = Type
                 geom_bar(binwidth = 0.013, alpha = 1, col = "black", fill = "skyblue") +
                 #     geom_vline(aes(xintercept = c(one)), col = "red") +
                 #     geom_vline(aes(xintercept = c(two)), col = "blue") +
                 theme_bw() +
                 RIVURplate::theme_science() +
                 facet_grid(plate_ID ~ .) -> p1


               # Annotatation - preassumed copy number
#                data_m_su$CopyNum_preassumed = ifelse(data_m_su$RelQ.mean<0.75, 1, 2)
#                data_m_su$CopyNum_preassumed = ifelse(data_m_su$zeroCopy, 0, data_m_su$CopyNum_preassumed)

               # p2 <- plt_prop(data_m_su, x = "CopyNum_preassumed", fill = "Type")
               p2 <- plt_prop(data_m_su, x = "Copy.Nbs.round", fill = "Type")
               # p2 <- plt_prop(data_m_su, x = "CopyNum_preassumed", fill = "Type")

               grid.arrange(p1, p2, nrow = 1)
             }

             )

    },

    qc = function(d = self$getData()){
      # graphics::plot(d$Copy.Nbs.round, d$diff)
      graphics::hist(d$Copy.Nbs.esti, breaks = 50, col = "skyblue", border = "black")
      # hist(d$diff)
    },

    makeTable = function(table){
  mutate(table,
         "hasSigP" = ifelse(p_value <0.05 | p_value <0.05 | p_value < 0.05, 1, 0),
         "p_value" = ifelse(p_value < 0.05, paste0(round(p_value,3), "*"), round(p_value,3))
  )%>%
    arrange(-hasSigP, p_value) %>%
    filter(hasSigP == 1) %>%
    select(Var, p_value, N, Test) %>%
    as.data.frame %>%
    xtable(digits = 3)%>%
    print(type = "html", digits = 3)

  var = mutate(table,
               "hasSigP" = ifelse(p_value <0.05 | p_value <0.05 | p_value < 0.05, 1, 0)) %>%
    arrange(p_value) %>%
    filter(hasSigP == 1)

  var = var$Var
  table  = list()
  table$va = var
  table$n = ceiling(length(var)/3)
  table$N = length(var)
  return(table)
},

    summaryF = function(filter = "", var){
      cohort = self$grand.table(filter = filter)
      cht = CNstat$new(cohort, var)
      gridExtra::grid.arrange(cht$plot(), digitalPCR::make_table3(cht$summary()), ncol = 2)
    },

    summaryP = function(filter = "", chisq = F){
      cohort = self$grand.table(filter = filter)
      va = c(names(cohort)[grepl("LABEL", names(cohort))][2:29],
             "PRIORUTI0101", "NUM_UTI01", "UTI_LIFETIME", "BBD_CHRCONS_1302", "BBD_CHRCONS_0102",
             "CHR_CONST0102", "CHR_CONST1302", "UTI_TYPE0103")

      data_frame(Var = va,
                 stat = lapply(va, function(v) CNstat$new(cohort, v)$stat(chisq = chisq) )) %>%
        mutate(N = sapply(stat, function(s) s$n),
               p_value = sapply(stat, function(s) s$p),
               Test = sapply(stat, function(s) s$method)) %>%
        arrange(p_value) -> p_values

      p_values

      },

    assoc = function(filter = "", var = "dtBBD_LABEL2"){
      cutie = self$grand.table(filter = filter)
      CNstat$new(cutie, var)
    }

  ),

  private = list(
    correct_isZero = function(target, diag){
      plate_ID_targ = unique(target$plate_ID)
      target_m = merge(target[target$Type == "Target Unknown",], diag[toupper(diag$plate_ID) == toupper(plate_ID_targ),], by.x = "Sample.Name", by.y = "Name")
      isZero = is.na(target_m$Cp.y) & target_m$Cp.x != 0
      transform(target_m, zero = ifelse(isZero, TRUE, FALSE))
    },

    txt2csv = function(dir){
      tab = read.delim(file = dir, skip = 1)

      write.table(tab, file="/Users/DL/Documents/R project/Github-DongL/Copy Number/RIVUR plates/RIVURplate/inst/RIVUR - BTNL3-ex1/Raw data/CORIELL plate/061615.csv",sep=",",col.names=T,row.names=FALSE)
    }

  )

)

CN.cutie = R6::R6Class(
  "CN.cutie",
  inherit = CN,

  public = list(

    grand.table = function(filter = ""){
      cutie = RIVURplate::cutie_with_annotation

      cutie$COPY.NBS.ESTI = NULL
      cutie$COPY.NBS.ROUND = NULL

      if(filter != ""){
        cutie = eval(parse(text = paste0( "cutie %>% ", "subset", "(", filter , ")"   ) ))
      }


      data_m_su = self$getData()
      cutie = merge(data_m_su[data_m_su$Type == "Plate",], cutie, by = "RUID")
      cutie$CopyNum = cutie$Copy.Nbs.round
      cutie
    },

#     dataBatch = function(data_raw = self$rawData$data_raw, keyMap_std = self$keyMap_std, zeroCopy = self$zeroCopy, calib = self$calib){
#       RIVURplate::makeData_cutie(data_raw = data_raw, keyMap_std = keyMap_std, zeroCopy = zeroCopy, calib = calib)
#     },


#     cohort.summary = function(filter = ""){
#       cht = self$grand.table(filter = filter)
#       par(mfrow = c(2,2))
#       hist(cht$AGE0101, xlab = "Age",breaks = 20, main = "")
#       plot(as.factor(cht$RACE0102_LABEL), breaks = 7, xlab = "Race", main = "")
#       plot(factor(cht$SEX0101), xlab = "Sex", main = "")
#       hist(cht$Copy.Nbs.round, xlab = "copy.number.round", breaks = 12, main = "")
#       # plot(cht$TXGROUP_LABEL, xlab = "Treatment group assignment")
#       # plot(as.factor(cht$TF02_LABEL), xlab = "Treatment failure")
#     }

    cohort.summary = function(filter = ""){
      cht = self$grand.table(filter = filter)

      ggplot(cht, aes(x = AGE0101)) +   # , fill = Type
        geom_bar(binwidth = 2, alpha = 1, col = "white", fill = "skyblue") +
        theme_bw() +
        RIVURplate::theme_science() -> p1

      ggplot(cht, aes(x = RACE0102_LABEL)) +   # , fill = Type
        geom_bar(binwidth = 1, alpha = 1, col = "white", aes(fill = RACE0102_LABEL)) +
        scale_x_discrete() +
        scale_fill_brewer(palette = "Accent") +
        theme_bw() +
        RIVURplate::theme_science() +
        theme(legend.position = "none") -> p2

      ggplot(cht, aes(x = SEX0101)) +   # , fill = Type
        geom_bar(binwidth = 1, alpha = 1, col = "white", aes(fill = SEX0101)) +
        theme_bw() +
        scale_fill_brewer(palette = "Set2") +
        RIVURplate::theme_science() +
        theme(legend.position = "none") -> p3

      plt_prop(cht, x = "Copy.Nbs.round", fill = "Type", col = "white", legend.position = "none" ) -> p4


      grid.arrange(p1, p2, p3, p4, nrow = 2)
    }

  )
)



CN.rivur = R6::R6Class(
  "CN.rivur",
  inherit = CN,

  public = list(

    getData = function(plateID = ""){
      r = self$dataBatch()$data_m_su
      r = r[!is.na(r$RUID), ]

      if(plateID != ""){
        r = r[r$plate_ID == plateID, ]
      }
      r
    },

    dataBatch = function(data_raw = self$rawData$data_raw, keyMap_std = self$keyMap_std, zeroCopy = self$zeroCopy, calib = self$calib, method = "std"){
      if(method != "std") {
        dB = RIVURplate::makeData_dCt(data_raw = data_raw, keyMap_std = keyMap_std, zeroCopy = zeroCopy, calib = calib)
      }else {
        dB = RIVURplate::makeData_RIVUR(data_raw = data_raw, keyMap_std = keyMap_std, zeroCopy = zeroCopy, calib = calib, comment = self$comment, mode = self$mode)
      }
      dB
    },

    grand.table = function(filter = ""){
      rv = RIVURplate::rv_with_annotation

      rv$COPY.NBS.ESTI = NULL
      rv$COPY.NBS.ROUND = NULL

#       if(filter == T){
#         rv= rv[rv$RACE0102_LABEL == "White" & rv$SEX0101 == "F", ]
#       }
      if(filter != ""){
        rv = eval(parse(text = paste0( "rv %>% ", "subset", "(", filter , ")"   ) ))
      }

      # cutie %>% filter()


      data_m_su = self$getData()
      rv = merge(data_m_su[data_m_su$Type == "RIVUR",], rv, by = "RUID")
      rv$Copy.Nbs.round = ifelse(rv$Copy.Nbs.round < 0, 0, rv$Copy.Nbs.round)
      rv$CopyNum = rv$Copy.Nbs.round
      rv
    },

    cohort.summary = function(filter = ""){
      cht = self$grand.table(filter = filter)

      ggplot(cht, aes(x = AGE0101)) +   # , fill = Type
        geom_bar(binwidth = 2, alpha = 1, col = "white", fill = "skyblue") +
        theme_bw() +
        RIVURplate::theme_science() -> p1

      ggplot(cht, aes(x = RACE0102_LABEL)) +   # , fill = Type
        geom_bar(binwidth = 1, alpha = 1, col = "white", aes(fill = RACE0102_LABEL)) +
        scale_x_discrete() +
        scale_fill_brewer(palette = "Accent") +
        theme_bw() +
        RIVURplate::theme_science() +
        theme(legend.position = "none") -> p2

      ggplot(cht, aes(x = SEX0101)) +   # , fill = Type
        geom_bar(binwidth = 1, alpha = 1, col = "white", aes(fill = SEX0101)) +
        theme_bw() +
        scale_fill_brewer(palette = "Set2") +
        RIVURplate::theme_science() +
        theme(legend.position = "none") -> p3

      plt_prop(cht, x = "Copy.Nbs.round", fill = "Type", col = "white", legend.position = "none" ) -> p4


      ggplot(cht, aes(x = TXGROUP_LABEL)) +   # , fill = Type
        geom_bar(binwidth = 1, alpha = 1, col = "white", aes(fill = TXGROUP_LABEL)) +
        theme_bw() +
        scale_fill_brewer(palette = "Set1") +
        RIVURplate::theme_science() +
        theme(legend.position = "none") -> p5

      ggplot(cht, aes(x = TF02_LABEL)) +   # , fill = Type
        geom_bar(binwidth = 1, alpha = 1, col = "white", aes(fill = TF02_LABEL)) +
        theme_bw() +
        scale_fill_brewer(palette = "Set3") +
        RIVURplate::theme_science() +
        theme(legend.position = "none") -> p6
            # plot(cht$TXGROUP_LABEL, xlab = "Treatment group assignment")
            # plot(as.factor(cht$TF02_LABEL), xlab = "Treatment failure")


      grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)
},

    cohort.summary2 = function(filter = ""){
      cht = self$grand.table(filter = filter)
      va = c(names(cht)[!grepl("LABEL", names(cht))][2:29],
             "PRIORUTI0101", "NUM_UTI01", "UTI_LIFETIME", "BBD_CHRCONS_1302", "BBD_CHRCONS_0102",
             "CHR_CONST0102", "CHR_CONST1302", "UTI_TYPE0103")

      v = c(names(cht)[!grepl("LABEL", names(cht))])
      cht_v = cht[,v]
      cht_va_mean = apply(cht_va, 2, as.numeric())


    },

    summaryF = function(filter = "", var){
      # print("its right!")

      if(filter != ""){
        # filter = filter
        filter_A = paste0(filter, "& TXGROUP_LABEL == 'TMP/SMZ'")
        filter_P = paste0(filter, "& TXGROUP_LABEL == 'Placebo'")

        cohort_A = self$grand.table(filter = filter_A)
        cohort_P = self$grand.table(filter = filter_P)
        cohort = self$grand.table(filter = filter)
      } else{
        filter = ""
        filter_A = "TXGROUP_LABEL == 'TMP/SMZ'"
        filter_P = "TXGROUP_LABEL == 'Placebo'"

        cohort_A = self$grand.table(filter = filter_A)
        cohort_P = self$grand.table(filter = filter_P)
        cohort = self$grand.table(filter = filter)
      }


      chta = CNstat$new(cohort_A, var)
      chtp = CNstat$new(cohort_P, var)
      cht = CNstat$new(cohort, var)

      # chta = CNstat$new(cr$grand.table(filter = ""), var = "BBD1302_LABEL")
#
#       gridExtra::grid.arrange(cht$plot(), digitalPCR::make_table3(cht$summary()), ncol = 2)

      gridExtra::grid.arrange(cht$plot(ggtitle = "Both"),
                              chta$plot(ggtitle ="TMP/SMZ"),
                              chtp$plot(ggtitle ="Placebo"),
                              digitalPCR::make_table3(cht$summary()),
                              digitalPCR::make_table3(chta$summary(), ylab =""),
                              digitalPCR::make_table3(chtp$summary(), ylab =""),
                              ncol = 3)


    },

    summaryP = function(filter = "", chisq = F){
      cohort = self$grand.table(filter = filter)
      cohort_P = cohort[cohort$TXGROUP == "P", ]
      cohort_A = cohort[cohort$TXGROUP == "A", ]

#       va = c(names(cohort)[grepl("LABEL", names(cohort))][2:29],
#              "PRIORUTI0101", "NUM_UTI01", "UTI_LIFETIME", "BBD_CHRCONS_1302", "BBD_CHRCONS_0102",
#              "CHR_CONST0102", "CHR_CONST1302", "UTI_TYPE0103")


      va = c(names(cohort)[grepl("LABEL", names(cohort))][-2],
             "PRIORUTI0101", "NUM_UTI01", "UTI_LIFETIME",
             "CHR_CONST0102", "CHR_CONST1302", "UTI_TYPE0103") %>% unique()

      va[grepl("REFLUX", va)]


      getPvalue <- function(cohort = cohort){
        data_frame(Var = va,
                   stat = lapply(va, function(v) CNstat$new(cohort, v)$stat(chisq = chisq) )) %>%
          mutate(N = sapply(stat, function(s) s$n),
                 p_value = sapply(stat, function(s) s$p),
                 Test = sapply(stat, function(s) s$method)) %>%
          arrange(p_value)
      }


      A_P = full_join(getPvalue(cohort_A), getPvalue(cohort_P), by = "Var") %>%
        rename(Test = Test.x, stat_A = stat.x, N_A = N.x,
               p_value_A = p_value.x, stat_P = stat.y, Test_P = Test.y,
               N_P = N.y, p_value_P = p_value.y)


      B_A_P = full_join(A_P, getPvalue(cohort), by = "Var") %>%
        rename(Test = Test.x,
               p_value_B = p_value, stat_B = stat,
               N_B = N, p_value_B = p_value) %>%
        mutate(
               "hasSigP" = ifelse(p_value_A <0.1 | p_value_P <0.1 | p_value_B < 0.1, 1, 0),
               "p_value_A" = ifelse(p_value_A < 0.05, paste0(round(p_value_A,3), "*"), round(p_value_A,3)),
               "p_value_P" = ifelse(p_value_P < 0.05, paste0(round(p_value_P,3), "*"), round(p_value_P,3)),
               "p_value_B" = ifelse(p_value_B < 0.05, paste0(round(p_value_B,3), "*"), round(p_value_B,3)))%>%
        arrange(-hasSigP, p_value_P, p_value_A) %>%
        select(Var, p_value_A, p_value_P, p_value_B, N_A, N_P, N_P,  N_B, Test)

      B_A_P

    }

  )
)





CN.coriell = R6::R6Class(
  "CN.coriell",
  inherit = CN,

  public = list(

    getData = function(){
      r = self$dataBatch()$data_m_su
      r
    },

    dataBatch = function(data_raw = self$rawData$data_raw, keyMap_std = self$keyMap_std, zeroCopy = self$zeroCopy, calib = self$calib){
      RIVURplate::makeData_Coriell(data_raw = data_raw, keyMap_std = keyMap_std, zeroCopy = zeroCopy, calib = calib)
    },

    grand.table = function(filter = ""){
#       rv = RIVURplate::rv_with_annotation
#
#       rv$COPY.NBS.ESTI = NULL
#       rv$COPY.NBS.ROUND = NULL
#
#       #       if(filter == T){
#       #         rv= rv[rv$RACE0102_LABEL == "White" & rv$SEX0101 == "F", ]
#       #       }
#       if(filter != ""){
#         rv = eval(parse(text = paste0( "rv %>% ", "subset", "(", filter , ")"   ) ))
#       }
#
      # cutie %>% filter()


      data_m_su = self$getData()
      # rv = merge(data_m_su[data_m_su$Type == "RIVUR",], rv, by = "RUID")
      data_m_su$Copy.Nbs.round = ifelse(data_m_su$Copy.Nbs.round < 0, 0, data_m_su$Copy.Nbs.round)
      data_m_su$CopyNum = data_m_su$Copy.Nbs.round
      data_m_su
    }

  )
)










#
# #
# # # #
# cc= CN.cutie$new(fileloc = fileloc, fileloc_diag = fileloc_diag, keyMap_std = stdInfo)
# cc$plot(type = "std")
# cc$plot(type = "stat")
# cc$qc()
# cc$assoc(var = "dtBBD_LABEL2")$summary()
# cc$assoc(var = "dtBBD_LABEL2")$plot()
# library(RIVURplate)
