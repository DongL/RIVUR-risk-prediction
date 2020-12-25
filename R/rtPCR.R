#' rtPCR object
#'
#' Layout of dPCR samples
#'
#' @param NULL
#' @importFrom  magrittr %>%
#' @importFrom  R6 R6Class
#' @import dplyr
#' @export rtPCR


rtPCR = R6::R6Class(
  "rtPCR",
  public = list(
    dataFile = NA,
    dataFile_diag = NA,
    data = NA,
    data_diag = NA,

    initialize = function(dataFile, dataFile_diag = dataFile_diag){   # dataFile, dataFile_diag
      self$dataFile = dataFile
      self$dataFile_diag = dataFile_diag
      d = read.csv(dataFile, skip = 1, header = T, na.strings = "NA")
      names(d) = make.names(names(d))
      self$data = d
      d$Cp = ifelse(d$Cp == 0, 40, d$Cp)
      self$data_diag = read.csv(dataFile_diag, skip = 1, header = T, na.strings = "NA")
    },

    dataSummary = function(data = self$data, data_diag = self$data_diag){
      data = dplyr::left_join(data, data_diag, by = "Pos")
      data$Cp.x = ifelse(is.na(data$Cp.y), NA, data$Cp.x)

      data %>%
        group_by(Sample.Name, Type, Target.Name) %>%
        summarize(Ct.mean = mean(Cp.x, na.rm = T), Ct.sem = sd(Cp.x, na.rm = T)/sqrt(length(Cp.x)))  %>%
        transform(Note = ifelse(Ct.sem > 0.5, "High SEM(>0.5)", "")) %>%
        as.data.frame %>%
        arrange(Type)
    },

    expression = function(){
      self$dataSummary() %>%
        group_by(Sample.Name) %>%
        summarize(dCt = Ct.mean[2] - Ct.mean[1]) %>%
        mutate(RelQ = 2^-dCt) %>%
        as.data.frame()
    },

    plot = function(xlab = "", ylab = "Relative Expression", main = ""){
      expression = self$expression()

      ggplot(expression, aes(x=Sample.Name, y = RelQ))+
        geom_bar(stat = "identity",
                 # fill = "skyblue",
                 fill = "#552683",
                 col = "black")+
        coord_polar() +
        scale_y_continuous(trans ='sqrt') +
        ylab(ylab) +
        ggtitle(main) +
        xlab(xlab) +
        theme_bw()+
        iMap::theme_science() +
        iMap::theme_poster() +
        theme(axis.text.x = element_text(angle = 0))
    }

  )
)

#
# dataFile = "/Users/DL/Documents/R project/Github-DongL/Copy Number/RIVUR plates/RIVURplate/inst/RIVUR - BTNL3-ex1/Raw data/Expression/061215-KeithPierce-TISSUEPANEL-BTNL3ex1-GAPDH.csv"
#
# dataFile_diag = "/Users/DL/Documents/R project/Github-DongL/Copy Number/RIVUR plates/RIVURplate/inst/RIVUR - BTNL3-ex1/Raw data/Expression/Diag/061215-KeithPierce-TISSUEPANEL-BTNL3ex1-GAPDH-diag.csv"
#
# exp = rtPCR$new(dataFile = dataFile, dataFile_diag = dataFile_diag)
# #
# # exp$data
# # exp$dataSummary()
# # exp$plot(main = "BTNL3_exon1")
# ex = exp$expression() %>% as.data.frame()
#
# ex1 = exp$expression() %>% as.data.frame()
# ex2 = exp$expression() %>% as.data.frame()
#
# ex = merge(ex1, ex2, by = "Sample.Name", all = T)
#
# ex = hm.matrix(ex)
# ex = log(ex)
#
# z = hm.zscore(ex)
# z = z[complete.cases(z),][,c(2,4)]
# # z[,"dCt"] = z[,"RelQ"]
#
#
# library("RColorBrewer")
# col.pal <- brewer.pal(9,"Blues")
# pheatmap(z, cellwidth = 20, cellheight = 20, cluster_rows = F,
#          color = col.pal )
#          # colorRampPalette(c("navy", "white", "firebrick3"))(12)) # , color = col.pal <- brewer.pal(9,"Blues")
