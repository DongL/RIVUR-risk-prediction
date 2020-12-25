#' pcrE object
#'
#' Layout of dPCR samples
#'
#' @param NULL
#' @importFrom  magrittr %>%
#' @importFrom  R6 R6Class
#' @import dplyr
#' @export rtPCR




pcrE = R6::R6Class(
  "pcrE",
  public = list(
    dataFile = NA,
    data = NA,

    initialize = function(dataFile, dataFile_diag = dataFile_diag){   # dataFile, dataFile_diag
      self$dataFile = dataFile
      self$data = read.csv(dataFile, skip = 1, header = T, na.strings = "NA")
    },

    dataSummary = function(data = self$data, data_diag = self$data_diag){
      data %>%
        group_by(Class, Target.Name, Dilution) %>%
        summarize(Ct.mean = mean(Cp, na.rm = T), Ct.sem = sd(Cp, na.rm = T)/sqrt(length(Cp)))  %>%
        filter(!is.na(Dilution)) %>%
        as.data.frame() # ->pl

    },

    plot = function(xlab = "", ylab = "Relative Expression", main = ""){
      pl = self$dataSummary()

      ggplot(pl, aes(x = Dilution, y = Ct.mean))+
        geom_point(aes(col = Class), size = 4) +
        geom_line(aes(col = Class)) +
        geom_errorbar(aes(x=Dilution, color= Class, ymin = Ct.mean -Ct.sem, ymax= Ct.mean + Ct.sem), width = 0.05, size = 0.5) +
        geom_smooth(method = "lm", col = "black", aes(group = Class)) +
        facet_grid(Target.Name ~ .) +
        scale_x_log10() +
        theme_bw()+
        iMap::theme_science()
    },

    slope = function(){
      data = self$dataSummary()
      library(data.table)
      setDT(data)
      data[, Slope := lm(Ct.mean ~ Dilution, data = .SD)$coefficients[2],
           by = c("Class", "Target.Name")
           ][, c("Class", "Target.Name", "Slope"), with = F] %>%
        unique %>%
        mutate(E = 10^(-1/Slope))

    }

  )
)


# dataFile = "/Users/DL/Documents/R project/Github-DongL/Copy Number/RIVUR plates/RIVURplate/inst/PCR efficiency/061115-Dong-DefA1A3-PCR_efficiency.csv"
#
# t = pcrE$new(dataFile = dataFile)
# t$dataSummary()
# t$plot()
# t$slope()
