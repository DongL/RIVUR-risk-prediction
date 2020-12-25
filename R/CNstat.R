#' CNstat object
#'
#'
#'
#' @param data
#' @importFrom  magrittr %>%
#' @importFrom  R6 R6Class
#' @import dplyr
#' @export CNstat




CNstat = R6::R6Class(
  "CNstat",
  public = list(
    data = NA,
    var = NA,

    initialize = function(data, var){   # dataFile, dataFile_diag
      self$data = data[!grepl("NA", data[, var]), ]
      self$var = var
    },

    summary = function(data = self$data, variable = self$var){
      table(data[, "Copy.Nbs.round"], data[,variable])
      # eval(parse(text = paste0("table(",  "data$Copy.Nbs.round,", "data$", var, ")" )))
      # table(data[,var])
    },

    stat = function(data = self$data, gene = "Copy.Nbs.round", fill = self$var, chisq = FALSE){
      RIVURplate::stat(data = data, gene = gene, fill = fill, chisq = chisq)
    },

    plot = function(data = self$data, fill = self$var, ggtitle = ""){
      RIVURplate::plt_prop(data, x = "Copy.Nbs.round", fill = fill, ggtitle = ggtitle)
    },


    poisson = function(data1, data2){
      d1 = data1
      d2 = data2

      d = do.call(rbind, list(d1[,c("Copy.Nbs.round", "Type")], d2[,c("Copy.Nbs.round", "Type")]))
      d$UTI = ifelse(d$Type == "Coriell", 0, 1)

      g.poisson = glm(UTI ~ Copy.Nbs.round, data = d,  family = poisson(link = log))

      summary(g.poisson)
    }


  )
)


#
# table(cutie_bbd_const$Copy.Nbs.round, cutie_bbd_const$BBD_CHRCONS)
# n = CNstat$new(cutie, "BBD_CHRCONS_1302")
# n$stat(chisq = T)
# print.table(n$summary(), zero.print = ".")
# cutie_bbd_const[,"Copy.Nbs.round"]
# data = cutie_bbd_const




