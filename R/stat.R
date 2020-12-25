#' This funciton perform statistical analysis.
#'
#' @param data a dataframe of the data that is pased
#' @param gene a character input stating the dependent variable to be used for statistical analysis
#' @param fill a character input stating the independent variable to be used for statistical analysis
#' @return a list of statistical test results including p values, number of obs, method used, and test itself as a class.
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
#' @importFrom stats kruskal.test
#' @importFrom stats wilcox.test
#' @import ggplot2
#' @importFrom xtable xtable
#' @importFrom knitr kable
#' @import reshape2
#' @export

# stat <- function(data, gene, fill){
#   Num_of_Var = unique(data[fill][!is.na(data[fill])]) %>% length
#   data = data[!is.na(data[fill]),]
#   n = nrow(data)
#   if(Num_of_Var > 2){
#     k = kruskal.test(as.formula(paste0(gene, "~", fill)), data = data)
#     #         cat("p_kruskal =", k$p.value)
#     p = k$p.value
#     method = "Kruskal Wallis test"
#     test = k
#   } else if(Num_of_Var == 2) {
#     w = wilcox.test(as.formula(paste0(gene, "~", fill)), data = data)
#     #         cat("p_wilcox =", w$p.value)
#     p = w$p.value
#     method = "Wilcoxon test"
#     test = w
#   } else {
#     #         print("Only one group!")
#     p = NA
#     method = NA
#     test = NA
#   }
#   stat = list()
#   stat$p = p
#   stat$method = method
#   stat$n = n
#   stat$test = test
#   stat$Num_of_Var = Num_of_Var
#   stat
# }





stat <- function(data, gene, fill, chisq = F){
  Num_of_Var = unique(data[fill][!is.na(data[fill])]) %>% length
  data = data[!is.na(data[fill]),]
  n = nrow(data)
  if(Num_of_Var > 2 & !chisq){
    k = kruskal.test(as.formula(paste0(gene, "~", fill)), data = data)
    #         cat("p_kruskal =", k$p.value)
    p = k$p.value
    method = "Kruskal Wallis test"
    test = k
  } else if(Num_of_Var == 2 & !chisq) {
    w = wilcox.test(as.formula(paste0(gene, "~", fill)), data = data)
    #         cat("p_wilcox =", w$p.value)
    p = w$p.value
    method = "Wilcoxon test"
    test = w
  } else {
    #         print("Only one group!")
    p = NA
    method = NA
    test = NA
  }

  if(Num_of_Var > 1 & chisq){
    # chi = chisq.test(as.formula(paste0(gene, ",", fill)), data = data)
    # chi = chisq.test(eval(parse(text = paste0("table", "(", "data$", gene, ",", "data$", fill, ")") )))
    chi = chisq.test(xtabs(as.formula(paste0("~", gene, "+", fill)), data = data))
    # chi = chisq.test(eval(parse(text = paste0("xtabs", "(", "~", gene, "+", fill,",", "data = data)"))))
    p = chi$p.value
    method = "Chisq test"
    test = chi
  }

  stat = list()
  stat$p = p
  stat$method = method
  stat$n = n
  stat$test = test
  stat$Num_of_Var = Num_of_Var
  stat
}
