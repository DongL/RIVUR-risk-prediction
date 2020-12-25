#' Plot for RIVUR assocaition study
#'
#' @param data
#' @param x
#' @param fill
#' @param va
#' @export plt_prop
#' @export multi_plt
#' @import dplyr
#'
#'
#'


plt_prop <- function (data, x, fill, legend.position = "top", col = "white", ggtitle = "") {
  data = data[!is.na(data[fill]),]
  TABLE_FUN = eval(parse(text = paste0("table(data$", x, ", data$", fill, ")")))
  data_prop = prop.table(TABLE_FUN,2) %>% data.frame
  data_prop = data_prop[complete.cases(data_prop),]
  names(data_prop) = c(x, fill, "Proportion")

  ggplot(data_prop, aes_string(x = x, y = "Proportion",  fill = fill))+
    geom_bar(stat = "identity", col = col, position = "dodge")+
    scale_x_discrete()+
    scale_fill_brewer(palette = "Accent")+
    theme_bw() +
    RIVURplate::theme_science() +
    theme(legend.position = legend.position) +
    ggtitle(ggtitle) +
    theme(legend.key.size = unit(0.8, "cm")) +
    theme(text=element_text(family="Impact"))
}

multi_plt <- function(data = data, va = va, x_axis = "CopyNum") {
  p <- lapply(data, function(d) plt_prop(data = d, x = x_axis, fill = va))
  do.call(grid.arrange, c(p, list(ncol =3)))
}
