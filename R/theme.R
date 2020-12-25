#' theme_science object
#'
#'
#'
#' @param data
#' @importFrom  magrittr %>%
#' @importFrom  R6 R6Class
#' @import dplyr
#' @export theme_science
#' @export theme_science2


theme_science = function(){
  theme(
        # panel.grid.major = element_line(size = 0.5, color = "grey"),
        axis.line = element_line(size = 0.7, color = "black"),

        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),

        legend.position = "bottom",
        #                       legend.position = c(0.15, 0.85),
        text = element_text(size = 12),
        plot.title = element_text(size = rel(1.5))
        #                       axis.text.x = element_text(angle = 90, hjust = 1)
  )
}

theme_science2 = function(){
  theme(
    plot.title = element_text(face = "bold",
                              size = rel(1.2), hjust = 0.5),
    text = element_text(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(face = "bold",size = rel(1)),
    axis.title.y = element_text(angle=90,vjust =2),
    axis.title.x = element_text(vjust = -0.2),
    axis.text = element_text(),
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(),
    panel.grid.major = element_blank(), #  element_line(colour="#f0f0f0"),
    panel.grid.minor = element_blank(),
    legend.key = element_rect(colour = NA),
    legend.position = c(.85, .85),  # "bottom",
    legend.margin = unit(0, "cm"),
    legend.title = element_blank(), # element_text(face="italic"),
    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
    strip.text = element_text(face="bold"),
    NULL
    # legend.direction = "horizontal",
    # legend.key.size= unit(0.2, "cm"),
    # plot.margin=unit(c(10,5,5,5),"mm"),
    # panel.background = element_rect(colour = NA),
    # plot.background = element_rect(colour = NA),
    # panel.border = element_rect(colour = NA),
  )


}
