suppressMessages(library("ggplot2"))

#' Consistent theming across plots
jrh_theme <- function() {
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(
            colour = "#000000",
            linetype = "solid"
        ),
        axis.text.y = element_text(
            colour = "#000000"
        ),
        axis.text.x = element_text(
            angle = 90,
            vjust = 0.5,
            hjust = 1,
            face = "italic",
            colour = "#000000"
        ),
        axis.line = element_line(
            colour = "#000000",
            linetype = "solid"
        )
    )
}
