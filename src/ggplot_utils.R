set_ggplot_theme <- function(map=FALSE) {
    require(ggplot2)
    require(grid)  # For unit
    theme_set(theme_bw())
    theme_update(axis.title.x=element_text(size=rel(1.2), margin=ggplot2::margin(10, 0, 0, 0)))
    theme_update(axis.title.y=element_text(angle=90, size=rel(1.2), margin=ggplot2::margin(0, 20, 0, 0)))
    theme_update(plot.margin=unit(c(1, 1, 2, 2), "lines"))
    theme_update(plot.title=element_text(size=rel(1.4)))
    theme_update(strip.background=element_blank())
    theme_update(legend.key=element_blank())
    theme_update(panel.border=element_blank())
    if(map) {
        theme_update(axis.text=element_blank())
        theme_update(axis.ticks=element_blank())
        theme_update(legend.background=element_blank())
        ## Warning: `legend.margin` must be specified using `margin()`. For the old behavior use legend.spacing
        ## theme_update(legend.margin=unit(0, "lines"))
        theme_update(legend.position="bottom")
        theme_update(panel.border=element_blank())
        theme_update(panel.grid.major=element_blank())
        theme_update(panel.grid.minor=element_blank())
    }
    return(invisible())
}