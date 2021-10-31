# Load necessary libraries, install using 'install.packages("LIBRARY")' if needed
library(tidyverse)
library(DescTools)
library(ggplot2)

# Set up function necessary for raincloud plots
geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(
                ymin = min(y),
                ymax = max(y),
                xmin = x,
                xmax = x + width / 2
              )
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x)
            )
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(
              plyr::arrange(transform(data, x = xminv), y),
              plyr::arrange(transform(data, x = xmaxv), -y)
            )
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1, ])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(
            weight = 1, colour = "grey20", fill = "white", size = 0.5,
            alpha = NA, linetype = "solid"
          ),
          
          required_aes = c("x", "y")
  )


# Load data. First command will prompt a window to locate the file with subjects demographics "subject_correlations.csv" on your machine
data.path <- file.choose()
df.correlations <- read_csv(data.path)

# normalize Pearson's r using Fisher r-to-z transformation
df.correlations$cor_z <- FisherZ(df.correlations$cor_r)

# this rather complex command constructs panel B of Figure
# note: due to the randomness of the "jitter" option data points ("rain") will be randomly scattered along the y-axis

df.correlations %>% 
  filter(roi_1 == "vmPFC" & roi_2 == "pm_ins") %>%
  ggplot(aes(x = Group, y = cor_z, fill = Group)) + 
  scale_fill_manual(name = "Group", values = c("HC" = "#004B9C", "MA" = "#BB0021FF"), labels = c("MA" = "GAD")) +
  geom_flat_violin() + 
  geom_boxplot(aes(x = Group, y = cor_z),outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK", show.legend = FALSE) + 
  scale_y_continuous(breaks = seq(-0.25, 0.9, by = 0.25)) +
  scale_x_discrete(labels= c("MA" = "GAD", "HC" = "HC"), limits = c("MA", "HC")) +
  labs(x = element_blank(), y = "z-transformed correlation coefficient",
       fill = "Group") + 
  geom_point(position = position_jitter(width = .23), size = 1, show.legend = FALSE) + 
  coord_flip() + 
  # remove shape from legend
  theme_classic() + 
  theme(axis.line = element_line(colour = 'black', size = 1.1),
        text = element_text(size = 14),
        axis.text = element_text(face = "bold", color = "black"),
        axis.ticks.x = element_line(size = 1.3, color = "black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = c(0.9, 0.55),
        legend.key=element_blank())