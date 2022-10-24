#' Theme for publication
#'
#' @param base_size 
#' @param base_family 
#' @param axis 
#' @param grid 
#' @param legend 
#'
#' @export
#'
theme_publication <- function(base_size = 12, 
                              base_family = 'Arial', 
                              axis = TRUE, 
                              grid = FALSE, 
                              legend.position = 'none',
                              rotate_text = 'none') {
  
  if(axis) {
    axis_element <- element_line(color = 'black')
  } else {
    axis_element <- element_blank()
  }
  
  if(grid) {
    grid_element <- element_line(color = 'grey90')
  } else {
    grid_element <- element_blank()
  }
  
  if(rotate_text == 'x') {
    x_axis_text <- element_text(angle = 90, hjust = 1)
  } else {
    x_axis_text <- element_text()
  }
  
  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) + 
    theme(plot.title = element_text(face = "plain", size = 14, hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(color = NA),
          plot.background = element_rect(color = NA),
          panel.border = element_blank(),
          axis.title = element_text(face = "plain", size = rel(1)),
          axis.title.y = element_text(angle = 90, vjust = 0.5),
          axis.title.x = element_text(vjust = 0),
          axis.text = element_text(color = 'black'),
          axis.text.x = x_axis_text,
          axis.line = axis_element,
          panel.grid.major = grid_element,
          panel.grid.minor = element_blank(),
          legend.key = element_rect(color = NA),
          legend.key.size = unit(0.4, "cm"),
          legend.background = element_rect(fill = NA),
          legend.spacing = unit(0, "cm"),
          legend.title = element_blank(), #element_text(face="italic"),
          legend.position = legend.position,
          strip.background = element_rect(color="#f0f0f0",fill="#f0f0f0"),
          strip.text = element_text(face="bold"))
}

#' Custom theme
#'
#' @param base_size 
#' @param base_family 
#' @param border 
#' @param axis 
#' @param grid 
#' @param legend.position 
#' @param facet_background 
#' @param facet_border 
#' @param facet_color 
#'
#' @export
#'
theme_dwu <- function(base_size = 12, 
                      base_family = 'Arial', 
                      border = TRUE,
                      axis = TRUE, 
                      grid = FALSE, 
                      legend.position = 'none',
                      facet_background = 'grey80',
                      facet_border = 'black',
                      facet_color = 'black') {
  if(border) {
    border_element <- element_rect(color = 'black', fill = NA)
  } else {
    border_element <- element_rect(color = NA, fill = NA)
  }
  
  if(axis) {
    axis_element <- element_line(color = 'black')
  } else {
    axis_element <- element_blank()
  }
  
  if(grid) {
    grid_element <- element_line(color = 'grey90')
  } else {
    grid_element <- element_blank()
  }
  
  
  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) + 
    theme(plot.title = element_text(face = "plain", size = 14, hjust = 0.5),
          text = element_text(),
          panel.background = element_rect(color = 'white'),
          plot.background = element_rect(color = 'white'),
          panel.border = border_element,
          axis.title = element_text(face = "plain", size = rel(1)),
          axis.title.y = element_text(angle = 90, vjust = 0.5),
          axis.title.x = element_text(vjust = 0),
          axis.text = element_text(color = 'black'),
          axis.line = axis_element,
          panel.grid.major = grid_element,
          panel.grid.minor = element_blank(),
          legend.key = element_rect(color = NA),
          legend.key.size = unit(0.4, "cm"),
          legend.background = element_rect(fill = NA),
          legend.spacing = unit(0, "cm"),
          legend.title = element_blank(), #element_text(face="italic"),
          legend.position = legend.position,
          strip.background = element_rect(color = facet_border, fill = facet_background),
          strip.text = element_text(face = "plain", color = facet_color))
}

#' Title
#'
#' @param palette 
#'
#' @export
#'
#' @examples
color_palette2 <- function(palette = 1) {
  
  if(palette == 1 | palette == 'RdBu') {
    colors <- RColorBrewer::brewer.pal(n = 11, name = 'RdBu')[c(3,9)]
  }
  
  colors
}