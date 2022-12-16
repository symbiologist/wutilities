#' Statsplot
#'
#' @param seuratobj 
#' @param x 
#' @param y 
#' @param marginal.type 
#' @param size 
#' @param alpha 
#' @param title 
#' @param cells 
#' @param logx 
#' @param logy 
#' @param vline 
#' @param hline 
#' @param median_lines 
#' @param xlim 
#' @param ylim 
#' @param ... 
#'
#' @export
#' @import tidyverse
#'
seurat_statsplot <- function(seuratobj, 
                             x = 'nCount_RNA', 
                             y = 'nFeature_RNA', 
                             marginal.type = 'density', 
                             size = 0.5, 
                             alpha = 0.05,
                             title = seuratobj@project.name, 
                             cells = NULL, 
                             logx = FALSE, 
                             logy = FALSE, 
                             vline = NULL, 
                             hline = NULL, 
                             xlab = NULL,
                             ylab = NULL,
                             median_lines = FALSE,
                             xlim = 'auto', 
                             ylim = 'auto',
                             ...){
  
  df <- Seurat::FetchData(seuratobj, vars = c(x, y), cells = cells, ...)
  
  x_min <- min(df[,x])
  x_max <- max(df[,x])
  x_med <- median(df[,x])
  
  y_min <- min(df[,y])
  y_max <- max(df[,y])
  y_med <- median(df[,y])
  
  suppressWarnings(
    if(xlim == 'auto') {
      if(logx) {
        lim_start <- 1
      } else {
        lim_start <- 0
      }
      xlim <- c(lim_start, x_max * 1.1)
    }
  )
  
  suppressWarnings(
    if(ylim == 'auto') {
      if(logy) {
        lim_start <- 1
      } else {
        lim_start <- 0
      }
      
      ylim <- c(lim_start, y_max * 1.1)
    }  
  )
  
  p <- df %>% 
    ggplot(aes(x = df[,1], y = df[,2])) + 
    geom_point(alpha = alpha, color = 'dodgerblue4', size = size) + 
    xlab(x) + 
    ylab(y) + 
    theme_dwu(axis = TRUE,
              grid = FALSE,
              border = TRUE,
              legend.position = 'none')
  
  if(logx) {
    p <- p + scale_x_log10(limits = xlim)
  } else {
    p <- p + xlim(xlim)
  }
  if(logy) {
    p <- p + scale_y_log10(limits = ylim)
  } else {
    p <- p + ylim(ylim)
    
  }
  
  if(!is.null(vline)) {
    p <- p + geom_vline(xintercept = vline, linetype = 'dashed')
  }
  if(!is.null(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = 'dashed')
  }
  if(median_lines) {
    p <- p + geom_vline(xintercept = x_med, linetype = 'dashed', alpha = 0.5) +
      geom_hline(yintercept = y_med, linetype = 'dashed', alpha = 0.5) +
      annotate('text', x = x_med, y = y_max, label = x_med, hjust = 1.5) +
      annotate('text', x = x_max, y = y_med, label = y_med, vjust = -0.5)
  }
  if(!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  if(!is.null(xlab)) {
    p <- p + xlab(xlab)
  }
  if(!is.null(ylab)) {
    p <- p + ylab(ylab)
  }
  
  ggExtra::ggMarginal(p, type = marginal.type, margins = "both", size = 4, fill = c('#D55E00'))
  
}

#' Feature plot
#'
#' @param seuratobj 
#' @param features 
#' @param cells 
#' @param facets 
#' @param nrow 
#' @param continuous 
#' @param show 
#' @param alpha 
#' @param size 
#' @param reduction 
#' @param dims 
#' @param title 
#' @param title_params 
#' @param color_package 
#' @param color_palette 
#' @param borders 
#' @param na_break 
#' @param normalize 
#' @param scaled 
#' @param label 
#' @param label_color 
#' @param axis_title_position
#' @param legend_position 
#' @param legend_size 
#' @param facet_background 
#' @param facet_color 
#' @param facet_size 
#' @param facet_hide
#' @param verbose 
#'
#' @export
#' @import tidyverse
seurat_feature <- function(seuratobj, 
                           features = 'ident', 
                           cells = NULL,
                           facets = NULL,
                           nrow = NULL,
                           continuous = 'auto',
                           show = TRUE, 
                           alpha = 'auto', 
                           size = 'auto', 
                           reduction = 'umap', 
                           dims = 1:2,
                           title = NULL,
                           title_params = 'italics', 
                           color_package = c('custom', 'generate', 'carto', 'ggplot'),
                           color_palette = 1, 
                           borders = TRUE,
                           na_break = TRUE,
                           normalize = TRUE,
                           scaled = FALSE, 
                           label = TRUE,
                           label_color = c('white', 'black'),
                           label_size = 5,
                           drop_points = FALSE,
                           axis_title_position = 0,
                           legend_position = 'auto',
                           legend_size = 10,
                           facet_background = 'dodgerblue4',
                           facet_color = 'white',
                           facet_size = 10,
                           facet_hide = FALSE,
                           rasterize = 'auto',
                           rasterize_scale = 0.75,
                           rasterize_dpi = 300,
                           rasterize_threshold = 1e4,
                           verbose = TRUE) {
  
  # coordinates
  if(!is.null(cells)) {
    seuratobj <- subset(seuratobj, cells = cells)
  }
  
  plot_input <- Seurat::Embeddings(seuratobj, reduction = reduction)[,dims] %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
  colnames(plot_input) <- c('rowname', 'dim1', 'dim2')
  
  # Obtain features and facets
  fetched_table <- Seurat::FetchData(object = seuratobj, 
                                     vars = c(features, setdiff(facets, 'auto'))) 
  
  # Determine continuity of variables
  if(continuous == 'auto') {
    
    continuous_var <- fetched_table %>% dplyr::select(where(is.numeric)) %>% colnames()
    discrete_var <- colnames(fetched_table) %>% setdiff(continuous_var)
    
    continuous <- ifelse(length(continuous_var > 0), TRUE, FALSE)
    discrete <- ifelse(length(discrete_var > 0), TRUE, FALSE)
  }
  
  if(continuous) {
    continuous_table <- fetched_table %>% 
      select(all_of(continuous_var)) %>% 
      rownames_to_column() %>% 
      as_tibble() %>% 
      pivot_longer(cols = -1,
                   names_to = 'feature',
                   values_to = 'value') 
  } 
  if(discrete) {
    discrete_table <- fetched_table %>% 
      select(all_of(discrete_var)) %>% 
      rownames_to_column() %>% 
      as_tibble()
  }
  
  # normalization and scaling
  if(continuous) {
    if(scaled) {
      continuous_table <- continuous_table %>% 
        group_by(feature) %>% 
        mutate(value = as.numeric(scale(value, center = TRUE, scale = TRUE))) %>% 
        ungroup()
      
      na_break <- FALSE
    }
    
    if(normalize) {
      continuous_table <- continuous_table %>% 
        group_by(feature) %>% 
        mutate(value = value/max(value)) %>% 
        ungroup()
    }
    
    # Set 0's to a separate color
    if(na_break) {
      continuous_table <- continuous_table %>% 
        mutate(value = replace(value, value == 0, NA))
    }
  }
  
  # merge table
  if(discrete) {
    plot_input <- plot_input %>% 
      left_join(discrete_table, by = c('rowname' = 'rowname'))
  }
  if(continuous) {
    plot_input <- plot_input %>% 
      left_join(continuous_table, by = c('rowname' = 'rowname')) %>% 
      dplyr::arrange(!is.na(value), value) %>% 
      dplyr::mutate(feature = factor(feature, levels = features))
  } else {
    plot_input[['value']] <- plot_input[[discrete_var[[1]]]]
  }
  
  # auto facet multiple continuous variables
  if(length(continuous_var) > 1) {
    facets <- 'feature'
  }
  
  # Color scales
  # continuous
  color_package <- color_package[1]
  
  if(continuous) {
    if (color_package == 'custom' & length(color_palette) == 1) {
      if(color_palette == 1) {
        color_palette = c('grey90', 'grey90', 'darkorchid4')
      } else if (color_palette == 2) {
        color_palette <- c('grey90','wheat3', 'darkorchid4')
      } else if (color_palette == 3) {
        color_palette <- c('grey90', 'dodgerblue4', 'orangered')
      } 
    }
  } else {
    # discrete colors
    n_colors <- dplyr::n_distinct(plot_input$value)
    
    if(n_colors > 12 & color_palette == 1) {
      color_package <- 'ggplot'
    }
    
    if (color_package == 'custom' & length(color_palette) == 1) {
      if (color_palette == 1) {
        color_package <- 'carto'
        color_palette <- 'Prism'
      }
    } 
    if (color_package == 'generate') {
      color_palette <- generate_colors(n = n_colors, palette = color_palette)
    } else if (color_package == 'carto') {
      if(n_colors > 12) {
        color_palette <- rep(rcartocolor::carto_pal(n = 12, name = color_palette), ceiling(n_colors/12))
      } else {
        color_palette <- rcartocolor::carto_pal(n = n_colors, name = color_palette)
      }
    } 
  }
  
  # Styling
  if (reduction == 'pca') {
    xlabel <- paste0('PC', dims[1])
    ylabel <- paste0('PC', dims[2])
  } else if (reduction == 'tsne') {
    xlabel <- bquote(italic('t')*'-SNE'*.(dims[1]))
    ylabel <- bquote(italic('t')*'-SNE'*.(dims[2]))
  } else {
    xlabel <- paste0(toupper(reduction), dims[1])
    ylabel <- paste0(toupper(reduction), dims[2])
  }
  
  if(title_params == 'italics') {
    title_params <- element_text(face = 'italic', size = rel(1), hjust = 0.5)
  } else if(title_params == 'bold') {
    title_params <- element_text(face = "bold", size = rel(1.2), hjust = 0.5)
  }
  
  if(is.null(facets)) {
    if(!continuous) {
      plot_input <- plot_input %>% mutate(feature = discrete_var[1])
    }
    facets <- 'feature'
  }
  
  n_cells <- ncol(seuratobj)
  
  if(alpha == 'auto') {
    alpha <- max(30/sqrt(n_cells), 0.2)
  }
  if(size == 'auto') {
    size <- max(20/sqrt(n_cells), 0.05)
  }
  if(verbose) {
    print2('Plotting with alpha ', round(alpha, 2), ' and size ', round(size, 2))
  }
  
  # plot
  p <- ggplot(plot_input, 
              aes(x = dim1,
                  y = dim2,
                  color = value)) +
    facet_wrap(reformulate(facets), nrow = nrow) +
    theme_minimal() +
    labs(x = xlabel,
         y = ylabel) +
    theme(plot.background = element_rect(fill = 'white', color = NA),
          panel.grid = element_blank(),
          text = element_text(family = 'Arial'),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_blank(),
          axis.title = element_text(hjust = axis_title_position),
          legend.title = element_blank(),
          legend.text = element_text(size = legend_size)) 
  

  # rasterize points if greater than rasterize threshold
  if(rasterize == 'auto') {
    if(n_cells > rasterize_threshold) {
      rasterize <- TRUE
    } else {
      rasterize <- FALSE
    }
  } 
  
  if(drop_points) {
    # do not add points
  } else if(rasterize) {
    p <- p + 
      ggrastr::rasterize(geom_point(size = size, 
                                    alpha = alpha), 
                         dpi = rasterize_dpi, 
                         scale = rasterize_scale)
  } else {
    p <- p + geom_point(size = size, 
                        alpha = alpha)
  }
  
  if(!continuous & label) {
    centers <- plot_input %>% dplyr::group_by(value) %>% dplyr::summarize(x = median(dim1), y = median(dim2))
    p <- p + ggrepel::geom_text_repel(data = centers, 
                                      aes(x = x, y = y, label = value), 
                                      size = label_size, 
                                      color = label_color[1], 
                                      bg.color = label_color[2],
                                      max.overlaps = Inf,
                                      point.size = NA)
  }
  
  if(borders) {
    p <- p + theme(panel.background = element_rect(fill = 'white', color = 'black'))
  } else {
    p <- p + theme(panel.background = element_rect(fill = 'white', color = NA))
  }
  
  if(continuous) {
    if (scaled) {
      p <- p +
        scale_color_gradient2(low = color_palette[2], mid = color_palette[1], high = color_palette[3], labels = NULL)
    } else {
      p <- p +
        scale_color_gradient(na.value = color_palette[1], low = color_palette[2], high = color_palette[3], labels = NULL)
    }
  } else {
    if(color_package == 'ggplot') {
      p <- p + scale_color_discrete()
    } else {
      p <- p + scale_color_manual(values = color_palette)  
    } 
  }
  
  if(legend_position == 'auto') {
    legend_position <- case_when(
      length(continuous_var) > 0 ~ 'none',
      n_distinct(plot_input$value) < 6 ~ 'bottom',
      TRUE ~ 'right'
    )
    
  }
  p <- p + theme(legend.position = legend_position) + guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))
  
  if(!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  if(!facet_hide) {
    p <- p + theme(strip.text = element_text(color = facet_color, size = facet_size),
                   strip.background = element_rect(color = 'black', fill = facet_background))
  } else {
    p <- p + theme(strip.text = element_blank(),
                   strip.background = element_blank())
    
  }
  
  p
}

#' Calculate density
#'
#' @param coordinates 
#' @param x 
#' @param y 
#' @param xlims 
#' @param ylims 
#' @param bins 
#' @param bw 
#' @param custom 
#' @param rel_threshold 
#' @param abs_threshold 
#'
#' @export
#'
density_calculate <- function(coordinates,
                              x = 'UMAP_1',
                              y = 'UMAP_2',
                              xlims = range(coordinates[[x]]), 
                              ylims = range(coordinates[[y]]), 
                              bins = 100,
                              bw = 'bcv',
                              custom = NULL,
                              rel_threshold = 0, # remove densities below this relative threshold (0 - 1)
                              abs_threshold = 0) { # remove densities below this absolute threshold 
  
  coord_df <- tibble(x = coordinates[[x]],
                     y = coordinates[[y]])
  
  # bandwidth selector, or manually provide custom bandwidths for x and y
  if(bw == 'custom') {
    bw.x <- custom[1]
    bw.y <- custom[2]
  } else if(bw == 'bcv') {
    bw.x <- bw.bcv(coord_df$x)
    bw.y <- bw.bcv(coord_df$y)
  } else if(bw == 'ucv') {
    bw.x <- bw.ucv(coord_df$x)
    bw.y <- bw.ucv(coord_df$y)
  } else if (bw == 'SJ') {
    bw.x <- bw.SJ(coord_df$x)
    bw.y <- bw.SJ(coord_df$y)
  } else if (bw == 'nrd0') {
    bw.x <- bw.nrd0(coord_df$x)
    bw.y <- bw.nrd0(coord_df$y)
  } else {
    bw.x <- bw.nrd(coord_df$x)
    bw.y <- bw.nrd(coord_df$y)
  }
  
  density_computation <- KernSmooth::bkde2D(as.matrix(coord_df),
                                            bandwidth = c(bw.x, bw.y),
                                            gridsize = c(bins, bins),
                                            range.x = list(xlims, ylims))
  
  names(density_computation) <- c('x', 'y', 'z')
  density_table <- crossing(x = density_computation$x,
                            y = density_computation$y) %>% 
    arrange(y, x) %>% 
    mutate(z = as.numeric(density_computation$z),
           z_scaled = z/max(z), # scale values between 0 and 1
           z_norm = z/norm(z, '2')) %>% # unit vector normalization
    mutate(z_norm = ifelse(z_scaled < rel_threshold, 0, z_norm), # threshold on densities above relative threshold
           z_norm = ifelse(z_norm < abs_threshold, 0, z_norm)) # threshold on densities above absolute threshold
  
  density_table
} 

#' Plot normalized density
#'
#' @param input 
#' @param x 
#' @param y 
#' @param xlims 
#' @param ylims 
#' @param group 
#' @param target 
#' @param background 
#' @param bins 
#' @param bw 
#' @param custom 
#' @param expansion 
#' @param plot 
#' @param rel_threshold 
#' @param abs_threshold 
#' @param downsample_data
#' @param downsample_points 
#' @param point_alpha 
#' @param point_size 
#' @param point_color 
#' @param clip_zeros 
#' @param colors 
#' @param background_color 
#' @param borders 
#' @param table_only 
#'
#' @export
#'

density_plot <- function(input,
                         x = 'UMAP_1',
                         y = 'UMAP_2',
                         xlims = NULL, 
                         ylims = NULL,
                         group = NULL,
                         target = NULL, 
                         background = 'all',
                         bins = 100, 
                         bw = 'custom', 
                         custom = NULL,
                         expansion = 0.05,
                         plot = TRUE,
                         rel_threshold = 0, # remove densities below this relative threshold (0 - 1)
                         abs_threshold = 0, # remove densities below this absolute threshold
                         downsample_data = 1, # downsampling all data for computational efficiency
                         downsample_points = 1, # if plotting points, downsample to this fraction for computational efficiency
                         point_alpha = 0.1, 
                         point_size = 0.1,
                         point_color = 'grey90',
                         clip_zeros = TRUE, # in background normalization, clip densities below 0
                         colors = 'custom1',
                         background_color = 'transparent',
                         borders = TRUE,
                         table_only = FALSE) {
  
  # Determine input format
  suppressWarnings(
    if(class(input) == 'Seurat') {
      coordinates <- FetchData(input, vars = c(x, y, group)) %>% as_tibble()
    } else {
      coordinates <- input %>% select(any_of(c(x, y, group)))
    })
  
  # Set column names
  colnames(coordinates) <- c('x', 'y', 'group')
  
  # Set limits
  if(is.null(xlims)) {
    xlims = range(coordinates$x)
  }
  
  if(is.null(ylims)) {
    ylims = range(coordinates$y)
  }
  
  if(downsample_data < 1) {
    coordinates <- coordinates %>% dplyr::sample_frac(size = downsample_data)
  } else if (downsample_data > 1) {
    coordinates <- coordinates %>% dplyr::sample_n(size = downsample_data)
  }
  
  if(is.null(target)) { 
    coordinates_subset <- coordinates # use all cells if null
  } else {
    coordinates_subset <- coordinates %>% filter(group %in% target) 
  }
  
  # expand x and y limits if desired
  xlims <- xlims %>% range() %>% scales::expand_range(mul = expansion) 
  ylims <- ylims %>% range() %>% scales::expand_range(mul = expansion) 
  
  if(bw == 'custom' & is.null(custom)) {
    custom <-  c(bw.bcv(coordinates$x) * 8,
                 bw.bcv(coordinates$y) * 8)
  }
  
  if(is.null(group)) {
    
    density_all <- density_calculate(coordinates,
                                     x = 'x',
                                     y = 'y',
                                     bins = bins,
                                     xlims = xlims,
                                     ylims = ylims,
                                     bw = bw,
                                     custom = custom,
                                     rel_threshold = rel_threshold,
                                     abs_threshold = abs_threshold) 
    
  } else {
    
    # target density 
    density_target <- density_calculate(coordinates_subset,
                                        x = 'x',
                                        y = 'y',
                                        bins = bins,
                                        xlims = xlims,
                                        ylims = ylims,
                                        bw = bw,
                                        custom = custom) 
    
    # background
    if(!is.null(background)) {
      
      #  background
      if(background == 'all') {
        coordinates_background <- coordinates
      } else {
        coordinates_background <- coordinates %>% filter(group %in% background)
      }
      
      
      density_background <- density_calculate(coordinates_background,
                                              x = 'x',
                                              y = 'y',
                                              bins = bins,
                                              xlims = xlims,
                                              ylims = ylims,
                                              bw = bw,
                                              custom = custom)
      
      # background subtraction
      density_subtraction <- density_target %>% 
        mutate(z_diff = z_norm - density_background$z_norm,
               z_scaled = z_diff/max(z_diff),
               z_diff = ifelse(z_scaled < rel_threshold, 0, z_diff),
               z_diff = ifelse(z_diff < abs_threshold, 0, z_diff))
      
      # set negative density differences to zero
      if(clip_zeros) {
        density_subtraction <- density_subtraction %>% 
          mutate(z_diff = ifelse(z_diff < 0, 0, z_diff),
                 z_diff = ifelse(z_diff < abs_threshold, 0, z_diff),
                 z_diff_scaled = z_diff / max(z_diff))
      }
    } 
  }
  
  # set up plotting parameters
  
  if(length(colors) == 1) {
    colors <- case_when(
      colors == 'custom1' ~ rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu"))[7:11],
      colors == 'custom2' ~ viridis::viridis(n = 5)
    )
  }
  
  if(!is.null(background_color)) {
    colors <- c(background_color, colors) 
  }
  
  color_scale <- colorRampPalette(colors)(bins) # color scale 
  coordinates_downsample <- coordinates %>% sample_frac(downsample_points) # downsampled points
  
  plot_layers <- list(geom_point(inherit.aes = F,
                                 data = coordinates_downsample, 
                                 aes(x = x, y = y),
                                 alpha = point_alpha,
                                 color = point_color,
                                 size = point_size),
                      geom_contour_filled(bins = bins, aes(alpha = ..level..)),
                      xlim(xlims),
                      ylim(ylims),
                      scale_fill_manual(values = color_scale),
                      theme_minimal(),
                      theme(text = element_text(family = 'Arial'),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            plot.title = element_text(face = 'plain', family = 'Arial', size = 10, hjust = 0.5),
                            panel.grid = element_blank(),
                            legend.position = 'none',
                            axis.title = element_blank()))
  
  if(borders) {
    theme_border <- theme(plot.background = element_rect(fill = NA, color = 'black', size = 0.5))  
  }
  
  
  if(!is.null(group)) {
    
    p_output <- density_subtraction %>% 
      ggplot(aes(x = x, 
                 y = y, 
                 z = z_diff)) +
      plot_layers 
    
  } else {
    
    p_output <- density_all %>% 
      ggplot(aes(x = x, 
                 y = y, 
                 z = z_norm)) +
      plot_layers 
    
  }
  
  if(borders) {
    p_output <- p_output + theme_border
  }
  
  
  if(table_only) {
    density_subtraction
  } else {
    p_output
  }
}


#' Title
#'
#' @param seuratobj 
#' @param cells 
#' @param genes 
#' @param total 
#' @param normalize 
#' @param log 
#' @param assay 
#' @param slot 
#'
#' @return
#' @export
#'
#' @examples
extract_expression <- function(seuratobj,
                               cells = NULL,
                               genes,
                               total = TRUE,
                               normalize = TRUE,
                               log = FALSE,
                               assay = 'RNA',
                               slot = 'counts') {
  
  data_subset <- FetchData(object = seuratobj, 
                           cells = cells,
                           vars = genes,
                           assay = assay,
                           slot = slot) %>% 
    rownames_to_column() %>% 
    pivot_longer(cols = -rowname,
                 names_to = 'gene',
                 values_to = 'counts')
  
  if(total) {
    total_count <- colSums(GetAssayData(seuratobj, 
                                        slot = slot, 
                                        assay = assay))  
    if(!is.null(cells)) {
      total_count <- total_count[cells]
    }
    
    data_subset <- data_subset %>% 
      left_join(data.frame('total' = total_count) %>% rownames_to_column(), by = 'rowname')
  }
  
  if(normalize) {
    data_subset <- data_subset %>% 
      mutate(cpm = 1e6 * counts / total) 
  }
  
  data_subset
  
}

#' Title
#'
#' @param seuratobj 
#' @param genes 
#' @param cells 
#' @param group_by 
#' @param pseudocount 
#' @param rounding 
#' @param assay 
#' @param slot 
#'
#' @return
#' @export
#'
#' @examples
pseudobulk <- function(seuratobj,
                       genes,
                       cells = NULL,
                       group_by = 'ident',
                       pseudocount = 1,
                       rounding = 3,
                       assay = 'RNA',
                       slot = 'count') {
  
  metadata_groups <- FetchData(seuratobj, vars = group_by, cells = cells) %>% rownames_to_column()
  
  expression <- extract_expression(seuratobj,
                                   genes = genes,
                                   cells = cells,
                                   normalize = FALSE,
                                   log = FALSE,
                                   assay = assay,
                                   slot = slot) %>% 
    left_join(metadata_groups, by = 'rowname')
  
  totals <- expression %>% 
    select(rowname, total, all_of(group_by)) %>% 
    unique() %>% 
    group_by(across(all_of(group_by))) %>% 
    summarize(total = sum(total))
  
  merged <- expression %>% 
    group_by(across(all_of(c('gene', group_by)))) %>% 
    summarize(counts = sum(counts)) %>% 
    ungroup() %>% 
    left_join(totals) %>% 
    mutate(cpm = counts/total * 1e6,
           log10cpm = log10(cpm + pseudocount))
  
  merged %>% mutate(across(contains('cpm'), round, rounding))
  
}

#' Title
#'
#' @param seuratobj 
#' @param cells 
#' @param group_by 
#'
#' @return
#' @export
#'
#' @examples
pseudobulk_matrix <- function(seuratobj,
                              cells = NULL,
                              group_by = 'ident') {
  
  metadata_groups <- FetchData(seuratobj, vars = group_by, cells = cells) 
  metadata_groups$identifier <- apply( metadata_groups[ , group_by ] , 1 , paste , collapse = ':')
  
  if(is.null(cells)) {
    cells <- colnames(seuratobj)
  }
  
  # Ensure identical cell ordering
  metadata_groups <- metadata_groups[cells, ]
  expression_matrix <- t(seuratobj@assays$RNA@counts[, cells])
  aggregated_matrix <- sapply(by(expression_matrix, metadata_groups$identifier, colSums), identity)
  
  # slower for unknown reason:
  # tic()
  # expression_matrix <- seuratobj@assays$RNA@counts[, cells]
  # aggregated_matrix <- aggregate_matrix(expression_matrix, cols = metadata_groups$identifier)
  # toc()
  
  aggregated_matrix
} 