#' Statsplot
#'
#' @param seurat_obj 
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
seurat_statsplot <- function(seurat_obj, 
                             x = 'nCount_RNA', 
                             y = 'nFeature_RNA', 
                             marginal.type = 'density', 
                             size = 0.5, 
                             alpha = 0.05,
                             title = seurat_obj@project.name, 
                             cells = NULL, 
                             logx = FALSE, 
                             logy = FALSE, 
                             vline = NULL, 
                             hline = NULL, 
                             median_lines = FALSE,
                             xlim = 'auto', 
                             ylim = 'auto',
                             ...){
  
  df <- Seurat::FetchData(seurat_obj, vars = c(x, y), cells = cells, ...)
  
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
  
  ggExtra::ggMarginal(p, type = marginal.type, margins = "both", size = 4, fill = c('#D55E00'))
  
}

#' Feature plot
#'
#' @param seurat_obj 
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
seurat_feature <- function(seurat_obj, 
                           features = 'seurat_clusters', 
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
                           legend_position = 'auto',
                           legend_size = 10,
                           facet_background = 'dodgerblue4',
                           facet_color = 'white',
                           facet_size = 10,
                           facet_hide = FALSE,
                           verbose = TRUE) {
  
  # coordinates
  plot_input <- Seurat::Embeddings(seurat_obj, reduction = reduction)[,dims] %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
  colnames(plot_input) <- c('rowname', 'dim1', 'dim2')
  
  # Obtain features and facets
  fetched_table <- Seurat::FetchData(object = seurat_obj, 
                                     vars = c(features, setdiff(facets, 'auto')),
                                     cells = cells) 
  
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
    
    if (color_package == 'custom' & length(color_palette) == 1) {
      if (color_palette == 1) {
        color_package <- 'carto'
        color_palette <- 'Prism'
      }
    } 
    if (color_package == 'generate') {
      color_palette <- generate_colors(n = n_colors, palette = color_palette)
    } else if (color_package == 'carto') {
      color_palette <- rcartocolor::carto_pal(n = n_colors, name = color_palette)
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
  
  n_cells <- ncol(seurat_obj)
  
  if(alpha == 'auto') {
    alpha <- 30/sqrt(n_cells)
  }
  if(size == 'auto') {
    size <- 20/sqrt(n_cells)
  }
  if(verbose) {
    print2('Plotting with alpha ', round(alpha, 2), ' and size ', round(size, 2))
  }
  
  # plot
  p <- ggplot(plot_input, 
              aes(x = dim1,
                  y = dim2,
                  color = value)) +
    geom_point(size = size, 
               alpha = alpha) +
    facet_wrap(reformulate(facets), nrow = nrow)
  
  if(!continuous & label) {
    centers <- plot_input %>% dplyr::group_by(value) %>% dplyr::summarize(x = median(dim1), y = median(dim2))
    p <- p + shadowtext::geom_shadowtext(data = centers, aes(x = x, y = y, label = value), color = label_color[1], size = 5, bg.color = label_color[2])
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
  
  p <- p + theme_minimal() +
    labs(x = xlabel,
         y = ylabel) +
    theme(plot.background = element_rect(fill = 'white', color = NA),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = legend_size)) 
  
  if(!facet_hide) {
    p <- p + theme(strip.text = element_text(color = facet_color, size = facet_size),
                   strip.background = element_rect(color = 'black', fill = facet_background))
  } else {
    p <- p + theme(strip.text = element_blank(),
                   strip.background = element_blank())
    
  }
  
  p
}