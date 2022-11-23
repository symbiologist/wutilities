# library(tidyverse)
# library(patchwork)
# library(ggrepel)
# library(tictoc)
# library(crayon)
# library(mixtools)
# library(parallel)

# Scatter plot with optional trendline and marginal plots

scatter_stats <- function(input_table,
                          x,
                          y,
                          color = NULL,
                          xlab = x,
                          ylab = y,
                          log_transform = F,
                          pseudocount = 1,
                          log_scale = F,
                          marginal = T,
                          marginal_type = 'density',
                          marginal_fill = 'grey80',
                          marginal_color = 'black',
                          marginal_alpha = 0.8,
                          title = '',
                          point_size = 0.5,
                          point_alpha = 'auto',
                          point_color = '#1D6996',
                          trendline = T,
                          trendline_se = F,
                          line_size = 1,
                          line_type = 'solid',
                          line_color = '#EDAD08',
                          xlims = NULL,
                          ylims = NULL,
                          xbreaks = NULL,
                          ybreaks = NULL,
                          correlation = T,
                          coord_equal = T) {
  
  table_subset <- input_table %>% select(one_of(x, y, color))
  
  if(is.null(color)) {
    colnames(table_subset) <- c('X', 'Y')
  } else {
    colnames(table_subset) <- c('X', 'Y', 'Color')
  }
  
  
  if(log_transform) {
    table_subset <- table_subset %>% mutate_all(function(i) {log10(i + pseudocount)}) %>% filter(is.finite(X), is.finite((Y)))
  }
  
  if (point_alpha == 'auto') {
    point_alpha <- 100*sqrt(point_size)/sqrt(nrow(table_subset)) 
  }
  
  if(correlation) {
    if(log_scale) {
      correlation_subset <- table_subset %>% mutate_all(function(i) {log10(i + pseudocount)}) %>% filter(is.finite(X), is.finite((Y)))
    } else {
      correlation_subset <- table_subset %>% drop_na()
    }
    r <- round(cor(correlation_subset$X, correlation_subset$Y, use = 'complete.obs'), 2)
    
  }
  
  if(is.null(color)) {
    p <- table_subset %>% 
      ggplot(aes(x = X,
                 y = Y)) + 
      geom_point(size = point_size,
                 alpha = point_alpha,
                 color = point_color) 
    
  } else {
    p <- table_subset %>% 
      ggplot(aes(x = X,
                 y = Y,
                 color = Color)) +
      geom_point(size = point_size,
                 alpha = point_alpha) 
  }
  
  
  if(coord_equal) {
    p <- p + coord_equal()  
  }
  
  if(log_scale) {
    if(is.null(xbreaks)) {
      p <- p + 
        scale_x_log10(limits = xlims) 
    } else {
      p <- p + 
        scale_x_log10(limits = xlims, breaks = xbreaks) 
    }
    if(is.null(ybreaks)) {
      p <- p + 
        scale_y_log10(limits = ylims) 
    } else {
      p <- p + 
        scale_y_log10(limits = ylims, breaks = ybreaks) 
    }
  } else {
    if(is.null(xbreaks)) {
      p <- p + 
        scale_x_continuous(limits = xlims) 
    } else {
      p <- p + 
        scale_x_continuous(limits = xlims, breaks = xbreaks) 
    }
    if(is.null(ybreaks)) {
      p <- p + 
        scale_y_continuous(limits = ylims) 
    } else {
      p <- p + 
        scale_y_continuous(limits = ylims, breaks = ybreaks) 
    }
  }
  if(trendline) {
    p <- p + geom_smooth(method = 'lm', color = line_color, size = line_size, se = trendline_se)
  }
  
  p <- p + theme_publication(grid = F) +
    labs(title = title,
         subtitle = paste('Pearson r =', r),
         x = xlab,
         y = ylab) +
    theme(plot.title = element_text(face = 'plain', size = rel(1.1)),
          plot.subtitle = element_text(hjust = 0.5))
  
  if(marginal) {
    ggExtra::ggMarginal(p, 
                        type = marginal_type, 
                        fill = marginal_fill, 
                        color = marginal_color,
                        alpha = marginal_alpha)
    
  } else {
    p
  }
}

# run an Fisher exact test by setting up a 2 x 2 contingency table
# Enrichment testing

# color scales



### mixture model
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

plot_mix_model <- function(mixmdl,
                           alpha = 0.5,
                           bins = 1000,
                           colors = c('orangered', 'dodgerblue4'),
                           xlab = '',
                           ylab = 'Density',
                           intercept = NULL) {
  
  data.frame(x = mixmdl$x) %>%
    ggplot() +
    geom_histogram(aes(x, ..density..),
                   alpha = alpha,
                   bins = bins) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl$mu[1], mixmdl$sigma[1], lam = mixmdl$lambda[1]),
                  color = colors[1], lwd = 1.5) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(mixmdl$mu[2], mixmdl$sigma[2], lam = mixmdl$lambda[2]),
                  color = colors[2], lwd = 1.5) +
    ylab(ylab) +
    xlab(xlab) +
    geom_vline(xintercept = intercept, linetype = 'dashed') +
    theme_publication()
  
}

mixture_model <- function(input_vector, 
                          log_transform = FALSE,
                          k = 2,
                          alpha = 0.5,
                          bins = 1000,
                          colors = c('orangered', 'dodgerblue4'),
                          xlab = '',
                          ylab = 'Density',
                          intercept = 'auto') {
  
  mixmdl <- normalmixEM(input_vector, k = k)
  
  if(intercept == 'auto') {
    intercept <- mixture_cutoff(mixmdl)
  }
  
  # plot
  mix_plot <- plot_mix_model(mixmdl,
                             alpha = alpha,
                             bins = bins,
                             colors = colors,
                             xlab = xlab,
                             ylab = ylab,
                             intercept = intercept)
  
  list('model' = mixmdl,
       'plot' = mix_plot,
       'cutoff' = intercept)
}

mixture_cutoff <- function(mixmod) {
  uniroot(f = function(x) 
    plot_mix_comps(x,
                   mixmod$mu[1],
                   mixmod$sigma[1],
                   mixmod$lambda[1]) -
      plot_mix_comps(x,
                     mixmod$mu[2],
                     mixmod$sigma[2],
                     mixmod$lambda[2]), mixmod$mu)$root
}

save_figure <- function(plot = last_plot(), 
                        filename, 
                        device = c('png', 'pdf'), 
                        dpi = 600, 
                        w = 5, 
                        h = 5, 
                        units = 'in', 
                        directory = 'figures', 
                        size = 'custom', 
                        gg = TRUE,
                        overwrite = TRUE) {
  
  # Some default sizes
  if(size == 'small') {
    w = 5
    h = 5
  } else if(size == 'medium') {
    w = 7
    h = 7
  } else if(size == 'large'){
    w = 10
    h = 10
  } else if(size == 'wide') {
    h = 5
    w = 8
  } else if(size == 'xwide') {
    h = 5
    w = 12
  } else if(size == 'long') {
    h = 12
    w = 8
  } else if(size == 'xlong') {
    h = 12
    w = 5
  }
  
  # prevent overwrites
  i = 0
  prefix <- filename
  map(device, function(x) {
    if(!overwrite){
      while(file.exists(paste0(directory, filename,'.', x)))
      {
        i = i + 1
        ifilename = paste0(prefix,i)
      }
    }
    cat('Saving as', paste0(filename, '.', x), '\n')
    
    filename = paste0(filename, '.', x)
    if(gg) { # for ggplot objects
      if(x == 'pdf') {
        x <- cairo_pdf
      }
      ggsave(filename = filename, plot = plot, dpi = dpi, path = directory, width = w, height = h, units = units)}
    else { # for others
      dev.copy(png, file = paste0(directory, filename), width = w, height = h, units = units, res = dpi)
      dev.off() }
    
  })
}

### enrichment
tabulate_frequency <- function(metadata_table,
                               group1 = guide_target,
                               group2 = seurat_clusters,
                               rename = TRUE) {
  frequency_table <- 
    metadata_table %>% 
    mutate_all(as.factor) %>% 
    mutate_all(droplevels) %>% 
    dplyr::select({{group1}}, 
                  {{group2}}) 
  
  if(rename) {
    colnames(frequency_table) <- c('target', 'cluster')
  }
  
  frequency_table %>% table()
  
}

calculate_proportions <- function(frequency_table) {
  
  target_proportion <- prop.table(frequency_table, margin = 1) %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(target_proportion = round(Freq, 3)) %>% 
    dplyr::select(-Freq)
  
  cluster_proportion <- prop.table(frequency_table, margin = 2) %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(cluster_proportion = round(Freq, 3)) %>% 
    dplyr::select(-Freq)
  
  frequency_long <- frequency_table %>% 
    as.data.frame() %>% 
    as_tibble()
  
  colnames(frequency_long) <- c('target', 'cluster', 'n_overlap')
  
  cluster_totals <- colSums(frequency_table) %>% 
    as.data.frame() %>% 
    rownames_to_column()
  
  colnames(cluster_totals) <- c('cluster', 'n_cluster')
  
  target_totals <- rowSums(frequency_table) %>% 
    as.data.frame() %>% 
    rownames_to_column()
  colnames(target_totals) <- c('target', 'n_target')
  
  proportion_table <- inner_join(cluster_proportion,
                                 target_proportion) %>% 
    inner_join(frequency_long) %>% 
    inner_join(target_totals) %>% 
    inner_join(cluster_totals) %>% 
    mutate(max_proportion = ifelse(cluster_proportion > target_proportion, cluster_proportion, target_proportion),
           score = round(cluster_proportion * target_proportion, 3)) %>% 
    arrange(-max_proportion)
  
  proportion_table
  
}

proportion_histograms <- function(proportion_table) {
  list(proportion_table %>% ggplot(aes(x = cluster_proportion)) + geom_histogram() + scale_x_log10(),
       proportion_table %>% ggplot(aes(x = target_proportion)) + geom_histogram() + scale_x_log10(),
       proportion_table %>% ggplot(aes(x = score)) + geom_histogram() + scale_x_log10()) %>% 
    wrap_plots()
}

create_proportion_matrix <- function(frequency_table) {
  target_proportion <- prop.table(frequency_table, margin = 1)
  
  cluster_proportion <- prop.table(frequency_table, margin = 2)
  
  list('cluster' = as.matrix(cluster_proportion),
       'target' = as.matrix(target_proportion))
}

proportion_analysis <- function(metadata_table,
                                group1 = guide_target,
                                group2 = seurat_clusters,
                                rename = TRUE,
                                heatmap_names = FALSE) {
  
  condition <- metadata_table$condition %>% unique()
  object <- metadata_table$object %>% unique()
  
  freq_table <- metadata_table %>% 
    tabulate_frequency(group1 = {{group1}},
                       group2 = {{group2}},
                       rename = rename)
  
  prop_table <- freq_table %>% calculate_proportions()
  prop_matrix <- freq_table %>% create_proportion_matrix()
  prop_hist <- prop_table %>% proportion_histograms()
  
  cluster_heatmap <- pheatmap(t(prop_matrix$cluster), 
                              scale = 'column', 
                              clustering_distance_rows = 'correlation', 
                              clustering_distance_cols = 'correlation',
                              main = paste(condition, object, 'cluster'),
                              show_colnames = heatmap_names)
  
  target_heatmap <- pheatmap(t(prop_matrix$target), 
                             scale = 'column', 
                             clustering_distance_rows = 'correlation', 
                             clustering_distance_cols = 'correlation',
                             main = paste(condition, object, 'target'),
                             show_colnames = heatmap_names)
  
  heatmaps <- list('cluster' = cluster_heatmap,
                   'target' = target_heatmap)
  
  list('frequency' = freq_table,
       'proportions' = prop_table,
       'matrix' = prop_matrix,
       'histogram' = prop_hist,
       'heatmap' = heatmaps)
}

cluster_enrichment <- function(metadata_table,
                               group1 = guide_target,
                               group2 = seurat_clusters,
                               rename = TRUE,
                               print = FALSE,
                               method = 'Fisher') {
  
  frequency_table <- metadata_table %>% 
    dplyr::select(barcode,
                  {{group1}}, 
                  {{group2}})
  
  if(rename) {
    colnames(frequency_table) <- c('barcode', 'target', 'cluster')
  }
  
  targets <- frequency_table$target %>% unique() %>% sort()
  clusters <- frequency_table$cluster %>% as.character() %>% unique() %>% sort()
  all_barcodes <- frequency_table$barcode
  
  enrich <- map(targets, function(i) {
    cat(red(paste('Running', i, '\n')))
    
    map(clusters, function(j) {
      
      target_subset <- frequency_table %>% filter(target == i) %>% pull(barcode) 
      cluster_subset <- frequency_table %>% filter(cluster == j) %>% pull(barcode)
      overlap_subset <- intersect(target_subset, cluster_subset)
      
      # run enrichment test
      enrichment_res <- enrichment_test(target_subset,
                                        cluster_subset,
                                        background = all_barcodes,
                                        print = print)
      
      odds <- enrichment_res$estimate
      pval <- enrichment_res$p.value
      
      output <- tibble('target' = i,
                       'cluster' = j,
                       'n_overlap' = length(overlap_subset),
                       'n_target' = length(target_subset),
                       'n_cluster' = length(cluster_subset),
                       'f_target' = n_overlap/n_target,
                       'f_cluster' = n_overlap/n_cluster,
                       'odds' = enrichment_res$fisher$estimate,
                       'log2odds' = log2(odds),
                       'pval' = enrichment_res$fisher$p.value)
      output
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  enrich %>% 
    arrange(pval)
}

normalized_cluster_enrichment <- function(metadata_table,
                                          group1 = guide_target,
                                          group2 = seurat_clusters,
                                          rename = TRUE,
                                          print = FALSE,
                                          method = 'Fisher') {
  
  frequency_table <- metadata_table %>% 
    dplyr::select(barcode,
                  {{group1}}, 
                  {{group2}})
  
  if(rename) {
    colnames(frequency_table) <- c('barcode', 'target', 'cluster')
  }
  
  targets <- frequency_table$target %>% unique() %>% sort() %>% setdiff('Non-Targeting')
  clusters <- frequency_table$cluster %>% as.character() %>% unique() %>% sort()
  
  enrich <- map(targets, function(i) {
    cat(red(paste('Running', i, '\n')))
    
    map(clusters, function(j) {
      
      frequency_subset <- frequency_table %>% filter(target %in% c(i, 'Non-Targeting'))
      target_subset <- frequency_subset %>% filter(target == i) %>% pull(barcode) 
      cluster_subset <- frequency_subset %>% filter(cluster == j) %>% pull(barcode)
      overlap_subset <- intersect(target_subset, cluster_subset)
      barcodes_subset <- frequency_subset$barcode
      
      # run enrichment test
      enrichment_res <- enrichment_test(target_subset,
                                        cluster_subset,
                                        background = barcodes_subset,
                                        print = print)
      
      odds <- enrichment_res$estimate
      pval <- enrichment_res$p.value
      
      output <- tibble('target' = i,
                       'cluster' = j,
                       'n_overlap' = length(overlap_subset),
                       'n_target' = length(target_subset),
                       'n_cluster' = length(cluster_subset),
                       'f_target' = n_overlap/n_target,
                       'f_cluster' = n_overlap/n_cluster,
                       'odds' = enrichment_res$fisher$estimate,
                       'log2odds' = log2(odds),
                       'pval' = enrichment_res$fisher$p.value)
      output
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  enrich %>% 
    arrange(pval)
}