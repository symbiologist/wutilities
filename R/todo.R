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