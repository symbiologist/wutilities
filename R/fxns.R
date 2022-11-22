#' cat-based printing
#'
#' @param ... 
#'
#' @export
#'
print2 <- function(...) {
  cat(crayon::magenta(paste0(..., '\n')))
}

#' Title
#'
#' @param n 
#' @param palette 
#' @param ... 
#'
#' @export
#'
generate_colors <- function(n, 
                            palette = 'default', 
                            ...) {
  
  colors <- dplyr::case_when(
    palette == 'default' ~ c(0, 360, 0, 180, 0, 100),
    palette == 'pastel' ~ c(0, 360, 0, 54, 67, 100),
    palette == 'pimp' ~ c(0, 360, 54, 180, 27, 67),
    palette == 'intense' ~ c(0, 360, 36, 180, 13, 73),
    palette == 'custom1' ~ c(0, 360, 30, 100, 25, 100),
    palette == 'custom2' ~ c(0, 360, 30, 70, 60, 100))
  
  unname(hues::iwanthue(n = n,
                        hmin = colors[1],
                        hmax = colors[2],
                        cmin = colors[3],
                        cmax = colors[4],
                        lmin = colors[5],
                        lmax = colors[6],
                        ...))
}

#' Enrichment testing
#'
#' @param list1 
#' @param list2 
#' @param cat1 
#' @param cat2 
#' @param background 
#' @param print 
#'
#' @export
#' @import tidyverse
enrichment_test <- function(list1, # list of positives in category 1
                            list2, # list of positives in category 2
                            cat1 = c(TRUE, FALSE), # category 1 names
                            cat2 = c(TRUE, FALSE), # category 2 names
                            background, # total background list of all elements
                            print = TRUE,
                            method = 'fisher') {
  
  # Make sure lists are within the full background set
  list1 <- intersect(list1, background)
  list2 <- intersect(list2, background)
  full_list <- unique(c(list1, list2, background))
  
  # contingency table
  df <- tibble(element = full_list) %>% 
    mutate(list1 = ifelse(element %in% list1, cat1[1], cat1[2]) %>% factor(levels = cat1), # assign list1 to category 1 and set as ordered factor
           list2 = ifelse(element %in% list2, cat2[1], cat2[2]) %>% factor(levels = cat2)) # assign list2 to category 2 and set as ordered factor
  
  contingency_table <- table(df$list1, df$list2)
  
  # Test
  if('fisher' %in% method) {
    fisher <- fisher.test(contingency_table)
  } 
  if('chisq' %in% method) {
    chisq <- chisq.test(contingency_table)
  }
  
  if(print) {
    print(contingency_table)
    if('fisher' %in% method) {
      print(fisher)  
    } 
    if('chisq' %in% method) {
      print(chisq)
    }
  }
  
  # Return fisher results as well as table
  output <- list('contingency' = contingency_table,
                 'table' = df)
  
  if('fisher' %in% method) {
    output[['fisher']] <- fisher
  } 
  if('chisq' %in% method) {
    output[['chisq']] <- chisq
  }
  
  output
}

#' Title
#'
#' @param input 
#' @param filters 
#' @param verbose
#' @export
#'
#'
apply_filters <- function(input, 
                          filters = list('column' = c('value1', 'value2')),
                          verbose = FALSE) {
  
  filter_vars <- names(filters)
  
  if(all(filter_vars %in% colnames(input))) {
    if(verbose) {
      print2('All variables present in input table') 
    }
  } else {
    missing_var <- setdiff(filter_vars, colnames(input))
    print2('Missing variables! Please check input for these variables: ', missing_var)
    stop()
  }
  
  index_to_retain <- map(filter_vars, function(i) {
    input[[i]] %in% filters[[i]]
  }) %>% Reduce("&", .)
  
  input[index_to_retain, ]
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

### sankey