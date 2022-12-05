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
#' @param colors 
#'
#' @return
#' @export
#'
#' @examples
view_colors <- function(colors) {
  barplot(rep(1, length(colors)), col = colors, names.arg = colors)
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

#' Title
#'
#' @param group1 
#' @param group2 
#' @param overlap 
#' @param total 
#' @param clean 
#'
#' @return
#' @export
#'
#' @examples
odds_quick <- function(group1,
                       group2,
                       overlap,
                       total,
                       clean = TRUE) {
  
  if(length(group1) + length(group2) + length(overlap) == 3) {
    
    a <- overlap
    b <- group1 - overlap
    c <- group2 - overlap
    d <- total - group1 - c
    
  }
  
  contingency <- matrix(c(a, b,
                          c, d),
                        nrow = 2)
  
  fisher <- fisher.test(contingency)
  
  if(clean) {
    tibble('odds' = round(fisher$estimate, 2),
           'log2odds' = round(log2(odds), 2),
           'pval' = fisher$p.value)
    
  } else {
    
    tibble('n_group1', 
           'n_group2',
           'n_overlap',
           'n_total',
           'odds' = round(as.numeric(fisher$estimate), 2),
           'log2odds' = round(log2(odds), 2),
           'pval' = fisher$p.value)
  }
}

#' Title
#'
#' @param df 
#' @param group1 
#' @param group2 
#'
#' @return
#' @export
#'
#' @examples
odds_table <- function(df,
                       group1,
                       group2) {
  
  df %>% 
    select(all_of(c(group1, group2))) %>% 
    set_names(c('group1', 'group2')) %>% 
    add_count(name = 'n_total') %>% 
    group_by(group1) %>% 
    add_count(name = 'n_group1') %>% 
    group_by(group2) %>% 
    add_count(name = 'n_group2') %>% 
    group_by_all() %>% 
    tally(name = 'n_overlap') %>% 
    rowwise() %>% 
    mutate(
      f_group1 = n_overlap/n_group1,
      f_group2 = n_overlap/n_group2,
      f_max = max(f_group1, f_group2),
      f_score = f_group1 * f_group2,
      jaccard = n_overlap/(n_group1 + n_group2 - n_overlap),
      overlap_coef = n_overlap/min(n_group1, n_group2),
      odds_quick(group1 = n_group1,
                 group2 = n_group2,
                 overlap = n_overlap,
                 total = n_total,
                 clean = TRUE)) %>% 
    ungroup() %>% 
    mutate(across(c(starts_with('f_'), jaccard, overlap_coef), round, 3)) %>% 
    as_tibble()
}

#' Title
#'
#' @param df 
#' @param nodes 
#' @param names 
#' @param reverse 
#' @param total 
#' @param prefix 
#' @param arrange 
#' @param title 
#' @param drop_na 
#'
#' @return
#' @export
#'
#' @examples
plot_sankey <- function(df,
                        nodes,
                        names = NULL,
                        reverse = TRUE,
                        total = TRUE,
                        prefix = 'Interactions:',
                        arrange = 1,
                        title = NULL,
                        drop_na = FALSE) {
  
  sankey_input <- df %>% select(all_of(nodes)) %>% mutate(across(.fns = as.character)) %>% arrange(across(nodes[arrange]))
  
  if(drop_na) {
    sankey_input <- sankey_input %>% drop_na()
  }
  
  if(is.null(names)) {
    names <- nodes
  } else {
    colnames(sankey_input) <- names
  }
  
  if(reverse) {
    sankey_input <- sankey_input %>% make_long(length(nodes):1)
  } else {
    sankey_input <- sankey_input %>% make_long(1:length(nodes))
  }
  
  
  if(total) {
    total <- paste(prefix, nrow(df))
  } else {
    total <- NULL
  }
  
  # Construct plot
  p <- ggplot(sankey_input %>% mutate(node = factor(node)), 
              aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = node,
                  label = node)) +
    geom_sankey(flow.alpha = 0.6,
                flow.color = 'black',
                node.color = 'black') +
    theme_sankey(base_size = 10, base_family = 'Arial') +
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 10, color = 'black'))+
    geom_sankey_label(size = 3, color = "black", fill = "white", family = 'Arial') + 
    scale_fill_viridis_d() + 
    labs(x = '',
         title = title,
         subtitle = total)
  
}

### sankey
#' Title
#'
#' @param df 
#' @param nodes 
#' @param names 
#' @param reverse 
#' @param top_n 
#' @param percentages 
#' @param total 
#' @param prefix 
#' @param title 
#'
#' @return
#' @export
#'
#' @examples
plot_sankey1 <- function(df,
                         nodes,
                         names = NULL,
                         reverse = TRUE,
                         top_n = 3, # only applicable for single category breakdown
                         percentages = TRUE,
                         total = TRUE,
                         prefix = 'Interactions:',
                         title = NULL,
                         drop_na = FALSE,
                         ordered = FALSE) {
  
  sankey_input <- df %>% select(all_of(nodes)) %>% mutate(across(.fns = as.character))
  
  if(drop_na) {
    sankey_input <- sankey_input %>% drop_na()
  }
  
  if(is.null(names)) {
    names <- nodes
  } else {
    colnames(sankey_input) <- names
  }
  
  if(!is.null(top_n) | percentages) {
    
    sankey_input <- sankey_input %>% 
      add_count(across(all_of(names))) %>% 
      mutate(pct = round(100 * n/nrow(.), 1))
    
  }
  
  if(!is.null(top_n)) {
    top_categories <- sankey_input %>% unique() %>% arrange(-pct) %>% mutate(row = row_number())
    
    bottom_pct <- top_categories %>% 
      filter(row > top_n) %>% 
      pull(pct) %>% 
      sum()
    
    sankey_input <- sankey_input %>% 
      left_join(top_categories %>% select(-n)) %>% 
      mutate(pct = ifelse(row > top_n, bottom_pct, pct),
             '{names[2]}' := ifelse(row > top_n, 'Other', .[[names[2]]])) 
    
    if(ordered) {
      
      sankey_input <- sankey_input %>% arrange(pct)
    }
    
  }
  
  if(percentages) {
    sankey_input <- sankey_input %>% mutate('{names[2]}' := paste0(.[[names[2]]], ', ', pct, '%')) 
  }
  
  plot_order <- map(names, function(i) {
    unique(sankey_input[[i]])
  }) %>% unlist() %>% unique()
  
  if(total) {
    total <- paste(prefix, nrow(df))
  } else {
    total <- NULL
  }
  
  if(reverse) {
    sankey_input <- sankey_input %>% make_long(2:1)
  } else {
    sankey_input <- sankey_input %>% make_long(1:2)
  }
  
  
  # Construct plot
  p <- ggplot(sankey_input %>% mutate(node = factor(node, levels = plot_order)), 
              aes(x = x, 
                  next_x = next_x, 
                  node = node, 
                  next_node = next_node,
                  fill = node,
                  label = node)) +
    geom_sankey(flow.alpha = 0.6,
                flow.color = 'black',
                node.color = 'black') +
    theme_sankey(base_size = 10, base_family = 'Arial') +
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 10, color = 'black'))+
    geom_sankey_label(size = 3, color = "black", fill = "white", family = 'Arial') + 
    scale_fill_viridis_d() + 
    labs(x = '',
         title = title,
         subtitle = total)
  
}

#' Title
#'
#' @param df 
#' @param nodes 
#' @param names 
#' @param top_n 
#' @param percentages 
#' @param total 
#' @param prefix 
#' @param title 
#'
#' @return
#' @export
#'
#' @examples
plot_sankey2 <- function(df,
                         nodes,
                         names = NULL,
                         top_n = 3, # only applicable for single category breakdown
                         percentages = TRUE,
                         total = TRUE,
                         prefix = 'Interactions:',
                         title = NULL,
                         drop_na = FALSE)
{
  if(drop_na) {
    sankey_input <- sankey_input %>% drop_na()
  }
  
  sankey_input <- df %>% select(all_of(nodes)) %>% mutate(across(.fns = as.character))
  
  if(is.null(names)) {
    names <- nodes
  } else {
    colnames(sankey_input) <- names
  }
  
  if(!is.null(top_n)) {
    
    sankey_input_left <- sankey_input %>% 
      select(all_of(names[-3])) %>% 
      group_by_all() %>% 
      tally() %>% 
      ungroup() %>% 
      mutate(pct = round(100 * n/nrow(sankey_input), 1))
    
    sankey_input_right <- sankey_input %>% 
      select(all_of(names[-1])) %>% 
      group_by_all() %>% 
      tally() %>% 
      ungroup() %>% 
      mutate(pct = round(100 * n/nrow(sankey_input), 1)) 
    
  }
  
  if(!is.null(top_n)) {
    
    top_categories_left <- sankey_input_left %>% unique() %>% arrange(-pct) %>% mutate(row = row_number())
    
    bottom_pct_left <- top_categories_left %>% 
      filter(row > top_n) %>% 
      pull(pct) %>% 
      sum()
    
    top_categories_right <- sankey_input_right %>% unique() %>% arrange(-pct) %>% mutate(row = row_number())
    
    bottom_pct_right <- top_categories_right %>% 
      filter(row > top_n) %>% 
      pull(pct) %>% 
      sum()
    
    sankey_input <- sankey_input %>% 
      left_join(top_categories_left %>% select(-n)) %>% 
      mutate(pct = ifelse(row > top_n, bottom_pct_left, pct),
             '{names[1]}' := ifelse(row > top_n, 'Other', .[[names[1]]])) %>% 
      mutate('{names[1]}' := paste0(.[[names[1]]], ', ', pct, '%')) %>% 
      arrange(pct) %>% 
      select(-row, -pct) %>% 
      left_join(top_categories_right %>% select(-n)) %>% 
      mutate(pct = ifelse(row > top_n, bottom_pct_right, pct),
             '{names[3]}' := ifelse(row > top_n, 'Other', .[[names[3]]])) %>% 
      mutate('{names[3]}' := paste0(.[[names[3]]], ', ', pct, '%')) %>% 
      select(-row, -pct) 
    
    # plot_order <- map(1:ncol(sankey_input), function(i) {
    #   unique(sankey_input[[i]])
    # }) %>% unlist()
    
    
  }
  
  if(total) {
    total <- paste(prefix, nrow(df))
  } else {
    total <- NULL
  }
  
  # Construct plot
  ggplot(sankey_input %>% make_long(1:3), 
         aes(x = x, 
             next_x = next_x, 
             node = node, 
             next_node = next_node,
             fill = node,
             label = node)) +
    geom_sankey(flow.alpha = 0.6,
                flow.color = 'black',
                node.color = 'black') +
    theme_sankey(base_size = 10, base_family = 'Arial') +
    theme(legend.position = 'none',
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 10, color = 'black'))+
    geom_sankey_label(size = 3, color = "black", fill = "white", family = 'Arial') + 
    scale_fill_viridis_d() + 
    labs(x = '',
         title = title,
         subtitle = total)
}


#' Title
#'
#' @param group1 
#' @param group2 
#'
#' @return
#' @export
#'
#' @examples
calculate_agreement <- function(group1,
                                group2) {
  
  tibble('AMI' = aricode::AMI(group1, group2),
         'ARI' = aricode::ARI(group1, group2))
  
}