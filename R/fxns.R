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
                            print = TRUE) {
  
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
  fisher <- fisher.test(contingency_table)
  chisq <- chisq.test(contingency_table)
  
  if(print) {
    print(contingency_table)
    print(fisher)  
    print(chisq)
  }
  
  # Return fisher results as well as table
  list('contingency' = contingency_table,
       'table' = df,
       'chisq' = chisq,
       'fisher' = fisher)
}

#' Title
#'
#' @param input 
#' @param filters 
#'
#' @export
#'
#'
apply_filters <- function(input, 
                          filters # list(column = c(value1, value2))
                          ) {
  
  filter_vars <- names(filters)
  
  if(all(filter_vars %in% colnames(input))) {
    print2('All variables present in input table')
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
