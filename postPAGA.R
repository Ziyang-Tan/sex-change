library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(tibble)
library(wesanderson)

figure_dir <- '/Users/tan/sex-change/figures_2022'
sub_names = dir(figure_dir)
for (sub_name in sub_names){
  # sub_name <- sub_names[1] # for debug
  file_path = file.path('PAGA_result_data_2022',
                        paste0(sub_name, '_sample5000.h5ad'))
  
  datah5 = readH5AD(file = file_path)
  data = t(assay(datah5)) %>% as_tibble()
  meta = colData(datah5) %>% 
    as_tibble() %>% 
    select(Visit, SubjectID, Subject, leiden)
  
  # subpop freq change
  data_freq = meta %>% 
    group_by(Visit, leiden) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(Visit) %>%
    mutate(freq = n / sum(n))
  
  g_list <- lapply(sort(unique(meta$leiden)), function(x){
    ggplot(data = data_freq %>% filter(leiden == x), 
           aes(x = Visit, y = freq, group = leiden)) +
      geom_line() +
      geom_point() +
      ylim(0, max(data_freq$freq) + 0.1) +
      labs(title = paste0(sub_name, '_cluster_', x))
  })
  ggarrange(plotlist = g_list, ncol = 4, nrow = 3) %>% 
    ggexport(filename = file.path(getwd(), 'figures_2022', sub_name, paste0(sub_name, '_cluster_changes_timepoint.pdf')),
             width = 20,
             height = 10)
  
  # subpop freq change by subjectID (individual)
  data_freq = meta %>% 
    group_by(Visit, leiden, SubjectID) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(Visit, SubjectID) %>%
    mutate(freq = n / sum(n))
  
  g_list <- lapply(sort(unique(meta$leiden)), function(x){
    ggplot(data = data_freq %>% filter(leiden == x), 
           aes(x = Visit, y = freq, group = SubjectID, color=SubjectID)) +
      geom_line() +
      geom_point() +
      ylim(0, max(data_freq$freq) + 0.1) +
      geom_smooth(aes(x = as.numeric(Visit), y = freq, group = leiden), method='loess') +
      labs(title = paste0(sub_name, '_cluster_', x))
  })
  ggarrange(plotlist = g_list, ncol = 4, nrow = 3) %>% 
    ggexport(filename = file.path(getwd(), 'figures_2022', sub_name, paste0(sub_name, '_cluster_changes_timepoint_by_SubjectID.pdf')),
             width = 20,
             height = 10)
  
  
  ################################
  
  # subpop marker expression
  data_anno <- data %>% add_column(leiden = meta$leiden)
  
  g_list <- lapply(colnames(data), function(x){
    d <- data_anno %>% select(x, leiden) %>% rename(marker_expression = x)
    ggplot(d, 
           aes(x = marker_expression, y = leiden, fill = leiden)) +
      geom_density_ridges(scale = 4, rel_min_height=.01) +
      scale_fill_manual(values = wes_palette("Rushmore1", nlevels(meta$leiden), type = "continuous")) +
      theme_ridges() + 
      theme(legend.position = "none") +
      xlim(quantile(d$marker_expression, 0.01),
           quantile(d$marker_expression, 0.99)) +
      labs(title = x)
  })
  ggarrange(plotlist = g_list, ncol = 5, nrow = 9) %>% 
    ggexport(filename = file.path(getwd(), 'figures_2022', sub_name, paste0(sub_name, '_cluster_marker_expressions.pdf')),
             width = 20,
             height = 40)
}
