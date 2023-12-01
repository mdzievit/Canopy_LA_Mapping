library(tidyverse)


data <- read_tsv("QTL_Analysis/QTL_Additive_Classes.txt")


(out <- data %>%
    mutate(`QTL Class` = factor(`QTL Class`),
           QTL = paste0(QTL," (",Source,")"),
           QTL = factor(QTL) %>% fct_reorder(Order)) %>% 
    ggplot(aes(x = Coded,y = Additive,group = QTL)) +
    geom_point(aes(color = QTL,
                   shape = Significant),
               size = 3) +
    geom_line(aes(linetype = `QTL Class`,
                  color = QTL),
              linewidth = 1) +
    theme_bw() +
    scale_x_continuous(breaks = c(-2,0,1,3,5),
                    labels = c("Two below ear",
                               "Ear leaf",
                               "One above ear",
                               "Three above ear",
                               "One below flag"),
                    name = "Canopy position") +
    # scale_x_reverse(breaks = c(-2,0,1,3,5),
    #                    labels = c("Two below ear",
    #                               "Ear leaf",
    #                               "One above ear",
    #                               "Three above ear",
    #                               "One below flag"),
    #                    name = "Canopy position") +
    coord_flip() +
    ylab("Additive effect") +
    scale_y_continuous(limits = c(-4.5,4.5),n.breaks = 10) +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size = 8),
          legend.title = element_text(size = 9)))


ggsave(plot = out,
       filename = "QTL_Analysis/Fig5-LA_QTL_Classes.png",
       width = 7,
       height = 6,
       dpi = 1200)
