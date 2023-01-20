# clear space
rm(list=ls()) 

# clear graphics
graphics.off() 

library(tidyverse)
library(cowplot)
library(ggpubr)
library(ggpattern)

filepath_project_group = "~/Desktop/COVID_AMR_GitHub/"

# upload different scenarios
source(paste0(filepath_project_group, "S0_no_cov.R"))
source(paste0(filepath_project_group, "S0_Snc.R"))
source(paste0(filepath_project_group, "S1_SAi.R"))
source(paste0(filepath_project_group, "S2_Sad.R"))
source(paste0(filepath_project_group, "S3_SL.R"))
source(paste0(filepath_project_group, "S4_SAiL.R"))
source(paste0(filepath_project_group, "S5_SadL.R"))


panel_A <- plot_grid(plotS0_a , plotS1_a, plotS2_a, plotS3_a, plotS4_a, plotS5_a,
                  labels = c(""),ncol = 6, nrow = 1, label_size = 24)
panA <- annotate_figure(panel_A,
                top = text_grob("Bacterial carriage and SARS-CoV-2 outbreak", color = "black", face = "bold", size = 28),
                bottom = text_grob("Time (days)", color = "black", size = 18),
                left = text_grob("Prevalance (proportion)", color = "black", rot = 90, size = 18),
)



panel_B <- plot_grid(plotS0_b , plotS1_b, plotS2_b, plotS3_b, plotS4_b, plotS5_b,
                  labels = c(""), ncol = 6, nrow = 1, label_size = 24)
panB <- annotate_figure(panel_B,
                top = text_grob("Antibiotic resistance rate", color = "black", face = "bold", size = 28),
                bottom = text_grob("Time (days)", color = "black", size = 18),
                left = text_grob("Prevalance (proportion)", color = "black", rot = 90, size = 18),
)



##### Data frame for panel C - incidence of the IBD total across all 6 scenarios compared to the baseline
total_short <- c(total_short_nc, total_short_S0, total_short_S1, total_short_S2, total_short_S3, total_short_S4, total_short_S5)
total_long <- c(total_long_nc, total_long_S0, total_long_S1, total_long_S2, total_long_S3, total_long_S4, total_long_S5)

prop_change <- c(0, round((total_short_S0-total_short_nc)/total_short_nc, digits = 3), 
                 round((total_short_S1-total_short_nc)/total_short_nc, digits = 3), 
                 round((total_short_S2-total_short_nc)/total_short_nc, digits = 3), 
                 round((total_short_S3-total_short_nc)/total_short_nc, digits = 3), 
                 round((total_short_S4-total_short_nc)/total_short_nc, digits = 3), 
                 round((total_short_S5-total_short_nc)/total_short_nc, digits = 3), 
                 0, round((total_long_S0-total_long_nc)/total_long_nc, digits = 3), 
                 round((total_long_S1-total_long_nc)/total_long_nc, digits = 3), 
                 round((total_long_S2-total_long_nc)/total_long_nc, digits = 3), 
                 round((total_long_S3-total_long_nc)/total_long_nc, digits = 3), 
                 round((total_long_S4-total_long_nc)/total_long_nc, digits = 3), 
                 round((total_long_S5-total_long_nc)/total_long_nc, digits = 3))

time_frame <- c('1short','1short','1short','1short','1short','1short','1short',
                '2long','2long','2long','2long', '2long','2long', '2long')

incidence <- c(total_short_nc, total_short_S0, total_short_S1, total_short_S2, total_short_S3, total_short_S4, total_short_S5,
               total_long_nc, total_long_S0, total_long_S1, total_long_S2, total_long_S3, total_long_S4, total_long_S5) # incidence of IBD added

scenario <- c('no_covid', 'scen_0', 'scen_1', 'scen_2', 'scen_3', 'scen_4', 'scen_5',
              'no_covid', 'scen_0', 'scen_1', 'scen_2', 'scen_3', 'scen_4', 'scen_5')

status <- c('total', 'total', 'total', 'total', 'total', 'total','total', 'total', 'total', 'total', 'total', 'total', 'total', 'total')

ipdat_total <- data.frame(status, incidence, prop_change, time_frame, scenario)  

ipdat_total$percent_change <- 100*ipdat_total$prop_change
ipdat_total$time_frame <- as.factor(ipdat_total$time_frame)
ipdat_total$scenario <- as.factor(ipdat_total$scenario)


ipdat_total <- ipdat_total %>% 
  ## add percentage label with `sprintf()`
  dplyr::mutate(perc = paste0(sprintf("%4.1f", percent_change), "%"))

# add to erase % in the baseline column
ipdat_total$new <- ifelse(ipdat_total$percent_change == "0",
                       ipdat_total$perc,ipdat_total$perc)

# name short and long term facet plots
names <- c(
  `1short` = "90-DAY",
  `2long` = "ANNUAL"
)

scales_y <- list(
  `1short` = scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1)),
  `2long` = scale_y_continuous(limits = c(0, 14), breaks = seq(0, 14, 2))
)

ipdat_total    # IBD incidence data frame
str(ipdat_total)





########## TOTAL incidence
TOTAL <- ggplot(data=ipdat_total, aes(x=scenario, y=incidence, group = status)) +
  geom_bar(stat="identity", color = "black", position = position_dodge(), fill = "white") +
  theme_classic() +
  labs(x="", y="") +
  #ggtitle("Relative change in total cumulative invasive bacterial disease incidence (%)") +
  theme(plot.title = element_text(hjust = 0.5), # center plot title
        legend.position = "bottom",
        axis.text.x = element_text(face="bold", size=16),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        strip.text = element_text(size=14)) +  # facet_wrap title size
  geom_text(aes(label = new), vjust = -0.5, position = position_dodge(1.2), size=8, fontface = "bold") +
  #facet_wrap(. ~ time_frame, ncol = 1, labeller = as_labeller(names), scales = "free_x")
  facet_grid_sc(rows = vars(time_frame), labeller = as_labeller(names), scales = list(y = scales_y)) + 
  scale_x_discrete(labels=c("no_covid"="Pre-pandemic", 
                            "scen_0" = "S0", 
                            "scen_1" = "S1", 
                            "scen_2" = "S2",
                            "scen_3" = "S3", 
                            "scen_4" = "S4", 
                            "scen_5" = "S5")) 
  






##### Data frame for plot D - incidence of the IBD S across all 6 scenarios compared to the baseline
S_short <- c(S_short_nc, S_short_S0, S_short_S1, S_short_S2, S_short_S3, S_short_S4, S_short_S5)
S_long <- c(S_long_nc, S_long_S0, S_long_S1, S_long_S2, S_long_S3, S_long_S4, S_long_S5)


prop_change <- c(0, round((S_short_S0-S_short_nc)/S_short_nc, digits = 3), 
                 round((S_short_S1-S_short_nc)/S_short_nc, digits = 3), 
                 round((S_short_S2-S_short_nc)/S_short_nc, digits = 3), 
                 round((S_short_S3-S_short_nc)/S_short_nc, digits = 3), 
                 round((S_short_S4-S_short_nc)/S_short_nc, digits = 3), 
                 round((S_short_S5-S_short_nc)/S_short_nc, digits = 3), 
                 0, round((S_long_S0-S_long_nc)/S_long_nc, digits = 3), 
                 round((S_long_S1-S_long_nc)/S_long_nc, digits = 3), 
                 round((S_long_S2-S_long_nc)/S_long_nc, digits = 3), 
                 round((S_long_S3-S_long_nc)/S_long_nc, digits = 3), 
                 round((S_long_S4-S_long_nc)/S_long_nc, digits = 3), 
                 round((S_long_S5-S_long_nc)/S_long_nc, digits = 3))

time_frame <- c('1short','1short','1short','1short','1short','1short','1short',
                '2long','2long','2long','2long', '2long','2long', '2long')

incidence <- c(S_short_nc, S_short_S0, S_short_S1, S_short_S2, S_short_S3, S_short_S4, S_short_S5,
               S_long_nc, S_long_S0, S_long_S1, S_long_S2, S_long_S3, S_long_S4, S_long_S5) # incidence of IBD added

scenario <- c('no_covid', 'scen_0', 'scen_1', 'scen_2', 'scen_3', 'scen_4', 'scen_5',
              'no_covid', 'scen_0', 'scen_1', 'scen_2', 'scen_3', 'scen_4', 'scen_5')

status <- c('S', 'S', 'S', 'S', 'S', 'S','S', 'S', 'S', 'S', 'S', 'S', 'S', 'S')

ipdat_S <- data.frame(status, incidence, prop_change, time_frame, scenario)  

ipdat_S$percent_change <- 100*ipdat_S$prop_change
ipdat_S$time_frame <- as.factor(ipdat_S$time_frame)
ipdat_S$scenario <- as.factor(ipdat_S$scenario)


ipdat_S <- ipdat_S %>% 
  ## add percentage label with `sprintf()`
  dplyr::mutate(perc = paste0(sprintf("%4.1f", percent_change), "%"))

# add to erase % in the baseline column
ipdat_S$new <- ifelse(ipdat_S$percent_change == "0",
                          ipdat_S$perc,ipdat_S$perc)




########## S IPD incidence
S <- ggplot(data=ipdat_S, aes(x=scenario, y=incidence, group = status)) +
  geom_bar(stat="identity", color = "black", position = position_dodge(), fill = "white") +
  theme_classic() +
  labs(x="", y="") +
  #ggtitle("Relative change in total cumulative invasive bacterial disease incidence (%)") +
  theme(plot.title = element_text(hjust = 0.5), # center plot title
        legend.position = "bottom",
        axis.text.x = element_text(face="bold", size=16),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        strip.text = element_text(size=14)) +  # facet_wrap title size
  geom_text(aes(label = new), vjust = -0.5, position = position_dodge(1.2), size=8, fontface = "bold") +
  #facet_wrap(. ~ time_frame, ncol = 1, labeller = as_labeller(names), scales = "free_x")
  facet_grid_sc(rows = vars(time_frame), labeller = as_labeller(names), scales = list(y = scales_y)) + 
  scale_x_discrete(labels=c("no_covid"="Pre-pandemic", 
                            "scen_0" = "S0", 
                            "scen_1" = "S1", 
                            "scen_2" = "S2",
                            "scen_3" = "S3", 
                            "scen_4" = "S4", 
                            "scen_5" = "S5")) 



##### Data frame for plot E - incidence of the IBD R across all 6 scenarios compared to the baseline
R_short <- c(R_short_nc, R_short_S0, R_short_S1, R_short_S2, R_short_S3, R_short_S4, R_short_S5)
R_long <- c(R_long_nc, R_long_S0, R_long_S1, R_long_S2, R_long_S3, R_long_S4, R_long_S5)


prop_change <- c(0, round((R_short_S0-R_short_nc)/R_short_nc, digits = 3), 
                 round((R_short_S1-R_short_nc)/R_short_nc, digits = 3), 
                 round((R_short_S2-R_short_nc)/R_short_nc, digits = 3), 
                 round((R_short_S3-R_short_nc)/R_short_nc, digits = 3), 
                 round((R_short_S4-R_short_nc)/R_short_nc, digits = 3), 
                 round((R_short_S5-R_short_nc)/R_short_nc, digits = 3), 
                 0, round((R_long_S0-R_long_nc)/R_long_nc, digits = 3), 
                 round((R_long_S1-R_long_nc)/R_long_nc, digits = 3), 
                 round((R_long_S2-R_long_nc)/R_long_nc, digits = 3), 
                 round((R_long_S3-R_long_nc)/R_long_nc, digits = 3), 
                 round((R_long_S4-R_long_nc)/R_long_nc, digits = 3), 
                 round((R_long_S5-R_long_nc)/R_long_nc, digits = 3))

time_frame <- c('1short','1short','1short','1short','1short','1short','1short',
                '2long','2long','2long','2long', '2long','2long', '2long')

incidence <- c(R_short_nc, R_short_S0, R_short_S1, R_short_S2, R_short_S3, R_short_S4, R_short_S5,
               R_long_nc, R_long_S0, R_long_S1, R_long_S2, R_long_S3, R_long_S4, R_long_S5) # incidence of IBD added

scenario <- c('no_covid', 'scen_0', 'scen_1', 'scen_2', 'scen_3', 'scen_4', 'scen_5',
              'no_covid', 'scen_0', 'scen_1', 'scen_2', 'scen_3', 'scen_4', 'scen_5')

status <- c('R', 'R', 'R', 'R', 'R', 'R','R', 'R', 'R', 'R', 'R', 'R', 'R', 'R')

ipdat_R <- data.frame(status, incidence, prop_change, time_frame, scenario)  

ipdat_R$percent_change <- 100*ipdat_R$prop_change
ipdat_R$time_frame <- as.factor(ipdat_R$time_frame)
ipdat_R$scenario <- as.factor(ipdat_R$scenario)


ipdat_R <- ipdat_R %>% 
  ## add percentage label with `sprintf()`
  dplyr::mutate(perc = paste0(sprintf("%4.1f", percent_change), "%"))

# add to erase % in the baseline column
ipdat_R$new <- ifelse(ipdat_R$percent_change == "0",
                      ipdat_R$perc,ipdat_R$perc)




########## R IPD incidence
R <- ggplot(data=ipdat_R, aes(x=scenario, y=incidence, group = status, color = scenario)) +
  geom_bar(stat="identity", position = position_dodge(), fill = "white") +
  theme_classic() +
  labs(x="", y="") +
  #ggtitle("Relative change in total cumulative invasive bacterial disease incidence (%)") +
  theme(plot.title = element_text(hjust = 0.5), # center plot title
        legend.position = "none",
        axis.text.x = element_text(face="bold", size=16),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size=12),
        strip.text = element_text(size=14)) +  # facet_wrap title size
  geom_text(aes(label = new), vjust = -0.5, position = position_dodge(1.2), size=8, fontface = "bold") +
  #facet_wrap(. ~ time_frame, ncol = 1, labeller = as_labeller(names), scales = "free") + 
  facet_grid_sc(rows = vars(time_frame), labeller = as_labeller(names), scales = list(y = scales_y)) + 
  scale_x_discrete(labels=c("no_covid"="Pre-pandemic", 
                            "scen_0" = "S0", 
                            "scen_1" = "S1", 
                            "scen_2" = "S2",
                            "scen_3" = "S3", 
                            "scen_4" = "S4", 
                            "scen_5" = "S5")) +
  scale_color_manual(values=c("no_covid"="black", 
                              "scen_0" = "black", 
                              "scen_1" = "black", 
                              "scen_2" = "black",
                              "scen_3" = "black", 
                              "scen_4" = "black", 
                              "scen_5" = "black")) 


panel_C <- plot_grid(TOTAL, S, R, 
                  labels = c("", "", ""), ncol = 3, nrow = 1, label_size = 24)
panC <- annotate_figure(panel_C,
                top = text_grob("Change in the cumulative IBD incidence compared to the pre-pandemic period", color = "black", face = "bold", size = 28),
                bottom = text_grob("", color = "black", size = 18),
                left = text_grob("Cumulative IBD incidence (per 100,000)", color = "black", rot = 90, size = 18),
)


final_plot = plot_grid(panA, panB, panC,
                       nrow = 3, rel_heights = c(1,1,1), align = "v", axis = "lr",
                       labels = c('A', 'B', 'C'), label_size = 32)

ggsave(final_plot, filename = paste0(filepath_project_group, "Figure_3.jpeg"), width = 75, height = 60, unit = 'cm')
