# Figure 4 A and B

# incidence

IPDr = as.data.frame(matrix(c(2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                              1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17, 1.17,
                              1.17, 1.17, 1.18, 1.19, 1.20, 1.21, 1.21, 1.23, 1.24, 1.26, 1.27,                     
                              1.17, 1.17, 1.18, 1.20, 1.21, 1.23, 1.24, 1.26, 1.29, 1.32, 1.35,
                              1.17, 1.17, 1.19, 1.20, 1.22, 1.24, 1.26, 1.29, 1.32, 1.36, 1.40,
                              1.17, 1.17, 1.19, 1.21, 1.23, 1.26, 1.28, 1.31, 1.35, 1.40, 1.44),
                              nrow = 11, ncol = 6))

#rownames(IPDr) <- c(2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0)
colnames(IPDr) <- c("R0", "0%", "05%", "10%", "15%", "20%")


IPDr_df <- IPDr %>% 
  pivot_longer(cols = "0%" : "20%",
               names_to = "azithromycin",
               values_to = "incidence")

IPDr_df






AR = as.data.frame(matrix(c(2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                            19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4, 19.4,
                            19.5, 19.6, 19.9, 20.4, 20.7, 20.9, 21.0, 21.0, 21.3, 21.8, 22.2,                    
                            19.5, 19.7, 20.3, 20.9, 21.5, 21.8, 22.0, 22.1, 22.6, 23.5, 24.2,
                            19.5, 19.7, 20.5, 21.3, 22.0, 22.4, 22.6, 22.8, 23.6, 24.7, 25.7,
                            19.5, 19.8, 20.6, 21.6, 22.3, 22.8, 23.1, 23.4, 24.2, 25.6, 26.8),
                          nrow = 11, ncol = 6))

#rownames(AR) <- c(2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0)
colnames(AR) <- c("R0", "0%", "05%", "10%", "15%", "20%")


AR_df <- AR %>% 
  pivot_longer(cols = "0%" : "20%",
               names_to = "azithromycin",
               values_to = "resistance")

AR_df



p_SA_model0 = IPDr_df %>%
  ggplot(aes(x = R0, y = azithromycin, fill = incidence))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = 'forestgreen',  mid = 'gold', high = '#e66101',  # #e66101 #5e3c99
                       limits = c(min(IPDr_df$incidence),
                                  max(IPDr_df$incidence)),
                       midpoint = (max(IPDr_df$incidence)+min(IPDr_df$incidence))/2, 
                       name = 'Annual\nantibiotic-resistant\nIPD incidence\nper 100,000\ninhabitants') +
  theme_classic() +
  theme(legend.position="right",
        legend.title.align = 0.5,
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14)) +
  xlab(expression(paste(R[0], " (SARS-CoV-2)")))+ylab("Percentage of COVID-19 infected\ntaking azithromycin")+
  scale_x_continuous(breaks = c(2,2.2,2.4,2.6,2.8,3,3.2,3.4,3.6,3.8,4.0)) +
  geom_text(aes(label = round((AR_df$resistance), digits=1)), color = "white", fontface='bold', size=6) #, fontface='bold'



heatmap <- p_SA_model0
heatmap ############# part A







library(tidyverse)

### Run all scenarios
# no lockdown implemented, just covid with different prophylactic use

filepath_project_group = "INSERT PATHWAY"

source(paste0(filepath_project_group, "1_GenPop_Baseline.R"))
source(paste0(filepath_project_group, "2_GenPop_0.R"))
source(paste0(filepath_project_group, "3_GenPop_5.R"))
source(paste0(filepath_project_group, "4_GenPop_10.R"))
source(paste0(filepath_project_group, "5_GenPop_15.R"))
source(paste0(filepath_project_group, "6_GenPop_20.R"))

### SAVE PLOT
#ggsave("heatmap.jpeg", width = 25, height = 10, units = "cm")


############################## PART B
# Single scenario S19, but with 5 different levels of prophylactic use (0-20%):
outputs = rbind(output_0%>%mutate(scenario = "0%")%>%dplyr::select(-c("p_abx_0", "beta_0")),
                output_5%>%mutate(scenario = "05%")%>%dplyr::select(-c("p_abx_5", "beta_5")),
                output_10%>%mutate(scenario = "10%")%>%dplyr::select(-c("p_abx_10", "beta_10")),
                output_15%>%mutate(scenario = "15%")%>%dplyr::select(-c("p_abx_15", "beta_15")),
                output_20%>%mutate(scenario = "20%")%>%dplyr::select(-c("p_abx_20", "beta_20")))%>%
  filter(time > 1000 & time < 1365)%>%    # 1 year, 364 days each if 1365 or 365 days is 1366 1000 and 1365 old
  mutate(time = time - 1000)


# write 2 functions for pwinter and psummer plus another for invasion risk to rotate 1 (2.5m), 0.2 (2m), 0.4 (7.5m)
times3 = seq(1:1820)


p <- as.data.frame(list(times = times3, pinv = times3)) 

# assign winter invasion rate
p$pinv <- ifelse((p$times > 89 & p$times < 271 | p$times > 453 & p$times < 635 |
                    p$times > 817 & p$times < 999 | p$times > 1181 & p$times < 1363 | 
                    p$times > 1545 & p$times < 1727), (psummer_0), (pwinter_0)) # winter invasion rate is 3x the summer rate (de Celles PNAS)


rate_IPD = p$pinv




risk_ncp <- as.data.frame(list(times = times3, risk = times3))

p$risk <- ifelse(p$times < 75 | p$times > 364 & p$times < 439 | p$times > 728 & p$times < 803 |
                   p$times > 1092 & p$times < 1167 | p$times > 1456 & p$times < 1531, 1,
                 ifelse(p$times > 74 & p$times < 136 | p$times > 438 & p$times < 500 | 
                          p$times > 802 & p$times < 864 | p$times > 1166 & p$times < 1228 | 
                          p$times > 1530 & p$times < 1592, 0.2, 0.4)) 

risk_IPD = p$risk



df_outputs = data.frame()
for(intx_ecol in c(1,40,100)){ # ecological interaction
  outputs_IPD = outputs%>%
    mutate(IPD = (Sr + Srr + 0.5*Ssr + Er + Err + 0.5*Esr + intx_ecol*(Ir + Irr + 0.5*Isr) + 
                    Rr + Rrr + 0.5*Rsr + Sra + Srra + Ssra + Era + Erra + Esra + intx_ecol*(Ira + Irra + Isra) + 
                    Rra + Rrra + Rsra + intx_ecol*(Irabx + Irrabx + Israbx) + 
                    Rrabx + Rrrabx + Rsrabx)*rate_IPD*risk_IPD) %>%  #
    mutate(intx_ecol = intx_ecol)%>%
    rowwise()%>%mutate(IPD_pois = rpois(1,IPD))
  
  df_outputs = data.frame(rbind(df_outputs, outputs_IPD))
}
#}



#baseline is based on pre-pandemic scenario only, no sars-cov-2 circulating in the population, what is the annual R incidence
dat_ecol = df_outputs%>%dplyr::filter(intx_ecol == as.numeric(intx_ecol))

### dynamics
ymin_dynamics = min(c(dat_ecol$IPD))
ymax_dynamics = max(c(dat_ecol$IPD))

values_breaks <- as.numeric(c("0", "60", "120", "180", "240", "300", "360", "420", "480"))
cols_lockdowns = c('forestgreen', 'gold', 'darkorange2', 'darkorange3', 'red')


cols_scenario = c('forestgreen', 'gold', 'darkorange2', 'darkorange3', 'red')


###################### Fig 4A
p_IPDs_ecol = ggplot(dat_ecol, aes(x = time, y = IPD, alpha = factor(intx_ecol), group = interaction(intx_ecol, scenario)))+
  geom_line(stat = 'identity', colour = 'darkgreen')+
  facet_wrap(facets = vars(scenario), ncol = 5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size=14, colour = 'black', face = 'bold'))+
  ylim(ymin_dynamics, ymax_dynamics) +
  ylab("Daily incidence per 100,000\n(resistant strain)")+
  xlab('Time (days)')+
  scale_alpha_manual(expression(atop(textstyle("Interaction"), 
                                     atop(textstyle("strength"),
                                          atop(textstyle((psi[c])))))), values = seq(0.2,1, length.out = 6))+
  scale_x_continuous(breaks = values_breaks) + 
  theme(legend.title=element_text(size=16), 
        legend.text=element_text(size=14))




IPD_baseline = 1.167697

### total IPDs
dat_ecol_totals = dat_ecol%>%
  group_by(scenario, intx_ecol)%>%
  summarise(IPD_annual = sum(IPD),
            IPD_pois_annual = sum(IPD_pois))%>%
  mutate(interaction = "ecol")

IPD_nointx = dat_ecol_totals%>%filter(intx_ecol == 1)%>%
  ungroup()%>%
  mutate(IPD_baseline = IPD_baseline)%>%
  dplyr::select(scenario, IPD_baseline)


dat_interactions = rbind(dat_ecol_totals)%>%
  left_join(IPD_nointx)%>%
  mutate(IPD_difference = IPD_annual - IPD_baseline)  


first_row <- plot_grid(p_IPDs_ecol+ggtitle(expression(paste(bold('Ecological interaction: Daily incidence of antibiotic-resistant IPDs')))) +
                         theme(plot.title = element_text(hjust = 0.5, size =18, color = 'darkgreen')),
                       labels = c('A'), label_size = 20)


dat_interactions$note <- round(dat_interactions$IPD_difference, digits = 2)


cols_scenario3 = c('forestgreen', 'gold', 'darkorange2', 'darkorange3', 'red')


p_excess_IPDs_ecol = ggplot(dat_interactions%>%filter(interaction == "ecol"), aes(x = scenario, y = IPD_difference, fill = scenario))+
  geom_bar(stat = 'identity', colour = 'black', alpha=0.3) +  
  theme_classic() +
  facet_wrap(facets = vars(intx_ecol), ncol = 6) +
  ylab("Excess cumulative incidence per 100,000\n(resistant strain)")+
  xlab("% of COVID-19 infected using azithromycin") +
  scale_y_continuous(breaks=seq(-1,6, 0.2), limits = c(-0.02,2.0)) +  
  theme(legend.position = 'none', 
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size=14, colour = 'white', face = 'bold'),
        strip.background = element_rect(colour="black", fill='darkgreen', linetype="solid"))+
  scale_fill_manual(values= cols_scenario3) +
  geom_text(aes(label = note), vjust=-0.5, position = position_dodge(1.2), size=4, fontface = "bold")







row1 <- plot_grid(heatmap + ggtitle(expression(paste(bold('Annual cumulative resistant IPD incidence and levels of antibiotic resistance (%)')))) +
                    theme(plot.title = element_text(hjust = 0.5, size =18, color = 'darkgreen')),
                  labels = c('A'), label_size = 20)


row2 <- plot_grid(p_excess_IPDs_ecol + ggtitle(expression(paste(bold('Ecological interaction: Annual excess incidence of antibiotic-resistant IPDs')))) +
                    theme(plot.title = element_text(hjust = 0.5, size =18, color = 'darkgreen')),
                    labels = c('B'), label_size = 20)

FIG4AB = plot_grid(row1, row2, 
                       nrow = 2, rel_heights = c(0.8,0.8,0.8,0.8), align = "v", axis = "lr")

FIG4AB

ggsave(FIG4AB, filename = paste0(filepath_project_group, "Figure4AB.jpeg"), width = 30, height = 25, unit = 'cm')


