# clear space
rm(list=ls()) 

# clear graphics
graphics.off() 

library(tidyverse)

### Run all four scenarios

filepath_project_group = "~/Desktop/COVID_AMR_GitHub/"

source(paste0(filepath_project_group, "S0_no_cov.R"))
source(paste0(filepath_project_group, "S0_Snc.R"))
source(paste0(filepath_project_group, "S1_SAi.R"))
source(paste0(filepath_project_group, "S2_Sad.R"))
source(paste0(filepath_project_group, "S3_SL.R"))
source(paste0(filepath_project_group, "S4_SAiL.R"))
source(paste0(filepath_project_group, "S5_SadL.R"))


outputs = rbind(output_0%>%mutate(scenario = "S0")%>%dplyr::select(-c("atb_0", "beta_0")),
                output_1%>%mutate(scenario = "S1")%>%dplyr::select(-c("atb_1", "beta_1")),
                output_2%>%mutate(scenario = "S2")%>%dplyr::select(-c("atb_2", "beta_2")),
                output_3%>%mutate(scenario = "S3")%>%dplyr::select(-c("atb_3", "beta_3")),
                output_4%>%mutate(scenario = "S4")%>%dplyr::select(-c("atb_4", "beta_4")),
                output_5%>%mutate(scenario = "S5")%>%dplyr::select(-c("atb_5", "beta_5")))%>%
  
  filter(time > 4000 & time < 4365)%>%    # 1 year 
  mutate(time = time - 4000)


rate_IPD = parameters["p"]


df_outputs = data.frame()
for(intx_ecol in c(1,10,25,85,90,100)){ # ecological interaction
  for(intx_abx in c(1,3,4,6,10,12)){ # antibiotic interaction
    outputs_IPD = outputs%>%
      mutate(IPD = (Sr + Ssr*0.5 + (Sar + Sasr*0.5)*intx_abx + 
                    Er + Esr*0.5 + (Ear + Easr*0.5)*intx_abx + 
                   (Ir + Isr*0.5 + (Iar + Iasr*0.5)*intx_abx)*intx_ecol + 
                    Rr + Rsr*0.5 + (Rar + Rasr*0.5)*intx_abx)*rate_IPD)%>%
      mutate(intx_ecol = intx_ecol, intx_abx = intx_abx)%>%
      rowwise()%>%mutate(IPD_pois = rpois(1,IPD))
    
    df_outputs = data.frame(rbind(df_outputs, outputs_IPD))
  }
}



#baseline is based on pre-pandemic scenario only, no sars-cov-2 circulating in the population, what is the annual R incidence
dat_ecol = df_outputs%>%dplyr::filter(intx_abx == 1)%>%mutate(intx_ecol = as.numeric(intx_ecol))
dat_abx = df_outputs%>%dplyr::filter(intx_ecol == 1)%>%mutate(intx_abx = as.numeric(intx_abx))

### dynamics
ymin_dynamics = min(c(dat_ecol$IPD, dat_abx$IPD))
ymax_dynamics = max(c(dat_ecol$IPD, dat_abx$IPD))

values_breaks <- as.numeric(c("0", "60", "120", "180", "240", "300", "360", "420", "480"))
cols_lockdowns = c('white', '#FED6D7', '#D8E9DF', 'grey', '#FFACAA',  '#AFD1BC')



p_IPDs_ecol = ggplot(dat_ecol, aes(x = time, y = IPD, alpha = factor(intx_ecol), group = interaction(intx_ecol, scenario)))+
  geom_rect(data=lock, aes(xmin=xmin-4000, xmax=xmax-4000, ymin=ymin, ymax=ymax),
            color="transparent", fill = cols_lockdowns,
            inherit.aes = FALSE)+
  geom_line(stat = 'identity', colour = 'darkgreen')+
  facet_wrap(facets = vars(scenario), ncol = 6)+
  theme_classic()+
  theme(strip.text = element_text(size=14, colour = 'black', face = 'bold'))+
  ylim(ymin_dynamics, ymax_dynamics) +
  ylab("Daily incidence per 100,000\n(resistant strain)")+
  xlab('Time (days)')+
  scale_alpha_manual(expression(atop(textstyle("Interaction"), 
                                     atop(textstyle("strength"),
                                          atop(textstyle((psi[c])))))), values = seq(0.2,1, length.out = 6))+
  scale_x_continuous(breaks = values_breaks) + 
  theme(legend.title=element_text(size=16), 
        legend.text=element_text(size=14))

p_IPDs_abx = ggplot(dat_abx, aes(x = time, y = IPD, alpha = factor(intx_abx), group = interaction(intx_abx, scenario)))+
  geom_rect(data=lock, aes(xmin=xmin-4000, xmax=xmax-4000, ymin=ymin, ymax=ymax),
            color="transparent", fill = cols_lockdowns,
            inherit.aes = FALSE)+
  geom_line(stat = 'identity', colour = 'navyblue')+
  facet_wrap(facets = vars(scenario), ncol = 6)+
  theme_classic()+
  theme(strip.text = element_text(size=14, colour = 'black', face = 'bold'))+
  ylim(ymin_dynamics, 0.04)+
  ylab("Daily incidence per 100,000\n(resistant strain)")+
  xlab('Time (days)')+
  scale_alpha_manual(expression(atop(textstyle("Interaction"), 
                                     atop(textstyle("strength"),
                                          atop(textstyle((psi[a])))))), values = seq(0.2,1, length.out = 6))+
  scale_x_continuous(breaks = values_breaks) +
  theme(legend.title = element_text(hjust = 0)) + 
  theme(legend.title=element_text(size=16), 
        legend.text=element_text(size=14))


IPD_baseline = 6.263677   # incidence of resistant disease per 100,000 (annual)

### total IPDs
dat_ecol_totals = dat_ecol%>%
  group_by(scenario, intx_ecol, intx_abx)%>%
  summarise(IPD_annual = sum(IPD),
            IPD_pois_annual = sum(IPD_pois))%>%
  mutate(interaction = "ecol")

IPD_nointx = dat_ecol_totals%>%filter(intx_ecol == 1, intx_abx == 1)%>%
  ungroup()%>%
  mutate(IPD_baseline = IPD_baseline)%>%
  dplyr::select(scenario, IPD_baseline)

dat_abx_totals = dat_abx%>%
  group_by(scenario, intx_ecol, intx_abx)%>%
  summarise(IPD_annual = sum(IPD),
            IPD_pois_annual = sum(IPD_pois))%>%
  mutate(interaction = "abx")

dat_interactions = rbind(dat_ecol_totals, dat_abx_totals)%>%
  left_join(IPD_nointx)%>%
  mutate(IPD_difference = IPD_annual - IPD_baseline)  


cols_scenario = c('orange','white', '#FED6D7', '#D8E9DF', 'grey', '#FFACAA',  '#AFD1BC')
IPD_baseline = 6.263677


first_row <- plot_grid(p_IPDs_ecol+ggtitle(expression(paste(bold('Ecological interaction: Daily incidence of antibiotic-resistant IBDs')))) +
                       theme(plot.title = element_text(hjust = 0.5, size =18, color = 'darkgreen')),
                       labels = c('A'), label_size = 20)

second_row <- plot_grid(p_IPDs_abx+ggtitle(expression(paste(bold('Antibiotic interaction: Daily incidence of antibiotic-resistant IBDs'))))+
                        theme(plot.title = element_text(hjust = 0.5, size =18, color = 'navyblue')),
                        labels = c('B'), label_size = 20)



dat_interactions$note <- round(dat_interactions$IPD_difference, digits = 2)

#write.csv(dat_interactions, file = "dat_interactions.csv")

######## difference only
cols_scenario3 = c('white', '#FED6D7', '#D8E9DF', 'grey', '#FFACAA',  '#AFD1BC')


p_excess_IPDs_ecol = ggplot(dat_interactions%>%filter(interaction == "ecol"), aes(x = scenario, y = IPD_difference, fill = scenario))+
  geom_bar(stat = 'identity', colour = 'black') +  
  theme_bw() +
  facet_wrap(facets = vars(intx_ecol), ncol = 6) +
  ylab("Excess cumulative incidence per 100,000\n(resistant strain)")+
  xlab("") +
  scale_y_continuous(breaks=seq(-4,16), limits = c(-4,16,1)) +  
  theme(legend.position = 'none', 
        strip.text = element_text(size=14, colour = 'white', face = 'bold'),
        strip.background = element_rect(colour="black", fill='darkgreen', linetype="solid"))+
  scale_fill_manual(values= cols_scenario3) +
  geom_text(aes(label = note), vjust=-0.5, position = position_dodge(1.2), size=4, fontface = "bold")

  


p_excess_IPDs_abx = ggplot(dat_interactions%>%filter(interaction == "abx"), aes(x = scenario, y = IPD_difference, fill = scenario))+
  geom_bar(stat = 'identity', colour = 'black') +  
  theme_bw()+
  facet_wrap(facets = vars(intx_abx), ncol = 6)+
  ylab("Excess cumulative incidence per 100,000\n(resistant strain)")+
  xlab("")+
  scale_y_continuous(breaks=seq(-4,16), limits = c(-4,16,1)) +  
  theme(legend.position = 'none', 
        strip.text = element_text(size=14, colour = 'white', face = 'bold'),
        strip.background = element_rect(colour="black", fill="navyblue", linetype="solid"))+
  scale_fill_manual(values= cols_scenario3) +
  geom_text(aes(label = note), vjust=-0.5, position = position_dodge(1.2), size=4, fontface = "bold")


row3 <- plot_grid(p_excess_IPDs_ecol + ggtitle(expression(paste(bold('Ecological interaction: Annual changes in the incidence of antibiotic-resistant IBDs')))) +
                         theme(plot.title = element_text(hjust = 0.5, size =18, color = 'darkgreen')),
                         labels = c('C'), label_size = 20)

row4 <- plot_grid(p_excess_IPDs_abx + ggtitle(expression(paste(bold('Antibiotic interaction: Annual changes in the incidence of antibiotic-resistant IBDs'))))+
                          theme(plot.title = element_text(hjust = 0.5, size=18, color = 'navyblue')),
                          labels = c('D'), label_size = 20)

final_plot = plot_grid(first_row, second_row, row3, row4,
                       nrow = 4, rel_heights = c(0.8,0.8,0.8,0.8), align = "v", axis = "lr")

ggsave(final_plot, filename = paste0(filepath_project_group, "Figure_4.jpeg"), width = 45, height = 40, unit = 'cm')
