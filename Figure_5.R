filepath_project_group = "~/Desktop/COVID_AMR_GitHub/"
load(paste0(filepath_project_group, "outputSA_2120_2485.Rdata"))


strain_i = "IPD_R"
outputSA_IPDtotals_model0 = filter(outputSA, model == "model_0", strain == strain_i)

df_VOC = data.frame(VOC = c("Wuhan", "Alpha", "Delta", "Omicron"),
           R0_est = c(2.5, 4.55, 5.94, 9.5),  #3.28
           Authors = c("\nLiu et al.\nSharma et al.", "Du et al.", "Du et al.", "Liu & Rocklöv"),
           IPD_total_12mo = c(0,0,0,0))


p_SA_model0 = outputSA_IPDtotals_model0%>%
  ggplot(aes(x = R0, y = R_init, fill = IPD_total_12mo))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = '#5e3c99',  mid = 'white', high = '#e66101',
                       limits = c(min(outputSA_IPDtotals_model0$IPD_total_12mo),
                                  max(outputSA_IPDtotals_model0$IPD_total_12mo)),
                       midpoint = (max(outputSA_IPDtotals_model0$IPD_total_12mo)+min(outputSA_IPDtotals_model0$IPD_total_12mo))/2,
                       name = 'Resistant\nIBD incidence\nper 100,000\n(12 months)')+
  theme_classic()+
  xlab(expression(paste(R[0], " (SARS-CoV-2)")))+ylab("Proportion immunized against SARS-CoV-2")+
  geom_vline(data=df_VOC, mapping = aes(xintercept = R0_est), linetype = 2, colour = "#525252")+
  coord_cartesian(ylim = c(0,1.3))+
  geom_text(data=df_VOC, mapping = aes(x = R0_est-0.22, label = VOC), y = 1.3, angle = 90, hjust = 1, size = 3, colour = "#525252")+
  geom_text(data=df_VOC, mapping = aes(x = R0_est+0.2, label = Authors), y = 1.3, angle = 90, hjust = 1, size = 2, colour = "#525252")+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1))

### SAVE PLOT
ggsave(p_SA_model0, filename = paste0(filepath_project_group, "Figure_5.jpeg"), width = 15, height = 12, unit = 'cm')
