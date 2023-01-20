# ECDC data

rm(list=ls())  # clear space
graphics.off() # clear graphics

library(readxl)
library(scales)
library(tidyverse)
library(ggtext)
library(glue)

data1 <- read_excel("AMR_trends.xlsx")
data1$Year <- as.character(data1$Year)

# subset
penicillin_data <- subset(data1, Antibiotic_class == "penicillin")
macrolides_data <- subset(data1, Antibiotic_class == "macrolides")

penicillin_data$Country <- as.factor(penicillin_data$Country)
macrolides_data$Country <- as.factor(macrolides_data$Country)

# margin of error and 95% CIs - penicillin data
penicillin_data$margin <- qnorm(0.975)*sqrt(penicillin_data$p*(1-penicillin_data$p)/penicillin_data$n)
penicillin_data$lower <- penicillin_data$p - penicillin_data$margin
penicillin_data$upper <- penicillin_data$p + penicillin_data$margin

# margin of error and 95% CIs - macrolide data
macrolides_data$margin <- qnorm(0.975)*sqrt(macrolides_data$p*(1-macrolides_data$p)/macrolides_data$n)
macrolides_data$lower <- macrolides_data$p - macrolides_data$margin
macrolides_data$upper <- macrolides_data$p + macrolides_data$margin





A <- ggplot(data=penicillin_data, aes(x=Country, y=p, ymin = lower, ymax = upper, fill=Year)) +
  geom_bar(stat="identity", position=position_dodge()) + geom_errorbar(position=position_dodge(0.9), size = 0.2, width=0.3) + 
  #scale_fill_brewer(palette="Paired", labels=c("2019", "2020")) + 
  scale_fill_manual(values=c("lightskyblue1", "cornflowerblue"), labels=c("2019", "2020")) +
  theme_classic(base_family='Arial') +
  ggtitle("Penicillin resistance") +
  theme(legend.position = c(0.5,0.9), 
        legend.direction = "horizontal", 
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + #face = bold.labels, 
  labs(x="", y="Resistant isolates\n (proportion)") +  
  #scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 5)) + 
  scale_x_discrete(limits = c("Malta","Romania","France","Iceland","Croatia","Spain","Ireland",
                              "Latvia", "Belgium", "Lithuania",
                              "Italy", "Slovenia", "Hungary", "Sweden",
                              "Norway", "United Kingdom", "Denmark", "Germany", "Netherlands",
                              "Poland", "Portugal", "Finland", "Austria", "Czechia")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), 
        legend.title = element_blank(),
        legend.text = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        axis.title.x = element_text(size=10),
        axis.text.x = element_text())



B <- ggplot(data=macrolides_data, aes(x=Country, y=p, ymin = lower, ymax = upper, fill=Year)) +
  geom_bar(stat="identity", position=position_dodge()) + geom_errorbar(position=position_dodge(0.9), size = 0.2, width=0.3) +
  #scale_fill_brewer(palette="Paired", labels=c("2019", "2020")) + 
  scale_fill_manual(values=c("aquamarine3", "forestgreen"), labels=c("2019", "2020")) +
  theme_classic(base_family='Arial') +
  ggtitle("Macrolide resistance") +
  theme(legend.position = c(0.5,0.9), 
        legend.direction = "horizontal", 
        legend.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) + #face = bold.labels, 
  labs(x="", y="Resistant isolates\n (proportion)") +  
  #scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 5)) + 
  scale_x_discrete(limits = c("Malta","Croatia","Iceland","Romania","Italy","Spain","France",
                              "Belgium", "Hungary","Portugal","Lithuania","Slovenia",
                              "Ireland","Finland","Estonia",
                              "Sweden","United Kingdom","Denmark", 
                              "Poland","Austria", "Czechia","Germany", "Norway","Netherlands")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), 
        legend.title = element_blank(),
        legend.text = element_text(size=10), 
        axis.title.y = element_text(size=10), 
        axis.title.x = element_text(size=10),
        axis.text.x = element_text())


library("cowplot")
fig <- plot_grid(A, B,labels = c("A", "B"), label_size = 16,
                           ncol = 1, nrow = 2)

ggsave("Figure_1.jpeg", width = 20, height = 15, units = "cm")

