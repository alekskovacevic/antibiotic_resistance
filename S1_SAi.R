# Community model:
################################## S1 - increase in atb use, no lockdown, S1_SAi a = 2 ###############

# Total population size
N <- 100000

initial_state_values <- c(S = N - 3*0.1*N, Ss = 0.1*N, Sr = 0.1*N, Ssr = 0.1*N, 
                          Sa = 0, Sas = 0, Sar = 0, Sasr = 0,
                          E = 0, Es = 0, Er = 0, Esr = 0, 
                          Ea = 0, Eas = 0, Ear = 0, Easr = 0,
                          I = 0, Is = 0, Ir = 0, Isr = 0, 
                          Ia = 0, Ias = 0, Iar = 0, Iasr = 0,
                          R = 0, Rs = 0, Rr = 0, Rsr = 0, 
                          Ra = 0, Ras = 0, Rar = 0, Rasr = 0,
                          Total = 0, TotalS = 0, TotalR = 0)


parameters <- c(betaS = 0.043,        # transmission rate (bacteria)
                c = 0.05,             # fitness cost of resistance 
                gamma = 1/30,         # natural clearance rate of bacterial carriage
                k = 0.5,              # probability of acquiring secondary bacterial carriage  
                q = 0.90,             # relative bacterial carriage transmission risk under quarantine (in addition to quarantine)
                w = 1/4,              # rate of antibiotic action
                r = 1/8,              # rate of return to antibiotic unexposed compartment
                tau = 0.002,          # average daily ppc (prescriptions per capita) in the community
                gammaC = 1/10,        # rate of recovery (COVID-19)
                A = 20,               # an increase in antibiotic use in covid-19 infected individuals
                betaC = 0.25,         # transmission rate of sars-cov-2, R0 = 2.5                                                                                                                                                                                                                         ,        # covid transmission rate
                alpha = 1/6,          # 1/incubation rate for sars-cov-2
                p = 1.4e-06)          # bacterial invasion rate (pathogenicity)  

# VARY ANTIBIOTIC USE IN THE COMMUNITY ACCORDING TO LOCKDOWN: a
times <- seq(0, 8000, by = 1)

# assume lockdown lasts from day 4120 to 4210
atb_1 <- as.data.frame(list(times = times, a_1 = rep(0, length(times))))
atb_1$a_1 <- ifelse((atb_1$times > 4119 & atb_1$times < 4211), 2, 1) 
atb_1
#Create the interpolating function
input_1 <- approxfun(atb_1)
input_1(seq(from = 0, to = 8000, by = 1))


# VARY RELATIVE BACTERIAL TRANSMISSION RISK IN THE COMMUNITY ACCORDING TO LOCKDOWN: thetabeta
# assume lockdown lasts from day 4120 to 4210
beta_1 <- as.data.frame(list(times = times, thetabeta_1 = rep(0, length(times))))
beta_1$thetabeta_1 <- ifelse((beta_1$times > 4119 & beta_1$times < 4211), 1, 1) 
beta_1
#Create the interpolating function
inputbeta_1 <- approxfun(beta_1)
inputbeta_1(seq(from = 0, to = 8000, by = 1))


################################### MODEL FUNCTION: ######################################
times_model <- seq(0, 8000, by = 1) # model time duration

# Model 1:
model_1 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {  
    covid = if(floor(time)==4000) { covid = 2e-5 } else {covid = 0}  # introduce covid on day 4000 
    
    a_1 <- input_1(time)
    thetabeta_1 <- inputbeta_1(time)

    
    N = S + Ss + Sr + Ssr + Sa + Sas + Sar + Sasr + 
      E + Es + Er + Esr + Ea + Eas + Ear + Easr + 
      I + Is + Ir + Isr + Ia + Ias + Iar + Iasr + 
      R + Rs + Rr + Rsr + Ra + Ras + Rar + Rasr  # sum all compartments
    
    #print(N)
    
    lambdaS <- betaS*thetabeta_1*(Ss/N + 1/2*Ssr/N + Sas/N + 1/2*Sasr/N + 
                                  Es/N + 1/2*Esr/N + Eas/N + 1/2*Easr/N + 
                                  q*(Is/N + 1/2*Isr/N + Ias/N + 1/2*Iasr/N) + 
                                  Rs/N + 1/2*Rsr/N + Ras/N + 1/2*Rasr/N)
    
    # betaR = betaS*(1-c)
    lambdaR <- betaS*(1-c)*thetabeta_1*(Sr/N + 1/2*Ssr/N + Sar/N + 1/2*Sasr/N + 
                                        Er/N + 1/2*Esr/N + Ear/N + 1/2*Easr/N + 
                                        q*(Ir/N + 1/2*Isr/N + Iar/N + 1/2*Iasr/N) + 
                                        Rr/N + 1/2*Rsr/N + Rar/N + 1/2*Rasr/N)
    
    lambdaC <- betaC*thetabeta_1*q*(I/N + Is/N + Ir/N + Isr/N + Ia/N + Ias/N + Iar/N + Iasr/N)   # covid force of infection
    
    
    
    # Equations:    
    
    # Susceptible to SARS-CoV-2:
    
    dS <- -(lambdaS + lambdaR)*S + gamma*(Ss+Sr) - tau*a_1*S + r*Sa - lambdaC*S - covid*N
    dSs <- lambdaS*S - gamma*Ss - k*lambdaR*Ss + gamma*1/2*Ssr - tau*a_1*Ss + r*Sas - lambdaC*Ss
    dSr <- lambdaR*S - gamma*Sr - k*lambdaS*Sr + gamma*1/2*Ssr - tau*a_1*Sr + r*Sar - lambdaC*Sr
    dSsr <- k*lambdaS*Sr + k*lambdaR*Ss - gamma*Ssr - tau*a_1*Ssr + r*Sasr - lambdaC*Ssr
    
    
    dSa <- -lambdaR*Sa + (gamma+w)*Sas + gamma*Sar + tau*a_1*S - r*Sa - lambdaC*Sa
    dSas <- -k*lambdaR*Sas - (gamma+w)*Sas + (gamma*1/2)*Sasr + tau*a_1*Ss - r*Sas - lambdaC*Sas 
    dSar <- lambdaR*Sa - gamma*Sar + (gamma*1/2+w)*Sasr + tau*a_1*Sr - r*Sar - lambdaC*Sar
    dSasr <- k*lambdaR*Sas - (gamma*1/2)*Sasr - (gamma*1/2+w)*Sasr + tau*a_1*Ssr - r*Sasr - lambdaC*Sasr
    
    
    # Exposed to SARS-CoV-2:
    
    dE <- -(lambdaS + lambdaR)*E + gamma*(Es+Er) - tau*a_1*E + r*Ea - alpha*E + lambdaC*S
    dEs <- lambdaS*E - gamma*Es - k*lambdaR*Es + gamma*1/2*Esr - tau*a_1*Es + r*Eas - alpha*Es + lambdaC*Ss
    dEr <- lambdaR*E - gamma*Er - k*lambdaS*Er + gamma*1/2*Esr - tau*a_1*Er + r*Ear - alpha*Er + lambdaC*Sr
    dEsr <- k*lambdaS*Er + k*lambdaR*Es - gamma*Esr - tau*a_1*Esr + r*Easr - alpha*Esr + lambdaC*Ssr
    
    
    dEa <- -lambdaR*Ea + (gamma+w)*Eas + gamma*Ear + tau*a_1*E - r*Ea + lambdaC*Sa - alpha*Ea
    dEas <- -k*lambdaR*Eas - (gamma+w)*Eas  + (gamma*1/2)*Easr + tau*a_1*Es - r*Eas + lambdaC*Sas - alpha*Eas
    dEar <- lambdaR*Ea - gamma*Ear + (gamma*1/2+w)*Easr + tau*a_1*Er - r*Ear + lambdaC*Sar - alpha*Ear
    dEasr <- k*lambdaR*Eas - (gamma*1/2)*Easr - (gamma*1/2+w)*Easr + tau*a_1*Esr - r*Easr + lambdaC*Sasr - alpha*Easr
    
    
    # Infected with SARS-CoV-2:
    
    dI <- -(lambdaS + lambdaR)*I + gamma*(Is+Ir) - tau*a_1*A*I + r*Ia + alpha*E - gammaC*I + covid*N
    dIs <- lambdaS*I - gamma*Is - k*lambdaR*Is + gamma*1/2*Isr - tau*a_1*A*Is + r*Ias + alpha*Es - gammaC*Is
    dIr <- lambdaR*I - gamma*Ir - k*lambdaS*Ir + gamma*1/2*Isr - tau*a_1*A*Ir + r*Iar + alpha*Er - gammaC*Ir
    dIsr <- k*lambdaS*Ir + k*lambdaR*Is - gamma*Isr - tau*a_1*A*Isr + r*Iasr + alpha*Esr - gammaC*Isr
    
    
    dIa <- -lambdaR*Ia + (gamma+w)*Ias + gamma*Iar + tau*a_1*A*I - r*Ia + alpha*Ea - gammaC*Ia
    dIas <-  -k*lambdaR*Ias - (gamma+w)*Ias + (gamma*1/2)*Iasr + tau*a_1*A*Is - r*Ias + alpha*Eas - gammaC*Ias
    dIar <- lambdaR*Ia - gamma*Iar + (gamma*1/2+w)*Iasr + tau*a_1*A*Ir - r*Iar + alpha*Ear - gammaC*Iar
    dIasr <- k*lambdaR*Ias - (gamma*1/2)*Iasr - (gamma*1/2+w)*Iasr + tau*a_1*A*Isr - r*Iasr + alpha*Easr - gammaC*Iasr
    
    # Recovered from SARS-CoV-2:
    
    dR <- -(lambdaS + lambdaR)*R + gamma*(Rs+Rr) - tau*a_1*R + r*Ra + gammaC*I
    dRs <- lambdaS*R - gamma*Rs - k*lambdaR*Rs + gamma*1/2*Rsr - tau*a_1*Rs + r*Ras + gammaC*Is
    dRr <- lambdaR*R - gamma*Rr - k*lambdaS*Rr + gamma*1/2*Rsr - tau*a_1*Rr + r*Rar + gammaC*Ir
    dRsr <- k*lambdaS*Rr + k*lambdaR*Rs - gamma*Rsr - tau*a_1*Rsr + r*Rasr + gammaC*Isr
    
    
    dRa <- -lambdaR*Ra + (gamma+w)*Ras + gamma*Rar + tau*a_1*R - r*Ra + gammaC*Ia
    dRas <-  - k*lambdaR*Ras - (gamma+w)*Ras + (gamma*1/2)*Rasr + tau*a_1*Rs - r*Ras + gammaC*Ias
    dRar <- lambdaR*Ra - gamma*Rar + (gamma*1/2+w)*Rasr + tau*a_1*Rr - r*Rar + gammaC*Iar
    dRasr <- k*lambdaR*Ras - (gamma*1/2)*Rasr - (gamma*1/2+w)*Rasr + tau*a_1*Rsr - r*Rasr + gammaC*Iasr
    
    # IPDs
    
    dTotal <- p*(Ss+Sr+Sas+Sar+Es+Er+Eas+Ear+Is+Ir+Ias+Iar+Rs+Rr+Ras+Rar) + p*(Ssr+Sasr+Esr+Easr+Isr+Iasr+Rsr+Rasr)
    
    dTotalS <- p*(Ss+Sas+Es+Eas+Is+Ias+Rs+Ras) + (p/2)*(Ssr+Sasr+Esr+Easr+Isr+Iasr+Rsr+Rasr)
    
    dTotalR <- p*(Sr+Sar+Er+Ear+Ir+Iar+Rr+Rar) + (p/2)*(Ssr+Sasr+Esr+Easr+Isr+Iasr+Rsr+Rasr)
    
    
    
    result <- c(dS, dSs, dSr, dSsr, dSa, dSas, dSar, dSasr, 
                dE, dEs, dEr, dEsr, dEa, dEas, dEar, dEasr, 
                dI, dIs, dIr, dIsr, dIa, dIas, dIar, dIasr,
                dR, dRs, dRr, dRsr, dRa, dRas, dRar, dRasr,
                dTotal, dTotalS, dTotalR)
    
    list(result, atb_1 = a_1, beta_1 = thetabeta_1)
    
    
  })
}



output_1 <- as.data.frame(ode(y = initial_state_values, 
                            times = times_model, 
                            func = model_1,
                            parms = parameters,
                            method = "euler"))




############################################## PLOTS #################################
values1 <- c ("0", "60", "120", "180", "240", "300", "360")
lock <- data.frame(xmin=4120, xmax=4210, ymin=-Inf, ymax=Inf)

colors <- c("S" = "cornflowerblue", "R" = "red", "none" = "grey", "dual" = "purple", "covid" = "orange")

# Bacterial carriage dynamics:
plotS1_a <- ggplot(output_1, aes(x=time)) + 
  geom_line(aes(y=(Ss + Sas + Es + Eas + Is + Ias + Rs + Ras +
                   Sr + Sar + Er + Ear + Ir + Iar + Rr + Rar +
                   Ssr + Sasr + Esr + Easr + Isr + Iasr + Rsr + Rasr)/N, color="none"), size=2) +  # total carriage, black
  geom_line(aes(y=(Ss + Sas + Es + Eas + Is + Ias + Rs + Ras)/N, color="S"), size=2) +  # S, blue
  geom_line(aes(y=(Sr + Sar + Er + Ear + Ir + Iar + Rr + Rar)/N, color="R"), size=2) +  # R, red
  geom_line(aes(y=(Ssr + Sasr + Esr + Easr + Isr + Iasr + Rsr + Rasr)/N, color="dual"), size=2) +  # SR, purple
  geom_line(aes(y=(I + Is + Ir + Isr + Ia + Ias + Iar + Iasr)/N, color="covid"), size=2) +  # SR, purple
  scale_x_continuous(limits = c(3995, 4365), breaks=seq(4000, 4365, 60), labels = values1) +
  scale_y_continuous(limits = c(0, 0.27), breaks=seq(0, 1, 0.05), labels = number) +
  theme_classic() + 
  labs(x="", y="") +
  #ggtitle("Bacterial carriage and SARS-CoV-2 infection") + 
  theme(legend.position="none", #legend.position="right",#legend.position = c(.8,.82) 
        legend.background=element_blank(),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + 
  scale_color_manual(name = "", values = colors, 
                     labels = c("S", 
                                "R", 
                                "Total carriage",
                                "SR",
                                "SARS-CoV-2")) +
  geom_rect(data=lock, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            color="transparent", fill = "red",
            alpha=0.2,
            inherit.aes = FALSE)



colors2 <- c("S" = "cornflowerblue", "R" = "red")

Carriage_S1 <- output_1$Ss + output_1$Sas + output_1$Es + output_1$Eas + output_1$Is + output_1$Ias + output_1$Rs + output_1$Ras
Carriage_R1 <- output_1$Sr + output_1$Sar + output_1$Er + output_1$Ear + output_1$Ir + output_1$Iar + output_1$Rr + output_1$Rar
Carriage_SR1 <- output_1$Ssr + output_1$Sasr + output_1$Esr + output_1$Easr + output_1$Isr + output_1$Iasr + output_1$Rsr + output_1$Rasr


output_1$R_carriage <- (output_1$Sr + output_1$Sar + output_1$Er + output_1$Ear + output_1$Ir + output_1$Iar + output_1$Rr + output_1$Rar + 
                          0.5*(output_1$Ssr + output_1$Sasr + output_1$Esr + output_1$Easr + output_1$Isr + output_1$Iasr + output_1$Rsr + output_1$Rasr))/
  (output_1$Ss + output_1$Sas + output_1$Es + output_1$Eas + output_1$Is + output_1$Ias + output_1$Rs + output_1$Ras + 
     output_1$Ssr + output_1$Sasr + output_1$Esr + output_1$Easr + output_1$Isr + output_1$Iasr + output_1$Rsr + output_1$Rasr + 
     output_1$Sr + output_1$Sar + output_1$Er + output_1$Ear + output_1$Ir + output_1$Iar + output_1$Rr + output_1$Rar)

R_rateS1 <- (output_1$R_carriage[4210] - output_1$R_carriage[4000])/output_1$R_carriage[4000] 
R_rateS1*100 # % change in antibiotic resistance rate

rs1 <- round(R_rateS1*100, digits = 2)
rs1



# Antibiotic resistance rate:
plotS1_b <- ggplot(output_1, aes(x=time)) + 
  geom_line(aes(y=(Carriage_S1 + 0.5*(Carriage_SR1))/
                  (Carriage_S1 + Carriage_SR1 + Carriage_R1), color="S"), size=2) +  # S, blue
  geom_line(aes(y=(Carriage_R1 + 0.5*(Carriage_SR1))/
                  (Carriage_S1 + Carriage_SR1 + Carriage_R1), color="R"), size=2) +
  scale_x_continuous(limits = c(3995, 4365), breaks=seq(4000, 4365, 60), labels = values1) +
  scale_y_continuous(limits = c(0.3, 0.7), breaks=seq(0, 1, 0.05), labels = number) +
  theme_classic() + 
  labs(x="", y="") +
  #ggtitle("Change in antibiotic resistance rate") + 
  theme(legend.position="none", #, #legend.position = c(.75,.65) 
        legend.background=element_blank(),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + 
  scale_color_manual(name = "Bacterial carriage:", values = colors2, 
                     labels = c("S and 1/2(SR)", "R and 1/2(SR)")) +
  geom_rect(data=lock, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            color="transparent", fill = "red",
            alpha=0.2,
            inherit.aes = FALSE) +
  annotate(geom="text", x=4170, y=0.67, label=rs1, color="black", size = 8, fontface = 'bold') +
  annotate(geom="text", x=4230, y=0.67, label="%", color="black", size = 8, fontface = 'bold')
#########################################################################################################

#---------- Cumulative IPD incidence short and long term:------------------------------------

# Cumulative ipd incidence (S, R, SR):
total_short_S1 <- output_1$Total[4210] - output_1$Total[4120]
total_long_S1 <- output_1$Total[4365] - output_1$Total[4000]

# IPD incidence caused by susceptible strain only (S + 0.5*SR):
S_short_S1 <- output_1$TotalS[4210] - output_1$TotalS[4120]
S_long_S1 <- output_1$TotalS[4365] - output_1$TotalS[4000]

# IPD incidence caused by resistant strain only (R + 0.5*SR):
R_short_S1 <- output_1$TotalR[4210] - output_1$TotalR[4120]
R_long_S1 <- output_1$TotalR[4365] - output_1$TotalR[4000]
