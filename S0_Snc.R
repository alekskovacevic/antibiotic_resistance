# Community model:
################################## SO - no change in atb use, no lockdown, S0_Snc

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
atb_0 <- as.data.frame(list(times = times, a_0 = rep(0, length(times))))
atb_0$a_0 <- ifelse((atb_0$times > 4119 & atb_0$times < 4211), 1, 1) 
atb_0
#Create the interpolating function
input_0 <- approxfun(atb_0)
input_0(seq(from = 0, to = 8000, by = 1))


# VARY RELATIVE BACTERIAL TRANSMISSION RISK IN THE COMMUNITY ACCORDING TO LOCKDOWN: thetabeta
# assume lockdown lasts from day 4120 to 4210
beta_0 <- as.data.frame(list(times = times, thetabeta_0 = rep(0, length(times))))
beta_0$thetabeta_0 <- ifelse((beta_0$times > 4119 & beta_0$times < 4211), 1, 1) 
beta_0
#Create the interpolating function
inputbeta_0 <- approxfun(beta_0)
inputbeta_0(seq(from = 0, to = 8000, by = 1))


################################### MODEL FUNCTION: ######################################
times_model <- seq(0, 8000, by = 1) # model time duration

# Model 1:
model_0 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {  
    covid = if(floor(time)==4000) { covid = 2e-5 } else {covid = 0}  # introduce covid on day 4000 
    
    a_0 <- input_0(time)
    thetabeta_0 <- inputbeta_0(time)

    
    N = S + Ss + Sr + Ssr + Sa + Sas + Sar + Sasr + 
      E + Es + Er + Esr + Ea + Eas + Ear + Easr + 
      I + Is + Ir + Isr + Ia + Ias + Iar + Iasr + 
      R + Rs + Rr + Rsr + Ra + Ras + Rar + Rasr  # sum all compartments
    
    #print(N)
    
    lambdaS <- betaS*thetabeta_0*(Ss/N + 1/2*Ssr/N + Sas/N + 1/2*Sasr/N + 
                                  Es/N + 1/2*Esr/N + Eas/N + 1/2*Easr/N + 
                               q*(Is/N + 1/2*Isr/N + Ias/N + 1/2*Iasr/N) + 
                                  Rs/N + 1/2*Rsr/N + Ras/N + 1/2*Rasr/N)
    
    # betaR = betaS*(1-c)
    lambdaR <- betaS*(1-c)*thetabeta_0*(Sr/N + 1/2*Ssr/N + Sar/N + 1/2*Sasr/N + 
                                        Er/N + 1/2*Esr/N + Ear/N + 1/2*Easr/N + 
                                     q*(Ir/N + 1/2*Isr/N + Iar/N + 1/2*Iasr/N) + 
                                        Rr/N + 1/2*Rsr/N + Rar/N + 1/2*Rasr/N)
    
    lambdaC <- betaC*thetabeta_0*q*(I/N + Is/N + Ir/N + Isr/N + Ia/N + Ias/N + Iar/N + Iasr/N)   # covid force of infection
    
    
    
    # Equations:    
    
    # Susceptible to SARS-CoV-2:
    
    dS <- -(lambdaS + lambdaR)*S + gamma*(Ss+Sr) - tau*a_0*S + r*Sa - lambdaC*S - covid*N
    dSs <- lambdaS*S - gamma*Ss - k*lambdaR*Ss + gamma*1/2*Ssr - tau*a_0*Ss + r*Sas - lambdaC*Ss
    dSr <- lambdaR*S - gamma*Sr - k*lambdaS*Sr + gamma*1/2*Ssr - tau*a_0*Sr + r*Sar - lambdaC*Sr
    dSsr <- k*lambdaS*Sr + k*lambdaR*Ss - gamma*Ssr - tau*a_0*Ssr + r*Sasr - lambdaC*Ssr
    
    
    dSa <- -lambdaR*Sa + (gamma+w)*Sas + gamma*Sar + tau*a_0*S - r*Sa - lambdaC*Sa
    dSas <- -k*lambdaR*Sas - (gamma+w)*Sas + (gamma*1/2)*Sasr + tau*a_0*Ss - r*Sas - lambdaC*Sas 
    dSar <- lambdaR*Sa - gamma*Sar + (gamma*1/2+w)*Sasr + tau*a_0*Sr - r*Sar - lambdaC*Sar
    dSasr <- k*lambdaR*Sas - (gamma*1/2)*Sasr - (gamma*1/2+w)*Sasr + tau*a_0*Ssr - r*Sasr - lambdaC*Sasr
    
    
    # Exposed to SARS-CoV-2:
    
    dE <- -(lambdaS + lambdaR)*E + gamma*(Es+Er) - tau*a_0*E + r*Ea - alpha*E + lambdaC*S
    dEs <- lambdaS*E - gamma*Es - k*lambdaR*Es + gamma*1/2*Esr - tau*a_0*Es + r*Eas - alpha*Es + lambdaC*Ss
    dEr <- lambdaR*E - gamma*Er - k*lambdaS*Er + gamma*1/2*Esr - tau*a_0*Er + r*Ear - alpha*Er + lambdaC*Sr
    dEsr <- k*lambdaS*Er + k*lambdaR*Es - gamma*Esr - tau*a_0*Esr + r*Easr - alpha*Esr + lambdaC*Ssr
    
    
    dEa <- -lambdaR*Ea + (gamma+w)*Eas + gamma*Ear + tau*a_0*E - r*Ea + lambdaC*Sa - alpha*Ea
    dEas <- -k*lambdaR*Eas - (gamma+w)*Eas  + (gamma*1/2)*Easr + tau*a_0*Es - r*Eas + lambdaC*Sas - alpha*Eas
    dEar <- lambdaR*Ea - gamma*Ear + (gamma*1/2+w)*Easr + tau*a_0*Er - r*Ear + lambdaC*Sar - alpha*Ear
    dEasr <- k*lambdaR*Eas - (gamma*1/2)*Easr - (gamma*1/2+w)*Easr + tau*a_0*Esr - r*Easr + lambdaC*Sasr - alpha*Easr
    
    
    # Infected with SARS-CoV-2:
    
    dI <- -(lambdaS + lambdaR)*I + gamma*(Is+Ir) - tau*a_0*A*I + r*Ia + alpha*E - gammaC*I + covid*N
    dIs <- lambdaS*I - gamma*Is - k*lambdaR*Is + gamma*1/2*Isr - tau*a_0*A*Is + r*Ias + alpha*Es - gammaC*Is
    dIr <- lambdaR*I - gamma*Ir - k*lambdaS*Ir + gamma*1/2*Isr - tau*a_0*A*Ir + r*Iar + alpha*Er - gammaC*Ir
    dIsr <- k*lambdaS*Ir + k*lambdaR*Is - gamma*Isr - tau*a_0*A*Isr + r*Iasr + alpha*Esr - gammaC*Isr
    
    
    dIa <- -lambdaR*Ia + (gamma+w)*Ias + gamma*Iar + tau*a_0*A*I - r*Ia + alpha*Ea - gammaC*Ia
    dIas <-  -k*lambdaR*Ias - (gamma+w)*Ias + (gamma*1/2)*Iasr + tau*a_0*A*Is - r*Ias + alpha*Eas - gammaC*Ias
    dIar <- lambdaR*Ia - gamma*Iar + (gamma*1/2+w)*Iasr + tau*a_0*A*Ir - r*Iar + alpha*Ear - gammaC*Iar
    dIasr <- k*lambdaR*Ias - (gamma*1/2)*Iasr - (gamma*1/2+w)*Iasr + tau*a_0*A*Isr - r*Iasr + alpha*Easr - gammaC*Iasr
    
    # Recovered from SARS-CoV-2:
    
    dR <- -(lambdaS + lambdaR)*R + gamma*(Rs+Rr) - tau*a_0*R + r*Ra + gammaC*I
    dRs <- lambdaS*R - gamma*Rs - k*lambdaR*Rs + gamma*1/2*Rsr - tau*a_0*Rs + r*Ras + gammaC*Is
    dRr <- lambdaR*R - gamma*Rr - k*lambdaS*Rr + gamma*1/2*Rsr - tau*a_0*Rr + r*Rar + gammaC*Ir
    dRsr <- k*lambdaS*Rr + k*lambdaR*Rs - gamma*Rsr - tau*a_0*Rsr + r*Rasr + gammaC*Isr
    
    
    dRa <- -lambdaR*Ra + (gamma+w)*Ras + gamma*Rar + tau*a_0*R - r*Ra + gammaC*Ia
    dRas <-  - k*lambdaR*Ras - (gamma+w)*Ras + (gamma*1/2)*Rasr + tau*a_0*Rs - r*Ras + gammaC*Ias
    dRar <- lambdaR*Ra - gamma*Rar + (gamma*1/2+w)*Rasr + tau*a_0*Rr - r*Rar + gammaC*Iar
    dRasr <- k*lambdaR*Ras - (gamma*1/2)*Rasr - (gamma*1/2+w)*Rasr + tau*a_0*Rsr - r*Rasr + gammaC*Iasr
    
    # IPDs
    
    dTotal <- p*(Ss+Sr+Sas+Sar+Es+Er+Eas+Ear+Is+Ir+Ias+Iar+Rs+Rr+Ras+Rar) + p*(Ssr+Sasr+Esr+Easr+Isr+Iasr+Rsr+Rasr)
    
    dTotalS <- p*(Ss+Sas+Es+Eas+Is+Ias+Rs+Ras) + (p/2)*(Ssr+Sasr+Esr+Easr+Isr+Iasr+Rsr+Rasr)
    
    dTotalR <- p*(Sr+Sar+Er+Ear+Ir+Iar+Rr+Rar) + (p/2)*(Ssr+Sasr+Esr+Easr+Isr+Iasr+Rsr+Rasr)
    
  
  
    result <- c(dS, dSs, dSr, dSsr, dSa, dSas, dSar, dSasr, 
                dE, dEs, dEr, dEsr, dEa, dEas, dEar, dEasr, 
                dI, dIs, dIr, dIsr, dIa, dIas, dIar, dIasr,
                dR, dRs, dRr, dRsr, dRa, dRas, dRar, dRasr,
                dTotal, dTotalS, dTotalR)
    
    list(result, atb_0 = a_0, beta_0 = thetabeta_0)
    
    
  })
}



output_0 <- as.data.frame(ode(y = initial_state_values, 
                            times = times_model, 
                            func = model_0,
                            parms = parameters,
                            method = "euler"))




############################################## PLOTS #################################
values1 <- c ("0", "60", "120", "180", "240", "300", "360")
lock <- data.frame(xmin=4120, xmax=4210, ymin=-Inf, ymax=Inf) # assume lockdown lasts from day 4120 to 4210
colors <- c("S" = "cornflowerblue", "R" = "red", "none" = "grey", "dual" = "purple", "covid" = "orange")


# Bacterial carriage dynamics:
plotS0_a <- ggplot(output_0, aes(x=time)) + 
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
  #ggtitle("") + 
  theme(legend.position=c(.3,.7), #legend.position="right",#legend.position = c(.8,.88) 
        legend.background=element_blank(),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + 
  scale_color_manual(name = "", values = colors, 
                     labels = c("S", 
                                "R", 
                                "Total carriage",
                                "SR",
                                "SARS-CoV-2"))



colors2 <- c("S" = "cornflowerblue", "R" = "red")

Carriage_S0 <- output_0$Ss + output_0$Sas + output_0$Es + output_0$Eas + output_0$Is + output_0$Ias + output_0$Rs + output_0$Ras
Carriage_R0 <- output_0$Sr + output_0$Sar + output_0$Er + output_0$Ear + output_0$Ir + output_0$Iar + output_0$Rr + output_0$Rar
Carriage_SR0 <- output_0$Ssr + output_0$Sasr + output_0$Esr + output_0$Easr + output_0$Isr + output_0$Iasr + output_0$Rsr + output_0$Rasr


output_0$R_carriage <- (output_0$Sr + output_0$Sar + output_0$Er + output_0$Ear + output_0$Ir + output_0$Iar + output_0$Rr + output_0$Rar + 
                          0.5*(output_0$Ssr + output_0$Sasr + output_0$Esr + output_0$Easr + output_0$Isr + output_0$Iasr + output_0$Rsr + output_0$Rasr))/
  (output_0$Ss + output_0$Sas + output_0$Es + output_0$Eas + output_0$Is + output_0$Ias + output_0$Rs + output_0$Ras + 
     output_0$Ssr + output_0$Sasr + output_0$Esr + output_0$Easr + output_0$Isr + output_0$Iasr + output_0$Rsr + output_0$Rasr + 
     output_0$Sr + output_0$Sar + output_0$Er + output_0$Ear + output_0$Ir + output_0$Iar + output_0$Rr + output_0$Rar)

R_rateS0 <- (output_0$R_carriage[4210] - output_0$R_carriage[4000])/output_0$R_carriage[4000] 
R_rateS0*100 # % change in antibiotic resistance rate

rs0 <- round(R_rateS0*100, digits = 2)
rs0


# Antibiotic resistance rate:
plotS0_b <- ggplot(output_0, aes(x=time)) + 
  geom_line(aes(y=(Carriage_S0 + 0.5*(Carriage_SR0))/
                  (Carriage_S0 + Carriage_SR0 + Carriage_R0), color="S"), size=2) +  # S, blue
  geom_line(aes(y=(Carriage_R0 + 0.5*(Carriage_SR0))/
                  (Carriage_S0 + Carriage_SR0 + Carriage_R0), color="R"), size=2) +
  scale_x_continuous(limits = c(3995, 4365), breaks=seq(4000, 4365, 60), labels = values1) +
  scale_y_continuous(limits = c(0.3, 0.7), breaks=seq(0, 1, 0.05), labels = number) +
  theme_classic() + 
  labs(x="", y="") +
  #ggtitle("") + 
  theme(legend.position=c(.3,.3), #, #legend.position = c(.75,.65) 
        legend.background=element_blank(),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")) + 
  scale_color_manual(name = "", values = colors2, 
                     labels = c(expression(paste("S and ", frac(1,2), "SR  ")),
                                expression(paste("R and ", frac(1,2), "SR  ")))) +
  annotate(geom="text", x=4180, y=0.67, label=rs0, color="black", size = 8, fontface = 'bold') +
  annotate(geom="text", x=4220, y=0.67, label="%", color="black", size = 8, fontface = 'bold') 
#########################################################################################################

#---------- Cumulative IPD incidence short and long term:------------------------------------

# Cumulative ipd incidence (S, R, SR):
total_short_S0 <- output_0$Total[4210] - output_0$Total[4120]
total_long_S0 <- output_0$Total[4365] - output_0$Total[4000]

# IPD incidence caused by susceptible strain only (S + 0.5*SR):
S_short_S0 <- output_0$TotalS[4210] - output_0$TotalS[4120]
S_long_S0 <- output_0$TotalS[4365] - output_0$TotalS[4000]

# IPD incidence caused by resistant strain only (R + 0.5*SR):
R_short_S0 <- output_0$TotalR[4210] - output_0$TotalR[4120]
R_long_S0 <- output_0$TotalR[4365] - output_0$TotalR[4000]
