# Community model:
################################## S2 - decrease in atb use, no lockdown, S2_Sad a = 0.5 ###############

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
atb_2 <- as.data.frame(list(times = times, a_2 = rep(0, length(times))))
atb_2$a_2 <- ifelse((atb_2$times > 4119 & atb_2$times < 4211), 0.5, 1) 
atb_2
#Create the interpolating function
input_2 <- approxfun(atb_2)
input_2(seq(from = 0, to = 8000, by = 1))


# VARY RELATIVE BACTERIAL TRANSMISSION RISK IN THE COMMUNITY ACCORDING TO LOCKDOWN: thetabeta
# assume lockdown lasts from day 4120 to 4210
beta_2 <- as.data.frame(list(times = times, thetabeta_2 = rep(0, length(times))))
beta_2$thetabeta_2 <- ifelse((beta_2$times > 4119 & beta_2$times < 4211), 1, 1) 
beta_2
#Create the interpolating function
inputbeta_2 <- approxfun(beta_2)
inputbeta_2(seq(from = 0, to = 8000, by = 1))


################################### MODEL FUNCTION: ######################################
times_model <- seq(0, 8000, by = 1) # model time duration

# Model 1:
model_2 <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {  
    covid = if(floor(time)==4000) { covid = 2e-5 } else {covid = 0}  # introduce covid on day 4000
    
    a_2 <- input_2(time)
    thetabeta_2 <- inputbeta_2(time)

    
    N = S + Ss + Sr + Ssr + Sa + Sas + Sar + Sasr + 
      E + Es + Er + Esr + Ea + Eas + Ear + Easr + 
      I + Is + Ir + Isr + Ia + Ias + Iar + Iasr + 
      R + Rs + Rr + Rsr + Ra + Ras + Rar + Rasr  # sum all compartments
    
    #print(N)
    
    lambdaS <- betaS*thetabeta_2*(Ss/N + 1/2*Ssr/N + Sas/N + 1/2*Sasr/N + 
                                  Es/N + 1/2*Esr/N + Eas/N + 1/2*Easr/N + 
                                  q*(Is/N + 1/2*Isr/N + Ias/N + 1/2*Iasr/N) + 
                                  Rs/N + 1/2*Rsr/N + Ras/N + 1/2*Rasr/N)
    
    # betaR = betaS*(1-c)
    lambdaR <- betaS*(1-c)*thetabeta_2*(Sr/N + 1/2*Ssr/N + Sar/N + 1/2*Sasr/N + 
                                        Er/N + 1/2*Esr/N + Ear/N + 1/2*Easr/N + 
                                        q*(Ir/N + 1/2*Isr/N + Iar/N + 1/2*Iasr/N) + 
                                        Rr/N + 1/2*Rsr/N + Rar/N + 1/2*Rasr/N)
    
    lambdaC <- betaC*thetabeta_2*q*(I/N + Is/N + Ir/N + Isr/N + Ia/N + Ias/N + Iar/N + Iasr/N)   # covid force of infection
    
    
    
    # Equations:    
    
    # Susceptible to SARS-CoV-2:
    
    dS <- -(lambdaS + lambdaR)*S + gamma*(Ss+Sr) - tau*a_2*S + r*Sa - lambdaC*S - covid*N
    dSs <- lambdaS*S - gamma*Ss - k*lambdaR*Ss + gamma*1/2*Ssr - tau*a_2*Ss + r*Sas - lambdaC*Ss
    dSr <- lambdaR*S - gamma*Sr - k*lambdaS*Sr + gamma*1/2*Ssr - tau*a_2*Sr + r*Sar - lambdaC*Sr
    dSsr <- k*lambdaS*Sr + k*lambdaR*Ss - gamma*Ssr - tau*a_2*Ssr + r*Sasr - lambdaC*Ssr
    
    
    dSa <- -lambdaR*Sa + (gamma+w)*Sas + gamma*Sar + tau*a_2*S - r*Sa - lambdaC*Sa
    dSas <- -k*lambdaR*Sas - (gamma+w)*Sas + (gamma*1/2)*Sasr + tau*a_2*Ss - r*Sas - lambdaC*Sas 
    dSar <- lambdaR*Sa - gamma*Sar + (gamma*1/2+w)*Sasr + tau*a_2*Sr - r*Sar - lambdaC*Sar
    dSasr <- k*lambdaR*Sas - (gamma*1/2)*Sasr - (gamma*1/2+w)*Sasr + tau*a_2*Ssr - r*Sasr - lambdaC*Sasr
    
    
    # Exposed to SARS-CoV-2:
    
    dE <- -(lambdaS + lambdaR)*E + gamma*(Es+Er) - tau*a_2*E + r*Ea - alpha*E + lambdaC*S
    dEs <- lambdaS*E - gamma*Es - k*lambdaR*Es + gamma*1/2*Esr - tau*a_2*Es + r*Eas - alpha*Es + lambdaC*Ss
    dEr <- lambdaR*E - gamma*Er - k*lambdaS*Er + gamma*1/2*Esr - tau*a_2*Er + r*Ear - alpha*Er + lambdaC*Sr
    dEsr <- k*lambdaS*Er + k*lambdaR*Es - gamma*Esr - tau*a_2*Esr + r*Easr - alpha*Esr + lambdaC*Ssr
    
    
    dEa <- -lambdaR*Ea + (gamma+w)*Eas + gamma*Ear + tau*a_2*E - r*Ea + lambdaC*Sa - alpha*Ea
    dEas <- -k*lambdaR*Eas - (gamma+w)*Eas  + (gamma*1/2)*Easr + tau*a_2*Es - r*Eas + lambdaC*Sas - alpha*Eas
    dEar <- lambdaR*Ea - gamma*Ear + (gamma*1/2+w)*Easr + tau*a_2*Er - r*Ear + lambdaC*Sar - alpha*Ear
    dEasr <- k*lambdaR*Eas - (gamma*1/2)*Easr - (gamma*1/2+w)*Easr + tau*a_2*Esr - r*Easr + lambdaC*Sasr - alpha*Easr
    
    
    # Infected with SARS-CoV-2:
    
    dI <- -(lambdaS + lambdaR)*I + gamma*(Is+Ir) - tau*a_2*A*I + r*Ia + alpha*E - gammaC*I + covid*N
    dIs <- lambdaS*I - gamma*Is - k*lambdaR*Is + gamma*1/2*Isr - tau*a_2*A*Is + r*Ias + alpha*Es - gammaC*Is
    dIr <- lambdaR*I - gamma*Ir - k*lambdaS*Ir + gamma*1/2*Isr - tau*a_2*A*Ir + r*Iar + alpha*Er - gammaC*Ir
    dIsr <- k*lambdaS*Ir + k*lambdaR*Is - gamma*Isr - tau*a_2*A*Isr + r*Iasr + alpha*Esr - gammaC*Isr
    
    
    dIa <- -lambdaR*Ia + (gamma+w)*Ias + gamma*Iar + tau*a_2*A*I - r*Ia + alpha*Ea - gammaC*Ia
    dIas <-  -k*lambdaR*Ias - (gamma+w)*Ias + (gamma*1/2)*Iasr + tau*a_2*A*Is - r*Ias + alpha*Eas - gammaC*Ias
    dIar <- lambdaR*Ia - gamma*Iar + (gamma*1/2+w)*Iasr + tau*a_2*A*Ir - r*Iar + alpha*Ear - gammaC*Iar
    dIasr <- k*lambdaR*Ias - (gamma*1/2)*Iasr - (gamma*1/2+w)*Iasr + tau*a_2*A*Isr - r*Iasr + alpha*Easr - gammaC*Iasr
    
    # Recovered from SARS-CoV-2:
    
    dR <- -(lambdaS + lambdaR)*R + gamma*(Rs+Rr) - tau*a_2*R + r*Ra + gammaC*I
    dRs <- lambdaS*R - gamma*Rs - k*lambdaR*Rs + gamma*1/2*Rsr - tau*a_2*Rs + r*Ras + gammaC*Is
    dRr <- lambdaR*R - gamma*Rr - k*lambdaS*Rr + gamma*1/2*Rsr - tau*a_2*Rr + r*Rar + gammaC*Ir
    dRsr <- k*lambdaS*Rr + k*lambdaR*Rs - gamma*Rsr - tau*a_2*Rsr + r*Rasr + gammaC*Isr
    
    
    dRa <- -lambdaR*Ra + (gamma+w)*Ras + gamma*Rar + tau*a_2*R - r*Ra + gammaC*Ia
    dRas <-  - k*lambdaR*Ras - (gamma+w)*Ras + (gamma*1/2)*Rasr + tau*a_2*Rs - r*Ras + gammaC*Ias
    dRar <- lambdaR*Ra - gamma*Rar + (gamma*1/2+w)*Rasr + tau*a_2*Rr - r*Rar + gammaC*Iar
    dRasr <- k*lambdaR*Ras - (gamma*1/2)*Rasr - (gamma*1/2+w)*Rasr + tau*a_2*Rsr - r*Rasr + gammaC*Iasr
    
    # IPDs
    
    dTotal <- p*(Ss+Sr+Sas+Sar+Es+Er+Eas+Ear+Is+Ir+Ias+Iar+Rs+Rr+Ras+Rar) + p*(Ssr+Sasr+Esr+Easr+Isr+Iasr+Rsr+Rasr)
    
    dTotalS <- p*(Ss+Sas+Es+Eas+Is+Ias+Rs+Ras) + (p/2)*(Ssr+Sasr+Esr+Easr+Isr+Iasr+Rsr+Rasr)
    
    dTotalR <- p*(Sr+Sar+Er+Ear+Ir+Iar+Rr+Rar) + (p/2)*(Ssr+Sasr+Esr+Easr+Isr+Iasr+Rsr+Rasr)
    

    
    result <- c(dS, dSs, dSr, dSsr, dSa, dSas, dSar, dSasr, 
                dE, dEs, dEr, dEsr, dEa, dEas, dEar, dEasr, 
                dI, dIs, dIr, dIsr, dIa, dIas, dIar, dIasr,
                dR, dRs, dRr, dRsr, dRa, dRas, dRar, dRasr,
                dTotal, dTotalS, dTotalR)
    
    list(result, atb_2 = a_2, beta_2 = thetabeta_2)
    
    
  })
}



output_2 <- as.data.frame(ode(y = initial_state_values, 
                            times = times_model, 
                            func = model_2,
                            parms = parameters,
                            method = "euler"))


############################################## PLOTS #################################
values1 <- c ("0", "60", "120", "180", "240", "300", "360")
lock <- data.frame(xmin=4120, xmax=4210, ymin=-Inf, ymax=Inf)
colors <- c("S" = "cornflowerblue", "R" = "red", "none" = "grey", "dual" = "purple", "covid" = "orange")


Carriage_S2 <- output_2$Ss + output_2$Sas + output_2$Es + output_2$Eas + output_2$Is + output_2$Ias + output_2$Rs + output_2$Ras
Carriage_R2 <- output_2$Sr + output_2$Sar + output_2$Er + output_2$Ear + output_2$Ir + output_2$Iar + output_2$Rr + output_2$Rar
Carriage_SR2 <- output_2$Ssr + output_2$Sasr + output_2$Esr + output_2$Easr + output_2$Isr + output_2$Iasr + output_2$Rsr + output_2$Rasr


output_2$R_carriage <- (output_2$Sr + output_2$Sar + output_2$Er + output_2$Ear + output_2$Ir + output_2$Iar + output_2$Rr + output_2$Rar + 
                          0.5*(output_2$Ssr + output_2$Sasr + output_2$Esr + output_2$Easr + output_2$Isr + output_2$Iasr + output_2$Rsr + output_2$Rasr))/
  (output_2$Ss + output_2$Sas + output_2$Es + output_2$Eas + output_2$Is + output_2$Ias + output_2$Rs + output_2$Ras + 
     output_2$Ssr + output_2$Sasr + output_2$Esr + output_2$Easr + output_2$Isr + output_2$Iasr + output_2$Rsr + output_2$Rasr + 
     output_2$Sr + output_2$Sar + output_2$Er + output_2$Ear + output_2$Ir + output_2$Iar + output_2$Rr + output_2$Rar)

R_rateS2 <- (output_2$R_carriage[4210] - output_2$R_carriage[4000])/output_2$R_carriage[4000] 
R_rateS2*100 # % change in antibiotic resistance rate

rs2 <- round(R_rateS2*100, digits = 2)
rs2




# Bacterial carriage dynamics:
plotS2_a <- ggplot(output_2, aes(x=time)) + 
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
            color="transparent", fill = "forestgreen",
            alpha=0.2,
            inherit.aes = FALSE)


colors2 <- c("S" = "cornflowerblue", "R" = "red")

# Antibiotic resistance rate:
plotS2_b <- ggplot(output_2, aes(x=time)) + 
  geom_line(aes(y=(Carriage_S2 + 0.5*(Carriage_SR2))/
                  (Carriage_S2 + Carriage_SR2 + Carriage_R2), color="S"), size=2) +  # S, blue
  geom_line(aes(y=(Carriage_R2 + 0.5*(Carriage_SR2))/
                  (Carriage_S2 + Carriage_SR2 + Carriage_R2), color="R"), size=2) +
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
            color="transparent", fill = "forestgreen",
            alpha=0.2,
            inherit.aes = FALSE) +
  annotate(geom="text", x=4170, y=0.67, label=rs2, color="black", size = 8, fontface = 'bold') +
  annotate(geom="text", x=4230, y=0.67, label="%", color="black", size = 8, fontface = 'bold')

#---------- Cumulative IPD incidence short and long term:------------------------------------

# Cumulative ipd incidence (S, R, SR):
total_short_S2 <- output_2$Total[4210] - output_2$Total[4120]
total_long_S2 <- output_2$Total[4365] - output_2$Total[4000]

# IPD incidence caused by susceptible strain only (S + 0.5*SR):
S_short_S2 <- output_2$TotalS[4210] - output_2$TotalS[4120]
S_long_S2 <- output_2$TotalS[4365] - output_2$TotalS[4000]

# IPD incidence caused by resistant strain only (R + 0.5*SR):
R_short_S2 <- output_2$TotalR[4210] - output_2$TotalR[4120]
R_long_S2 <- output_2$TotalR[4365] - output_2$TotalR[4000]
