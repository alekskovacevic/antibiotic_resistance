# Author: Aleksandra Kovacevic 12/07/2023
# Code for co-circulation model of SARS-CoV-2 tranmission and infection, and pneumococcal colonization


# load the packages:
library(deSolve)
library(reshape2)
library(ggplot2)
library(scales)
library(readxl)
library(shiny)
library(xlsx)
#library(ggh4x)
library(dplyr)
library(cowplot)
library(ggpubr)
library(ggpattern)

library(facetscales)

# fixed parameters for pneumococcal carriage
q=1/2        # relative infectiousness with each strain for dually colonized
c=1/2        # fraction of dually colonized returning to single-colonized upon reinfection with S (R) 
k=1/2        # probability of acquiring secondary bacterial carriage
ps=1/2       # probability of transmitting S strain
psingle=1/2  # probability of a single infection


# Model inputs: 
# Total population size
N <- 100000

# pneumoccocal carriage assumed to be 20% where resistance represents 4/20 (20% AR)
initial_state <- c(S = N - 0.20*N, Ss = 0.16*N, Sr = 0.04*N, Ssr = 0, Sss = 0, Srr = 0,
                   Sa = 0, Ssa = 0, Sra = 0, Ssra = 0, Sssa = 0, Srra = 0,
                   E = 0, Es = 0, Er = 0, Esr = 0, Ess = 0, Err = 0,
                   Ea = 0, Esa = 0, Era = 0, Esra = 0, Essa = 0, Erra = 0,
                   I = 0, Is = 0, Ir = 0, Isr = 0, Iss = 0, Irr = 0,
                   Ia = 0, Isa = 0, Ira = 0, Isra = 0, Issa = 0, Irra = 0,
                   Iabx = 0, Isabx = 0, Irabx = 0, Israbx = 0, Issabx = 0, Irrabx = 0, # infected + prophylaxis
                   R = 0, Rs = 0, Rr = 0, Rsr = 0, Rss = 0, Rrr = 0,
                   Ra = 0, Rsa = 0, Rra = 0, Rsra = 0, Rssa = 0, Rrra = 0,
                   Rabx = 0, Rsabx = 0, Rrabx = 0, Rsrabx = 0, Rssabx = 0, Rrrabx = 0, # recovered + prophylaxis
                   IPDs = 0, IPDr = 0, IPD = 0, ATB = 0) 


parameters_15 <- c(betaS = 0.042,    # pneumococcal transmission rate 
                   gammaS = 1/30,    # pneumococcal recovery rate (sensitive strain) 
                   gammaR = 1/30,    # pneumococcal recovery rate (resistant strain) 
                   betaC = 0.46,     # covid transmission rate
                   alpha = 1/5,      # covid incubation rate 5 days, latent period
                   gammaC = 1/7,     # covid recovery rate 7 days
                   fit= 0.949,       # fitness cost to be at pre-pandemic equilibrium 
                   w = 1/3,          # rate of antibiotic action; takes 3 days for antibiotics to clear susceptible carriage
                   tau = 0.0014,     # antibiotic prescription rate (per capita), daily rate
                   r = 1/7,          # baseline antibiotic treatment duration
                   r_abx = 1/11.5    # 11.5 (if 3-day treatment) 
)

times <- seq(0, 3650, by = 1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters to change according to the scenario:
p_abx_15 = 0.15                  # proportion of infected individuals getting abx prophylaxis
psummer_15 = 1e-06              # summer invasion rate
pwinter_15 = 2.5e-06            # winter invasion rate 
a_lockdown_15 = 0.77            # reduction in the nb of abx prescriptions
covid_transm_lockdown_15 = 0.23  # reduction of the covid transmission during lockdown 
IPD_risk_lockdown_15 = 0.2      # risk of IPD during lockdown 
IPD_risk_rest_15 = 0.4          # risk of IPD during the rest of the year 
covid.present_15 = 1            # presence (1)/absence (0) of covid (pre-pandemic and pandemic period)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Interpolating functions:

#Pneumococcal invasion rate (summer vs winter invasion rate) (p):
p_15 <- as.data.frame(list(times = times, p_inv_15 = rep(0, length(times))))
p_15$p_inv_15 <- ifelse((p_15$times > 1089 & p_15$times < 1271), (psummer_15), (pwinter_15)) 
invasion_15 <- approxfun(p_15)       #Create the interpolating function
invasion_15(seq(from = 1, to = 3650, by = 1))

# EFFECT OF LOCKDOWN on antibiotic use in the community - during lockdown and over the entire year (a)
atb_15 <- as.data.frame(list(times = times, a_15 = rep(0, length(times))))   
atb_15$a_15 <- ifelse(atb_15$times > 1074 & atb_15$times < 1365, a_lockdown_15, 1)
input_15 <- approxfun(atb_15)
input_15(seq(from = 0, to = 3650, by = 1)) ##################### 

# SARS-CoV-2 transmission risk (during lockdown vs without)
beta_15 <- as.data.frame(list(times = times, thetabeta_15 = rep(0, length(times))))
beta_15$thetabeta_15 <- ifelse(beta_15$times < 1075, 1,
                               ifelse((beta_15$times > 1074 & beta_15$times < 1136), 
                                      covid_transm_lockdown_15, 
                                      ifelse(beta_15$times > 1135 & beta_15$times < 1231, 0.45, 0.85))) 
inputbeta_15 <- approxfun(beta_15)
inputbeta_15(seq(from = 0, to = 3650, by = 1))

# Pneumococcal transmission risk during lockdown and the rest of the year 
beta_15p <- as.data.frame(list(times = times, thetab_15p = rep(0, length(times))))
beta_15p$thetab_15p <- ifelse((beta_15p$times > 1074 & beta_15p$times < 1136|beta_15p$times > 1304 & beta_15p$times < 1366), 1, 1) 
inputbeta_15p <- approxfun(beta_15p)
inputbeta_15p(seq(from = 0, to = 3650, by = 1))

# IPD disease risk reduction during NPI implementation before and during lockdown 
risk_15 <- as.data.frame(list(times = times, IPDrisk_15 = rep(0, length(times))))
risk_15$IPDrisk_15 <- ifelse(risk_15$times < 1075, 1,
                             ifelse(risk_15$times > 1074 & risk_15$times < 1136, IPD_risk_lockdown_15, IPD_risk_rest_15)) 
inputrisk_15 <- approxfun(risk_15)
inputrisk_15(seq(from = 1, to = 3650, by = 1))




sim_time <- seq(from = 0, to = 3650, by = 1)  # model simulation time


################################### MODEL FUNCTION: ######################################

model_15 <- function(time, state, parameters) {  
  with(as.list(c(state, parameters)), { 
    covid = if(floor(time) == 1000) { covid = covid.present_15 } else {covid = 0}  # COVID
    
    inv_rate_15 <- invasion_15(time)       # invasion rate varies between summer and winter
    a_15 <- input_15(time)                   # antibiotic use in the community
    thetabeta_15 <- inputbeta_15(time)       # covid transmission risk 
    thetab_15 <- inputbeta_15p(time)            # carriage transmission risk
    IPDrisk_15 <- inputrisk_15(time)           # IPD risk (varies with the present or absence of ILIs)
    
    N <- S+Ss+Sr+Ssr+Sss+Srr+Sa+Ssa+Sra+Ssra+Sssa+Srra+
      E+Es+Er+Esr+Ess+Err+Ea+Esa+Era+Esra+Essa+Erra+ 
      I+Is+Ir+Isr+Iss+Irr+Ia+Isa+Ira+Isra+Issa+Irra+Iabx+Isabx+Irabx+Israbx+Issabx+Irrabx+
      R+Rs+Rr+Rsr+Rss+Rrr+Ra+Rsa+Rra+Rsra+Rssa+Rrra+Rabx+Rsabx+Rrabx+Rsrabx+Rssabx+Rrrabx
    
    #print(sum(N))
    
    # Forces of infection:
    # Forces of infection (competition between strains for dual carriage)
    lambdaS <- betaS*thetab_15*(Ss + 2*q*ps*Ssr + 2*q*Sss + Ssa + 2*q*ps*Ssra + 2*q*Sssa + 
                                  Es + 2*q*ps*Esr + 2*q*Ess + Esa + 2*q*ps*Esra + 2*q*Essa + 
                                  Is + 2*q*ps*Isr + 2*q*Iss + Isa + 2*q*ps*Isra + 2*q*Issa +
                                  Isabx + 2*q*ps*Israbx + 2*q*Issabx +
                                  Rs + 2*q*ps*Rsr + 2*q*Rss + Rsa + 2*q*ps*Rsra + 2*q*Rssa +
                                  Rsabx + 2*q*ps*Rsrabx + 2*q*Rssabx)/N
    
    
    lambdaR <- betaS*fit*thetab_15*(Sr + 2*q*(1-ps)*Ssr + 2*q*Srr + Sra + 2*q*(1-ps)*Ssra + 2*q*Srra +
                                      Er + 2*q*(1-ps)*Esr + 2*q*Err + Era + 2*q*(1-ps)*Esra + 2*q*Erra +
                                      Ir + 2*q*(1-ps)*Isr + 2*q*Irr + Ira + 2*q*(1-ps)*Isra + 2*q*Irra +
                                      Irabx + 2*q*(1-ps)*Israbx + 2*q*Irrabx +
                                      Rr + 2*q*(1-ps)*Rsr + 2*q*Rrr + Rra + 2*q*(1-ps)*Rsra + 2*q*Rrra + 
                                      Rrabx + 2*q*(1-ps)*Rsrabx + 2*q*Rrrabx)/N
    
    
    lambdaXS <- betaS*thetab_15*(Ss + psingle*(2*q*ps*Ssr + 2*q*Sss) + Ssa + psingle*(2*q*ps*Ssra + 2*q*Sssa) + 
                                   Es + psingle*(2*q*ps*Esr + 2*q*Ess) + Esa + psingle*(2*q*ps*Esra + 2*q*Essa) + 
                                   Is + psingle*(2*q*ps*Isr + 2*q*Iss) + Isa + psingle*(2*q*ps*Isra + 2*q*Issa) + 
                                   Isabx + psingle*(2*q*ps*Israbx + 2*q*Issabx) +
                                   Rs + psingle*(2*q*ps*Rsr + 2*q*Rss) + Rsa + psingle*(2*q*ps*Rsra + 2*q*Rssa) + 
                                   Rsabx + psingle*(2*q*ps*Rsrabx + 2*q*Rssabx))/N
    
    
    lambdaXR <- betaS*fit*thetab_15*(Sr + psingle*(2*q*(1-ps)*Ssr + 2*q*Srr) + Sra + psingle*(2*q*(1-ps)*Ssra + 2*q*Srra) + 
                                       Er + psingle*(2*q*(1-ps)*Esr + 2*q*Err) + Era + psingle*(2*q*(1-ps)*Esra + 2*q*Erra) + 
                                       Ir + psingle*(2*q*(1-ps)*Isr + 2*q*Irr) + Ira + psingle*(2*q*(1-ps)*Isra + 2*q*Irra) + 
                                       Irabx + psingle*(2*q*(1-ps)*Israbx + 2*q*Irrabx) + 
                                       Rr + psingle*(2*q*(1-ps)*Rsr + 2*q*Rrr) + Rra + psingle*(2*q*(1-ps)*Rsra + 2*q*Rrra) + 
                                       Rrabx + psingle*(2*q*(1-ps)*Rsrabx + 2*q*Rrrabx))/N
    
    
    lambdaXSR <- ((1-psingle)*(betaS*thetab_15*ps + betaS*thetab_15*fit*(1-ps))*2*q*(Ssr+Ssra) + 
                    (1-psingle)*(betaS*thetab_15*ps + betaS*thetab_15*fit*(1-ps))*2*q*(Esr+Esra) + 
                    (1-psingle)*(betaS*thetab_15*ps + betaS*thetab_15*fit*(1-ps))*2*q*(Isr+Isra+Israbx) +
                    (1-psingle)*(betaS*thetab_15*ps + betaS*thetab_15*fit*(1-ps))*2*q*(Rsr+Rsra+Rsrabx))/N
    
    lambdaXSS <- ((1-psingle)*betaS*thetab_15*2*q*(Sss+Sssa) + 
                    (1-psingle)*betaS*thetab_15*2*q*(Ess+Essa) + 
                    (1-psingle)*betaS*thetab_15*2*q*(Iss+Issa+Issabx) + 
                    (1-psingle)*betaS*thetab_15*2*q*(Rss+Rssa+Rssabx))/N
    
    lambdaXRR <- ((1-psingle)*betaS*thetab_15*fit*2*q*(Srr+Srra) + 
                    (1-psingle)*betaS*thetab_15*fit*2*q*(Err+Erra) + 
                    (1-psingle)*betaS*thetab_15*fit*2*q*(Irr+Irra+Irrabx) + 
                    (1-psingle)*betaS*thetab_15*fit*2*q*(Rrr+Rrra+Rrrabx))/N
    
    # Force of infection (SARS-CoV-2)
    lambdaC <- betaC*thetabeta_15*((I+Is+Ir+Isr+Iss+Irr+Ia+Isa+Ira+Isra+Issa+Irra+
                                      Iabx+Isabx+Irabx+Israbx+Issabx+Irrabx)/N)
    
    ########################################## Model equations: #############################################
    
    # Susceptible to SARS-CoV-2: 
    dS <- -lambdaXS*S - lambdaXR*S - lambdaXSR*S - lambdaXSS*S - lambdaXRR*S + 
      gammaS*Ss + gammaS*Sss + gammaR*Sr + gammaR*Ssr + gammaR*Srr + 
      gammaR*Sra + gammaR*Srra + (gammaS+w)*Ssa + (gammaS+w)*Sssa - 
      tau*a_15*S + r*Sa - lambdaC*S - covid  
    
    dSs <- lambdaXS*S - gammaS*Ss - k*lambdaS*Ss - k*lambdaR*Ss - tau*a_15*Ss - lambdaC*Ss  
    dSr <- lambdaXR*S - gammaR*Sr - k*lambdaS*Sr - k*lambdaR*Sr - tau*a_15*Sr + r*Sra + (gammaS+w)*Ssra - lambdaC*Sr 
    
    dSsr <- lambdaXSR*S - gammaR*Ssr + k*lambdaR*Ss + k*lambdaS*Sr - 
      k*c*lambdaS*Ssr - k*c*lambdaR*Ssr + 2*k*c*lambdaR*Sss + 2*k*c*lambdaS*Srr - 
      tau*a_15*Ssr - lambdaC*Ssr
    
    dSss <- lambdaXSS*S - gammaS*Sss + k*lambdaS*Ss + k*c*lambdaS*Ssr - 2*k*c*lambdaR*Sss - 
      tau*a_15*Sss - lambdaC*Sss
    
    dSrr <- lambdaXRR*S - gammaR*Srr + k*lambdaR*Sr + k*c*lambdaR*Ssr - 2*k*c*lambdaS*Srr - 
      tau*a_15*Srr + r*Srra - lambdaC*Srr
    
    # Antibiotic exposed compartments:
    dSa <- -lambdaXR*Sa - lambdaXRR*Sa + tau*a_15*S - r*Sa - lambdaC*Sa
    dSsa <- -(gammaS+w)*Ssa + tau*a_15*Ss - k*lambdaR*Ssa - lambdaC*Ssa
    dSra <- - gammaR*Sra + lambdaXR*Sa + tau*a_15*Sr - r*Sra - lambdaC*Sra
    dSsra <- -(gammaS+w)*Ssra + k*lambdaR*Ssa + 2*k*c*lambdaR*Sssa + tau*a_15*Ssr - lambdaC*Ssra
    dSssa <- -(gammaS+w)*Sssa - 2*k*c*lambdaR*Sssa + tau*a_15*Sss - lambdaC*Sssa
    dSrra <-  lambdaXRR*Sa + tau*a_15*Srr - r*Srra - gammaR*Srra - lambdaC*Srra
    
    
    # Exposed to SARS-CoV-2: 
    dE <- -lambdaXS*E - lambdaXR*E - lambdaXSR*E - lambdaXSS*E - lambdaXRR*E + 
      gammaS*Es + gammaS*Ess + gammaR*Er + gammaR*Esr + gammaR*Err +
      (gammaS+w)*Esa + (gammaS+w)*Essa + gammaR*Era + gammaR*Erra - 
      tau*a_15*E + r*Ea + lambdaC*S - alpha*E   
    
    dEs <- lambdaXS*E - gammaS*Es - k*lambdaS*Es - k*lambdaR*Es - tau*a_15*Es + lambdaC*Ss - alpha*Es
    dEr <- lambdaXR*E - gammaR*Er - k*lambdaS*Er - k*lambdaR*Er - tau*a_15*Er + r*Era + lambdaC*Sr - alpha*Er + (gammaS+w)*Esra
    
    dEsr <- lambdaXSR*E - gammaR*Esr + k*lambdaR*Es + k*lambdaS*Er - 
      k*c*lambdaS*Esr - k*c*lambdaR*Esr + 2*k*c*lambdaR*Ess + 2*k*c*lambdaS*Err - 
      tau*a_15*Esr + lambdaC*Ssr - alpha*Esr 
    
    dEss <- lambdaXSS*E - gammaS*Ess + k*lambdaS*Es + k*c*lambdaS*Esr - 2*k*c*lambdaR*Ess - 
      tau*a_15*Ess + lambdaC*Sss - alpha*Ess
    
    dErr <- lambdaXRR*E - gammaR*Err + k*lambdaR*Er + k*c*lambdaR*Esr - 2*k*c*lambdaS*Err - 
      tau*a_15*Err + r*Erra + lambdaC*Srr - alpha*Err
    
    # Antibiotic exposed compartments:
    dEa <- -lambdaXR*Ea - lambdaXRR*Ea + tau*a_15*E - r*Ea + lambdaC*Sa - alpha*Ea
    dEsa <- -(gammaS+w)*Esa - k*lambdaR*Esa + tau*a_15*Es + lambdaC*Ssa - alpha*Esa
    dEra <- lambdaXR*Ea - gammaR*Era + tau*a_15*Er - r*Era + lambdaC*Sra - alpha*Era
    dEsra <- -(gammaS+w)*Esra + k*lambdaR*Esa + 2*k*c*lambdaR*Essa + tau*a_15*Esr + lambdaC*Ssra - alpha*Esra
    dEssa <- -(gammaS+w)*Essa - 2*k*c*lambdaR*Essa + tau*a_15*Ess + lambdaC*Sssa - alpha*Essa
    dErra <- lambdaXRR*Ea + tau*a_15*Err - r*Erra - gammaR*Erra + lambdaC*Srra - alpha*Erra
    
    
    # Infected with SARS-CoV-2: 
    dI <- -lambdaXS*I - lambdaXR*I - lambdaXSR*I - lambdaXSS*I - lambdaXRR*I + 
      gammaS*Is + gammaS*Iss + gammaR*Ir + gammaR*Isr + gammaR*Irr + 
      (gammaS+w)*Isa + (gammaS+w)*Issa + gammaR*Ira + gammaR*Irra - 
      (1-p_abx_15)*tau*a_15*I + r*Ia + alpha*E + covid - gammaC*I - p_abx_15*I + 
      gammaR*Irabx + gammaR*Irrabx
    
    dIs <- lambdaXS*I - gammaS*Is - k*lambdaS*Is - k*lambdaR*Is - (1-p_abx_15)*tau*a_15*Is + alpha*Es - gammaC*Is - p_abx_15*Is
    dIr <- lambdaXR*I - gammaR*Ir - k*lambdaS*Ir - k*lambdaR*Ir + (gammaS+w)*Isra - (1-p_abx_15)*tau*a_15*Ir + r*Ira + alpha*Er - gammaC*Ir - p_abx_15*Ir 
    
    dIsr <- lambdaXSR*I - gammaR*Isr + k*lambdaR*Is + k*lambdaS*Ir - 
      k*c*lambdaS*Isr - k*c*lambdaR*Isr + 2*k*c*lambdaR*Iss + 2*k*c*lambdaS*Irr - 
      (1-p_abx_15)*tau*a_15*Isr + alpha*Esr - gammaC*Isr - p_abx_15*Isr 
    
    dIss <- lambdaXSS*I - gammaS*Iss + k*lambdaS*Is + k*c*lambdaS*Isr - 2*k*c*lambdaR*Iss - 
      (1-p_abx_15)*tau*a_15*Iss + alpha*Ess - gammaC*Iss - p_abx_15*Iss 
    
    dIrr <- lambdaXRR*I - gammaR*Irr + k*lambdaR*Ir + k*c*lambdaR*Isr - 2*k*c*lambdaS*Irr - 
      (1-p_abx_15)*tau*a_15*Irr + r*Irra + alpha*Err - gammaC*Irr - p_abx_15*Irr 
    
    # Antibiotic exposed compartments:
    dIa <- -lambdaXR*Ia - lambdaXRR*Ia + (1-p_abx_15)*tau*a_15*I - r*Ia + alpha*Ea - gammaC*Ia 
    dIsa <- -(gammaS+w)*Isa - k*lambdaR*Isa + (1-p_abx_15)*tau*a_15*Is + alpha*Esa - gammaC*Isa
    dIra <- lambdaXR*Ia - gammaR*Ira + (1-p_abx_15)*tau*a_15*Ir - r*Ira + alpha*Era - gammaC*Ira
    dIsra <- -(gammaS+w)*Isra + k*lambdaR*Isa + 2*k*c*lambdaR*Issa + (1-p_abx_15)*tau*a_15*Isr + alpha*Esra - gammaC*Isra 
    dIssa <- -(gammaS+w)*Issa - 2*k*c*lambdaR*Issa + (1-p_abx_15)*tau*a_15*Iss + alpha*Essa - gammaC*Issa  
    dIrra <-  lambdaXRR*Ia + (1-p_abx_15)*tau*a_15*Irr - r*Irra - gammaR*Irra + alpha*Erra - gammaC*Irra 
    
    # Additional antibiotic (abx) exposed compartments:
    dIabx <- p_abx_15*I - gammaC*Iabx - lambdaXR*Iabx - lambdaXRR*Iabx + (gammaS+w)*Isabx + (gammaS+w)*Issabx 
    dIsabx <- p_abx_15*Is - gammaC*Isabx - (gammaS+w)*Isabx - k*lambdaR*Isabx   
    dIrabx <- p_abx_15*Ir - gammaC*Irabx + lambdaXR*Iabx - gammaR*Irabx + (gammaS+w)*Israbx 
    dIsrabx <- p_abx_15*Isr - gammaC*Israbx - (gammaS+w)*Israbx + k*lambdaR*Isabx + 2*k*c*lambdaR*Issabx  
    dIssabx <- p_abx_15*Iss - gammaC*Issabx - (gammaS+w)*Issabx - 2*k*c*lambdaR*Issabx   
    dIrrabx <- p_abx_15*Irr - gammaC*Irrabx + lambdaXRR*Iabx - gammaR*Irrabx 
    
    
    # Recovered from SARS-CoV-2: 
    dR <- -lambdaXS*R - lambdaXR*R - lambdaXSR*R - lambdaXSS*R - lambdaXRR*R + 
      gammaS*Rs + gammaS*Rss + gammaR*Rr + gammaR*Rsr + gammaR*Rrr + 
      gammaR*Rra + gammaR*Rrra + (gammaS+w)*Rsa + (gammaS+w)*Rssa - 
      tau*a_15*R + r*Ra + gammaC*I + r_abx*Rabx + 
      gammaR*Rrabx + gammaR*Rrrabx
    
    dRs <- lambdaXS*R - gammaS*Rs - k*lambdaS*Rs - k*lambdaR*Rs - tau*a_15*Rs + gammaC*Is + r_abx*Rsabx
    dRr <- lambdaXR*R - gammaR*Rr - k*lambdaS*Rr - k*lambdaR*Rr + (gammaS+w)*Rsra - tau*a_15*Rr + r*Rra + gammaC*Ir + r_abx*Rrabx
    
    dRsr <- lambdaXSR*R - gammaR*Rsr + k*lambdaR*Rs + k*lambdaS*Rr - 
      k*c*lambdaS*Rsr - k*c*lambdaR*Rsr + 2*k*c*lambdaR*Rss + 2*k*c*lambdaS*Rrr - 
      tau*a_15*Rsr + gammaC*Isr + r_abx*Rsrabx
    
    dRss <- lambdaXSS*R - gammaS*Rss + k*lambdaS*Rs + k*c*lambdaS*Rsr - 2*k*c*lambdaR*Rss - 
      tau*a_15*Rss + gammaC*Iss + r_abx*Rssabx
    
    dRrr <- lambdaXRR*R - gammaR*Rrr + k*lambdaR*Rr + k*c*lambdaR*Rsr - 2*k*c*lambdaS*Rrr - 
      tau*a_15*Rrr + r*Rrra + gammaC*Irr + r_abx*Rrrabx
    
    # Antibiotic exposed compartments:
    dRa <- gammaC*Ia + tau*a_15*R - r*Ra - lambdaXR*Ra - lambdaXRR*Ra 
    dRsa <- gammaC*Isa + tau*a_15*Rs - (gammaS+w)*Rsa - k*lambdaR*Rsa
    dRra <- gammaC*Ira + tau*a_15*Rr - r*Rra + lambdaXR*Ra - gammaR*Rra
    dRsra <- gammaC*Isra + tau*a_15*Rsr - (gammaS+w)*Rsra + k*lambdaR*Rsa + 2*k*c*lambdaR*Rssa 
    dRssa <- gammaC*Issa + tau*a_15*Rss - (gammaS+w)*Rssa - 2*k*c*lambdaR*Rssa
    dRrra <- gammaC*Irra + tau*a_15*Rrr - r*Rrra + lambdaXRR*Ra - gammaR*Rrra 
    
    # Additional antibiotic (abx) exposed compartments:
    dRabx <- gammaC*Iabx - r_abx*Rabx - lambdaXR*Rabx - lambdaXRR*Rabx + (gammaS+w)*Rsabx + (gammaS+w)*Rssabx
    dRsabx <- gammaC*Isabx - r_abx*Rsabx - (gammaS+w)*Rsabx - k*lambdaR*Rsabx
    dRrabx <- gammaC*Irabx - r_abx*Rrabx + lambdaXR*Rabx - gammaR*Rrabx + (gammaS+w)*Rsrabx
    dRsrabx <- gammaC*Israbx - r_abx*Rsrabx - (gammaS+w)*Rsrabx + k*lambdaR*Rsabx + 2*k*c*lambdaR*Rssabx 
    dRssabx <- gammaC*Issabx - r_abx*Rssabx - (gammaS+w)*Rssabx - 2*k*c*lambdaR*Rssabx
    dRrrabx <- gammaC*Irrabx - r_abx*Rrrabx + lambdaXRR*Rabx - gammaR*Rrrabx  
    
    # Cumulative incidence of invasive pneumococcal disease:
    
    # IPD caused by S cannot be caused by individuals receiving antibiotic treatment:
    dIPDs <- inv_rate_15*IPDrisk_15*(Ss+Sss+0.5*Ssr+Es+Ess+0.5*Esr+Is+Iss+0.5*Isr+Rs+Rss+0.5*Rsr)
    
    # IPD caused by R can be caused by all individuals carrying R:
    dIPDr <- inv_rate_15*IPDrisk_15*(Sr+Srr+0.5*Ssr+Er+Err+0.5*Esr+Ir+Irr+0.5*Isr+Rr+Rrr+0.5*Rsr+
                                       Sra+Srra+Ssra+Era+Erra+Esra+Ira+Irra+Isra+Rra+Rrra+Rsra+
                                       Irabx+Irrabx+Israbx+Rrabx+Rrrabx+Rsrabx)
    
    dIPD <- inv_rate_15*IPDrisk_15*(Ss+Sss+0.5*Ssr+Es+Ess+0.5*Esr+Is+Iss+0.5*Isr+Rs+Rss+0.5*Rsr)+
      inv_rate_15*IPDrisk_15*(Sr+Srr+0.5*Ssr+Er+Err+0.5*Esr+Ir+Irr+0.5*Isr+Rr+Rrr+0.5*Rsr+
                                Sra+Srra+Ssra+Era+Erra+Esra+Ira+Irra+Isra+Rra+Rrra+Rsra+
                                Irabx+Irrabx+Israbx+Rrabx+Rrrabx+Rsrabx)
    
    # Cumulative number of antibiotic exposed individuals:
    dATB <- tau*a_15*(S+Ss+Ssr+Ssr+Sss+Srr+E+Es+Er+Esr+Ess+Err) +
      (1-p_abx_15)*tau*a_15*(I+Is+Ir+Isr+Iss+Irr) + 
      p_abx_15*tau*a_15*(I+Is+Ir+Isr+Iss+Irr) + 
      tau*a_15*(R+Rs+Rr+Rsr+Rss+Rrr) 
    
    return(list(c(dS, dSs, dSr, dSsr, dSss, dSrr, dSa, dSsa, dSra, dSsra, dSssa, dSrra,  
                  dE, dEs, dEr, dEsr, dEss, dErr, dEa, dEsa, dEra, dEsra, dEssa, dErra,
                  dI, dIs, dIr, dIsr, dIss, dIrr, dIa, dIsa, dIra, dIsra, dIssa, dIrra,
                  dIabx, dIsabx, dIrabx, dIsrabx, dIssabx, dIrrabx,
                  dR, dRs, dRr, dRsr, dRss, dRrr, dRa, dRsa, dRra, dRsra, dRssa, dRrra, 
                  dRabx, dRsabx, dRrabx, dRsrabx, dRssabx, dRrrabx, 
                  dIPDs, dIPDr, dIPD, dATB), p_abx_15=p_abx_15,beta_15=thetabeta_15)) 
  })
  
}

# Solving the differential equations using the ode integration algorithm
output_15 <- as.data.frame(ode(y = initial_state, 
                               times = sim_time, 
                               func = model_15, 
                               parms = parameters_15, 
                               method = "euler"))
output_15
#####################################################################################################################








##################### Pneumococcal carriage prevalence over time:
# Total pneumococcal carriage
out15 <- dplyr::mutate(output_15, PC_15 = Ss + Sr + Ssr + Sss + Srr + Ssa + Sra + Ssra + Sssa + Srra + 
                        Es + Er + Esr + Ess + Err + Esa + Era + Esra + Essa + Erra + 
                        Is + Ir + Isr + Iss + Irr + Isa + Ira + Isra + Issa + Irra + 
                        Isabx + Irabx + Israbx + Issabx + Irrabx +
                        Rs + Rr + Rsr + Rss + Rrr + Rsa + Rra + Rsra + Rssa + Rrra +
                        Rsabx + Rrabx + Rsrabx + Rssabx + Rrrabx,
                      R_15 = Sr + 0.5*Ssr + Srr + Sra + 0.5*Ssra + Srra + 
                        Er + 0.5*Esr + Err + Era + 0.5*Esra + Erra + 
                        Ir + 0.5*Isr + Irr + Ira + 0.5*Isra + Irra + 
                        Irabx + 0.5*Israbx + Irrabx +
                        Rr + 0.5*Rsr + Rrr + Rra + 0.5*Rsra + Rrra + 
                        Rrabx + 0.5*Rsrabx + Rrrabx,
                      S_15 = Ss + 0.5*Ssr + Sss + Ssa + 0.5*Ssra + Sssa + 
                        Es + 0.5*Esr + Ess + Esa + 0.5*Esra + Essa + 
                        Is + 0.5*Isr + Iss + Isa + 0.5*Isra + Issa + 
                        Isabx + 0.5*Israbx + Issabx + 
                        Rs + 0.5*Rsr + Rss + Rsa + 0.5*Rsra + Rssa + 
                        Rsabx + 0.5*Rsrabx + Rssabx)

# prevalence of the resistant carriage at the beginning of the year 1000 (%)
R_prev_1000_15 <- round((out15$R_15[1000]/out15$PC_15[1000])*100, digits = 2)
R_prev_1000_15
# prevalence of the resistant carriage at the end of lockdown 1135 (%)
R_prev_1135_15 <- round((out15$R_15[1135]/out15$PC_15[1135])*100, digits = 2)
R_prev_1135_15
# prevalence of the resistant carriage at the end of the year 1365 (%)
R_prev_1365_15 <- round((out15$R_15[1365]/out15$PC_15[1365])*100, digits = 2)
R_prev_1365_15



plot(out15$R_15/out15$PC_15, type = "l")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~ Change in the resistance prevalence ~~~~~~~~~~~~~~~~~~~~~~~~~~~
resistant_percent_change_lockdown_15 <- round(R_prev_1135_15 - R_prev_1000_15, digits = 2)
print('% change in the prevalence of resistant carriage at the end of the lockdown')
resistant_percent_change_lockdown_15

# change in the prevalence of RESISTANT bacterial carriage at the end of the YEAR
resistant_percent_change_15 <- round(R_prev_1365_15 - R_prev_1000_15, digits = 2)
print('% change in the prevalence of resistant carriage at the end of the year')
resistant_percent_change_15 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


############################################## PLOTS #############################################
values1 <- c ("0", "60", "120", "180", "240", "300", "360")
lock <- data.frame(xmin=1075, xmax=1135, ymin=-Inf, ymax=Inf) # assume lockdown lasts from day 4120 to 4210
colors <- c("S" = "cornflowerblue", "R" = "red", "Carriers" = "grey", "SR" = "purple", "SS" = "blue", "RR" = "brown")

p1 <- ggplot(output_15, aes(x=time)) + 
  #geom_line(aes(y=(S+E+I+R+Sa+Ea+Ia+Ra+Iabx+Rabx)/N, color="Non-carriers"), size=1) +  # S, blue
  geom_line(aes(y=(Ss + Sr + Ssr + Sss + Srr + Ssa + Sra + Ssra + Sssa + Srra + 
                     Es + Er + Esr + Ess + Err + Esa + Era + Esra + Essa + Erra +
                     Is + Ir + Isr + Iss + Irr + Isa + Ira + Isra + Issa + Irra + 
                     Isabx + Irabx + Israbx + Issabx + Irrabx + 
                     Rs + Rr + Rsr + Rss + Rrr + Rsa + Rra + Rsra + Rssa + Rrra + 
                     Rsabx + Rrabx + Rsrabx + Rssabx + Rrrabx)/N, color="Carriers"), size=1) +  # total carriage, black
  geom_line(aes(y=(Ss+Es+Is+Rs+Ssa+Esa+Isa+Rsa+Isabx+Rsabx)/N, color="S"), size=1) +  # S, blue
  geom_line(aes(y=(Sr+Er+Ir+Rr+Sra+Era+Ira+Rra+Irabx+Rrabx)/N, color="R"), size=1) +  # R, red
  geom_line(aes(y=(Ssr+Esr+Isr+Rsr+Ssra+Esra+Isra+Rsra+Israbx+Rsrabx)/N, color="SR"), size=1) +  # SR, purple
  geom_line(aes(y=(Sss+Ess+Iss+Rss+Sssa+Essa+Issa+Rssa+Issabx+Rssabx)/N, color="SS"), size=1) +  # SS, blue
  geom_line(aes(y=(Srr+Err+Irr+Rrr+Srra+Erra+Irra+Rrra+Irrabx+Rrrabx)/N, color="RR"), size=1) +  # RR, red
  #geom_line(aes(y=(I+Is+Ir+Isr+Iss+Irr+Ia+Isa+Ira+Isra+Issa+Irra+Iabx+Isabx+Irabx+Israbx+Issabx+Irrabx)/N, color="covid"), size=1) +
  theme_classic() + 
  scale_x_continuous(limits = c(995, 1365), breaks=seq(1000, 1365, 60), labels = values1) +
  scale_y_continuous(limits = c(0, 0.20), breaks=seq(0, 1, 0.02), labels = number) +
  labs(x="Time (days)", y="Proportion") +
  ggtitle("Prevalance of pneumococcal carriage") + 
  theme(legend.position="top", #legend.position="right",#legend.position = c(.8,.88) c(.8,.8)
        legend.background=element_blank(),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 14, face = "plain")) + 
  scale_color_manual(name = "Pneumococcal carriage:", values = colors, 
                     labels = c("S","R","Total","SR", "SS", "RR")) +
  geom_rect(data=lock, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            color="transparent", fill = "orange",alpha=0.3,inherit.aes = FALSE) 
##################################################################################################################




colors2 <- c("R" = "red", "covid" = "orange")
# Antibiotic resistance rate:
p2 <- ggplot(out15, aes(x=time)) + 
  geom_line(aes(y=R_15/PC_15, color="R"), size=1) +
  #geom_line(aes(y=(I+Is+Ir+Isr+Iss+Irr+Ia+Isa+Ira+Isra+Issa+Irra+Iabx+Isabx+Irabx+Israbx+Issabx+Irrabx)/N, color="covid"), size=1) +
  scale_x_continuous(limits = c(995, 1365), breaks=seq(1000, 1365, 60), labels = values1) +
  scale_y_continuous(limits = c(0, 0.5), breaks=seq(0, 1, 0.05), labels = number) +
  theme_classic() + 
  labs(x="Time (days)", y="Proportion") +
  ggtitle("Change in the prevalance of antibiotic-resistant carriage") + 
  theme(legend.position="none", #, #legend.position = c(.75,.65) c(.2,.5)
        legend.background=element_blank(),
        legend.title=element_text(size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 11),
        plot.title = element_text(hjust = 0.5, size = 14, face = "plain")) + 
  scale_color_manual(name = "Prevalance of antibiotic-resistant carriage", values = colors2, 
                     labels = c(expression(paste("R, RR and ", frac(1,2), "SR  ")))) +
  annotate(geom="text", x=1135, y=0.15, label=resistant_percent_change_lockdown_15, color="black", size = 8, fontface = 'bold') +
  annotate(geom="text", x=1185, y=0.15, label="%", color="black", size = 8, fontface = 'bold') +
  annotate(geom="text", x=1305, y=0.15, label=resistant_percent_change_15, color="black", size = 8, fontface = 'bold') +
  annotate(geom="text", x=1355, y=0.15, label="%", color="black", size = 8, fontface = 'bold') +
  geom_rect(data=lock, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            color="transparent", fill = "orange",alpha=0.3,inherit.aes = FALSE) + 
  annotate("segment", x = 1135, y = 0.05, xend = 1135, yend = 0,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("segment", x = 1365, y = 0.05, xend = 1365, yend = 0,
           arrow = arrow(type = "closed", length = unit(0.02, "npc"))) + 
  annotate(geom="text", x=1145, y=0.06, label="End of the lockdown") + 
  annotate(geom="text", x=1350, y=0.06, label="End of the year")




# Cumulative IPD incidence short and long term:------------------------------------
# Cumulative ipd incidence (S, R, SR):
# short term and long term
total_short_15 <- output_15$IPD[1135] - output_15$IPD[1075] # per 100,000 people cumulative ipd incidence during the 60-day period
total_long_15 <- output_15$IPD[1365] - output_15$IPD[1000] # per 100,000 people cumulative ipd incidence during one year since the outbreak

# IPD incidence caused by susceptible strain only (S + 0.5*SR):
# short term
S_short_15 <- output_15$IPDs[1135] - output_15$IPDs[1075] # per 100,000 people cumulative ipd incidence during the 60-day period
S_long_15 <- output_15$IPDs[1365] - output_15$IPDs[1000] # per 100,000 people cumulative ipd incidence during one year since the outbreak

# IPD incidence caused by resistant strain only (R + 0.5*SR):
# short term
R_short_15 <- output_15$IPDr[1135] - output_15$IPDr[1075] # per 100,000 people cumulative ipd incidence during the 60-day period
R_long_15 <- output_15$IPDr[1365] - output_15$IPDr[1000] # per 100,000 people cumulative ipd incidence during one year since the outbreak


# change in the prevalence of TOTAL bacterial carriage at the end of the 60-day lockdown period
relative_percent_change_15 <- ((out15$PC_15[1135] - out15$PC_15[1075])/out15$PC_15[1075])*100
print('% change in the prevalence of TOTAL pneumococcal carriage at the end of the 60-day period')
relative_percent_change_15

# Daily values:
output_15$daily_IPDr <- diff(c(0, output_15$IPDr))               # true daily incidence of Virus 1


values1 <- c ("0", "60", "120", "180", "240", "300", "360")
# baseline resistant IPD incidence over one year
p3 <- ggplot(data=output_15, aes(x=time, y=daily_IPDr)) + 
  geom_bar(stat = "identity", position = "dodge", width = 1, colour = "darkgrey") +
  #scale_x_continuous(limits = c(1000, 1365), breaks = seq(1000, 1365, 30)) +
  scale_y_continuous(limits = c(0, 0.025), breaks=seq(0, 0.07, 0.002), labels = number) +
  theme_classic() +
  labs(x="Time (days)", y="Daily cases of resistant IPDs") +
  ggtitle("Daily incidence of resistant IPDs over a one year period") +
  scale_x_continuous(limits = c(995, 1365), breaks=seq(1000, 1365, 60), labels = values1)


################################################################################
total_short_15                 # 60-day IPD incidence
total_long_15                  # Annual IPD incidence
R_long_15                      # Annual resistant IPD incidence
AR_a15 <- (R_long_15/total_long_15)*100  
AR_a15 # ANTIBIOTIC RESISTANCE (%)
################################################################################
output_15$ATB[1365]-output_15$ATB[1000] 



library(cowplot)
Fig_5 <- plot_grid(p1, p2, p3, 
                   #labels = c("a", "b", "c"), label_size = 12,
                   ncol = 3, nrow = 1)
Fig_5
#ggsave("Fig_1.jpeg", width = 55, height = 12, units = "cm")
