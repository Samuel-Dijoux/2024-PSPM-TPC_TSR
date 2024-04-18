####################################################################################################
#                                 PSPM_TPC-TSR_PSPMequi_simulations
#
# This script contains all the simulations using the PSPMequi functions we performed in our analysis
# The data of each simulations are stored in .... .RData
#
# Physiologically-structured population model describing a tri-trophic chain:
# Resource-Roach-Perch along gradients of temperature and resource productivity
#
# The model account for:
#   Temperature-size rule (TSR) & Temperature performance curve (TPC) in both consumer and predator species.
#   TSR implemented in Consumer size thresholds (maturation size 'Lj', asymptotic size 'Lm'),
#       in Predator foraging size range, through the upper size threshold of consumer exposition to predator 'Lv'
#
#   TPC implemented in consumer birth, growth, ingestion and mortality rates, 
#       in Predator functional response  and biomass loss rate
####################################################################################################
# library(PSPManalysis);
load("PSPM_TPC-TSR_PSPMequi_data.Rdata")

modelname = "PSPM_TPC-TSR_model.R" # Baseline model for scenarios 1-10 (Sections 1-3)
modelname_2 = "PSPM_TPC-TSR_Thermal-Mismatch_model.R" # Model to implement Species thermal mismatch in Scenarios 7-10 (Section 4)
########## Section 1: TPC only (S1) ====
# 9 scenario of TPC implementation:
# S1A) TPC in consumer ingestion rate;
# S1B) TPC in consumer birth rate;
# S1C) TPC in consumer growth rate;
# S1D) TPC in consumer mortality rate;
# S1E) TPC in all consumer rates *****
# S1F) TPC in predator functional response;
# S1G) TPC in predator mortality rate
# S1H) TPC in all predator rates *****
# S1I) TPC in all consumer & predator rates *****

## Community transition along productivity gradient ====
R_K <- PSPMequi(modelname, "EQ", c(1.0E-06, 1.0E-06), 0.8, c(1, 0, 5E-5), NULL,
                options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_K <- PSPMequi(modelname, "EQ", R_K$bifpoints[c(1, 2, 5)], 2, c(1, 0, 1E-3), NULL,
                 options = c("envZE", "1", "envZE", "2"));
PCR_K <- PSPMequi(modelname, "EQ", CR_K$bifpoints[,c(1, 2, 3, 7, 5)], -2, c(1, 0, 1E-4), NULL, NULL);

R_K$bifpoints; R_K$biftypes;
CR_K$bifpoints; CR_K$biftypes; 
PCR_K$bifpoints; PCR_K$biftypes;

## Data and points saved along K productivity gradient for each panel ====
# Dataset along productivity gradient and temperature = 20C
# S1A_K_R <- R_K$curvepoints; S1A_K_CR <- CR_K$curvepoints; S1A_K_PCR <- PCR_K$curvepoints; # for 15?C
# S1B_K_R <- R_K$curvepoints; S1B_K_CR <- CR_K$curvepoints; S1B_K_PCR <- PCR_K$curvepoints; # for 20?C
# S1C_K_R <- R_K$curvepoints; S1C_K_CR <- CR_K$curvepoints; S1C_K_PCR <- PCR_K$curvepoints; # for 20?C
# S1D_K_R <- R_K$curvepoints; S1D_K_CR <- CR_K$curvepoints; S1D_K_PCR <- PCR_K$curvepoints; # for 20?C
# S1E_K_R <- R_K$curvepoints; S1E_K_CR <- CR_K$curvepoints; S1E_K_PCR <- PCR_K$curvepoints; # for 20?C
# S1F_K_R <- R_K$curvepoints; S1F_K_CR <- CR_K$curvepoints; S1F_K_PCR <- PCR_K$curvepoints; # for 20?C
# S1G_K_R <- R_K$curvepoints; S1G_K_CR <- CR_K$curvepoints; S1G_K_PCR <- PCR_K$curvepoints; # for 20?C
# S1H_K_R <- R_K$curvepoints; S1H_K_CR <- CR_K$curvepoints; S1H_K_PCR <- PCR_K$curvepoints; # for 20?C
# S1I_K_R <- R_K$curvepoints; S1I_K_CR <- CR_K$curvepoints; S1I_K_PCR <- PCR_K$curvepoints; # for 20?C

# Bifurcation points
# S1A_K_BP <- R_K$bifpoints[1,]; S1A_K_BPE <- CR_K$bifpoints[1,]; S1A_K_LP <- PCR_K$bifpoints[1,];
# S1B_K_BP <- R_K$bifpoints[1,]; S1B_K_BPE <- CR_K$bifpoints[1,]; S1B_K_LP <- PCR_K$bifpoints[1,];
# S1C_K_BP <- R_K$bifpoints[1,]; S1C_K_BPE <- CR_K$bifpoints[1,]; S1C_K_LP <- PCR_K$bifpoints[1,];
# S1D_K_BP <- R_K$bifpoints[1,]; S1D_K_BPE <- CR_K$bifpoints[1,]; S1D_K_LP <- PCR_K$bifpoints[1,];
# S1E_K_BP <- R_K$bifpoints[1,]; S1E_K_BPE <- CR_K$bifpoints[1,]; S1E_K_LP <- PCR_K$bifpoints[1,];
# S1F_K_BP <- R_K$bifpoints[1,]; S1F_K_BPE <- CR_K$bifpoints[1,]; S1F_K_LP <- PCR_K$bifpoints[1,];
# S1G_K_BP <- R_K$bifpoints[1,]; S1G_K_BPE <- CR_K$bifpoints[1,]; S1G_K_LP <- PCR_K$bifpoints[1,];
# S1H_K_BP <- R_K$bifpoints[1,]; S1H_K_BPE <- CR_K$bifpoints[1,]; S1H_K_LP <- PCR_K$bifpoints[1,];
# S1I_K_BP <- R_K$bifpoints[1,]; S1I_K_BPE <- CR_K$bifpoints[1,]; S1I_K_LP <- PCR_K$bifpoints[1,];

## Community transition along productivity & Temperature gradients ====
# S1A: TPC in Consumer ingestion rate ====
S1A_KT_BPa <- PSPMequi(modelname, "BP", c(S1A_K_BP[c(1, 2, 5)], 15), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S1A_KT_BPb <- PSPMequi(modelname, "BP", c(S1A_K_BP[c(1, 2, 5)], 15), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S1A_KT_BP<-rbind(S1A_KT_BPa$curvepoints[dim(S1A_KT_BPa$curvepoints)[1]:1,], S1A_KT_BPb$curvepoints[2:dim(S1A_KT_BPb$curvepoints)[1],]);

S1A_KT_BPEa <- PSPMequi(modelname, "BPE", c(S1A_K_BPE[c(1,2,4,5)], 15), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1A_KT_BPEb <- PSPMequi(modelname, "BPE", c(S1A_K_BPE[c(1,2,4,5)], 15), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1A_KT_BPE<-rbind(S1A_KT_BPEa$curvepoints[dim(S1A_KT_BPEa$curvepoints)[1]:1,], S1A_KT_BPEb$curvepoints[2:dim(S1A_KT_BPEb$curvepoints)[1],]);

S1A_KT_LPa <- PSPMequi(modelname, "LP", c(S1A_K_LP[c(1:5)], 15), -0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
#S1A_KT_LPb <- PSPMequi(modelname, "LP", c(S1A_KT_LPa$curvepoints[3,c(1:6)]), 0.8, c(1, 1E-6, 1E-3, 5.01601260, 1E-2, 30), NULL, options = NULL);
#S1A_KT_LP<-rbind( S1A_KT_LPa$curvepoints[dim(S1A_KT_LPa$curvepoints)[1]:1,], S1A_KT_LPb$curvepoints[2:dim(S1A_KT_LPb$curvepoints)[1],]);
## LP points to detect manually
S1A_KT_LP16 <- PCR_K$bifpoints[1,];# for 16°C
S1A_KT_LP17 <- PCR_K$bifpoints[1,];# for 17°C
S1A_KT_LP18 <- PCR_K$bifpoints[1,];# for 18°C
S1A_KT_LP19 <- PCR_K$bifpoints[1,];# for 19°C
S1A_KT_LP20 <- PCR_K$bifpoints[1,];# for 20°C
S1A_KT_LP21 <- PCR_K$bifpoints[1,];# for 21°C
S1A_KT_LP22 <- PCR_K$bifpoints[1,];# for 22°C
S1A_KT_LP23 <- PCR_K$bifpoints[1,];# for 23°C
S1A_KT_LP24 <- PCR_K$bifpoints[1,];# for 24°C
S1A_KT_LP245 <- PCR_K$bifpoints[1,];# for 24.5°C
S1A_KT_LP247 <- PCR_K$bifpoints[1,];# for 24.7°C
S1A_KT_LP249 <- PCR_K$bifpoints[1,];# for 24.9°C

S1A_KT_LP_DF <- rbind(as.numeric(S1A_KT_LP16), as.numeric(S1A_KT_LP17), as.numeric(S1A_KT_LP18), as.numeric(S1A_KT_LP19),
                                 as.numeric(S1A_KT_LP20), as.numeric(S1A_KT_LP21), as.numeric(S1A_KT_LP22), as.numeric(S1A_KT_LP23),
                                 as.numeric(S1A_KT_LP24), as.numeric(S1A_KT_LP245), as.numeric(S1A_KT_LP247), as.numeric(S1A_KT_LP249));
S1A_KT_LP_DFb <- cbind(S1A_KT_LP_DF[,1:5], as.vector(c(seq(16,24,by=1), 24.5, 24.7, 24.9)), S1A_KT_LP_DF[,6:12]);
S1A_KT_LP<-rbind( S1A_KT_LPa$curvepoints[dim(S1A_KT_LPa$curvepoints)[1]:1,], S1A_KT_LP_DFb);

#graph test
BPcurve <- S1A_KT_BP;
BPEcurve <- S1A_KT_BPE;
LPcurve <- S1A_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S1B: TPC in Consumer birth rate ====
S1B_KT_BPa <- PSPMequi(modelname, "BP", c(S1B_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S1B_KT_BPb <- PSPMequi(modelname, "BP", c(S1B_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S1B_KT_BP<-rbind(S1B_KT_BPa$curvepoints[dim(S1B_KT_BPa$curvepoints)[1]:1,], S1B_KT_BPb$curvepoints[2:dim(S1B_KT_BPb$curvepoints)[1],]);

S1B_KT_BPEa <- PSPMequi(modelname, "BPE", c(S1B_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1B_KT_BPEb <- PSPMequi(modelname, "BPE", c(S1B_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1B_KT_BPE<-rbind(S1B_KT_BPEa$curvepoints[dim(S1B_KT_BPEa$curvepoints)[1]:1,], S1B_KT_BPEb$curvepoints[2:dim(S1B_KT_BPEb$curvepoints)[1],]);

S1B_KT_LPa <- PSPMequi(modelname, "LP", c(S1B_K_LP[c(1:5)], 20), -0.1, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S1B_KT_LPb <- PSPMequi(modelname, "LP", c(S1B_K_LP[c(1:5)], 20), 0.2, c(1, 1E-6, 1E-3, 0, 1E-2, 24.995), NULL, options = NULL);
S1B_KT_LP<-rbind(S1B_KT_LPa$curvepoints[dim(S1B_KT_LPa$curvepoints)[1]:1,], S1B_KT_LPb$curvepoints[2:dim(S1B_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve <- S1B_KT_BP;
BPEcurve <- S1B_KT_BPE;
LPcurve <- S1B_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S1C: TPC in Consumer growth rate ====
S1C_KT_BPa <- PSPMequi(modelname, "BP", c(S1C_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S1C_KT_BPb <- PSPMequi(modelname, "BP", c(S1C_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S1C_KT_BP <-rbind(S1C_KT_BPa$curvepoints[dim(S1C_KT_BPa$curvepoints)[1]:1,], S1C_KT_BPb$curvepoints[2:dim(S1C_KT_BPb$curvepoints)[1],]);

S1C_KT_BPEa <- PSPMequi(modelname, "BPE", c(S1C_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1C_KT_BPEb <- PSPMequi(modelname, "BPE", c(S1C_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1C_KT_BPE <-rbind(S1C_KT_BPEa$curvepoints[dim(S1C_KT_BPEa$curvepoints)[1]:1,], S1C_KT_BPEb$curvepoints[2:dim(S1C_KT_BPEb$curvepoints)[1],]);

S1C_KT_LPa <- PSPMequi(modelname, "LP", c(S1C_K_LP[c(1:5)], 20), -0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S1C_KT_LPb <- PSPMequi(modelname, "LP", c(S1C_K_LP[c(1:5)], 20), 0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S1C_KT_LPc <- PSPMequi(modelname, "LP", S1C_KT_LPb$curvepoints[dim(S1C_KT_LPb$curvepoints)[1],1:6], 0.8, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S1C_KT_LP <-rbind(S1C_KT_LPa$curvepoints[dim(S1C_KT_LPa$curvepoints)[1]:1,], S1C_KT_LPb$curvepoints[2:dim(S1C_KT_LPb$curvepoints)[1],], S1C_KT_LPc$curvepoints[2:dim(S1C_KT_LPc$curvepoints)[1],]);

#graph test
BPcurve <- S1C_KT_BP;
BPEcurve <- S1C_KT_BPE;
LPcurve <- S1C_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S1D) TPC in Consumer mortality rate ====
S1D_KT_BPa <- PSPMequi(modelname, "BP", c(S1D_K_BP[c(1, 2, 5)], 20), -0.8, c(1, 0, 6E-3, 10, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S1D_KT_BPb <- PSPMequi(modelname, "BP", c(S1D_K_BP[c(1, 2, 5)], 20), 0.8, c(1, 0, 6E-3, 10, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S1D_KT_BP<-rbind(S1D_KT_BPa$curvepoints[dim(S1D_KT_BPa$curvepoints)[1]:1,], S1D_KT_BPb$curvepoints[2:dim(S1D_KT_BPb$curvepoints)[1],]);

S1D_KT_BPEa <- PSPMequi(modelname, "BPE", c(S1D_K_BPE[c(1,2,4,5)], 20), -0.8, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1D_KT_BPEb <- PSPMequi(modelname, "BPE", c(S1D_K_BPE[c(1,2,4,5)], 20), 0.8, c(1, 1E-6, 6E-3, 0.1, 1E-2, 30), NULL, options = c("envBP", "1"));
S1D_KT_BPE<-rbind(S1D_KT_BPEa$curvepoints[dim(S1D_KT_BPEa$curvepoints)[1]:1,], S1D_KT_BPEb$curvepoints[2:dim(S1D_KT_BPEb$curvepoints)[1],]);

S1D_KT_LPa <- PSPMequi(modelname, "LP", c(S1D_K_LP[c(1:5)], 20), -0.8, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S1D_KT_LPb <- PSPMequi(modelname, "LP", c(S1D_K_LP[c(1:5)], 20), 0.8, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S1D_KT_LP <- rbind(S1D_KT_LPa$curvepoints[dim(S1D_KT_LPa$curvepoints)[1]:1,], S1D_KT_LPb$curvepoints[2:dim(S1D_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve <- S1D_KT_BP;
BPEcurve <- S1D_KT_BPE;
LPcurve <- S1D_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);
abline(v=max(BPEcurve[,6]), lty=4)

# S1E) TPC in all Consumer rates ***** ==== 
S1E_KT_BPa <- PSPMequi(modelname, "BP", c(S1E_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S1E_KT_BPb <- PSPMequi(modelname, "BP", c(S1E_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S1E_KT_BP<-rbind(S1E_KT_BPa$curvepoints[dim(S1E_KT_BPa$curvepoints)[1]:1,], S1E_KT_BPb$curvepoints[2:dim(S1E_KT_BPb$curvepoints)[1],]);

S1E_KT_BPEa <- PSPMequi(modelname, "BPE", c(S1E_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1E_KT_BPEb <- PSPMequi(modelname, "BPE", c(S1E_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1E_KT_BPE<-rbind(S1E_KT_BPEa$curvepoints[dim(S1E_KT_BPEa$curvepoints)[1]:1,], S1E_KT_BPEb$curvepoints[2:dim(S1E_KT_BPEb$curvepoints)[1],]);

S1E_KT_LPa <- PSPMequi(modelname, "LP", c(S1E_K_LP[c(1:5)], 20), -0.2, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S1E_KT_LPb <- PSPMequi(modelname, "LP", c(S1E_K_LP[c(1:5)], 20), 0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S1E_KT_LP<-rbind(S1E_KT_LPa$curvepoints[dim(S1E_KT_LPa$curvepoints)[1]:1,], S1E_KT_LPb$curvepoints[2:dim(S1E_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve <- S1E_KT_BP;
BPEcurve <- S1E_KT_BPE;
LPcurve <- S1E_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S1F: TPC in Predator functional response ====
# no influence on consumer invasion threshold
S1F_KT_BP <- cbind.data.frame(rep(S1F_K_BP[1], 7), rep(S1F_K_BP[2], 7), rep(S1F_K_BP[5], 7), seq(0, 30, by=5));
colnames(S1F_KT_BP) <- c("Rmax", "E[0]", "b[0]", "Temperature")

S1F_KT_BPEa <- PSPMequi(modelname, "BPE", c(S1F_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1F_KT_BPEb <- PSPMequi(modelname, "BPE", c(S1F_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1F_KT_BPE<-rbind(S1F_KT_BPEa$curvepoints[dim(S1F_KT_BPEa$curvepoints)[1]:1,], S1F_KT_BPEb$curvepoints[2:dim(S1F_KT_BPEb$curvepoints)[1],]);

S1F_KT_LPa <- PSPMequi(modelname, "LP", c(S1F_K_LP[c(1:5)], 20), -0.2, c(1, 1E-6, 1E-2, 0, 1E-2, 30), NULL, options = NULL);
S1F_KT_LPb <- PSPMequi(modelname, "LP", c(S1F_K_LP[c(1:5)], 20), 0.2, c(1, 1E-6, 1E-2, 0, 1E-2, 30), NULL, options = NULL);
S1F_KT_LP<-rbind(S1F_KT_LPa$curvepoints[dim(S1F_KT_LPa$curvepoints)[1]:1,], S1F_KT_LPb$curvepoints[2:dim(S1F_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve <- S1F_KT_BP;
BPEcurve <- S1F_KT_BPE;
LPcurve <- S1F_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,4], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);
# S1G: TPC in Predator Biomass loss rate ====
# no influence on consumer invasion threshold
S1G_KT_BP <- cbind.data.frame(rep(S1G_K_BP[1], 7), rep(S1G_K_BP[2], 7), rep(S1G_K_BP[5], 7), seq(0, 30, by=5));
colnames(S1G_KT_BP) <- c("Rmax", "E[0]", "b[0]", "Temperature")

S1G_KT_BPEa <- PSPMequi(modelname, "BPE", c(S1G_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1G_KT_BPEb <- PSPMequi(modelname, "BPE", c(S1G_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1G_KT_BPE<-rbind(S1G_KT_BPEa$curvepoints[dim(S1G_KT_BPEa$curvepoints)[1]:1,], S1G_KT_BPEb$curvepoints[2:dim(S1G_KT_BPEb$curvepoints)[1],]);

S1G_KT_LPa <- PSPMequi(modelname, "LP", c(S1G_K_LP[c(1:5)], 20), -0.2, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S1G_KT_LPb <- PSPMequi(modelname, "LP", c(S1G_K_LP[c(1:5)], 20), 0.2, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S1G_KT_LP<-rbind(S1G_KT_LPa$curvepoints[dim(S1G_KT_LPa$curvepoints)[1]:1,], S1G_KT_LPb$curvepoints[2:dim(S1G_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve <- S1G_KT_BP;
BPEcurve <- S1G_KT_BPE;
LPcurve <- S1G_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,4], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);
# S1H) TPC in all Predator rates ***** ====
# no influence on consumer invasion threshold
S1H_KT_BP <- cbind.data.frame(rep(S1H_K_BP[1], 7), rep(S1H_K_BP[2], 7), rep(S1H_K_BP[5], 7), seq(0, 30, by=5));
colnames(S1H_KT_BP) <- c("Rmax", "E[0]", "b[0]", "Temperature")

S1H_KT_BPEa <- PSPMequi(modelname, "BPE", c(S1H_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1H_KT_BPEb <- PSPMequi(modelname, "BPE", c(S1H_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1H_KT_BPE<-rbind(S1H_KT_BPEa$curvepoints[dim(S1H_KT_BPEa$curvepoints)[1]:1,], S1H_KT_BPEb$curvepoints[2:dim(S1H_KT_BPEb$curvepoints)[1],]);

S1H_KT_LPa <- PSPMequi(modelname, "LP", c(S1H_K_LP[c(1:5)], 20), -0.2, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S1H_KT_LPb <- PSPMequi(modelname, "LP", c(S1H_K_LP[c(1:5)], 20), 0.2, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S1H_KT_LP<-rbind(S1H_KT_LPa$curvepoints[dim(S1H_KT_LPa$curvepoints)[1]:1,], S1H_KT_LPb$curvepoints[2:dim(S1H_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve <- S1H_KT_BP;
BPEcurve <- S1H_KT_BPE;
LPcurve <- S1H_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,4], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S1I) TPC in all Consumer & Predator rates ***** ====
S1I_KT_BPa <- PSPMequi(modelname, "BP", c(S1I_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S1I_KT_BPb <- PSPMequi(modelname, "BP", c(S1I_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S1I_KT_BP<-rbind(S1I_KT_BPa$curvepoints[dim(S1I_KT_BPa$curvepoints)[1]:1,], S1I_KT_BPb$curvepoints[2:dim(S1I_KT_BPb$curvepoints)[1],]);

S1I_KT_BPEa <- PSPMequi(modelname, "BPE", c(S1I_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1I_KT_BPEb <- PSPMequi(modelname, "BPE", c(S1I_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S1I_KT_BPE<-rbind(S1I_KT_BPEa$curvepoints[dim(S1I_KT_BPEa$curvepoints)[1]:1,], S1I_KT_BPEb$curvepoints[2:dim(S1I_KT_BPEb$curvepoints)[1],]);

S1I_KT_LPa <- PSPMequi(modelname, "LP", c(S1I_K_LP[c(1:5)], 20), -0.2, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S1I_KT_LPb <- PSPMequi(modelname, "LP", c(S1I_K_LP[c(1:5)], 20), 0.2, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S1I_KT_LP<-rbind(S1I_KT_LPa$curvepoints[dim(S1I_KT_LPa$curvepoints)[1]:1,], S1I_KT_LPb$curvepoints[2:dim(S1I_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve <- S1I_KT_BP;
BPEcurve <- S1I_KT_BPE;
LPcurve <- S1I_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

############ Section 2: TSR only (S2) ====
# 7 scenario of TSR implementation:
# S2A) TSR in Predator max foraging size Lv;  ******
# S2B) TSR in Consumer maturation size Lj;
# S2C) TSR in Consumer asymptotic size Lm;
# S2D) TSR in Lj & Lm; ******
# S1E) TSR in Lv & Lj;
# S2F) TSR in Lv & Lm
# S2G) TSR in Lv, Lj & Lm.  ******
## Community transition along productivity gradient ====
R_K <- PSPMequi(modelname, "EQ", c(1.0E-06, 1.0E-06), 0.1, c(1, 0, 5E-5), NULL,
                options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_K <- PSPMequi(modelname, "EQ", R_K$bifpoints[c(1, 2, 5)], 2, c(1, 0, 5E-4), NULL,
                 options = c("envZE", "1", "envZE", "2"));
PCR_K <- PSPMequi(modelname, "EQ", CR_K$bifpoints[,c(1, 2, 3, 7, 5)], -2, c(1, 0, 1E-3), NULL, NULL);

R_K$bifpoints; R_K$biftypes;
CR_K$bifpoints; CR_K$biftypes; 
PCR_K$bifpoints; PCR_K$biftypes;

## Data and points saved along K productivity gradient for each panel ====
# Dataset along productivity gradient and temperature = 20C
# S2A_K_R <- R_K$curvepoints; S2A_K_CR <- CR_K$curvepoints; S2A_K_PCR <- PCR_K$curvepoints;
# S2B_K_R <- R_K$curvepoints; S2B_K_CR <- CR_K$curvepoints; S2B_K_PCR <- PCR_K$curvepoints;
# S2C_K_R <- R_K$curvepoints; S2C_K_CR <- CR_K$curvepoints; S2C_K_PCR <- PCR_K$curvepoints;
# S2D_K_R <- R_K$curvepoints; S2D_K_CR <- CR_K$curvepoints; S2D_K_PCR <- PCR_K$curvepoints;
# S2E_K_R <- R_K$curvepoints; S2E_K_CR <- CR_K$curvepoints; S2E_K_PCR <- PCR_K$curvepoints;
# S2F_K_R <- R_K$curvepoints; S2F_K_CR <- CR_K$curvepoints; S2F_K_PCR <- PCR_K$curvepoints;
# S2G_K_R <- R_K$curvepoints; S2G_K_CR <- CR_K$curvepoints; S2G_K_PCR <- PCR_K$curvepoints;

# Bifurcation points
# S2A_K_BP <- R_K$bifpoints[1,]; S2A_K_BPE <- CR_K$bifpoints[1,]; S2A_K_LP <- PCR_K$bifpoints[1,];
# S2B_K_BP <- R_K$bifpoints[1,]; S2B_K_BPE <- CR_K$bifpoints[1,]; S2B_K_LP <- PCR_K$bifpoints[1,];
# S2C_K_BP <- R_K$bifpoints[1,]; S2C_K_BPE <- CR_K$bifpoints[1,]; S2C_K_LP <- PCR_K$bifpoints[1,];
# S2D_K_BP <- R_K$bifpoints[1,]; S2D_K_BPE <- CR_K$bifpoints[1,]; S2D_K_LP <- PCR_K$bifpoints[1,];
# S2E_K_BP <- R_K$bifpoints[1,]; S2E_K_BPE <- CR_K$bifpoints[1,]; S2E_K_LP <- PCR_K$bifpoints[1,];
# S2F_K_BP <- R_K$bifpoints[1,]; S2F_K_BPE <- CR_K$bifpoints[1,]; S2F_K_LP <- PCR_K$bifpoints[1,];
# S2G_K_BP <- R_K$bifpoints[1,]; S2G_K_BPE <- CR_K$bifpoints[1,]; S2G_K_LP <- PCR_K$bifpoints[1,];

## Community transition along productivity & Temperature gradients ====
# S2A: TSR in Predator mxm foraging size threshold Lv ***** ====
S2A_KT_BPa <- PSPMequi(modelname, "BP", c(S2A_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2A_KT_BPb <- PSPMequi(modelname, "BP", c(S2A_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2A_KT_BP<-rbind(S2A_KT_BPa$curvepoints[dim(S2A_KT_BPa$curvepoints)[1]:1,], S2A_KT_BPb$curvepoints[2:dim(S2A_KT_BPb$curvepoints)[1],]);

S2A_KT_BPEa <- PSPMequi(modelname, "BPE", c(S2A_K_BPE[c(1,2,4,5)], 15), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2A_KT_BPEb <- PSPMequi(modelname, "BPE", c(S2A_K_BPE[c(1,2,4,5)], 15), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2A_KT_BPE<-rbind(S2A_KT_BPEa$curvepoints[dim(S2A_KT_BPEa$curvepoints)[1]:1,], S2A_KT_BPEb$curvepoints[2:dim(S2A_KT_BPEb$curvepoints)[1],]);

S2A_KT_LPa <- PSPMequi(modelname, "LP", c(S2A_K_LP[c(1:5)], 20), -0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S2A_KT_LPb <- PSPMequi(modelname, "LP", c(S2A_K_LP[c(1:5)], 20), 0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S2A_KT_LP<-rbind(S2A_KT_LPa$curvepoints[dim(S2A_KT_LPa$curvepoints)[1]:1,], S2A_KT_LPb$curvepoints[2:dim(S2A_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S2A_KT_BP;
BPEcurve <- S2A_KT_BPE;
LPcurve  <- S2A_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S2B: TSR in Consumer maturation size Lj ====
S2B_KT_BPa <- PSPMequi(modelname, "BP", c(S2B_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2B_KT_BPb <- PSPMequi(modelname, "BP", c(S2B_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2B_KT_BP<-rbind(S2B_KT_BPa$curvepoints[dim(S2B_KT_BPa$curvepoints)[1]:1,], S2B_KT_BPb$curvepoints[2:dim(S2B_KT_BPb$curvepoints)[1],]);

S2B_KT_BPEa <- PSPMequi(modelname, "BPE", c(S2B_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2B_KT_BPEb <- PSPMequi(modelname, "BPE", c(S2B_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2B_KT_BPE<-rbind(S2B_KT_BPEa$curvepoints[dim(S2B_KT_BPEa$curvepoints)[1]:1,], S2B_KT_BPEb$curvepoints[2:dim(S2B_KT_BPEb$curvepoints)[1],]);

S2B_KT_LPa <- PSPMequi(modelname, "LP", c(S2B_K_LP[c(1:5)], 20), -0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S2B_KT_LPb <- PSPMequi(modelname, "LP", c(S2B_K_LP[c(1:5)], 20), 0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S2B_KT_LP<-rbind(S2B_KT_LPa$curvepoints[dim(S2B_KT_LPa$curvepoints)[1]:1,], S2B_KT_LPb$curvepoints[2:dim(S2B_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S2B_KT_BP;
BPEcurve <- S2B_KT_BPE;
LPcurve  <- S2B_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S2C: TSR in Consumer asymptotic size Lm ====
S2C_KT_BPa <- PSPMequi(modelname, "BP", c(S2C_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2C_KT_BPb <- PSPMequi(modelname, "BP", c(S2C_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2C_KT_BP<-rbind(S2C_KT_BPa$curvepoints[dim(S2C_KT_BPa$curvepoints)[1]:1,], S2C_KT_BPb$curvepoints[2:dim(S2C_KT_BPb$curvepoints)[1],]);

S2C_KT_BPEa <- PSPMequi(modelname, "BPE", c(S2C_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2C_KT_BPEb <- PSPMequi(modelname, "BPE", c(S2C_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2C_KT_BPE<-rbind(S2C_KT_BPEa$curvepoints[dim(S2C_KT_BPEa$curvepoints)[1]:1,], S2C_KT_BPEb$curvepoints[2:dim(S2C_KT_BPEb$curvepoints)[1],]);

S2C_KT_LPa <- PSPMequi(modelname, "LP", c(S2C_K_LP[c(1:5)], 20), -0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S2C_KT_LPb <- PSPMequi(modelname, "LP", c(S2C_K_LP[c(1:5)], 20), 0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S2C_KT_LP<-rbind(S2C_KT_LPa$curvepoints[dim(S2C_KT_LPa$curvepoints)[1]:1,], S2C_KT_LPb$curvepoints[2:dim(S2C_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S2C_KT_BP;
BPEcurve <- S2C_KT_BPE;
LPcurve  <- S2C_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S2D: TSR in lj & Lm ***** ====
S2D_KT_BPa <- PSPMequi(modelname, "BP", c(S2D_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2D_KT_BPb <- PSPMequi(modelname, "BP", c(S2D_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2D_KT_BP<-rbind(S2D_KT_BPa$curvepoints[dim(S2D_KT_BPa$curvepoints)[1]:1,], S2D_KT_BPb$curvepoints[2:dim(S2D_KT_BPb$curvepoints)[1],]);

S2D_KT_BPEa <- PSPMequi(modelname, "BPE", c(S2D_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2D_KT_BPEb <- PSPMequi(modelname, "BPE", c(S2D_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2D_KT_BPE<-rbind(S2D_KT_BPEa$curvepoints[dim(S2D_KT_BPEa$curvepoints)[1]:1,], S2D_KT_BPEb$curvepoints[2:dim(S2D_KT_BPEb$curvepoints)[1],]);

S2D_KT_LPa <- PSPMequi(modelname, "LP", c(S2D_K_LP[c(1:5)], 20), -0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S2D_KT_LPb <- PSPMequi(modelname, "LP", c(S2D_K_LP[c(1:5)], 20), 0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S2D_KT_LP<-rbind(S2D_KT_LPa$curvepoints[dim(S2D_KT_LPa$curvepoints)[1]:1,], S2D_KT_LPb$curvepoints[2:dim(S2D_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S2D_KT_BP;
BPEcurve <- S2D_KT_BPE;
LPcurve  <- S2D_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S2E: TSR in Lv & Lj ====
S2E_KT_BPa <- PSPMequi(modelname, "BP", c(S2E_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2E_KT_BPb <- PSPMequi(modelname, "BP", c(S2E_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2E_KT_BP<-rbind(S2E_KT_BPa$curvepoints[dim(S2E_KT_BPa$curvepoints)[1]:1,], S2E_KT_BPb$curvepoints[2:dim(S2E_KT_BPb$curvepoints)[1],]);

S2E_KT_BPEa <- PSPMequi(modelname, "BPE", c(S2E_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2E_KT_BPEb <- PSPMequi(modelname, "BPE", c(S2E_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2E_KT_BPE<-rbind(S2E_KT_BPEa$curvepoints[dim(S2E_KT_BPEa$curvepoints)[1]:1,], S2E_KT_BPEb$curvepoints[2:dim(S2E_KT_BPEb$curvepoints)[1],]);

S2E_KT_LPa <- PSPMequi(modelname, "LP", c(S2E_K_LP[c(1:5)], 20), -0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S2E_KT_LPb <- PSPMequi(modelname, "LP", c(S2E_K_LP[c(1:5)], 20), 0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S2E_KT_LP<-rbind(S2E_KT_LPa$curvepoints[dim(S2E_KT_LPa$curvepoints)[1]:1,], S2E_KT_LPb$curvepoints[2:dim(S2E_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S2E_KT_BP;
BPEcurve <- S2E_KT_BPE;
LPcurve  <- S2E_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S2F: TSR in Lv & Lm ====
S2F_KT_BPa <- PSPMequi(modelname, "BP", c(S2F_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2F_KT_BPb <- PSPMequi(modelname, "BP", c(S2F_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2F_KT_BP<-rbind(S2F_KT_BPa$curvepoints[dim(S2F_KT_BPa$curvepoints)[1]:1,], S2F_KT_BPb$curvepoints[2:dim(S2F_KT_BPb$curvepoints)[1],]);

S2F_KT_BPEa <- PSPMequi(modelname, "BPE", c(S2F_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2F_KT_BPEb <- PSPMequi(modelname, "BPE", c(S2F_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2F_KT_BPE<-rbind(S2F_KT_BPEa$curvepoints[dim(S2F_KT_BPEa$curvepoints)[1]:1,], S2F_KT_BPEb$curvepoints[2:dim(S2F_KT_BPEb$curvepoints)[1],]);

S2F_KT_LPa <- PSPMequi(modelname, "LP", c(S2F_K_LP[c(1:5)], 20), -0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S2F_KT_LPb <- PSPMequi(modelname, "LP", c(S2F_K_LP[c(1:5)], 20), 0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S2F_KT_LP<-rbind(S2F_KT_LPa$curvepoints[dim(S2F_KT_LPa$curvepoints)[1]:1,], S2F_KT_LPb$curvepoints[2:dim(S2F_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S2F_KT_BP;
BPEcurve <- S2F_KT_BPE;
LPcurve  <- S2F_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S2G: TSR in Lv, Lj & Lm ***** ====
S2G_KT_BPa <- PSPMequi(modelname, "BP", c(S2G_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2G_KT_BPb <- PSPMequi(modelname, "BP", c(S2G_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S2G_KT_BP<-rbind(S2G_KT_BPa$curvepoints[dim(S2G_KT_BPa$curvepoints)[1]:1,], S2G_KT_BPb$curvepoints[2:dim(S2G_KT_BPb$curvepoints)[1],]);

S2G_KT_BPEa <- PSPMequi(modelname, "BPE", c(S2G_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2G_KT_BPEb <- PSPMequi(modelname, "BPE", c(S2G_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S2G_KT_BPE<-rbind(S2G_KT_BPEa$curvepoints[dim(S2G_KT_BPEa$curvepoints)[1]:1,], S2G_KT_BPEb$curvepoints[2:dim(S2G_KT_BPEb$curvepoints)[1],]);

S2G_KT_LPa <- PSPMequi(modelname, "LP", c(S2G_K_LP[c(1:5)], 20), -0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S2G_KT_LPb <- PSPMequi(modelname, "LP", c(S2G_K_LP[c(1:5)], 20), 0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S2G_KT_LP<-rbind(S2G_KT_LPa$curvepoints[dim(S2G_KT_LPa$curvepoints)[1]:1,], S2G_KT_LPb$curvepoints[2:dim(S2G_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S2G_KT_BP;
BPEcurve <- S2G_KT_BPE;
LPcurve  <- S2G_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

############ Section 3: TPC+TSR only (S3) ====
# 3*3 combinations of TPC x TSR implemented in consumer only, predator only, and both consumer & predator:
# S3A-C) TPC(C)+ A) TSR(C), B) TSR(P), C) TSR(C+P)
# S3D-F) TPC(P)+ A) TSR(C), B) TSR(P), C) TSR(C+P)
# S3G-I) TPC(C+P)+ A) TSR(C), B) TSR(P), C) TSR(C+P)
## Community transition along productivity gradient ====
R_K <- PSPMequi(modelname, "EQ", c(1.0E-06, 1.0E-06), 0.1, c(1, 0, 5E-5), NULL,
                options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_K <- PSPMequi(modelname, "EQ", R_K$bifpoints[c(1, 2, 5)], 2, c(1, 0, 5E-4), NULL,
                 options = c("envZE", "1", "envZE", "2"));
PCR_K <- PSPMequi(modelname, "EQ", CR_K$bifpoints[,c(1, 2, 3, 7, 5)], -2, c(1, 0, 1E-3), NULL, NULL);

R_K$bifpoints; R_K$biftypes;
CR_K$bifpoints; CR_K$biftypes; 
PCR_K$bifpoints; PCR_K$biftypes;
## Data and points saved along K productivity gradient for each panel ====
# Dataset along productivity gradient and temperature = 20C
# S3A_K_R <- R_K$curvepoints; S3A_K_CR <- CR_K$curvepoints; S3A_K_PCR <- PCR_K$curvepoints;
# S3B_K_R <- R_K$curvepoints; S3B_K_CR <- CR_K$curvepoints; S3B_K_PCR <- PCR_K$curvepoints;
# S3C_K_R <- R_K$curvepoints; S3C_K_CR <- CR_K$curvepoints; S3C_K_PCR <- PCR_K$curvepoints;
# S3D_K_R <- R_K$curvepoints; S3D_K_CR <- CR_K$curvepoints; S3D_K_PCR <- PCR_K$curvepoints;
# S3E_K_R <- R_K$curvepoints; S3E_K_CR <- CR_K$curvepoints; S3E_K_PCR <- PCR_K$curvepoints;
# S3F_K_R <- R_K$curvepoints; S3F_K_CR <- CR_K$curvepoints; S3F_K_PCR <- PCR_K$curvepoints;
# S3G_K_R <- R_K$curvepoints; S3G_K_CR <- CR_K$curvepoints; S3G_K_PCR <- PCR_K$curvepoints;
# S3H_K_R <- R_K$curvepoints; S3H_K_CR <- CR_K$curvepoints; S3H_K_PCR <- PCR_K$curvepoints;
# S3I_K_R <- R_K$curvepoints; S3I_K_CR <- CR_K$curvepoints; S3I_K_PCR <- PCR_K$curvepoints;

# S3I_K_R_20 <- R_K$curvepoints; S3I_K_CR_20 <- CR_K$curvepoints; S3I_K_PCR_20 <- PCR_K$curvepoints; #T=20
# S3I_K_R_13 <- R_K$curvepoints; S3I_K_CR_13 <- CR_K$curvepoints; S3I_K_PCR_13 <- PCR_K$curvepoints; #T=13

# Bifurcation points
# S3A_K_BP <- R_K$bifpoints[1,]; S3A_K_BPE <- CR_K$bifpoints[1,]; S3A_K_LP <- PCR_K$bifpoints[1,];
# S3B_K_BP <- R_K$bifpoints[1,]; S3B_K_BPE <- CR_K$bifpoints[1,]; S3B_K_LP <- PCR_K$bifpoints[1,];
# S3C_K_BP <- R_K$bifpoints[1,]; S3C_K_BPE <- CR_K$bifpoints[1,]; S3C_K_LP <- PCR_K$bifpoints[1,];
# S3D_K_BP <- R_K$bifpoints[1,]; S3D_K_BPE <- CR_K$bifpoints[1,]; S3D_K_LP <- PCR_K$bifpoints[1,];
# S3E_K_BP <- R_K$bifpoints[1,]; S3E_K_BPE <- CR_K$bifpoints[1,]; S3E_K_LP <- PCR_K$bifpoints[1,];
# S3F_K_BP <- R_K$bifpoints[1,]; S3F_K_BPE <- CR_K$bifpoints[1,]; S3F_K_LP <- PCR_K$bifpoints[1,];
# S3G_K_BP <- R_K$bifpoints[1,]; S3G_K_BPE <- CR_K$bifpoints[1,]; S3G_K_LP <- PCR_K$bifpoints[1,];
# S3H_K_BP <- R_K$bifpoints[1,]; S3H_K_BPE <- CR_K$bifpoints[1,]; S3H_K_LP <- PCR_K$bifpoints[1,];
# S3I_K_BP <- R_K$bifpoints[1,]; S3I_K_BPE <- CR_K$bifpoints[1,]; S3I_K_LP <- PCR_K$bifpoints[1,];

# S3I_K_BP_20 <- R_K$bifpoints[1,]; S3I_K_BPE_20 <- CR_K$bifpoints[1,]; S3I_K_LP_20 <- PCR_K$bifpoints[1,]; #T=20
# S3I_K_BP_13 <- R_K$bifpoints[1,]; S3I_K_BPE_13 <- CR_K$bifpoints[1,]; S3I_K_LP_13 <- PCR_K$bifpoints[1,]; #T=13

## Community transition along Temperature gradient at K~1E-4g/L (only for Scenario TSR*TPC) ====
# for Scenario S3A (TSR*TPC in C): ====
R_T <- PSPMequi(modelname, "EQ", c(1, 1E-4), 0.1, c(0, 0, 30), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_T <- PSPMequi(modelname, "EQ", R_T$bifpoints[c(1, 2, 5)], 0.1, c(0, 0, 30), NULL, options = c("envZE", "1", "envZE", "2"));
PCR_T_a <- PSPMequi(modelname, "EQ", CR_T$bifpoints[1,c(1, 2, 3, 7, 5)], 0.1, c(0, 0, 30), NULL, NULL);
PCR_T_b <- PSPMequi(modelname, "EQ", CR_T$bifpoints[2,c(1, 2, 3, 7, 5)], 0.1, c(0, 0, 30), NULL, NULL);

S3A_T_R <- R_T$curvepoints;
S3A_T_CR <- CR_T$curvepoints; 
S3A_T_PCR_a <- PCR_T_a$curvepoints;
S3A_T_PCR_b <- PCR_T_b$curvepoints;
S3A_T_PCR_b <- PCR_T_b$curvepoints;

S3A_T_BP_a <- R_T$bifpoints[1,];
S3A_T_BP_b <- CR_T$bifpoints[5,];
S3A_T_BPE_a <- CR_T$bifpoints[1,];
S3A_T_BPE_b <- CR_T$bifpoints[2,];
S3A_T_BPE_c <- CR_T$bifpoints[3,];
S3A_T_BPE_d <- CR_T$bifpoints[4,];

# for Scenario S3E (TSR*TPC in P): ====
R_K <- PSPMequi(modelname, "EQ", c(1.0E-06, 1.0E-06), 0.1, c(1, 0, 9E-6), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_K <- PSPMequi(modelname, "EQ", R_K$bifpoints[c(1, 2, 5)], 2, c(1, 0, 1.1E-4), NULL, options = c("envZE", "1", "envZE", "2"));
#at 1.12#-4
CR_T <- PSPMequi(modelname, "EQ", c(1, CR_K$curvepoints[8, c(2, 5)]), 0.1, c(0, 0, 30), NULL, options = c("envZE", "1", "envZE", "2"));
CR_T_b <- PSPMequi(modelname, "EQ", c(1, CR_K$curvepoints[8, c(2, 5)]), -0.1, c(0, 0, 30), NULL, options = c("envZE", "1", "envZE", "2"));
PCR_T_a <- PSPMequi(modelname, "EQ", CR_T$bifpoints[1,c(1, 2, 3, 7, 5)], -0.1, c(0, 0, 30), NULL, NULL);

S3E_T_CR <- rbind(CR_T_b$curvepoints[dim(CR_T_b$curvepoints)[1]:1,], CR_T$curvepoints[2:dim(CR_T$curvepoints)[1],]);
S3E_T_PCR <- PCR_T_a$curvepoints;

S3E_T_BPE_a <- CR_T$bifpoints[1,];
S3E_T_BPE_b <- CR_T$bifpoints[2,];
S3E_T_BPE_c <- PCR_T_a$bifpoints[3,];
S3E_T_LP_a <- PCR_T_a$bifpoints[1,];
S3E_T_LP_b <- PCR_T_a$bifpoints[2,];

# for Scenario S3I (TSR*TPC in C+P): =====
R_T <- PSPMequi(modelname, "EQ", c(1, 1E-4), 0.1, c(0, 0, 30), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_T <- PSPMequi(modelname, "EQ", R_T$bifpoints[c(1, 2, 5)], 0.1, c(0, 0, 30), NULL, options = c("envZE", "1", "envZE", "2"));
CR_T_b <- PSPMequi(modelname, "EQ", R_T$bifpoints[c(1, 2, 5)], -0.1, c(0, 0, 30), NULL, options = c("envZE", "1", "envZE", "2"));

PCR_T_a <- PSPMequi(modelname, "EQ", CR_T$bifpoints[1,c(1, 2, 3, 7, 5)], 0.1, c(0, 0, 30), NULL, NULL);
PCR_T_b <- PSPMequi(modelname, "EQ", CR_T$bifpoints[2,c(1, 2, 3, 7, 5)], 0.1, c(0, 0, 30), NULL, NULL);

S3I_T_R <- R_T$curvepoints;
S3I_T_CR <- rbind(CR_T_b$curvepoints[dim(CR_T_b$curvepoints)[1]:1,], CR_T$curvepoints[2:dim(CR_T$curvepoints)[1],]);
S3I_T_PCR_a <- PCR_T_a$curvepoints;
S3I_T_PCR_b <- PCR_T_b$curvepoints;

S3I_T_BP_a <- R_T$bifpoints[1,];
S3I_T_BP_b <- CR_T$bifpoints[5,];
S3I_T_BPE_a <- CR_T$bifpoints[1,];
S3I_T_BPE_b <- CR_T$bifpoints[2,];
S3I_T_BPE_c <- CR_T$bifpoints[3,];
S3I_T_BPE_d <- CR_T$bifpoints[4,];
S3I_T_BPE_e <- PCR_T_a$bifpoints[1,];
S3I_T_BPE_f <- PCR_T_b$bifpoints[1,];

## Community transition along productivity & Temperature gradients ====
# S3A: TPC(C) + TSR(C) **** ====
S3A_KT_BPa <- PSPMequi(modelname, "BP", c(S3A_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3A_KT_BPb <- PSPMequi(modelname, "BP", c(S3A_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3A_KT_BP<-rbind(S3A_KT_BPa$curvepoints[dim(S3A_KT_BPa$curvepoints)[1]:1,], S3A_KT_BPb$curvepoints[2:dim(S3A_KT_BPb$curvepoints)[1],]);

S3A_KT_BPEa <- PSPMequi(modelname, "BPE", c(S3A_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3A_KT_BPEb <- PSPMequi(modelname, "BPE", c(S3A_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3A_KT_BPE<-rbind(S3A_KT_BPEa$curvepoints[dim(S3A_KT_BPEa$curvepoints)[1]:1,], S3A_KT_BPEb$curvepoints[2:dim(S3A_KT_BPEb$curvepoints)[1],]);

S3A_KT_LPa <- PSPMequi(modelname, "LP", c(S3A_K_LP[c(1:5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S3A_KT_LPb <- PSPMequi(modelname, "LP", c(S3A_K_LP[c(1:5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S3A_KT_LP<-rbind(S3A_KT_LPa$curvepoints[dim(S3A_KT_LPa$curvepoints)[1]:1,], S3A_KT_LPb$curvepoints[2:dim(S3A_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S3A_KT_BP;
BPEcurve <- S3A_KT_BPE;
LPcurve  <- S3A_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S3B: TPC(C) + TSR(P) ====
S3B_KT_BPa <- PSPMequi(modelname, "BP", c(S3B_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3B_KT_BPb <- PSPMequi(modelname, "BP", c(S3B_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3B_KT_BP<-rbind(S3B_KT_BPa$curvepoints[dim(S3B_KT_BPa$curvepoints)[1]:1,], S3B_KT_BPb$curvepoints[2:dim(S3B_KT_BPb$curvepoints)[1],]);

S3B_KT_BPEa <- PSPMequi(modelname, "BPE", c(S3B_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3B_KT_BPEb <- PSPMequi(modelname, "BPE", c(S3B_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3B_KT_BPE<-rbind(S3B_KT_BPEa$curvepoints[dim(S3B_KT_BPEa$curvepoints)[1]:1,], S3B_KT_BPEb$curvepoints[2:dim(S3B_KT_BPEb$curvepoints)[1],]);

S3B_KT_LPa <- PSPMequi(modelname, "LP", c(S3B_K_LP[c(1:5)], 20), -0.4, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S3B_KT_LPb <- PSPMequi(modelname, "LP", c(S3B_K_LP[c(1:5)], 20), 0.4, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S3B_KT_LP<-rbind(S3B_KT_LPa$curvepoints[dim(S3B_KT_LPa$curvepoints)[1]:1,], S3B_KT_LPb$curvepoints[2:dim(S3B_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S3B_KT_BP;
BPEcurve <- S3B_KT_BPE;
LPcurve  <- S3B_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S3C: TPC(C) + TSR(C+P) ====
S3C_KT_BPa <- PSPMequi(modelname, "BP", c(S3C_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3C_KT_BPb <- PSPMequi(modelname, "BP", c(S3C_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3C_KT_BP<-rbind(S3C_KT_BPa$curvepoints[dim(S3C_KT_BPa$curvepoints)[1]:1,], S3C_KT_BPb$curvepoints[2:dim(S3C_KT_BPb$curvepoints)[1],]);

S3C_KT_BPEa <- PSPMequi(modelname, "BPE", c(S3C_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3C_KT_BPEb <- PSPMequi(modelname, "BPE", c(S3C_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3C_KT_BPE<-rbind(S3C_KT_BPEa$curvepoints[dim(S3C_KT_BPEa$curvepoints)[1]:1,], S3C_KT_BPEb$curvepoints[2:dim(S3C_KT_BPEb$curvepoints)[1],]);

S3C_KT_LPa <- PSPMequi(modelname, "LP", c(S3C_K_LP[c(1:5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S3C_KT_LPb <- PSPMequi(modelname, "LP", c(S3C_K_LP[c(1:5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S3C_KT_LP<-rbind(S3C_KT_LPa$curvepoints[dim(S3C_KT_LPa$curvepoints)[1]:1,], S3C_KT_LPb$curvepoints[2:dim(S3C_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S3C_KT_BP;
BPEcurve <- S3C_KT_BPE;
LPcurve  <- S3C_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S3D: TPC(P) + TSR(C) ====
S3D_KT_BPa <- PSPMequi(modelname, "BP", c(S3D_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3D_KT_BPb <- PSPMequi(modelname, "BP", c(S3D_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3D_KT_BP<-rbind(S3D_KT_BPa$curvepoints[dim(S3D_KT_BPa$curvepoints)[1]:1,], S3D_KT_BPb$curvepoints[2:dim(S3D_KT_BPb$curvepoints)[1],]);

S3D_KT_BPEa <- PSPMequi(modelname, "BPE", c(S3D_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3D_KT_BPEb <- PSPMequi(modelname, "BPE", c(S3D_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3D_KT_BPE<-rbind(S3D_KT_BPEa$curvepoints[dim(S3D_KT_BPEa$curvepoints)[1]:1,], S3D_KT_BPEb$curvepoints[2:dim(S3D_KT_BPEb$curvepoints)[1],]);

S3D_KT_LPa <- PSPMequi(modelname, "LP", c(S3D_K_LP[c(1:5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S3D_KT_LPb <- PSPMequi(modelname, "LP", c(S3D_K_LP[c(1:5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S3D_KT_LP<-rbind(S3D_KT_LPa$curvepoints[dim(S3D_KT_LPa$curvepoints)[1]:1,], S3D_KT_LPb$curvepoints[2:dim(S3D_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S3D_KT_BP;
BPEcurve <- S3D_KT_BPE;
LPcurve  <- S3D_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S3E: TPC(P) + TSR(P) **** ====
S3E_KT_BPa <- PSPMequi(modelname, "BP", c(S3E_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3E_KT_BPb <- PSPMequi(modelname, "BP", c(S3E_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3E_KT_BP  <- rbind(S3E_KT_BPa$curvepoints[dim(S3E_KT_BPa$curvepoints)[1]:1,], S3E_KT_BPb$curvepoints[2:dim(S3E_KT_BPb$curvepoints)[1],]);

S3E_KT_BPEa <- PSPMequi(modelname, "BPE", c(S3E_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3E_KT_BPEb <- PSPMequi(modelname, "BPE", c(S3E_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3E_KT_BPE  <- rbind(S3E_KT_BPEa$curvepoints[dim(S3E_KT_BPEa$curvepoints)[1]:1,], S3E_KT_BPEb$curvepoints[2:dim(S3E_KT_BPEb$curvepoints)[1],]);

S3E_KT_LPa <- PSPMequi(modelname, "LP", c(S3E_K_LP[c(1:5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S3E_KT_LPb <- PSPMequi(modelname, "LP", c(S3E_K_LP[c(1:5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S3E_KT_LP  <- rbind(S3E_KT_LPa$curvepoints[dim(S3E_KT_LPa$curvepoints)[1]:1,], S3E_KT_LPb$curvepoints[2:dim(S3E_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S3E_KT_BP;
BPEcurve <- S3E_KT_BPE;
LPcurve  <- S3E_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S3F: TPC(P) + TSR(C+P) ====
S3F_KT_BPa <- PSPMequi(modelname, "BP", c(S3F_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3F_KT_BPb <- PSPMequi(modelname, "BP", c(S3F_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3F_KT_BP<-rbind(S3F_KT_BPa$curvepoints[dim(S3F_KT_BPa$curvepoints)[1]:1,], S3F_KT_BPb$curvepoints[2:dim(S3F_KT_BPb$curvepoints)[1],]);

S3F_KT_BPEa <- PSPMequi(modelname, "BPE", c(S3F_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3F_KT_BPEb <- PSPMequi(modelname, "BPE", c(S3F_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3F_KT_BPE<-rbind(S3F_KT_BPEa$curvepoints[dim(S3F_KT_BPEa$curvepoints)[1]:1,], S3F_KT_BPEb$curvepoints[2:dim(S3F_KT_BPEb$curvepoints)[1],]);

S3F_KT_LPa <- PSPMequi(modelname, "LP", c(S3F_K_LP[c(1:5)], 20), -0.4, c(1, 1E-6, 6E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S3F_KT_LPb <- PSPMequi(modelname, "LP", c(S3F_K_LP[c(1:5)], 20), 0.4, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S3F_KT_LP<-rbind(S3F_KT_LPa$curvepoints[dim(S3F_KT_LPa$curvepoints)[1]:1,], S3F_KT_LPb$curvepoints[2:dim(S3F_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S3F_KT_BP;
BPEcurve <- S3F_KT_BPE;
LPcurve  <- S3F_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S3G: TPC(C+P) + TSR(C) ====
S3G_KT_BPa <- PSPMequi(modelname, "BP", c(S3G_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3G_KT_BPb <- PSPMequi(modelname, "BP", c(S3G_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3G_KT_BP<-rbind(S3G_KT_BPa$curvepoints[dim(S3G_KT_BPa$curvepoints)[1]:1,], S3G_KT_BPb$curvepoints[2:dim(S3G_KT_BPb$curvepoints)[1],]);

S3G_KT_BPEa <- PSPMequi(modelname, "BPE", c(S3G_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3G_KT_BPEb <- PSPMequi(modelname, "BPE", c(S3G_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3G_KT_BPE<-rbind(S3G_KT_BPEa$curvepoints[dim(S3G_KT_BPEa$curvepoints)[1]:1,], S3G_KT_BPEb$curvepoints[2:dim(S3G_KT_BPEb$curvepoints)[1],]);

S3G_KT_LPa <- PSPMequi(modelname, "LP", c(S3G_K_LP[c(1:5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S3G_KT_LPb <- PSPMequi(modelname, "LP", c(S3G_K_LP[c(1:5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S3G_KT_LP<-rbind(S3G_KT_LPa$curvepoints[dim(S3G_KT_LPa$curvepoints)[1]:1,], S3G_KT_LPb$curvepoints[2:dim(S3G_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S3G_KT_BP;
BPEcurve <- S3G_KT_BPE;
LPcurve  <- S3G_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S3H: TPC(C+P) + TSR(P) ====
S3H_KT_BPa <- PSPMequi(modelname, "BP", c(S3H_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3H_KT_BPb <- PSPMequi(modelname, "BP", c(S3H_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3H_KT_BP<-rbind(S3H_KT_BPa$curvepoints[dim(S3H_KT_BPa$curvepoints)[1]:1,], S3H_KT_BPb$curvepoints[2:dim(S3H_KT_BPb$curvepoints)[1],]);

S3H_KT_BPEa <- PSPMequi(modelname, "BPE", c(S3H_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3H_KT_BPEb <- PSPMequi(modelname, "BPE", c(S3H_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3H_KT_BPE<-rbind(S3H_KT_BPEa$curvepoints[dim(S3H_KT_BPEa$curvepoints)[1]:1,], S3H_KT_BPEb$curvepoints[2:dim(S3H_KT_BPEb$curvepoints)[1],]);

S3H_KT_LPa <- PSPMequi(modelname, "LP", c(S3H_K_LP[c(1:5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S3H_KT_LPb <- PSPMequi(modelname, "LP", c(S3H_K_LP[c(1:5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = NULL);
S3H_KT_LP<-rbind(S3H_KT_LPa$curvepoints[dim(S3H_KT_LPa$curvepoints)[1]:1,], S3H_KT_LPb$curvepoints[2:dim(S3H_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S3H_KT_BP;
BPEcurve <- S3H_KT_BPE;
LPcurve  <- S3H_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

# S3I: TPC(C+P) + TSR(C+P) **** ====
S3I_KT_BPa <- PSPMequi(modelname, "BP", c(S3I_K_BP[c(1, 2, 5)], 20), -0.1, c(1, 0, 6E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3I_KT_BPb <- PSPMequi(modelname, "BP", c(S3I_K_BP[c(1, 2, 5)], 20), 0.1, c(1, 0, 1E-3, 0, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));
S3I_KT_BP<-rbind(S3I_KT_BPa$curvepoints[dim(S3I_KT_BPa$curvepoints)[1]:1,], S3I_KT_BPb$curvepoints[2:dim(S3I_KT_BPb$curvepoints)[1],]);

S3I_KT_BPEa <- PSPMequi(modelname, "BPE", c(S3I_K_BPE[c(1,2,4,5)], 20), -0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3I_KT_BPEb <- PSPMequi(modelname, "BPE", c(S3I_K_BPE[c(1,2,4,5)], 20), 0.1, c(1, 1E-6, 6E-3, 0, 1E-2, 30), NULL, options = c("envBP", "1"));
S3I_KT_BPE<-rbind(S3I_KT_BPEa$curvepoints[dim(S3I_KT_BPEa$curvepoints)[1]:1,], S3I_KT_BPEb$curvepoints[2:dim(S3I_KT_BPEb$curvepoints)[1],]);

S3I_KT_LPa <- PSPMequi(modelname, "LP", c(S3I_K_LP[c(1:5)], 20), -0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 24.9950), NULL, options = NULL);
S3I_KT_LPb <- PSPMequi(modelname, "LP", c(S3I_K_LP[c(1:5)], 20), 0.4, c(1, 1E-6, 1E-3, 0, 1E-2, 30), NULL, options = NULL);
S3I_KT_LP<-rbind(S3I_KT_LPa$curvepoints[dim(S3I_KT_LPa$curvepoints)[1]:1,], S3I_KT_LPb$curvepoints[2:dim(S3I_KT_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S3I_KT_BP;
BPEcurve <- S3I_KT_BPE;
LPcurve  <- S3I_KT_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(log10(1.7E-6), log10(5E-3)));
axis(1, at=seq(0,30, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1)), label=F, las=2, cex.axis=1.8);
lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

############ Section 4: TPC+TSR implementing species thermal (mis)match (S4) ====
# S4A-C) Thermal niche shift of consumer relative to predator; TSR(C+P) + A: TPC(C), B: TPC(P), C: TPC(C+P)
# S4D-F) Thermal niche shift of predator relative to consumer; TSR(C+P) + D: TPC(C), E: TPC(P), F: TPC(C+P)

# Note that the environmental temperature is kept at 20°C and the thermal niche of the focal species is kept identical at in the previous scenario (~5-25°C, optimum at 20°C)
# Delta TPC < 0 => Cold-adapted species relative to the focal species
# Delta TPC > 0 => Warm-adapted species relative to the focal species

## Community transition along productivity gradient ====
modelname_1 = "PSPM_TPC-TSR_model.R" 
modelname_2 = "PSPM_TPC-TSR_Thermal-Mismatch_model.R" 

R_K <- PSPMequi(modelname_2, "EQ", c(1.0E-06, 1.0E-06), 0.8, c(1, 0, 5E-5), NULL,
                options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_K <- PSPMequi(modelname_2, "EQ", R_K$bifpoints[c(1, 2, 5)], 2, c(1, 0, 1E-3), NULL,
                 options = c("envZE", "1", "envZE", "2"));
PCR_K <- PSPMequi(modelname_2, "EQ", CR_K$bifpoints[,c(1, 2, 3, 7, 5)], -4, c(1, 0, 1E-3), NULL, NULL);

R_K$bifpoints; R_K$biftypes;
CR_K$bifpoints; CR_K$biftypes; 
PCR_K$bifpoints; PCR_K$biftypes;

#S4B_Delta0_R
#CR_Delta <- PSPMequi(modelname_2, "EQ", c(7, CR_K$curvepoints[4, c(2, 5)]), 0.01, c(30, 5, 8), NULL, options = c("envZE", "1", "envZE", "2"));
#PCR_Delta <- PSPMequi(modelname_2, "EQ", CR_Delta$bifpoints[, c(1:3, 7, 5)], 0.5, c(30, 5, 8), NULL, NULL);

## Data and points saved along K productivity gradient for each panel ====
# S4A_K_R <- R_K$curvepoints; S4A_K_CR <- CR_K$curvepoints; S4A_K_PCR <- PCR_K$curvepoints; # for 20?C + Delta= 1
# S4B_K_R <- R_K$curvepoints; S4B_K_CR <- CR_K$curvepoints; S4B_K_PCR <- PCR_K$curvepoints; # for 20?C + Delta= 1
# S4C_K_R <- R_K$curvepoints; S4C_K_CR <- CR_K$curvepoints; S4C_K_PCR <- PCR_K$curvepoints; # for 13?C + Delta= 1
# S4D_K_R <- R_K$curvepoints; S4D_K_CR <- CR_K$curvepoints; S4D_K_PCR <- PCR_K$curvepoints; # for 13?C + Delta = -5

# Bifurcation points
# S4A_K_BP <- R_K$bifpoints[1,]; S4A_K_BPE <- CR_K$bifpoints[1,]; S4A_K_LP <- PCR_K$bifpoints[1,];
# S4B_K_BP <- R_K$bifpoints[1,]; S4B_K_BPE <- CR_K$bifpoints[1,]; S4B_K_LP <- PCR_K$bifpoints[1,];
# S4C_K_BP <- R_K$bifpoints[1,]; S4C_K_BPE <- CR_K$bifpoints[1,]; S4C_K_LP <- PCR_K$bifpoints[1,];
# S4D_K_BP <- R_K$bifpoints[1,]; S4D_K_BPE <- CR_K$bifpoints[1,]; S4D_K_LP <- PCR_K$bifpoints[1,];

## Community transition along productivity gradient & axis of thermal mismatch ====
# S4A: TSR+TPC(C+P): Thermal niche shift of consumer relative to predator at 20degC ====
S4A_TM_BPa <- PSPMequi(modelname_2, "BP", c(S4A_K_BP[c(1, 2, 5)], 1), -1, c(1, 0, 1E-3, 30, -10, 15), NULL, options = c("envZE", "1", "popBP", "0"));
S4A_TM_BPb <- PSPMequi(modelname_2, "BP", c(S4A_K_BP[c(1, 2, 5)], 1), 1, c(1, 0, 1E-3, 30, -10, 15), NULL, options = c("envZE", "1", "popBP", "0"));
S4A_TM_BP<-rbind(S4A_TM_BPa$curvepoints[dim(S4A_TM_BPa$curvepoints)[1]:1,], S4A_TM_BPb$curvepoints[2:dim(S4A_TM_BPb$curvepoints)[1],]);

S4A_TM_BPEa <- PSPMequi(modelname_2, "BPE", c(S4A_K_BPE[c(1,2,4,5)], 1), -1, c(1, 1E-6, 6E-3, 30, -10, 15), NULL, options = c("envBP", "1"));
S4A_TM_BPEb <- PSPMequi(modelname_2, "BPE", c(S4A_K_BPE[c(1,2,4,5)], 1), 1, c(1, 1E-6, 6E-3, 30, -10, 15), NULL, options = c("envBP", "1"));
S4A_TM_BPE<-rbind(S4A_TM_BPEa$curvepoints[dim(S4A_TM_BPEa$curvepoints)[1]:1,], S4A_TM_BPEb$curvepoints[2:dim(S4A_TM_BPEb$curvepoints)[1],]);

S4A_TM_LPa <- PSPMequi(modelname_2, "LP", c(S4A_K_LP[c(1:5)], 1), -1, c(1, 1E-6, 1E-3, 30, -10, 15), NULL, options = NULL);
S4A_TM_LPb <- PSPMequi(modelname_2, "LP", c(S4A_K_LP[c(1:5)], 1), 1, c(1, 1E-6, 1E-3, 30, -10, 15), NULL, options = NULL);
S4A_TM_LPc <- PSPMequi(modelname_2, "LP", c(S4A_TM_LPb$curvepoints[dim(S4A_TM_LPb$curvepoints)[1],1:6]), 1, c(1, 1E-6, 1E-4, 30, -4.89, 15), NULL, options = NULL);
S4A_TM_LPd <- PSPMequi(modelname_2, "LP", c(S4A_TM_LPc$curvepoints[dim(S4A_TM_LPc$curvepoints)[1],1:6]), 1, c(1, 1E-6, 1E-4, 30, -4.89, 15), NULL, options = NULL);
S4A_TM_LPe <- PSPMequi(modelname_2, "LP", c(S4A_TM_LPd$curvepoints[dim(S4A_TM_LPd$curvepoints)[1],1:6]), 1, c(1, 1E-6, 1E-4, 30, -4.89, 15), NULL, options = NULL);
S4A_TM_LPf <- PSPMequi(modelname_2, "LP", c(S4A_TM_LPe$curvepoints[dim(S4A_TM_LPe$curvepoints)[1],1:6]), -2, c(1, 1E-6, 1E-4, 30, -4.89, 15), NULL, options = NULL);

S4A_TM_LP<-rbind( 
  S4A_TM_LPf$curvepoints[dim(S4A_TM_LPf$curvepoints)[1]:2,],
  S4A_TM_LPe$curvepoints[dim(S4A_TM_LPe$curvepoints)[1]:2,], 
  S4A_TM_LPd$curvepoints[dim(S4A_TM_LPd$curvepoints)[1]:2,],
  S4A_TM_LPc$curvepoints[dim(S4A_TM_LPc$curvepoints)[1]:2,],
  S4A_TM_LPb$curvepoints[dim(S4A_TM_LPb$curvepoints)[1]:2,],
  S4A_TM_LPa$curvepoints[1:dim(S4A_TM_LPa$curvepoints)[1],]);

#graph test
BPcurve  <- S4A_TM_BP;
BPEcurve <- S4A_TM_BPE;
LPcurve  <- S4A_TM_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(-15,15), ylim=c(log10(1.7E-6), log10(1E-3)));
axis(1, at=seq(-10,20, by=5), label=T, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-5, 1E-4)), label=c(expression(10^-5), expression(10^-4)), las=2, cex.axis=1.8);
axis(3, at=seq(-10,20, by=5), label=seq(10,40, by=5), las=1, cex.axis=1.8);

lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(S4A_TM_LPa$curvepoints[,1])~S4A_TM_LPa$curvepoints[,6], type="l", lwd=2, lty=2);
lines(log10(S4A_TM_LPb$curvepoints[,1])~S4A_TM_LPb$curvepoints[,6], type="l", lwd=2, lty=2);
lines(log10(S4A_TM_LPc$curvepoints[,1])~S4A_TM_LPc$curvepoints[,6], type="l", lwd=2, lty=2);

lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

abline(v=0, lwd=2, lty=4, col="grey");

# S4B: TSR+TPC(C+P): Thermal niche shift of predator relative to consumer at 20degC ====
# When consumer is the focal species, its invasion threshold at 20 degC is not affected by shifts in predator thermal niche 
S4B_TM_BP = S4B_K_BP;

S4B_TM_BPEa <- PSPMequi(modelname_2, "BPE", c(S4B_K_BPE[c(1,2,4,5)], 1), -1, c(1, 1E-6, 6E-3, 30, -10, 15), NULL, options = c("envBP", "1"));
S4B_TM_BPEb <- PSPMequi(modelname_2, "BPE", c(S4B_K_BPE[c(1,2,4,5)], 1), 1, c(1, 1E-6, 6E-3, 30, -10, 15), NULL, options = c("envBP", "1"));
S4B_TM_BPE<-rbind(S4B_TM_BPEa$curvepoints[dim(S4B_TM_BPEa$curvepoints)[1]:1,], S4B_TM_BPEb$curvepoints[2:dim(S4B_TM_BPEb$curvepoints)[1],]);

S4B_TM_LPa <- PSPMequi(modelname_2, "LP", c(S4B_K_LP[c(1:5)], 1), -1, c(1, 1E-6, 1E-3, 30, -10, 15), NULL, options = NULL);
S4B_TM_LPb <- PSPMequi(modelname_2, "LP", S4B_TM_LPa$curvepoints[dim(S4B_TM_LPa$curvepoints)[1],c(1:6)], -1, c(1, 1E-6, 1E-3, 30, -10, 15), NULL, options = NULL);

S4B_TM_LPc <- PSPMequi(modelname_2, "LP", c(S4B_K_LP[c(1:5)], 1), 1, c(1, 1E-6, 1E-3, 30, -10, 15), NULL, options = NULL);
S4B_TM_LPd <- PSPMequi(modelname_2, "LP", S4B_TM_LPc$curvepoints[dim(S4B_TM_LPc$curvepoints)[1],c(1:6)], 2, c(1, 1E-6, 1E-3, 30, -10, 15), NULL, options = NULL);
S4B_TM_LPe <- PSPMequi(modelname_2, "LP", S4B_TM_LPd$curvepoints[dim(S4B_TM_LPd$curvepoints)[1],c(1:6)], 2, c(1, 1E-6, 1E-3, 30, -10, 15), NULL, options = NULL);
S4B_TM_LPf <- PSPMequi(modelname_2, "LP", S4B_TM_LPe$curvepoints[dim(S4B_TM_LPe$curvepoints)[1],c(1:6)], 1, c(1, 1E-6, 1E-3, 30, -10, 15), NULL, options = NULL);

S4B_TM_LP<-rbind(
  S4B_TM_LPf$curvepoints[dim(S4B_TM_LPf$curvepoints)[1]:2,],
  S4B_TM_LPe$curvepoints[dim(S4B_TM_LPe$curvepoints)[1]:2,],
  S4B_TM_LPd$curvepoints[dim(S4B_TM_LPd$curvepoints)[1]:2,],
  S4B_TM_LPc$curvepoints[dim(S4B_TM_LPc$curvepoints)[1]:2,],
  S4B_TM_LPa$curvepoints[1:dim(S4B_TM_LPa$curvepoints)[1],],
  S4B_TM_LPb$curvepoints[2:dim(S4B_TM_LPb$curvepoints)[1],]);

#graph test
BPcurve  <- S4B_TM_BP;
BPEcurve <- S4B_TM_BPE;
LPcurve  <- S4B_TM_LP;

plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(-15,15), ylim=c(log10(1.7E-6), log10(1E-3)));
axis(1, at=seq(-10,20, by=5), label=T, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-5, 1E-4)), label=c(expression(10^-5), expression(10^-4)), las=2, cex.axis=1.8);
abline(h=log10(S4B_TM_BP[1]), type="l", lwd=2, lty=3);
lines(log10(S4B_TM_BPE[,1])~S4B_TM_BPE[,6], type="l", lwd=2, lty=1);
lines(log10(S4B_TM_LP[,1])~S4B_TM_LP[,6], type="l", lwd=2, lty=2);

abline(v=0, lwd=2, lty=4, col="grey");
abline(v=min(S4B_TM_LP[,6]), lwd=2, lty=5);

abline(h=log10(5.47244780E-05)); abline(v=7)
abline(h=log10(1.323151e-04))
# S4C: TSR+TPC(C+P): Thermal niche shift of consumer relative to predator at 13degC ====
S4C_TM_BPa <- PSPMequi(modelname_2, "BP", c(S4C_K_BP[c(1, 2, 5)], 1), -1, c(1, 0, 1E-3, 30, -15, 15), NULL, options = c("envZE", "1", "popBP", "0"));
S4C_TM_BPb <- PSPMequi(modelname_2, "BP", c(S4C_K_BP[c(1, 2, 5)], 1), 1, c(1, 0, 1E-3, 30, -15, 15), NULL, options = c("envZE", "1", "popBP", "0"));
S4C_TM_BP<-rbind(S4C_TM_BPa$curvepoints[dim(S4C_TM_BPa$curvepoints)[1]:1,], S4C_TM_BPb$curvepoints[2:dim(S4C_TM_BPb$curvepoints)[1],]);

S4C_TM_BPEa <- PSPMequi(modelname_2, "BPE", c(S4C_K_BPE[c(1,2,4,5)], 1), -1, c(1, 1E-6, 6E-3, 30, -15, 10), NULL, options = c("envBP", "1"));
S4C_TM_BPEb <- PSPMequi(modelname_2, "BPE", c(S4C_K_BPE[c(1,2,4,5)], 1), 1, c(1, 1E-6, 6E-3, 30, -15, 10), NULL, options = c("envBP", "1"));
S4C_TM_BPE<-rbind(S4C_TM_BPEa$curvepoints[dim(S4C_TM_BPEa$curvepoints)[1]:1,], S4C_TM_BPEb$curvepoints[2:dim(S4C_TM_BPEb$curvepoints)[1],]);

S4C_TM_LPa <- PSPMequi(modelname_2, "LP", c(S4C_K_LP[c(1:5)], 1), -1, c(1, 1E-6, 1E-4, 30, -10, 10), NULL, options = NULL);
S4C_TM_LPb <- PSPMequi(modelname_2, "LP", c(S4C_K_LP[c(1:5)], 1), 2, c(1, 1E-6, 1E-4, 30, -10, 10), NULL, options = NULL);

S4C_TM_LP<-rbind(
  S4C_TM_LPb$curvepoints[dim(S4C_TM_LPb$curvepoints)[1]:2,],
  S4C_TM_LPa$curvepoints[dim(S4C_TM_LPa$curvepoints)[1]:1,]);

#graph test
BPcurve  <- S4C_TM_BP;
BPEcurve <- S4C_TM_BPE;
LPcurve  <- S4C_TM_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(-20,20), ylim=c(log10(1.7E-6), log10(1E-3)));
axis(1, at=seq(-20,20, by=5), label=T, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-5, 1E-4)), label=c(expression(10^-5), expression(10^-4)), las=2, cex.axis=1.8);
axis(3, at=seq(-20,20, by=5), label=seq(-7,33, by=5), las=1, cex.axis=1.8);

lines(log10(BPcurve[,1])~BPcurve[,6], type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);
abline(v=0, lwd=2, lty=4, col="grey");

# S4D: TSR+TPC(C+P): Thermal niche shift of predator relative to consumer at 13degC ====
# When consumer is the focal species, its invasion threshold at 13 degC is not affected by shifts in predator thermal niche 
S4D_TM_BP = S4D_K_BP;

S4D_TM_BPEa <- PSPMequi(modelname_2, "BPE", c(S4D_K_BPE[c(1,2,4,5)], -5), -1, c(1, 1E-6, 6E-3, 30, -15, 10), NULL, options = c("envBP", "1"));
S4D_TM_BPEb <- PSPMequi(modelname_2, "BPE", c(S4D_K_BPE[c(1,2,4,5)], -5), 1, c(1, 1E-6, 6E-3, 30, -15, 10), NULL, options = c("envBP", "1"));
S4D_TM_BPEc <- PSPMequi(modelname_2, "BPE", S4D_TM_BPEa$curvepoints[dim(S4D_TM_BPEa$curvepoints)[1],c(1:2,4:6)], -1, c(1, 1E-6, 6E-3, 30, -15, 10), NULL, options = c("envBP", "1"));

S4D_TM_BPE<-rbind(
  S4D_TM_BPEb$curvepoints[dim(S4D_TM_BPEb$curvepoints)[1]:2,],
  S4D_TM_BPEa$curvepoints[1:dim(S4D_TM_BPEa$curvepoints)[1],],
  S4D_TM_BPEc$curvepoints[2:dim(S4D_TM_BPEc$curvepoints)[1],]);

S4D_TM_LPa <- PSPMequi(modelname_2, "LP", c(S4D_K_LP[c(1:5)], -5), 1, c(1, 1E-6, 1E-4, 30, -10, 10), NULL, options = NULL);
S4D_TM_LPb <- PSPMequi(modelname_2, "LP", c(S4D_K_LP[c(1:5)], -5), -1, c(1, 1E-6, 1E-4, 30, -10, 10), NULL, options = NULL);

S4D_TM_LPc <- PSPMequi(modelname_2, "LP", S4D_TM_LPb$curvepoints[dim(S4D_TM_LPb$curvepoints)[1],c(1:6)], -2, c(1, 1E-6, 1E-3, 30, -20, 15), NULL, options = NULL);
S4D_TM_LPd <- PSPMequi(modelname_2, "LP", S4D_TM_LPb$curvepoints[dim(S4D_TM_LPb$curvepoints)[1],c(1:6)], 10, c(1, 1E-6, 1E-3, 30, -20, 15), NULL, options = NULL);

S4D_TM_LP<-rbind(
  S4D_TM_LPa$curvepoints[dim(S4D_TM_LPa$curvepoints)[1]:1,],
  S4D_TM_LPb$curvepoints[2:dim(S4D_TM_LPb$curvepoints)[1],],
  S4D_TM_LPc$curvepoints[2:dim(S4D_TM_LPc$curvepoints)[1],]);

#graph test
BPcurve  <- S4D_TM_BP;
BPEcurve <- S4D_TM_BPE;
LPcurve  <- S4D_TM_LP;
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(-15,15), ylim=c(log10(1.7E-6), log10(1E-3)));
axis(1, at=seq(-20,20, by=5), label=T, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-5, 1E-4)), label=c(expression(10^-5), expression(10^-4)), las=2, cex.axis=1.8);
#axis(3, at=seq(-20,20, by=5), label=seq(-7,33, by=5), las=1, cex.axis=1.8);
abline(h=log10(S4D_TM_BP[1]), type="l", lwd=2, lty=3);
lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);

abline(v=0, lwd=2, lty=4, col="grey");
abline(v=min(S4D_TM_BPE[,6]), lwd=2, lty=4, col="grey");

############ Section 5: Consumer life histories over temperature gradient under combinations of TPC+TSR (S5) ====
modelname_1 = "PSPM_TPC-TSR_model.R"
modelname_3 = "PSPM_TPC-TSR_Demo_model.R" 

# Community transition along productivity gradient for T = 1
R_K <- PSPMequi(modelname_1, "EQ", c(1.0E-06, 1.0E-06), 0.5, c(1, 0, 1E-4), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_K <- PSPMequi(modelname_1, "EQ", R_K$bifpoints[c(1, 2, 5)], 2, c(1, 0, 1E-3), NULL, options = c("envZE", "1", "envZE", "2"));
PCR_K <- PSPMequi(modelname_1, "EQ", CR_K$bifpoints[,c(1, 2, 3, 7, 5)], -2, c(1, 0, 2E-3), NULL, NULL);

R_K$bifpoints; R_K$biftypes;
CR_K$bifpoints; CR_K$biftypes;
PCR_K$bifpoints; PCR_K$biftypes;
# Community transitions along Temperature gradient under joint influences of TSR & TPC in consumer sizes ====
# S5A-S5D TSR in traits in absence of Lv
# S5E-S5G extend S5B-S5D to include TSR in Lv
# S5H-S5J TSR & TPC in absence of Lv
# S5K-S5M extend S5H-S5J to include TSR in Lv
# S5A: TSR in Lv ====
CR_T <- PSPMequi(modelname_1, "EQ", c(1, CR_K$curvepoints[9,c(2,5)]), 0.1, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"));
PCR_T <- PSPMequi(modelname_1, "EQ", CR_T$bifpoints[,c(1:3, 7, 5)], 0.4, c(0, 1, 30), NULL, NULL);

S5A_T_CR_II <- CR_T$curvepoints; S5A_T_PCR_II <- PCR_T$curvepoints;
S5A_T_BPE <- CR_T$bifpoints; S5A_T_LP <- PCR_T$bifpoints;

PGR <- PSPMdemo(modelname_3, c(0, 1, 0.8, 1, 30));
S5A_PGR_II <- PGR$curvepoints;

# S5B: TSR in Lmat ====
R_T <- PSPMequi(modelname, "EQ", c(1, 1E-5), 0.1, c(0, 1, 30), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_T <- PSPMequi(modelname, "EQ", R_T$bifpoints[,c(1,2,5)], 0.1, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);

CR_T_II <- PSPMequi(modelname_1, "EQ", c(1, CR_K$curvepoints[6,c(2,5)]), 0.1, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"));
# Obtained With T=30:
PCR_T_IIa <- PSPMequi(modelname_1, "EQ", c(30, PCR_K$curvepoints[4,c(2:3,7,5)]), -0.4, c(0, 1, 30), NULL, NULL);

S5B_T_R_I <- R_T$curvepoints;
S5B_T_CR_I <- CR_T$curvepoints;
S5B_T_CR_II <- CR_T_II$curvepoints;
S5B_T_PCR_IIa <- PCR_T_IIa$curvepoints;

S5B_T_BP <- R_T$bifpoints;
S5B_T_LPa <- PCR_T_IIa$bifpoints;

PGR <- PSPMdemo(modelname_3, c(0, 1, 0.8, 1, 30));
S5B_PGR_II <- PGR$curvepoints;

S1D_KT_BPb <- PSPMequi(modelname, "BP", c(S1D_K_BP[c(1, 2, 5)], 20), 0.8, c(1, 0, 6E-3, 10, 1E-2, 30), NULL, options = c("envZE", "1", "popBP", "0"));

# S5C: TSR in Linf ====
R_T <- PSPMequi(modelname, "EQ", c(1, 1E-5), 0.1, c(0, 1, 30), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_T <- PSPMequi(modelname, "EQ", R_T$bifpoints[,c(1,2,5)], -0.1, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);

CR_T_II <- PSPMequi(modelname_1, "EQ", c(1, CR_K$curvepoints[10,c(2,5)]), 0.1, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"));
PCR_T_IIa <- PSPMequi(modelname_1, "EQ", c(1, PCR_K$curvepoints[9,c(2:3,7,5)]), 0.4, c(0, 1, 30), NULL, NULL);
PCR_T_IIb <- PSPMequi(modelname_1, "EQ", c(1, PCR_K$curvepoints[20,c(2:3,7,5)]), 0.4, c(0, 1, 30), NULL, NULL);

S5C_T_R_I <- R_T$curvepoints;
S5C_T_CR_I <- CR_T$curvepoints;
S5C_T_CR_II <- CR_T_II$curvepoints;
S5C_T_PCR_IIa <- PCR_T_IIa$curvepoints;
S5C_T_PCR_IIb <- PCR_T_IIb$curvepoints;

S5C_T_BP <- R_T$bifpoints;
S5C_T_LPa <- PCR_T_IIa$bifpoints;
S5C_T_LPb <- PCR_T_IIb$bifpoints;

PGR <- PSPMdemo(modelname_3, c(0, 1, 0.8, 1, 30));
S5C_PGR_II <- PGR$curvepoints;

# S5D: TSR in Lmat + Linf ====
CR_T_II <- PSPMequi(modelname_1, "EQ", c(1, CR_K$curvepoints[8,c(2,5)]), 0.1, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"));
# Obtained With T=30:
PCR_T_IIa <- PSPMequi(modelname_1, "EQ", c(30, PCR_K$curvepoints[6,c(2:3,7,5)]), -0.4, c(0, 1, 30), NULL, NULL);
PCR_T_IIb <- PSPMequi(modelname_1, "EQ", c(30, PCR_K$curvepoints[17,c(2:3,7,5)]), -0.4, c(0, 1, 30), NULL, NULL);

S5D_T_CR_II <- CR_T_II$curvepoints;
S5D_T_PCR_IIa <- PCR_T_IIa$curvepoints;
S5D_T_PCR_IIb <- PCR_T_IIb$curvepoints;

S5D_T_LPa <- PCR_T_IIa$bifpoints;
S5D_T_LPb <- PCR_T_IIb$bifpoints;

PGR <- PSPMdemo(modelname_3, c(0, 1, 0.8, 1, 30));
S5D_PGR_II <- PGR$curvepoints;

# S5E: TSR in Lv + Lmat ====
R_T <- PSPMequi(modelname, "EQ", c(1, 1E-5), 0.1, c(0, 1, 30), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_T <- PSPMequi(modelname, "EQ", R_T$bifpoints[,c(1,2,5)], 0.1, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);

CR_T_II <- PSPMequi(modelname_1, "EQ", c(1, CR_K$curvepoints[6,c(2,5)]), 0.1, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"));
PCR_T_IIa <- PSPMequi(modelname_1, "EQ", c(1, PCR_K$curvepoints[11,c(2:3,7,5)]), 0.1, c(0, 1, 30), NULL, NULL);
PCR_T_IIb <- PSPMequi(modelname_1, "EQ", c(1, PCR_K$curvepoints[16,c(2:3,7,5)]), 0.1, c(0, 1, 30), NULL, NULL);

S5E_T_R_I <- R_T$curvepoints;
S5E_T_CR_I <- CR_T$curvepoints;
S5E_T_CR_II <- CR_T_II$curvepoints;
S5E_T_PCR_IIa <- PCR_T_IIa$curvepoints;
S5E_T_PCR_IIb <- PCR_T_IIb$curvepoints;
S5E_T_BP <- R_T$bifpoints;

PGR <- PSPMdemo(modelname_3, c(0, 1, 0.8, 1, 30));
S5E_PGR_II <- PGR$curvepoints;

# S5F: TSR in Lv + Linf ====
R_T <- PSPMequi(modelname, "EQ", c(1, 1E-5), 0.1, c(0, 1, 30), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_T <- PSPMequi(modelname, "EQ", R_T$bifpoints[,c(1,2,5)], -0.1, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);

CR_T_II <- PSPMequi(modelname_1, "EQ", c(1, CR_K$curvepoints[11,c(2,5)]), 0.1, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"));
PCR_T_II <- PSPMequi(modelname_1, "EQ", CR_T_II$bifpoints[,c(1:3, 7, 5)], 0.1, c(0, 1, 30), NULL, NULL);

S5F_T_R_I <- R_T$curvepoints; S5F_T_CR_I <- CR_T$curvepoints;
S5F_T_CR_II <- CR_T_II$curvepoints;
S5F_T_PCR_II <- PCR_T_II$curvepoints;

S5F_T_BP <- R_T$bifpoints;
S5F_T_BPE <- CR_T_II$bifpoints;
S5F_T_LP <- PCR_T_II$bifpoints;

PGR <- PSPMdemo(modelname_3, c(0, 1, 0.8, 1, 30));
S5F_PGR_II <- PGR$curvepoints;
# S5G: TSR in Lv + Lmat + Linf ====
CR_T_II <- PSPMequi(modelname_1, "EQ", c(1, CR_K$curvepoints[8,c(2,5)]), 0.1, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"));
PCR_T_IIa <- PSPMequi(modelname_1, "EQ", c(1, PCR_K$curvepoints[9,c(2:3,7,5)]), 0.4, c(0, 1, 30), NULL, NULL);
PCR_T_IIb <- PSPMequi(modelname_1, "EQ", c(1, PCR_K$curvepoints[20,c(2:3,7,5)]), 0.4, c(0, 1, 30), NULL, NULL);

S5G_T_CR_II <- CR_T_II$curvepoints;
S5G_T_PCR_IIa <- PCR_T_IIa$curvepoints;
S5G_T_PCR_IIb <- PCR_T_IIb$curvepoints;

S5G_T_LPa <- PCR_T_IIa$bifpoints;
S5G_T_LPb <- PCR_T_IIb$bifpoints;

PGR <- PSPMdemo(modelname_3, c(0, 1, 0.8, 1, 30));
S5G_PGR_II <- PGR$curvepoints;

# S5H: TPC+TSR (Lmat) ====
R_T <- PSPMequi(modelname, "EQ", c(1, 1E-5), 0.1, c(0, 1, 30), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_T <- PSPMequi(modelname, "EQ", R_T$bifpoints[,c(1,2,5)], 0.4, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);

PCR_Ta <- PSPMequi(modelname, "EQ", CR_T$bifpoints[1,c(1:3, 7, 5)], 0.4, c(0, 1, 30), NULL, NULL, clean = TRUE, force=TRUE);
PCR_Tb <- PSPMequi(modelname, "EQ", CR_T$bifpoints[2,c(1:3, 7, 5)], 0.4, c(0, 1, 30), NULL, NULL, clean = TRUE, force=TRUE);#unstable branch

S5H_T_R_II <- R_T$curvepoints;
S5H_T_CR_II <- CR_T$curvepoints; 
S5H_T_PCR_IIa <- PCR_Ta$curvepoints; 
S5H_T_PCR_IIb <- PCR_Tb$curvepoints;

S5H_T_BP_a <- R_T$bifpoints;
S5H_T_BP_b <- CR_T$bifpoints[5,];
S5H_T_BPE_a <- CR_T$bifpoints[1,]; 
S5H_T_BPE_b <- CR_T$bifpoints[2,];
S5H_T_BPE_c <- CR_T$bifpoints[3,]; 
S5H_T_BPE_d <- CR_T$bifpoints[4,];
S5H_T_BPE_e <- PCR_Ta$bifpoints[1,]; 
S5H_T_BPE_f <- PCR_Tb$bifpoints[1,];

PGR <- PSPMdemo(modelname_3, c(0, 5.5, 0.8, 5.5, 24.96));
S5H_PGR_II <- PGR$curvepoints;

S5H_KT_BPa <- PSPMequi(modelname, "BP", c(S5H_T_BP_a[c(1, 2, 5)], 1E-4), -2, c(0, 0.4, 30, 1, 0, 6E-3), NULL, options = c("envZE", "1", "popBP", "0"));
S5H_KT_BPb <- PSPMequi(modelname, "BP", c(S5H_T_BP_a[c(1, 2, 5)], 1E-4), 3, c(0, 0.8, 30, 1, 0, 6E-3), NULL, options = c("envZE", "1", "popBP", "0"));
S5H_KT_BP<-rbind(S5H_KT_BPa$curvepoints[dim(S5H_KT_BPa$curvepoints)[1]:1,], S5H_KT_BPb$curvepoints[2:dim(S5H_KT_BPb$curvepoints)[1],]);

# S5I: TPC+TSR (Linf) ====
R_T <- PSPMequi(modelname, "EQ", c(1, 1E-5), 0.1, c(0, 1, 30), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_T <- PSPMequi(modelname, "EQ", R_T$bifpoints[,c(1,2,5)], 0.4, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);

PCR_Ta <- PSPMequi(modelname, "EQ", CR_T$bifpoints[1,c(1:3, 7, 5)], 0.4, c(0, 1, 30), NULL, NULL, clean = TRUE, force=TRUE);
PCR_Tb <- PSPMequi(modelname, "EQ", CR_T$bifpoints[2,c(1:3, 7, 5)], 0.4, c(0, 1, 30), NULL, NULL, clean = TRUE, force=TRUE);#unstable branch

S5I_T_R_II <- R_T$curvepoints;
S5I_T_CR_II <- CR_T$curvepoints; 
S5I_T_PCR_IIa <- PCR_Ta$curvepoints;
S5I_T_PCR_IIb <- PCR_Tb$curvepoints;

S5I_T_BP_a <- R_T$bifpoints;
S5I_T_BP_b <- CR_T$bifpoints[5,];
S5I_T_BPE_a <- CR_T$bifpoints[1,];
S5I_T_BPE_b <- CR_T$bifpoints[2,];
S5I_T_BPE_c <- CR_T$bifpoints[3,];
S5I_T_BPE_d <- CR_T$bifpoints[4,];
S5I_T_BPE_e <- PCR_Ta$bifpoints[1,];
S5I_T_BPE_f <- PCR_Tb$bifpoints[1,];

PGR <- PSPMdemo(modelname_3, c(0, 5.5, 0.8, 5.5, 24.96));
S5I_PGR_II <- PGR$curvepoints;

S5I_KT_BPa <- PSPMequi(modelname, "BP", c(S5I_T_BP_a[c(1, 2, 5)], 1E-4), -2, c(0, 0.4, 30, 1, 0, 6E-3), NULL, options = c("envZE", "1", "popBP", "0"));
S5I_KT_BPb <- PSPMequi(modelname, "BP", c(S5I_T_BP_a[c(1, 2, 5)], 1E-4), 3, c(0, 0.8, 30, 1, 0, 6E-3), NULL, options = c("envZE", "1", "popBP", "0"));
S5I_KT_BP<-rbind(S5I_KT_BPa$curvepoints[dim(S5I_KT_BPa$curvepoints)[1]:1,], S5I_KT_BPb$curvepoints[2:dim(S5I_KT_BPb$curvepoints)[1],]);

# S5J: TPC+TSR (Lmat + Linf) ====
R_T <- PSPMequi(modelname, "EQ", c(1, 1E-5), 0.1, c(0, 1, 30), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_T <- PSPMequi(modelname, "EQ", R_T$bifpoints[,c(1,2,5)], 0.4, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);

PCR_Ta <- PSPMequi(modelname, "EQ", CR_T$bifpoints[1,c(1:3, 7, 5)], 0.8, c(0, 1, 30), NULL, NULL, clean = TRUE, force=TRUE);
PCR_Tb <- PSPMequi(modelname, "EQ", CR_T$bifpoints[2,c(1:3, 7, 5)], 0.4, c(0, 1, 30), NULL, NULL, clean = TRUE, force=TRUE);#unstable branch

S5J_T_R_II <- R_T$curvepoints;
S5J_T_CR_II <- CR_T$curvepoints; 
S5J_T_PCR_IIa <- PCR_Ta$curvepoints;
S5J_T_PCR_IIb <- PCR_Tb$curvepoints;

S5J_T_BP_a <- R_T$bifpoints;
S5J_T_BP_b <- CR_T$bifpoints[5,];
S5J_T_BPE_a <- CR_T$bifpoints[1,];
S5J_T_BPE_b <- CR_T$bifpoints[2,];
S5J_T_BPE_c <- CR_T$bifpoints[3,];
S5J_T_BPE_d <- CR_T$bifpoints[4,];
S5J_T_BPE_e <- PCR_Ta$bifpoints[1,];
S5J_T_BPE_f <- PCR_Tb$bifpoints[1,];

PGR <- PSPMdemo(modelname_3, c(0, 5.5, 0.8, 5.5, 24.96));
S5J_PGR_II <- PGR$curvepoints;

S5J_KT_BPa <- PSPMequi(modelname, "BP", c(S5J_T_BP_a[c(1, 2, 5)], 1E-4), -2, c(0, 0.4, 30, 1, 0, 6E-3), NULL, options = c("envZE", "1", "popBP", "0"));
S5J_KT_BPb <- PSPMequi(modelname, "BP", c(S5J_T_BP_a[c(1, 2, 5)], 1E-4), 3, c(0, 0.8, 30, 1, 0, 6E-3), NULL, options = c("envZE", "1", "popBP", "0"));
S5J_KT_BP<-rbind(S5J_KT_BPa$curvepoints[dim(S5J_KT_BPa$curvepoints)[1]:1,], S5J_KT_BPb$curvepoints[2:dim(S5J_KT_BPb$curvepoints)[1],]);

# S5K: TPC+TSR (Lmat) + TSR(LV) ====
R_T <- PSPMequi(modelname, "EQ", c(1, 1E-5), 0.1, c(0, 1, 30), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_T <- PSPMequi(modelname, "EQ", R_T$bifpoints[,c(1,2,5)], 0.4, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);

PCR_Ta <- PSPMequi(modelname, "EQ", CR_T$bifpoints[1,c(1:3, 7, 5)], 0.4, c(0, 1, 30), NULL, NULL, clean = TRUE, force=TRUE);
PCR_Tb <- PSPMequi(modelname, "EQ", CR_T$bifpoints[2,c(1:3, 7, 5)], 0.4, c(0, 1, 30), NULL, NULL, clean = TRUE, force=TRUE);#unstable branch

S5K_T_R_II <- R_T$curvepoints;
S5K_T_CR_II <- CR_T$curvepoints; 
S5K_T_PCR_IIa <- PCR_Ta$curvepoints; 
S5K_T_PCR_IIb <- PCR_Tb$curvepoints;

S5K_T_BP_a <- R_T$bifpoints;
S5K_T_BP_b <- CR_T$bifpoints[5,];
S5K_T_BPE_a <- CR_T$bifpoints[1,]; 
S5K_T_BPE_b <- CR_T$bifpoints[2,];
S5K_T_BPE_c <- CR_T$bifpoints[3,]; 
S5K_T_BPE_d <- CR_T$bifpoints[4,];
S5K_T_BPE_e <- PCR_Ta$bifpoints[1,]; 
S5K_T_BPE_f <- PCR_Tb$bifpoints[1,];
# S5L: TPC+TSR (Linf) + TSR(LV) ====
R_T <- PSPMequi(modelname, "EQ", c(1, 1E-5), 0.1, c(0, 1, 30), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_T <- PSPMequi(modelname, "EQ", R_T$bifpoints[,c(1,2,5)], 0.4, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);

PCR_Ta <- PSPMequi(modelname, "EQ", CR_T$bifpoints[1,c(1:3, 7, 5)], 0.4, c(0, 1, 30), NULL, NULL, clean = TRUE, force=TRUE);
PCR_Tb <- PSPMequi(modelname, "EQ", CR_T$bifpoints[2,c(1:3, 7, 5)], 0.4, c(0, 1, 30), NULL, NULL, clean = TRUE, force=TRUE);#unstable branch

S5L_T_R_II <- R_T$curvepoints;
S5L_T_CR_II <- CR_T$curvepoints; 
S5L_T_PCR_IIa <- PCR_Ta$curvepoints;
S5L_T_PCR_IIb <- PCR_Tb$curvepoints;

S5L_T_BP_a <- R_T$bifpoints;
S5L_T_BP_b <- CR_T$bifpoints[5,];
S5L_T_BPE_a <- CR_T$bifpoints[1,];
S5L_T_BPE_b <- CR_T$bifpoints[2,];
S5L_T_BPE_c <- CR_T$bifpoints[3,];
S5L_T_BPE_d <- CR_T$bifpoints[4,];
S5L_T_BPE_e <- PCR_Ta$bifpoints[1,];
S5L_T_BPE_f <- PCR_Tb$bifpoints[1,];
# S5M: TPC+TSR (Lmat + Linf) + TSR(LV) ====
R_T <- PSPMequi(modelname, "EQ", c(1, 1E-5), 0.1, c(0, 1, 30), NULL, options = c("popZE", "0", "envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);
CR_T <- PSPMequi(modelname, "EQ", R_T$bifpoints[,c(1,2,5)], 0.4, c(0, 1, 30), NULL, options = c("envZE", "1", "envZE", "2"), clean = TRUE, force=TRUE);

PCR_Ta <- PSPMequi(modelname, "EQ", CR_T$bifpoints[1,c(1:3, 7, 5)], 0.4, c(0, 1, 30), NULL, NULL, clean = TRUE, force=TRUE);
PCR_Tb <- PSPMequi(modelname, "EQ", CR_T$bifpoints[2,c(1:3, 7, 5)], 0.4, c(0, 1, 30), NULL, NULL, clean = TRUE, force=TRUE);#unstable branch

S5M_T_R_II <- R_T$curvepoints;
S5M_T_CR_II <- CR_T$curvepoints; 
S5M_T_PCR_IIa <- PCR_Ta$curvepoints;
S5M_T_PCR_IIb <- PCR_Tb$curvepoints;

S5M_T_BP_a <- R_T$bifpoints;
S5M_T_BP_b <- CR_T$bifpoints[5,];
S5M_T_BPE_a <- CR_T$bifpoints[1,];
S5M_T_BPE_b <- CR_T$bifpoints[2,];
S5M_T_BPE_c <- CR_T$bifpoints[3,];
S5M_T_BPE_d <- CR_T$bifpoints[4,];
S5M_T_BPE_e <- PCR_Ta$bifpoints[1,];
S5M_T_BPE_f <- PCR_Tb$bifpoints[1,];
