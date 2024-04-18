####################################################################################################
#                                 Figures
#
# Tri-trophic chain: Resource-Roach-Perch
# along gradients of temperature and resource productivity
#
# The model account for:
#   Temperature-size rule (TSR) & Temperature-dependent performance curve (TPC):
#   TSR in Consumer lengths (maturation size 'Lj' and/or asymptotic size 'Lm') & Predator foraging size range 'Lv'
#   TPC in consumer birth, growth and ingestion rates, and in Predator functional response and metabolic rate
####################################################################################################

load("PSPM_TPC-TSR_PSPMequi_data.RData") #load the generated data from the simulations
########## Fig.1 - Temperature-dependent functions ====
## Parameters ====
a  <- 9E-6;#  # Length weight allometric coefficient
b  <- 3;#     # Length weight allometric exponent
Lb <- 7.0;#   # Length at birth
Lv <- 27.0;#  # Length threshold of exposition to predation
Lj <- 110;#   # Length at maturation
Lm <- 300;#   # Length asymptotic
IR <- 1E-4;#  # Maximum ingestion rate
GR <- 0.006;# # Somatic growth rate
BR <- 0.003;# # Birth rate
Tlr <- 5;#    # Lower viable temperature
Tor <- 20;#   # Optimal temperature (average between lower and upper temperature)
Tur <- 25;#   # Upper viable temperature
Ap  <- 5000;# # Predator attack rate
Hp  <- 0.1;#  # Predator handling time
Eps <- 0.5;#  # Predator assimilation efficiency
Muo <- 0.005;## Scaling coeffocoamt
Mub <- 0.005;## Background mortality rate for consumer species
Mup <- 0.005;## Predator mortality rate
Eax <- 0.55;#   # Activation energy (eV)
TN  <- 293.15;# # Normalization temperature at 20?C (in K)
T0K <- 273.15;# # Conversion factor for Temperature from degC to K (=0 degC)
Boltz = 8.617*10^(-5);#  # Boltzmann constant (eV/K)
Tsr_lv = -0.05;## Slope of size reduction over temperature per degree celsius for Lv
Tsr_lj = -0.05;## Slope of size reduction over temperature per degree celsius for Lj
Tsr_lm = -0.05;## Slope of size reduction over temperature per degree celsius for Lm

Cv <- 1E-5; # Biomass of juvenile consumer vulnerable for predator

TT <- seq(0,30, by=0.1);# Temperature gradient
Tseq <- seq(0,30, by=5);

LinfT=round(Lm*exp((Tsr_lm/b)*(TT-20)));
LjT=round(Lj*exp((Tsr_lj/b)*(TT-20)));
LvT=round(Lv*exp((Tsr_lv/3)*(TT-20)));

FR <- (Ap*Cv/(1+Ap*Hp*Cv)); # Predator functional response
FRT <- FR*(((TT-Tlr)*(TT-Tur))/((TT-Tlr)*(TT-Tur) - (TT-Tor)^2) );
FRT[252:301] <- NA;
FRT[which(FRT<0)] <- NA
XT <- Mup + Muo*exp(-Eax*(TN-(TT+T0K))/(Boltz*TN*(TT+T0K)))

## Figure ====
par(mfrow=c(1,3)); #To save in landscape, 6*15

par(mfg=c(1,1), mar = c(3, 3, 3, 1.5));
plot("", xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(Lb, 450))

## Temperature-Size Rule function
par(mfg=c(1,2), mar = c(3, 1.5, 3, 1.5));
plot("", xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(Lb, 450));
axis(1, at=Tseq, label=T, las=1, cex.axis=1.5);
axis(2, at=c(Lb, max(LvT), max(LjT), max(LinfT)), label=F, las=2);
lines(LinfT~TT, lwd=2);
lines(LjT~TT, lty=2, lwd=2);
lines(LvT~TT, lty=3, lwd=2);
segments(20, Lb, 20, 300, lwd=1, lty=3);

## Temperature-dependent Performance Curves using predator functional response
par(mfg=c(1,3), mar = c(3, 1.5, 3, 3));

plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(0,30), ylim=c(0,0.06));
axis(1, at=Tseq, label=T, las=1, cex.axis=1.5);
axis(2, at=max(FRT, na.rm=T), label=F, las=1, cex.axis=1.5);
lines(FRT~TT, lwd=2);
lines(XT~TT, lwd=2);
segments(20, 0, 20, max(FRT, na.rm=T), lwd=1, lty=3); 

########## Fig.2 - TPC * TSR combinations on trophic chain structure ====
# a-c) Temperature-dependent functions in Consumer only: TPC, TSR, TPC+TSR
# d-f) Temperature-dependent functions in Predator only: TPC, TSR, TPC+TSR
# g-i) Temperature-dependent functions in Consumer & Predator: TPC, TSR, TPC+TSR
##
X_lim = c(0,30);
Y_lim = c(log10(1.7E-6), log10(5E-3));
X_seq = seq(0,30, by=5);
Y_seq = log10(c(1E-5, 1E-4, 1E-3));
Xexpr = c(0, NA, 10, NA, 20, NA, 30);
Yexpr = c(expression(10^{-5}), expression(10^{-4}), expression(10^{-3}));

# Function to generate the graphics
Graph <- function(BP, BPE, LP, XX, YY, ZZ, TCOL){
  BPcurve <- BP; BPEcurve <- BPE; LPcurve <- LP; X_lab <- XX; Y_lab <- YY; Z_lab <- ZZ;
  COL <- TCOL;
  plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=X_lim, ylim=Y_lim);
  if(X_lab==1){axis(1, at=X_seq, label=Xexpr, las=1, cex.axis=1.8)}else{axis(1, at=X_seq, label=F, las=1, cex.axis=1.8)}
  if(Y_lab==1){axis(2, at=Y_seq, label=Yexpr, las=1, cex.axis=1.8)}else{axis(2, at=Y_seq, label=F, las=1, cex.axis=1.8)}
  if(Z_lab==1){axis(4, at=Y_seq, label=Yexpr, las=1, cex.axis=1.8)}else{axis(4, at=Y_seq, label=F, las=1, cex.axis=1.8)}
  lines(log10(BPcurve[,1])~BPcurve[,COL], type="l", lwd=2, lty=3);
  lines(log10(BPEcurve[,1])~BPEcurve[,6], type="l", lwd=2, lty=1);
  lines(log10(LPcurve[,1])~LPcurve[,6], type="l", lwd=2, lty=2);
}

## Figure ====
par(mfrow=c(3,3)); #To save in portrait, 12*12 or 9*9

# S1E) TPC in all Consumer rates
par(mfg=c(1,1), mar = c(1.5, 3, 3, 2.5));
Graph(BP=S1E_KT_BP, BPE=S1E_KT_BPE, LP=S1E_KT_LP, XX=0, YY=0, ZZ=1, TCOL=6);
# S2D: TSR in lj & Lm
par(mfg=c(2,1), mar = c(1.5, 3, 1.5, 2.5));
Graph(BP=S2D_KT_BP, BPE=S2D_KT_BPE, LP=S2D_KT_LP, XX=0, YY=0, ZZ=1, TCOL=6);
# S3A: TPC(C) + TSR(C)
par(mfg=c(3,1), mar = c(3, 3, 1.5, 2.5));
Graph(BP=S3A_KT_BP, BPE=S3A_KT_BPE, LP=S3A_KT_LP, XX=1, YY=0, ZZ=1, TCOL=6);

# S1H) TPC in all Predator rates ****
par(mfg=c(1,2), mar = c(1.5, 2.5, 3, 2.5));
Graph(BP=S1H_KT_BP, BPE=S1H_KT_BPE, LP=S1H_KT_LP, XX=0, YY=0, ZZ=1, TCOL=4);
# S2A: TSR in Predator mxm foraging size threshold Lv
par(mfg=c(2,2), mar = c(1.5, 2.5, 1.5, 2.5));
Graph(BP=S2A_KT_BP, BPE=S2A_KT_BPE, LP=S2A_KT_LP, XX=0, YY=0, ZZ=1, TCOL=6);
# S3E: TPC(P) + TSR(P)
par(mfg=c(3,2), mar = c(3, 2.5, 1.5, 2.5));
Graph(BP=S3E_KT_BP, BPE=S3E_KT_BPE, LP=S3E_KT_LP, XX=1, YY=0, ZZ=1, TCOL=6);

# S1I) TPC in all Consumer & Predator rates *****
par(mfg=c(1,3), mar = c(1.5, 2.5, 3, 3));
Graph(BP=S1I_KT_BP, BPE=S1I_KT_BPE, LP=S1I_KT_LP, XX=0, YY=0, ZZ=0, TCOL=6);
# S2G: TSR in Lv, Lj & Lm
par(mfg=c(2,3), mar = c(1.5, 2.5, 1.5, 3));
Graph(BP=S2G_KT_BP, BPE=S2G_KT_BPE, LP=S2G_KT_LP, XX=0, YY=0, ZZ=0, TCOL=6);
# S3I: TPC(C+P) + TSR(C+P)
par(mfg=c(3,3), mar = c(3, 2.5, 1.5, 3));
Graph(BP=S3I_KT_BP, BPE=S3I_KT_BPE, LP=S3I_KT_LP, XX=1, YY=0, ZZ=0, TCOL=6);

########## Fig.3 - Species Thermal mismatch between species ====
Tvec <- c(rep("13", 3),rep("20",3));
Nich <- c(rep("belowNiche",3),rep("Niche",3),rep("aboveNiche",3));
par(mfg=c(1,1), mar = c(3, 1.5, 3, 2));
Cold <- c(2, 20, 8)
Warm <- c(10, 20, 0)
barplot(as.matrix( cbind.data.frame( Cold, Warm)), col=c("grey", "#67a9cf", "grey", "grey", "#d8b365", "grey"), horiz=T, axes=F)
# -15; -12; 8 # 
# -15; -5; 35 # 0-15-35 
# Focal species is Predator

par(mfrow=c(1,2));
### Focal species: Predator ====
par(mfg=c(1,1), mar = c(2.5, 2.5, 2.5, 2.5));
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(-15,15), ylim=c(log10(1.7E-6), log10(1E-3)));
axis(1, at=c(-10,0,10), label=T, las=1, cex.axis=1.8); axis(1, at=seq(-20,20, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-5, 1E-4)), label=F, las=2, cex.axis=1.8);
axis(2, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=log10(c(1E-5, 1E-4)), label=c(expression(10^-5), expression(10^-4)), las=2, cex.axis=1.8);
axis(3, at=c(-12, -5, 8, 15), label=F, las=2, cex.axis=1.8);
#axis(3, at=c(-5, 15), label=c(15,35), las=2, cex.axis=1.8, col.axis="#d8b365");
#axis(3, at=c(-12,8), label=c(1,21), las=2, cex.axis=1.8, col.axis="#67a9cf");

lines(log10(S4A_TM_BP[,1])~S4A_TM_BP[,6], type="l", lwd=2, lty=3, col="#d8b365");
lines(log10(S4A_TM_BPE[,1])~S4A_TM_BPE[,6], type="l", lwd=2, lty=1, col="#d8b365");
lines(log10(S4A_TM_LP[,1])~S4A_TM_LP[,6], type="l", lwd=2, lty=2, col="#d8b365");

lines(log10(S4C_TM_BP[,1])~S4C_TM_BP[,6], type="l", lwd=2, lty=3, col="#67a9cf");
lines(log10(S4C_TM_BPE[,1])~S4C_TM_BPE[,6], type="l", lwd=2, lty=1, col="#67a9cf");
lines(log10(S4C_TM_LP[,1])~S4C_TM_LP[,6], type="l", lwd=2, lty=2, col="#67a9cf");
abline(v=0, lwd=2, lty=4);

### Focal species: Consumer ====
par(mfg=c(1,2), mar = c(2.5, 2.5, 2.5, 2.5));
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(-15,15), ylim=c(log10(1.7E-6), log10(1E-3)));
axis(1, at=c(-10,0,10), label=T, las=1, cex.axis=1.8); axis(1, at=seq(-20,20, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(3, at=c(-12, -5, 8, 15), label=F, las=2, cex.axis=1.8);
#axis(3, at=c(-5, 15), label=c(15,35), las=2, cex.axis=1.8, col.axis="#d8b365");
#axis(3, at=c(-12,8), label=c(1,21), las=2, cex.axis=1.8, col.axis="#67a9cf");

abline(h=log10(S4B_TM_BP[1]), type="l", lwd=2, lty=3, col="#d8b365");
lines(log10(S4B_TM_BPE[,1])~S4B_TM_BPE[,6], type="l", lwd=2, lty=1, col="#d8b365");
lines(log10(S4B_TM_LP[,1])~S4B_TM_LP[,6], type="l", lwd=2, lty=2, col="#d8b365");

abline(h=log10(S4D_TM_BP[1]), type="l", lwd=2, lty=3, col="#67a9cf");
lines(log10(S4D_TM_BPE[,1])~S4D_TM_BPE[,6], type="l", lwd=2, lty=1, col="#67a9cf");
lines(log10(S4D_TM_LP[,1])~S4D_TM_LP[,6], type="l", lwd=2, lty=2, col="#67a9cf");
abline(v=0, lwd=2, lty=4);


### T=13degC, removed from the figure ====
par(mfg=c(2,1), mar = c(4, 4, 1, 2.5));
plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(-15,15), ylim=c(log10(1.7E-6), log10(1E-3)));
axis(1, at=c(-10,0,10), label=T, las=1, cex.axis=1.8); axis(1, at=seq(-20,20, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=log10(c(1E-5, 1E-4)), label=F, las=2, cex.axis=1.8);
axis(2, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=log10(c(1E-5, 1E-4)), label=c(expression(10^-5), expression(10^-4)), las=2, cex.axis=1.8);

axis(3, at=c(-12,8), label=F, las=2, cex.axis=1.8, col.axis="#67a9cf");

lines(log10(S4C_TM_BP[,1])~S4C_TM_BP[,6], type="l", lwd=2, lty=3, col="#67a9cf");
lines(log10(S4C_TM_BPE[,1])~S4C_TM_BPE[,6], type="l", lwd=2, lty=1, col="#67a9cf");
lines(log10(S4C_TM_LP[,1])~S4C_TM_LP[,6], type="l", lwd=2, lty=2, col="#67a9cf");
abline(v=0, lwd=2, lty=4);

abline(h=log10(S4D_TM_BP[1]), type="l", lwd=2, lty=3, col="darkblue");
lines(log10(S4D_TM_BPE[,1])~S4D_TM_BPE[,6], type="l", lwd=2, lty=1, col="darkblue");
lines(log10(S4D_TM_LP[,1])~S4D_TM_LP[,6], type="l", lwd=2, lty=2, col="darkblue");

### T=20degC, removed from the figure ====
par(mfg=c(2,2), mar = c(4, 2.5, 1, 4));

plot("", type="l", lwd=2, xlab="", ylab="", xaxs="i",yaxs="i", xaxt="n", yaxt="n", xlim=c(-15,15), ylim=c(log10(1.7E-6), log10(1E-3)));
axis(1, at=c(-10,0,10), label=T, las=1, cex.axis=1.8); axis(1, at=seq(-20,20, by=5), label=F, las=1, cex.axis=1.8);
axis(2, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(3, at=c(-5, 15), label=F, las=2, cex.axis=1.8, col.axis="#d8b365");

lines(log10(S4A_TM_BP[,1])~S4A_TM_BP[,6], type="l", lwd=2, lty=3, col="#d8b365");
lines(log10(S4A_TM_BPE[,1])~S4A_TM_BPE[,6], type="l", lwd=2, lty=1, col="#d8b365");
lines(log10(S4A_TM_LP[,1])~S4A_TM_LP[,6], type="l", lwd=2, lty=2, col="#d8b365");

abline(h=log10(S4B_TM_BP[1]), type="l", lwd=2, lty=3, col="darkorange");
lines(log10(S4B_TM_BPE[,1])~S4B_TM_BPE[,6], type="l", lwd=2, lty=1, col="darkorange");
lines(log10(S4B_TM_LP[,1])~S4B_TM_LP[,6], type="l", lwd=2, lty=2, col="darkorange");
abline(v=0, lwd=2, lty=4);
########## Fig. S1 - Bifurcation plots ~ Productivity Gradient ====
SEQ3<-c(1E-3, 2E-3, 3E-3, 4E-3, 5E-3, 6E-3, 7E-3, 8E-3, 9E-3)
SEQL3<-log10(SEQ3+1)
SEQ4<-c(1E-4, 2E-4, 3E-4, 4E-4, 5E-4, 6E-4, 7E-4, 8E-4, 9E-4)
SEQL4<-log10(SEQ4+1)
SEQ5<-c(1E-5, 2E-5, 3E-5, 4E-5, 5E-5, 6E-5, 7E-5, 8E-5, 9E-5)
SEQL5<-log10(SEQ5+1)
SEQ6<-c(1E-6, 2E-6, 3E-6, 4E-6, 5E-6, 6E-6, 7E-6, 8E-6, 9E-6)
SEQL6<-log10(SEQ6+1)
SEQ7<-c(1E-7, 2E-7, 3E-7, 4E-7, 5E-7, 6E-7, 7E-7, 8E-7, 9E-7)
SEQL7<-log10(SEQ7+1)
SEQ8<-c(1E-8, 2E-8, 3E-8, 4E-8, 5E-8, 6E-8, 7E-8, 8E-8, 9E-8)
SEQL8<-log10(SEQ8+1)

graphics.off()
plot.new()
par(mfrow = c(4, 2))
par(cex = 1.0, cex.lab=1.3)
par(oma = c(2.5, 1.5, 1.5, 1.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))
##### T = 20 ====
## R
par(mfg=c(4,1), mar = c(0, 4, 0, 0.5))

plot(1, 1, type="l", xlim=c(log10(1E-6), log10(1E-3)), xlab="", xaxt="n", xaxs="i", ylim=c(log10(1E-6), log10(1E-3)), ylab="", yaxt="n", yaxs="i");
mtext(~R, 2, line=1, cex=1.5);
axis(1, at=seq(-6,-3), labels=c(expression(10^{-6}), expression(10^{-5}),expression(10^{-4}), expression(10^{-3})))
axis(1, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

lines(log10(S3I_K_R_20[1:which(S3I_K_R_20[,1] == S3I_K_BP_20[1]),1]), log10(S3I_K_R_20[1:which(S3I_K_R_20[,1] == S3I_K_BP_20[1]),2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(log10(S3I_K_CR_20[1:which(S3I_K_CR_20[,1] == S3I_K_BPE_20[1]),1]), log10(S3I_K_CR_20[1:which(S3I_K_CR_20[,1] == S3I_K_BPE_20[1]),2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(log10(S3I_K_PCR_20[1:which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]),1]), log10(S3I_K_PCR_20[1:which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]),2]), type="l", col=rgb(0,0,0), lwd=2, lty=2);
lines(log10(S3I_K_PCR_20[which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]):dim(S3I_K_PCR_20)[1],1]), log10(S3I_K_PCR_20[which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]):dim(S3I_K_PCR_20)[1],2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);

points(log10(S3I_K_BP_20[1]), log10(S3I_K_BP_20[2]), col=rgb(0,0,.9), pch=8, lwd=2);
points(log10(S3I_K_BPE_20[1]), log10(S3I_K_BPE_20[2]), col=rgb(0,0,0), pch=8, lwd=2);
points(log10(S3I_K_LP_20[1]), log10(S3I_K_LP_20[2]), col=rgb(.9,0,0), pch=8, lwd=2);
abline(v=log10(S3I_K_BP_20[1]), lty=3); abline(v=log10(S3I_K_BPE_20[1]), lty=3); abline(v=log10(S3I_K_LP_20[1]), lty=3);

## C[juveniles]
par(mfg=c(3,1), mar = c(0, 4, 0, 0.5))

plot(1, 1, type="l", xlim=c(log10(1E-6), log10(1E-3)), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-6), log10(9E-4)), ylab="", yaxt="n", yaxs="i");
mtext(~C[juveniles], 2, line=1, cex=1.5);
axis(4, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

S3I_K_CR_20[1,7:9] <- 1E-9; #Artificially replacing the initial null biomass of consumer by a very low biomass to have a continuous line from the invasion threshold
lines(log10(S3I_K_CR_20[1:which(S3I_K_CR_20[,1] == S3I_K_BPE_20[1]),1]), log10(S3I_K_CR_20[1:which(S3I_K_CR_20[,1] == S3I_K_BPE_20[1]),8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(log10(S3I_K_PCR_20[1:which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]),1]), log10(S3I_K_PCR_20[1:which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]),8]), type="l", col=rgb(0,0,0), lwd=2, lty=2);
lines(log10(S3I_K_PCR_20[which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]):dim(S3I_K_PCR_20)[1],1]), log10(S3I_K_PCR_20[which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]):dim(S3I_K_PCR_20)[1],8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);

points(log10(S3I_K_BPE_20[1]), log10(S3I_K_BPE_20[8]), col=rgb(0,0,0), pch=8, lwd=2);
points(log10(S3I_K_LP_20[1]), log10(S3I_K_LP_20[8]), col=rgb(.9,0,0), pch=8, lwd=2);
abline(v=log10(S3I_K_BP_20[1]), lty=3); abline(v=log10(S3I_K_BPE_20[1]), lty=3); abline(v=log10(S3I_K_LP_20[1]), lty=3);

## C[adults]
par(mfg=c(2,1), mar = c(0, 4, 0, 0.5))

plot(1, 1, type="l", xlim=c(log10(1E-6), log10(1E-3)), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-8), log10(2E-4)), ylab="", yaxt="n", yaxs="i");
mtext(~C[adults], 2, line=1, cex=1.5);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

lines(log10(S3I_K_CR_20[1:which(S3I_K_CR_20[,1] == S3I_K_BPE_20[1]),1]), log10(S3I_K_CR_20[1:which(S3I_K_CR_20[,1] == S3I_K_BPE_20[1]),9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(log10(S3I_K_PCR_20[1:which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]),1]), log10(S3I_K_PCR_20[1:which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]),9]), type="l", col=rgb(0,0,0), lwd=2, lty=2);
lines(log10(S3I_K_PCR_20[which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]):dim(S3I_K_PCR_20)[1],1]), log10(S3I_K_PCR_20[which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]):dim(S3I_K_PCR_20)[1],9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);

points(log10(S3I_K_BPE_20[1]), log10(S3I_K_BPE_20[9]), col=rgb(0,0,0), pch=8, lwd=2);
points(log10(S3I_K_LP_20[1]), log10(S3I_K_LP_20[9]), col=rgb(.9,0,0), pch=8, lwd=2);
abline(v=log10(S3I_K_BP_20[1]), lty=3); abline(v=log10(S3I_K_BPE_20[1]), lty=3); abline(v=log10(S3I_K_LP_20[1]), lty=3);

## P
par(mfg=c(1,1), mar = c(0, 4, 0, 0.5))

plot(1, 1, type="l", xlim=c(log10(1E-6), log10(1E-3)), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-6), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
mtext("P", 2, line=1, cex=1.8);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

lines(log10(S3I_K_PCR_20[1:which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]),1]), log10(S3I_K_PCR_20[1:which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]),3]), type="l", col=rgb(0,0,0), lwd=2, lty=2);
lines(log10(S3I_K_PCR_20[which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]):dim(S3I_K_PCR_20)[1],1]), log10(S3I_K_PCR_20[which(S3I_K_PCR_20[,1] == S3I_K_LP_20[1]):dim(S3I_K_PCR_20)[1],3]), type="l", col=rgb(0,0,0), lwd=2, lty=1);

points(log10(S3I_K_LP_20[1]), log10(S3I_K_LP_20[3]), col=rgb(.9,0,0), pch=8, lwd=2);
abline(v=log10(S3I_K_BP_20[1]), lty=3); abline(v=log10(S3I_K_BPE_20[1]), lty=3); abline(v=log10(S3I_K_LP_20[1]), lty=3);

##### T = 13 ====
## R
par(mfg=c(4,2), mar = c(0, 2, 0, 0.5))

plot(1, 1, type="l", xlim=c(log10(1E-6), log10(1E-3)), xlab="", xaxt="n", xaxs="i", ylim=c(log10(1E-6), log10(1E-3)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(-6,-3), labels=c(expression(10^{-6}), expression(10^{-5}),expression(10^{-4}), expression(10^{-3})))
axis(2, at=seq(-6,-3), labels=c(expression(10^{-6}), expression(10^{-5}),expression(10^{-4}), expression(10^{-3})), las=2)
axis(1, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(2, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

lines(log10(S3I_K_R_13[1:which(S3I_K_R_13[,1] == S3I_K_BP_13[1]),1]), log10(S3I_K_R_13[1:which(S3I_K_R_13[,1] == S3I_K_BP_13[1]),2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(log10(S3I_K_CR_13[1:which(S3I_K_CR_13[,1] == S3I_K_BPE_13[1]),1]), log10(S3I_K_CR_13[1:which(S3I_K_CR_13[,1] == S3I_K_BPE_13[1]),2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(log10(S3I_K_PCR_13[1:which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]),1]), log10(S3I_K_PCR_13[1:which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]),2]), type="l", col=rgb(0,0,0), lwd=2, lty=2);
lines(log10(S3I_K_PCR_13[which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]):dim(S3I_K_PCR_13)[1],1]), log10(S3I_K_PCR_13[which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]):dim(S3I_K_PCR_13)[1],2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);

points(log10(S3I_K_BP_13[1]), log10(S3I_K_BP_13[2]), col=rgb(0,0,.9), pch=8, lwd=2);
points(log10(S3I_K_BPE_13[1]), log10(S3I_K_BPE_13[2]), col=rgb(0,0,0), pch=8, lwd=2);
points(log10(S3I_K_LP_13[1]), log10(S3I_K_LP_13[2]), col=rgb(.9,0,0), pch=8, lwd=2);
abline(v=log10(S3I_K_BP_13[1]), lty=3); abline(v=log10(S3I_K_BPE_13[1]), lty=3); abline(v=log10(S3I_K_LP_13[1]), lty=3);

## C[juveniles]
par(mfg=c(3,2), mar = c(0, 2, 0, 0.5))

plot(1, 1, type="l", xlim=c(log10(1E-6), log10(1E-3)), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-6), log10(9E-4)), ylab="", yaxt="n", yaxs="i");
axis(2, at=seq(-6,-3), labels=c(expression(10^{-6}), expression(10^{-5}),expression(10^{-4}), expression(10^{-3})), las=2);
axis(2, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

S3I_K_CR_13[1,7:9] <- 1E-9; #Artificially replacing the initial null biomass of consumer by a very low biomass to have a continuous line from the invasion threshold
lines(log10(S3I_K_CR_13[1:which(S3I_K_CR_13[,1] == S3I_K_BPE_13[1]),1]), log10(S3I_K_CR_13[1:which(S3I_K_CR_13[,1] == S3I_K_BPE_13[1]),8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(log10(S3I_K_PCR_13[1:which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]),1]), log10(S3I_K_PCR_13[1:which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]),8]), type="l", col=rgb(0,0,0), lwd=2, lty=2);
lines(log10(S3I_K_PCR_13[which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]):dim(S3I_K_PCR_13)[1],1]), log10(S3I_K_PCR_13[which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]):dim(S3I_K_PCR_13)[1],8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);

points(log10(S3I_K_BPE_13[1]), log10(S3I_K_BPE_13[8]), col=rgb(0,0,0), pch=8, lwd=2);
points(log10(S3I_K_LP_13[1]), log10(S3I_K_LP_13[8]), col=rgb(.9,0,0), pch=8, lwd=2);
abline(v=log10(S3I_K_BP_13[1]), lty=3); abline(v=log10(S3I_K_BPE_13[1]), lty=3); abline(v=log10(S3I_K_LP_13[1]), lty=3);

## C[adults]
par(mfg=c(2,2), mar = c(0, 2, 0, 0.5))

plot(1, 1, type="l", xlim=c(log10(1E-6), log10(1E-3)), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-8), log10(2E-4)), ylab="", yaxt="n", yaxs="i");
axis(2, at=seq(-7,-3), labels=c(expression(10^{-7}), expression(10^{-6}), expression(10^{-5}),expression(10^{-4}), expression(10^{-3})), las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

lines(log10(S3I_K_CR_13[1:which(S3I_K_CR_13[,1] == S3I_K_BPE_13[1]),1]), log10(S3I_K_CR_13[1:which(S3I_K_CR_13[,1] == S3I_K_BPE_13[1]),9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(log10(S3I_K_PCR_13[1:which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]),1]), log10(S3I_K_PCR_13[1:which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]),9]), type="l", col=rgb(0,0,0), lwd=2, lty=2);
lines(log10(S3I_K_PCR_13[which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]):dim(S3I_K_PCR_13)[1],1]), log10(S3I_K_PCR_13[which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]):dim(S3I_K_PCR_13)[1],9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);

points(log10(S3I_K_BPE_13[1]), log10(S3I_K_BPE_13[9]), col=rgb(0,0,0), pch=8, lwd=2);
points(log10(S3I_K_LP_13[1]), log10(S3I_K_LP_13[9]), col=rgb(.9,0,0), pch=8, lwd=2);
abline(v=log10(S3I_K_BP_13[1]), lty=3); abline(v=log10(S3I_K_BPE_13[1]), lty=3); abline(v=log10(S3I_K_LP_13[1]), lty=3);

## P
par(mfg=c(1,2), mar = c(0, 2, 0, 0.5))

plot(1, 1, type="l", xlim=c(log10(1E-6), log10(1E-3)), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-6), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(2, at=seq(-7,-3), labels=c(expression(10^{-7}), expression(10^{-6}), expression(10^{-5}),expression(10^{-4}), expression(10^{-3})), las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

lines(log10(S3I_K_PCR_13[1:which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]),1]), log10(S3I_K_PCR_13[1:which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]),3]), type="l", col=rgb(0,0,0), lwd=2, lty=2);
lines(log10(S3I_K_PCR_13[which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]):dim(S3I_K_PCR_13)[1],1]), log10(S3I_K_PCR_13[which(S3I_K_PCR_13[,1] == S3I_K_LP_13[1]):dim(S3I_K_PCR_13)[1],3]), type="l", col=rgb(0,0,0), lwd=2, lty=1);

points(log10(S3I_K_LP_13[1]), log10(S3I_K_LP_13[3]), col=rgb(.9,0,0), pch=8, lwd=2);
abline(v=log10(S3I_K_BP_13[1]), lty=3); abline(v=log10(S3I_K_BPE_13[1]), lty=3); abline(v=log10(S3I_K_LP_13[1]), lty=3);

########## Fig. S2 - Bifurcation plots ~ Temperature Gradient ====
SEQ3<-c(1E-3, 2E-3, 3E-3, 4E-3, 5E-3, 6E-3, 7E-3, 8E-3, 9E-3)
SEQL3<-log10(SEQ3+1)
SEQ4<-c(1E-4, 2E-4, 3E-4, 4E-4, 5E-4, 6E-4, 7E-4, 8E-4, 9E-4)
SEQL4<-log10(SEQ4+1)
SEQ5<-c(1E-5, 2E-5, 3E-5, 4E-5, 5E-5, 6E-5, 7E-5, 8E-5, 9E-5)
SEQL5<-log10(SEQ5+1)
SEQ6<-c(1E-6, 2E-6, 3E-6, 4E-6, 5E-6, 6E-6, 7E-6, 8E-6, 9E-6)
SEQL6<-log10(SEQ6+1)
SEQ7<-c(1E-7, 2E-7, 3E-7, 4E-7, 5E-7, 6E-7, 7E-7, 8E-7, 9E-7)
SEQL7<-log10(SEQ7+1)
SEQ8<-c(1E-8, 2E-8, 3E-8, 4E-8, 5E-8, 6E-8, 7E-8, 8E-8, 9E-8)
SEQL8<-log10(SEQ8+1)


graphics.off()
plot.new()
par(mfrow = c(4, 3))
par(cex = 1.0, cex.lab=1.3)
par(oma = c(2.5, 1.5, 1.5, 1.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))
##### S3A ====
## R
par(mfg=c(4,1), mar = c(0, 4, 0, 0.5))

plot(1, 1, type="l", xlim=c(0, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-6), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
mtext(~R, 2, line=1, cex=1.5);
axis(1, at=seq(0,30, by=5), labels=T)
axis(4, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

segments(0, log10(S3A_T_R[1,2]), S3A_T_BP_a[1], log10(S3A_T_R[1,2]), col=rgb(0,0,0), lwd=2);
segments(S3A_T_BP_b[1], log10(S3A_T_R[1,2]), 30, log10(S3A_T_R[1,2]), col=rgb(0,0,0), lwd=2);
lines(S3A_T_CR[1:which(S3A_T_CR[,1]==S3A_T_BPE_a[1]),1], log10(S3A_T_CR[1:which(S3A_T_CR[,1]==S3A_T_BPE_a[1]),2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3A_T_CR[271:dim(S3A_T_CR)[1],1], log10(S3A_T_CR[271:dim(S3A_T_CR)[1],2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3A_T_CR[,1], log10(S3A_T_CR[,2]), type="l", col=rgb(0,0,0), lwd=1, lty=1);
lines(S3A_T_PCR_a[,1], log10(S3A_T_PCR_a[,2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3A_T_PCR_b[,1], log10(S3A_T_PCR_b[,2]), type="l", col=rgb(0,0,0), lwd=2, lty=2);

points(S3A_T_BP_a[1], log10(S3A_T_BP_a[2]), col=rgb(0,0,.9), pch=8, lwd=2);
points(S3A_T_BP_b[1], log10(S3A_T_BP_b[2]), col=rgb(0,0,.9), pch=8, lwd=2);
points(S3A_T_BPE_a[1], log10(S3A_T_BPE_a[2]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3A_T_BPE_b[1], log10(S3A_T_BPE_b[2]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3A_T_BPE_c[1], log10(S3A_T_BPE_c[2]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3A_T_BPE_d[1], log10(S3A_T_BPE_d[2]), col=rgb(0,0,0), pch=8, lwd=2);
abline(v=S3A_T_BP_a[1], lty=3);
abline(v=S3A_T_BP_b[1], lty=3)
abline(v=S3A_T_BPE_b[1], lty=3)
abline(v=S3A_T_BPE_c[1], lty=3);

## C[juveniles] 
par(mfg=c(3,1), mar = c(0, 4, 0, 0.5))

plot(1, 1, type="l", xlim=c(0, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-6), log10(9E-4)), ylab="", yaxt="n", yaxs="i");
mtext(~C[juveniles], 2, line=1, cex=1.5);
axis(4, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

segments(S3A_T_CR[1,1], log10(1E-6), S3A_T_CR[2,1], log10(S3A_T_CR[2,8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3A_T_CR[1:which(S3A_T_CR[,1]==S3A_T_BPE_a[1]),1], log10(S3A_T_CR[1:which(S3A_T_CR[,1]==S3A_T_BPE_a[1]),8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3A_T_CR[271:dim(S3A_T_CR)[1],1], log10(S3A_T_CR[271:dim(S3A_T_CR)[1],8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3A_T_CR[,1], log10(S3A_T_CR[,8]), type="l", col=rgb(0,0,0), lwd=1, lty=1);
lines(S3A_T_PCR_a[,1], log10(S3A_T_PCR_a[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3A_T_PCR_b[,1], log10(S3A_T_PCR_b[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=2);

points(S3A_T_BPE_a[1], log10(S3A_T_BPE_a[8]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3A_T_BPE_b[1], log10(S3A_T_BPE_b[8]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3A_T_BPE_c[1], log10(S3A_T_BPE_c[8]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3A_T_BPE_d[1], log10(S3A_T_BPE_d[8]), col=rgb(0,0,0), pch=8, lwd=2);
abline(v=S3A_T_BP_a[1], lty=3);
abline(v=S3A_T_BP_b[1], lty=3);
abline(v=S3A_T_BPE_b[1], lty=3);
abline(v=S3A_T_BPE_c[1], lty=3);

## C[adults]
par(mfg=c(2,1), mar = c(0, 4, 0, 0.5))

plot(1, 1, type="l", xlim=c(0, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-7), log10(2E-4)), ylab="", yaxt="n", yaxs="i");
mtext(~C[adults], 2, line=1, cex=1.5);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

segments(S3A_T_CR[1,1], log10(1E-8), S3A_T_CR[2,1], log10(S3A_T_CR[2,9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
segments(S3A_T_CR[dim(S3A_T_CR)[1]-2,1], log10(1E-8), S3A_T_CR[dim(S3A_T_CR)[1]-2,1], log10(S3A_T_CR[dim(S3A_T_CR)[1]-2,9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3A_T_CR[1:which(S3A_T_CR==S3A_T_BPE_a[1]),1], log10(S3A_T_CR[1:which(S3A_T_CR==S3A_T_BPE_a[1]),9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3A_T_CR[271:295,1], log10(S3A_T_CR[271:295,9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3A_T_CR[,1], log10(S3A_T_CR[,9]), type="l", col=rgb(0,0,0), lwd=1, lty=1);
lines(S3A_T_PCR_a[,1], log10(S3A_T_PCR_a[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3A_T_PCR_b[,1], log10(S3A_T_PCR_b[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=2);

points(S3A_T_BPE_a[1], log10(S3A_T_BPE_a[9]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3A_T_BPE_b[1], log10(S3A_T_BPE_b[9]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3A_T_BPE_c[1], log10(S3A_T_BPE_c[9]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3A_T_BPE_d[1], log10(S3A_T_BPE_d[9]), col=rgb(0,0,0), pch=8, lwd=2);
abline(v=S3A_T_BP_a[1], lty=3);
abline(v=S3A_T_BP_b[1], lty=3);
abline(v=S3A_T_BPE_b[1], lty=3);
abline(v=S3A_T_BPE_c[1], lty=3);

## P
par(mfg=c(1,1), mar = c(0, 4, 0, 0.5))

plot(1, 1, type="l", xlim=c(0, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(1E-7), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
mtext("P", 2, line=1, cex=1.8);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

lines(S3A_T_PCR_a[,1], log10(S3A_T_PCR_a[,3]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3A_T_PCR_b[,1], log10(S3A_T_PCR_b[,3]), type="l", col=rgb(0,0,0), lwd=2, lty=2);

abline(v=S3A_T_BP_a[1], lty=3);
abline(v=S3A_T_BP_b[1], lty=3);
abline(v=S3A_T_BPE_b[1], lty=3);
abline(v=S3A_T_BPE_c[1], lty=3);

##### S3E ====
## R
par(mfg=c(4,2), mar = c(0, 2, 0, 0.5))

plot(1, 1, type="l", xlim=c(0, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-6), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30, by=5), labels=T)
axis(2, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(2, at=seq(-5,-4), labels=c(expression(10^{-5}),expression(10^{-4})), las=2)

lines(S3E_T_CR[1:which(S3E_T_CR[,1]==S3E_T_BPE_a[1]),1], log10(S3E_T_CR[1:which(S3E_T_CR[,1]==S3E_T_BPE_a[1]),2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
segments(S3E_T_BPE_b[1], log10(S3E_T_BPE_b[2]), 30, log10(S3E_T_BPE_b[2]), col=rgb(0,0,0), lwd=2, lty=1);
lines(S3E_T_PCR[1:which(S3E_T_PCR[,1]==S3E_T_LP_a[1]),1], log10(S3E_T_PCR[1:which(S3E_T_PCR[,1]==S3E_T_LP_a[1]),2]), type="l", col=rgb(0,0,0), lwd=2, lty=2);#dashed
lines(S3E_T_PCR[which(S3E_T_PCR[,1]==S3E_T_LP_a[1]):521,1], log10(S3E_T_PCR[which(S3E_T_PCR[,1]==S3E_T_LP_a[1]):521,2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3E_T_PCR[521:dim(S3E_T_PCR)[1],1], log10(S3E_T_PCR[521:dim(S3E_T_PCR)[1],2]), type="l", col=rgb(0,0,0), lwd=2, lty=2);#dashed

points(S3E_T_BPE_a[1], log10(S3E_T_BPE_a[2]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3E_T_BPE_b[1], log10(S3E_T_BPE_b[2]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3E_T_LP_a[1], log10(S3E_T_LP_a[2]), col=rgb(.9,0,0), pch=8, lwd=2);
points(S3E_T_LP_b[1], log10(S3E_T_LP_b[2]), col=rgb(.9,0,0), pch=8, lwd=2);
abline(v=S3E_T_BPE_a[1], lty=3); abline(v=S3E_T_BPE_b[1], lty=3);
abline(v=S3E_T_LP_a[1], lty=3); abline(v=S3E_T_LP_b[1], lty=3);

## C[juveniles]
par(mfg=c(3,2), mar = c(0, 2, 0, 0.5))

plot(1, 1, type="l", xlim=c(0, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-6), log10(9E-4)), ylab="", yaxt="n", yaxs="i");
axis(2, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(2, at=seq(-6,-4), labels=c(expression(10^{-6}), expression(10^{-5}),expression(10^{-4})), las=2)
axis(4, at=c(log10(SEQ6), log10(SEQ5),log10(SEQ4)), label=F, las=2);

lines(S3E_T_CR[1:which(S3E_T_CR[,1]==S3E_T_BPE_a[1]),1], log10(S3E_T_CR[1:which(S3E_T_CR[,1]==S3E_T_BPE_a[1]),8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3E_T_CR[508:dim(S3E_T_CR)[1],1], log10(S3E_T_CR[508:dim(S3E_T_CR)[1],8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3E_T_PCR[1:which(S3E_T_PCR[,1]==S3E_T_LP_a[1]),1], log10(S3E_T_PCR[1:which(S3E_T_PCR[,1]==S3E_T_LP_a[1]),8]), type="l", col=rgb(0,0,0), lwd=2, lty=2);#dashed
lines(S3E_T_PCR[which(S3E_T_PCR[,1]==S3E_T_LP_a[1]):521,1], log10(S3E_T_PCR[which(S3E_T_PCR[,1]==S3E_T_LP_a[1]):521,8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3E_T_PCR[521:dim(S3E_T_PCR)[1],1], log10(S3E_T_PCR[521:dim(S3E_T_PCR)[1],8]), type="l", col=rgb(0,0,0), lwd=2, lty=2);#dashed

points(S3E_T_BPE_a[1], log10(S3E_T_BPE_a[8]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3E_T_BPE_b[1], log10(S3E_T_BPE_b[8]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3E_T_LP_a[1], log10(S3E_T_LP_a[8]), col=rgb(.9,0,0), pch=8, lwd=2);
points(S3E_T_LP_b[1], log10(S3E_T_LP_b[8]), col=rgb(.9,0,0), pch=8, lwd=2);
abline(v=S3E_T_BPE_a[1], lty=3); abline(v=S3E_T_BPE_b[1], lty=3);
abline(v=S3E_T_LP_a[1], lty=3); abline(v=S3E_T_LP_b[1], lty=3);

## C[adults]
par(mfg=c(2,2), mar = c(0, 2, 0, 0.5))

plot(1, 1, type="l", xlim=c(0, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-7), log10(2E-4)), ylab="", yaxt="n", yaxs="i");
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(2, at=seq(-7,-4), labels=c(expression(10^{-7}), expression(10^{-6}), expression(10^{-5}),expression(10^{-4})), las=2)
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

lines(S3E_T_CR[1:which(S3E_T_CR[,1]==S3E_T_BPE_a[1]),1], log10(S3E_T_CR[1:which(S3E_T_CR[,1]==S3E_T_BPE_a[1]),9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3E_T_CR[508:dim(S3E_T_CR)[1],1], log10(S3E_T_CR[508:dim(S3E_T_CR)[1],9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3E_T_PCR[1:which(S3E_T_PCR[,1]==S3E_T_LP_a[1]),1], log10(S3E_T_PCR[1:which(S3E_T_PCR[,1]==S3E_T_LP_a[1]),9]), type="l", col=rgb(0,0,0), lwd=2, lty=2);#dashed
lines(S3E_T_PCR[which(S3E_T_PCR[,1]==S3E_T_LP_a[1]):521,1], log10(S3E_T_PCR[which(S3E_T_PCR[,1]==S3E_T_LP_a[1]):521,9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3E_T_PCR[521:dim(S3E_T_PCR)[1],1], log10(S3E_T_PCR[521:dim(S3E_T_PCR)[1],9]), type="l", col=rgb(0,0,0), lwd=2, lty=2);#dashed

points(S3E_T_BPE_a[1], log10(S3E_T_BPE_a[9]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3E_T_BPE_b[1], log10(S3E_T_BPE_b[9]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3E_T_LP_a[1], log10(S3E_T_LP_a[9]), col=rgb(.9,0,0), pch=8, lwd=2);
points(S3E_T_LP_b[1], log10(S3E_T_LP_b[9]), col=rgb(.9,0,0), pch=8, lwd=2);
abline(v=S3E_T_BPE_a[1], lty=3); abline(v=S3E_T_BPE_b[1], lty=3);
abline(v=S3E_T_LP_a[1], lty=3); abline(v=S3E_T_LP_b[1], lty=3);

## P
par(mfg=c(1,2), mar = c(0, 2, 0, 0.5))

plot(1, 1, type="l", xlim=c(0, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(1E-7), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(2, at=seq(-6,-3), labels=c(expression(10^{-6}), expression(10^{-5}),expression(10^{-4}), expression(10^{-3})), las=2)
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);

lines(S3E_T_PCR[1:which(S3E_T_PCR[,1]==S3E_T_LP_a[1]),1], log10(S3E_T_PCR[1:which(S3E_T_PCR[,1]==S3E_T_LP_a[1]),3]), type="l", col=rgb(0,0,0), lwd=2, lty=2);#dashed
lines(S3E_T_PCR[which(S3E_T_PCR[,1]==S3E_T_LP_a[1]):521,1], log10(S3E_T_PCR[which(S3E_T_PCR[,1]==S3E_T_LP_a[1]):521,3]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3E_T_PCR[521:dim(S3E_T_PCR)[1],1], log10(S3E_T_PCR[521:dim(S3E_T_PCR)[1],3]), type="l", col=rgb(0,0,0), lwd=2, lty=2);#dashed

points(S3E_T_LP_a[1], log10(S3E_T_LP_a[3]), col=rgb(.9,0,0), pch=8, lwd=2);
points(S3E_T_LP_b[1], log10(S3E_T_LP_b[3]), col=rgb(.9,0,0), pch=8, lwd=2);
abline(v=S3E_T_BPE_a[1], lty=3); abline(v=S3E_T_BPE_b[1], lty=3);
abline(v=S3E_T_LP_a[1], lty=3); abline(v=S3E_T_LP_b[1], lty=3);

##### S3I ====
## R
par(mfg=c(4,3), mar = c(0, 2, 0, 0.5))

plot(1, 1, type="l", xlim=c(0, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-6), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30, by=5), labels=T)
axis(2, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(2, at=seq(-5,-4), labels=c(expression(10^{-5}),expression(10^{-4})), las=2)

segments(0, log10(S3I_T_BP_a[2]), S3I_T_BP_a[1], log10(S3I_T_BP_a[2]), col=rgb(0,0,0), lwd=2);
segments(S3I_T_BP_b[1], log10(S3I_T_BP_b[2]), 30, log10(S3I_T_BP_b[2]), col=rgb(0,0,0), lwd=2);
lines(S3I_T_CR[21:dim(S3I_T_CR)[1],1], log10(S3I_T_CR[21:dim(S3I_T_CR)[1],2]), type="l", col=rgb(0,0,0), lwd=1, lty=1);
lines(S3I_T_PCR_a[,1], log10(S3I_T_PCR_a[,2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3I_T_PCR_b[,1], log10(S3I_T_PCR_b[,2]), type="l", col=rgb(0,0,0), lwd=2, lty=2);

points(S3I_T_BP_a[1], log10(S3I_T_BP_a[2]), col=rgb(0,0,.9), pch=8, lwd=2);
points(S3I_T_BP_b[1], log10(S3I_T_BP_b[2]), col=rgb(0,0,.9), pch=8, lwd=2);
points(S3I_T_BPE_a[1], log10(S3I_T_BPE_a[2]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3I_T_BPE_b[1], log10(S3I_T_BPE_b[2]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3I_T_BPE_c[1], log10(S3I_T_BPE_c[2]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3I_T_BPE_d[1], log10(S3I_T_BPE_d[2]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3I_T_BPE_e[1], log10(S3I_T_BPE_e[2]), col=rgb(0,0,0), pch=8, lwd=2);
abline(v=S3I_T_BP_a[1], lty=3);
abline(v=S3I_T_BP_b[1], lty=3);
abline(v=S3I_T_BPE_a[1], lty=3)
abline(v=S3I_T_BPE_b[1], lty=3);
abline(v=S3I_T_BPE_c[1], lty=3);
abline(v=S3I_T_BPE_d[1], lty=3);

## C[juveniles] 
par(mfg=c(3,3), mar = c(0, 2, 0, 0.5))

plot(1, 1, type="l", xlim=c(0, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-6), log10(9E-4)), ylab="", yaxt="n", yaxs="i");
axis(2, at=c(log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(2, at=seq(-6,-4), labels=c(expression(10^{-6}), expression(10^{-5}),expression(10^{-4})), las=2)

S3I_T_CR[22,7:9] <- 1E-9; #Artificially replacing the initial null biomass of consumer by a very low biomass to have a continuous line from the invasion threshold
lines(S3I_T_CR[22:dim(S3I_T_CR)[1],1], log10(S3I_T_CR[22:dim(S3I_T_CR)[1],8]), type="l", col=rgb(0,0,0), lwd=1, lty=1);
lines(S3I_T_PCR_a[,1], log10(S3I_T_PCR_a[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3I_T_PCR_b[,1], log10(S3I_T_PCR_b[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=2);

points(S3I_T_BPE_a[1], log10(S3I_T_BPE_a[8]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3I_T_BPE_b[1], log10(S3I_T_BPE_b[8]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3I_T_BPE_c[1], log10(S3I_T_BPE_c[8]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3I_T_BPE_d[1], log10(S3I_T_BPE_d[8]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3I_T_BPE_e[1], log10(S3I_T_BPE_e[8]), col=rgb(0,0,0), pch=8, lwd=2);
abline(v=S3I_T_BP_a[1], lty=3);
abline(v=S3I_T_BP_b[1], lty=3);
abline(v=S3I_T_BPE_a[1], lty=3)
abline(v=S3I_T_BPE_b[1], lty=3);
abline(v=S3I_T_BPE_c[1], lty=3);
abline(v=S3I_T_BPE_d[1], lty=3);

## C[adults]
par(mfg=c(2,3), mar = c(0, 2, 0, 0.5))

plot(1, 1, type="l", xlim=c(0, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-7), log10(2E-4)), ylab="", yaxt="n", yaxs="i");
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(2, at=seq(-7,-4), labels=c(expression(10^{-7}), expression(10^{-6}), expression(10^{-5}),expression(10^{-4})), las=2)

lines(S3I_T_CR[22:dim(S3I_T_CR)[1],1], log10(S3I_T_CR[22:dim(S3I_T_CR)[1],9]), type="l", col=rgb(0,0,0), lwd=1, lty=1);
lines(S3I_T_PCR_a[,1], log10(S3I_T_PCR_a[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3I_T_PCR_b[,1], log10(S3I_T_PCR_b[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=2);

points(S3I_T_BPE_a[1], log10(S3I_T_BPE_a[9]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3I_T_BPE_b[1], log10(S3I_T_BPE_b[9]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3I_T_BPE_c[1], log10(S3I_T_BPE_c[9]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3I_T_BPE_d[1], log10(S3I_T_BPE_d[9]), col=rgb(0,0,0), pch=8, lwd=2);
points(S3I_T_BPE_e[1], log10(S3I_T_BPE_e[9]), col=rgb(0,0,0), pch=8, lwd=2);
abline(v=S3I_T_BP_a[1], lty=3);
abline(v=S3I_T_BP_b[1], lty=3);
abline(v=S3I_T_BPE_a[1], lty=3)
abline(v=S3I_T_BPE_b[1], lty=3);
abline(v=S3I_T_BPE_c[1], lty=3);
abline(v=S3I_T_BPE_d[1], lty=3);

## P
par(mfg=c(1,3), mar = c(0, 2, 0, 0.5))

plot(1, 1, type="l", xlim=c(0, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(1E-7), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(2, at=seq(-6,-3), labels=c(expression(10^{-6}), expression(10^{-5}),expression(10^{-4}), expression(10^{-3})), las=2)

lines(S3I_T_PCR_a[,1], log10(S3I_T_PCR_a[,3]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S3I_T_PCR_b[,1], log10(S3I_T_PCR_b[,3]), type="l", col=rgb(0,0,0), lwd=2, lty=2);

abline(v=S3I_T_BP_a[1], lty=3);
abline(v=S3I_T_BP_b[1], lty=3);
abline(v=S3I_T_BPE_a[1], lty=3)
abline(v=S3I_T_BPE_b[1], lty=3);
abline(v=S3I_T_BPE_c[1], lty=3);
abline(v=S3I_T_BPE_d[1], lty=3);

########## Fig. S3 - Influence of TPC alone on the trophic chain structure when implemented in a single parameter ====
# a-d) consumer traits
# e-f) predator traits
## Figure ====
par(mfrow=c(2,3)); #To save in portrait, 8*12
# S1A: TPC in Consumer ingestion rate
par(mfg=c(1,1), mar = c(1.5, 3, 3, 2.5));
Graph(BP=S1A_KT_BP, BPE=S1A_KT_BPE, LP=S1A_KT_LP, XX=0, YY=0, ZZ=1, TCOL=6);
abline(v=5, lwd=2, lty=4)
abline(v=25, lwd=2, lty=4)
# S1B: TPC in Consumer birth rate
par(mfg=c(1,2), mar = c(1.5, 2.5, 3, 2.5));
Graph(BP=S1B_KT_BP, BPE=S1B_KT_BPE, LP=S1B_KT_LP, XX=0, YY=0, ZZ=1, TCOL=6);
# S1C: TPC in Consumer growth rate
par(mfg=c(1,3), mar = c(1.5, 2.5, 3, 3));
Graph(BP=S1C_KT_BP, BPE=S1C_KT_BPE, LP=S1C_KT_LP, XX=0, YY=0, ZZ=0, TCOL=6);
# S1D: TPC in Consumer mortality rate
par(mfg=c(2,1), mar = c(3, 3, 1.5, 2.5));
Graph(BP=S1D_KT_BP, BPE=S1D_KT_BPE, LP=S1D_KT_LP, XX=1, YY=0, ZZ=1, TCOL=6);
# S1F: TPC in Predator Functional response
par(mfg=c(2,2), mar = c(3, 2.5, 1.5, 2.5));
Graph(BP=S1F_KT_BP, BPE=S1F_KT_BPE, LP=S1F_KT_LP, XX=1, YY=0, ZZ=1, TCOL=4);
# S1G: TPC in Predator mortality rate
par(mfg=c(2,3), mar = c(3, 2.5, 1.5, 3));
Graph(BP=S1G_KT_BP, BPE=S1G_KT_BPE, LP=S1G_KT_LP, XX=1, YY=0, ZZ=0, TCOL=4);

########## Fig. S4 - Influence of TSR alone on the trophic chain structure when implemented in a single/more traits ====
# a-c) TSR in a single traits
# d-f) TSR in two of the three traits
# g) TSR in all traits
## Figure ====
par(mfrow=c(2,4)); #To save in portrait, 9*12
# S2A: TSR in Predator mxm foraging size threshold Lv
par(mfg=c(1,1), mar = c(1.5, 3, 3, 2.5));
Graph(BP=S2A_KT_BP, BPE=S2A_KT_BPE, LP=S2A_KT_LP, XX=0, YY=0, ZZ=1, TCOL=6);
# S2B: TSR in Consumer maturation size Lj
par(mfg=c(1,2), mar = c(1.5, 2.5, 3, 2.5));
Graph(BP=S2B_KT_BP, BPE=S2B_KT_BPE, LP=S2B_KT_LP, XX=0, YY=0, ZZ=1, TCOL=6);
# S2C: TSR in Consumer asymptotic size Lm
par(mfg=c(1,3), mar = c(1.5, 2.5, 3, 2.5));
Graph(BP=S2C_KT_BP, BPE=S2C_KT_BPE, LP=S2C_KT_LP, XX=0, YY=0, ZZ=0, TCOL=6);
# S2D: TSR in lj & Lm
par(mfg=c(2,1), mar = c(3, 3, 1.5, 2.5));
Graph(BP=S2D_KT_BP, BPE=S2D_KT_BPE, LP=S2D_KT_LP, XX=1, YY=0, ZZ=1, TCOL=6);
# S2E: TSR in Lv & Lj
par(mfg=c(2,2), mar = c(3, 2.5, 1.5, 2.5));
Graph(BP=S2E_KT_BP, BPE=S2E_KT_BPE, LP=S2E_KT_LP, XX=1, YY=0, ZZ=1, TCOL=6);
# S2F: TSR in Lv & Lm
par(mfg=c(2,3), mar = c(3, 2.5, 1.5, 2.5));
Graph(BP=S2F_KT_BP, BPE=S2F_KT_BPE, LP=S2F_KT_LP, XX=1, YY=0, ZZ=1, TCOL=6);
# S2G: TSR in Lv, Lj & Lm
par(mfg=c(2,4), mar = c(3, 2.5, 1.5, 3));
Graph(BP=S2G_KT_BP, BPE=S2G_KT_BPE, LP=S2G_KT_LP, XX=1, YY=0, ZZ=0, TCOL=6);

########## Fig. S5 - 3x3 Combinations of TSR & TPC in species (C alone, P alone, and both species) ====
par(mfrow=c(3,3)); #To save in portrait, 9*9

# S3A: TPC(C) + TSR(C)
par(mfg=c(1,1), mar = c(1.5, 3, 3, 2.5));
Graph(BP=S3A_KT_BP, BPE=S3A_KT_BPE, LP=S3A_KT_LP, XX=0, YY=0, ZZ=1, TCOL=6);
# S3B: TPC(C) + TSR(P)
par(mfg=c(1,2), mar = c(1.5, 2.5, 3, 2.5));
Graph(BP=S3B_KT_BP, BPE=S3B_KT_BPE, LP=S3B_KT_LP, XX=0, YY=0, ZZ=1, TCOL=6);
# S3C: TPC(C) + TSR(C+P)
par(mfg=c(1,3), mar = c(1.5, 2.5, 3, 3));
Graph(BP=S3C_KT_BP, BPE=S3C_KT_BPE, LP=S3C_KT_LP, XX=0, YY=0, ZZ=0, TCOL=6);

# S3D: TPC(P) + TSR(C)
par(mfg=c(2,1), mar = c(1.5, 3, 1.5, 2.5));
Graph(BP=S3D_KT_BP, BPE=S3D_KT_BPE, LP=S3D_KT_LP, XX=0, YY=0, ZZ=1, TCOL=6);
# S3E: TPC(P) + TSR(P)
par(mfg=c(2,2), mar = c(1.5, 2.5, 1.5, 2.5));
Graph(BP=S3E_KT_BP, BPE=S3E_KT_BPE, LP=S3E_KT_LP, XX=0, YY=0, ZZ=1, TCOL=6);
# S3F: TPC(P) + TSR(C+P)
par(mfg=c(2,3), mar = c(1.5, 2.5, 1.5, 3));
Graph(BP=S3F_KT_BP, BPE=S3F_KT_BPE, LP=S3F_KT_LP, XX=0, YY=0, ZZ=0, TCOL=6);

# S3G: TPC(C+P) + TSR(C)
par(mfg=c(3,1), mar = c(3, 3, 1.5, 2.5));
Graph(BP=S3G_KT_BP, BPE=S3G_KT_BPE, LP=S3G_KT_LP, XX=1, YY=0, ZZ=1, TCOL=6);
# S3H: TPC(C+P) + TSR(P)
par(mfg=c(3,2), mar = c(3, 2.5, 1.5, 2.5));
Graph(BP=S3H_KT_BP, BPE=S3H_KT_BPE, LP=S3H_KT_LP, XX=1, YY=0, ZZ=1, TCOL=6);
# S3I: TPC(C+P) + TSR(C+P)
par(mfg=c(3,3), mar = c(3, 2.5, 1.5, 3));
Graph(BP=S3I_KT_BP, BPE=S3I_KT_BPE, LP=S3I_KT_LP, XX=1, YY=0, ZZ=0, TCOL=6);


########## Fig. S6 - Consumer Life histories ====
SEQ1<-c(1E-1, 2E-1, 3E-1, 4E-1, 5E-1, 6E-1, 7E-1, 8E-1, 9E-1)
SEQL1<-log10(SEQ1+1)
SEQ2<-c(1E-2, 2E-2, 3E-2, 4E-2, 5E-2, 6E-2, 7E-2, 8E-2, 9E-2)
SEQL2<-log10(SEQ2+1)
SEQ4<-c(1E-4, 2E-4, 3E-4, 4E-4, 5E-4, 6E-4, 7E-4, 8E-4, 9E-4)
SEQL4<-log10(SEQ4+1)
SEQ5<-c(1E-5, 2E-5, 3E-5, 4E-5, 5E-5, 6E-5, 7E-5, 8E-5, 9E-5)
SEQL5<-log10(SEQ5+1)
SEQ6<-c(1E-6, 2E-6, 3E-6, 4E-6, 5E-6, 6E-6, 7E-6, 8E-6, 9E-6)
SEQL6<-log10(SEQ6+1)
SEQ7<-c(1E-7, 2E-7, 3E-7, 4E-7, 5E-7, 6E-7, 7E-7, 8E-7, 9E-7)
SEQL7<-log10(SEQ7+1)
SEQ8<-c(1E-8, 2E-8, 3E-8, 4E-8, 5E-8, 6E-8, 7E-8, 8E-8, 9E-8)
SEQL8<-log10(SEQ8+1)

graphics.off()
plot.new()
par(mfrow = c(3,2))
par(cex = 1.0, cex.lab=1.3)
par(oma = c(2.5, 1.5, 1.5, 1.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))

### Row 1: Generation time ====
par(mfg=c(1,1), mar = c(.5, 2, 1, .5))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(20), log10(2E3)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(1,2,3), labels=c(expression(10^{1}), expression(10^{2}), expression(10^{3})), las=2, cex.axis=1.5);
axis(4, at=c(1,2,3), labels=F, las=2);
lines(S5E_PGR_II[,1], log10(S5E_PGR_II[,3]), type="l", col=rgb(0,0,0), lwd=2, lty=1); # L mat
lines(S5C_PGR_II[,1], log10(S5C_PGR_II[,3]), type="l", col=rgb(0,0,0), lwd=2, lty=2); # L inf
lines(S5G_PGR_II[,1], log10(S5G_PGR_II[,3]), type="l", col=rgb(0,0,0), lwd=2, lty=3); # L mat + L inf

par(mfg=c(1,2), mar = c(.5, .5, 1, 2))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(20), log10(2E3)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(1,2,3), labels=F, las=2);
axis(4, at=c(1,2,3), labels=c(expression(10^{1}), expression(10^{2}), expression(10^{3})), las=2, cex.axis=1.5);
lines(S5H_PGR_II[,1], log10(S5H_PGR_II[,3]), type="l", col=rgb(0,0,0), lwd=2, lty=1); # L mat
lines(S5I_PGR_II[,1], log10(S5I_PGR_II[,3]), type="l", col=rgb(0,0,0), lwd=2, lty=2); # L inf
lines(S5J_PGR_II[,1], log10(S5J_PGR_II[,3]), type="l", col=rgb(0,0,0), lwd=2, lty=3);# L mat + L inf

### Row 2: Critical resource density ====
par(mfg=c(2,1), mar = c(.5, 2, .5, .5))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-6), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=seq(-5,-4), labels=c(expression(10^{-5}), expression(10^{-4})), las=2, cex.axis=1.5)
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
lines(S2E_KT_BP[,6], log10(S2E_KT_BP[,1]), type="l", col=rgb(0,0,0), lwd=2, lty=1); # L mat
lines(S2C_KT_BP[,6], log10(S2C_KT_BP[,1]), type="l", col=rgb(0,0,0), lwd=2, lty=2); # L inf
lines(S2G_KT_BP[,6], log10(S2G_KT_BP[,1]), type="l", col=rgb(0,0,0), lwd=2, lty=3); # L mat + L inf

par(mfg=c(2,2), mar = c(.5, .5, .5, 2))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-6), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(4, at=seq(-5,-4), labels=c(expression(10^{-5}), expression(10^{-4})), las=2, cex.axis=1.5)
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
lines(S5H_KT_BP[,1], log10(S5H_KT_BP[,6]), type="l", col=rgb(0,0,0), lwd=2, lty=1); # L mat
lines(S5I_KT_BP[,1], log10(S5I_KT_BP[,6]), type="l", col=rgb(0,0,0), lwd=2, lty=2); # L inf
lines(S5J_KT_BP[,1], log10(S5J_KT_BP[,6]), type="l", col=rgb(0,0,0), lwd=2, lty=3); # L mat + L inf

### Row 3: Population growth rate ====
par(mfg=c(3,1), mar = c(2, 2, .5, .5))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(4E-3), log10(1.5E-1)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(1, at=c(5,15,25), label=T, las=1, cex.axis=1.5);
axis(2, at=seq(-2,-1), labels=c(expression(10^{-3}), expression(10^{-2})), las=2, cex.axis=1.5)
axis(2, at=c(log10(SEQ3), log10(SEQ2), log10(SEQ1)), label=F, las=2);
axis(4, at=c(log10(SEQ3), log10(SEQ2), log10(SEQ1)), label=F, las=2);
lines(S5E_PGR_II[,1], log10(S5E_PGR_II[,2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);# L mat
lines(S5C_PGR_II[,1], log10(S5C_PGR_II[,2]), type="l", col=rgb(0,0,0), lwd=2, lty=2);# L inf
lines(S5G_PGR_II[,1], log10(S5G_PGR_II[,2]), type="l", col=rgb(0,0,0), lwd=2, lty=3);# L mat + L inf

par(mfg=c(3,2), mar = c(2, .5, .5, 2))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(4E-3), log10(1.5E-1)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(1, at=c(5,15,25), label=T, las=1, cex.axis=1.5);
axis(4, at=seq(-2,-1), labels=c(expression(10^{-3}), expression(10^{-2})), las=2, cex.axis=1.5)
axis(2, at=c(log10(SEQ3), log10(SEQ2), log10(SEQ1)), label=F, las=2);
axis(4, at=c(log10(SEQ3), log10(SEQ2), log10(SEQ1)), label=F, las=2);
lines(S5H_PGR_II[,1], log10(S5H_PGR_II[,2]), type="l", col=rgb(0,0,0), lwd=2, lty=1);# L mat
lines(S5I_PGR_II[,1], log10(S5I_PGR_II[,2]), type="l", col=rgb(0,0,0), lwd=2, lty=2);# L inf
lines(S5J_PGR_II[,1], log10(S5J_PGR_II[,2]), type="l", col=rgb(0,0,0), lwd=2, lty=3);# L mat + L inf

########## Fig. S7 - Consumer biomasses & Birth rate ====
SEQ1<-c(1E-1, 2E-1, 3E-1, 4E-1, 5E-1, 6E-1, 7E-1, 8E-1, 9E-1)
SEQL1<-log10(SEQ1+1)
SEQ2<-c(1E-2, 2E-2, 3E-2, 4E-2, 5E-2, 6E-2, 7E-2, 8E-2, 9E-2)
SEQL2<-log10(SEQ2+1)
SEQ4<-c(1E-4, 2E-4, 3E-4, 4E-4, 5E-4, 6E-4, 7E-4, 8E-4, 9E-4)
SEQL4<-log10(SEQ4+1)
SEQ5<-c(1E-5, 2E-5, 3E-5, 4E-5, 5E-5, 6E-5, 7E-5, 8E-5, 9E-5)
SEQL5<-log10(SEQ5+1)
SEQ6<-c(1E-6, 2E-6, 3E-6, 4E-6, 5E-6, 6E-6, 7E-6, 8E-6, 9E-6)
SEQL6<-log10(SEQ6+1)
SEQ7<-c(1E-7, 2E-7, 3E-7, 4E-7, 5E-7, 6E-7, 7E-7, 8E-7, 9E-7)
SEQL7<-log10(SEQ7+1)
SEQ8<-c(1E-8, 2E-8, 3E-8, 4E-8, 5E-8, 6E-8, 7E-8, 8E-8, 9E-8)
SEQL8<-log10(SEQ8+1)

graphics.off()
plot.new()
par(mfrow = c(4,6))
par(cex = 1.0, cex.lab=1.3)
par(oma = c(2.5, 1.5, 1.5, 1.0))
par(tcl = 0.4)
par(mgp = c(2, 0.6, 0))

### Row 1: Birth rate ====
## TSR(Lmat) alone, and +TSR(Lv) ====
par(mfg=c(1,1), mar = c(.5, 1, .5, .25))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(9E-8), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=seq(-7,-4), labels=c("", expression(10^{-6}), expression(10^{-5}), expression(10^{-4})), las=2, cex.axis=1.1)
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TSR(Lmat)
lines(S5B_T_CR_II[,1], log10(S5B_T_CR_II[,5]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5B_T_PCR_IIa[1:44,1], log10(S5B_T_PCR_IIa[1:44,5]), type="l", col=rgb(.9,0,0), lwd=2, lty=2); # unstable branch
lines(S5B_T_PCR_IIa[44:dim(S5B_T_PCR_IIa)[1],1], log10(S5B_T_PCR_IIa[44:dim(S5B_T_PCR_IIa)[1],5]), type="l", col=rgb(.9,0,0), lwd=2, lty=1); # stable branch
points(S5B_T_LPa[1], log10(S5B_T_LPa[5]), col=rgb(.9,0,0), pch=16, lwd=2);
# TSR(Lv+Lmat)
lines(S5E_T_CR_II[,1], log10(S5E_T_CR_II[,5]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5E_T_PCR_IIa[,1], log10(S5E_T_PCR_IIa[,5]), type="l", col=rgb(.6,0,0), lwd=2, lty=4); # unstable branch
lines(S5E_T_PCR_IIb[,1], log10(S5E_T_PCR_IIb[,5]), type="l", col=rgb(.6,0,0), lwd=2, lty=3); # stable branch

## TPC +TSR(Lmat + Lv) ====
par(mfg=c(1,2), mar = c(.5, .25, .5, .75))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(9E-8), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TPC + TSR(Lmat)
S5H_T_CR_II[1,c(3:9)] <- 1E-10; #Artificially replacing the initial null values of consumer by a very low number to draw a continuous line from the invasion threshold
lines(S5H_T_CR_II[,1], log10(S5H_T_CR_II[,5]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5H_T_PCR_IIb[,1], log10(S5H_T_PCR_IIb[,5]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5H_T_PCR_IIa[,1], log10(S5H_T_PCR_IIa[,5]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5H_T_BPE_a[1], log10(S5H_T_BPE_a[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5H_T_BPE_b[1], log10(S5H_T_BPE_b[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5H_T_BPE_c[1], log10(S5H_T_BPE_c[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5H_T_BPE_d[1], log10(S5H_T_BPE_d[5]), col=rgb(0,0,0), pch=16, lwd=2);

# TPC + TSR(Lmat + Lv)
S5K_T_CR_II[1,c(3:9)] <- 1E-10; #Artificially replacing the initial null values of consumer by a very low number to draw a continuous line from the invasion threshold
lines(S5K_T_CR_II[,1], log10(S5K_T_CR_II[,5]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5K_T_PCR_IIb[,1], log10(S5K_T_PCR_IIb[,5]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5K_T_PCR_IIa[,1], log10(S5K_T_PCR_IIa[,5]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5K_T_BPE_a[1], log10(S5K_T_BPE_a[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5K_T_BPE_b[1], log10(S5K_T_BPE_b[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5K_T_BPE_c[1], log10(S5K_T_BPE_c[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5K_T_BPE_d[1], log10(S5K_T_BPE_d[5]), col=rgb(0,0,0), pch=16, lwd=2);

## TSR(Linf) alone, and +TSR(Lv) ====
par(mfg=c(1,3), mar = c(.5, .75, .5, .25))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(9E-8), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TSR(Linf)
lines(S5C_T_CR_II[,1], log10(S5C_T_CR_II[,5]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5C_T_PCR_IIa[1:48,1], log10(S5C_T_PCR_IIa[1:48,5]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5C_T_PCR_IIa[48:dim(S5C_T_PCR_IIa)[1],1], log10(S5C_T_PCR_IIa[48:dim(S5C_T_PCR_IIa)[1],5]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5C_T_LPa[1], log10(S5C_T_LPa[5]), col=rgb(.9,0,0), pch=16, lwd=2);
# TSR(Lv+Linf)
lines(S5F_T_CR_II[,1], log10(S5F_T_CR_II[,5]), type="l", col=rgb(0,0,0), lwd=2, lty=3); 
lines(S5F_T_PCR_II[1:264,1], log10(S5F_T_PCR_II[1:264,5]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5F_T_PCR_II[264:dim(S5F_T_PCR_II)[1],1], log10(S5F_T_PCR_II[264:dim(S5F_T_PCR_II)[1],5]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5F_T_LP[1], log10(S5F_T_LP[5]), col=rgb(.6,0,0), pch=16, lwd=2);

## TPC + TSR(Linf+Lv) ====
par(mfg=c(1,4), mar = c(.5, .25, .5, .75))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(9E-8), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TPC + TSR(Linf)
S5I_T_CR_II[1,c(3:9)] <- 1E-10; #Artificially replacing the initial null values of consumer by a very low number to draw a continuous line from the invasion threshold
lines(S5I_T_CR_II[,1], log10(S5I_T_CR_II[,5]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5I_T_PCR_IIb[,1], log10(S5I_T_PCR_IIb[,5]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5I_T_PCR_IIa[,1], log10(S5I_T_PCR_IIa[,5]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5I_T_BPE_a[1], log10(S5I_T_BPE_a[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5I_T_BPE_b[1], log10(S5I_T_BPE_b[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5I_T_BPE_c[1], log10(S5I_T_BPE_c[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5I_T_BPE_d[1], log10(S5I_T_BPE_d[5]), col=rgb(0,0,0), pch=16, lwd=2);

# TSR(Lv) + TSR+TPC(Linf)
S5L_T_CR_II[1,c(3:9)] <- 1E-10; #Artificially replacing the initial null values of consumer by a very low number to draw a continuous line from the invasion threshold
lines(S5L_T_CR_II[,1], log10(S5L_T_CR_II[,5]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5L_T_PCR_IIb[,1], log10(S5L_T_PCR_IIb[,5]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5L_T_PCR_IIa[,1], log10(S5L_T_PCR_IIa[,5]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5L_T_BPE_a[1], log10(S5L_T_BPE_a[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5L_T_BPE_b[1], log10(S5L_T_BPE_b[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5L_T_BPE_c[1], log10(S5L_T_BPE_c[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5L_T_BPE_d[1], log10(S5L_T_BPE_d[5]), col=rgb(0,0,0), pch=16, lwd=1);

## TSR(Lmat+Linf) alone, and +TSR(Lv) ====
par(mfg=c(1,5), mar = c(.5, .75, .5, .25))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(9E-8), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TSR(Lmat+Linf)
lines(S5D_T_CR_II[,1], log10(S5D_T_CR_II[,5]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5D_T_PCR_IIa[1:32,1], log10(S5D_T_PCR_IIa[1:32,5]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);#unstable branch
lines(S5D_T_PCR_IIa[32:45,1], log10(S5D_T_PCR_IIa[32:45,5]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);#stable branch
points(S5D_T_LPa[1], log10(S5D_T_LPa[5]), col=rgb(.9,0,0), pch=16, lwd=2);

# TSR(Lv+Lmat+Linf)
lines(S5G_T_CR_II[,1], log10(S5G_T_CR_II[,5]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5G_T_PCR_IIa[1:51,1], log10(S5G_T_PCR_IIa[1:51,5]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);#unstable branch
lines(S5G_T_PCR_IIa[51:dim(S5G_T_PCR_IIa)[1],1], log10(S5G_T_PCR_IIa[51:dim(S5G_T_PCR_IIa)[1],5]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);#stable branch
points(S5G_T_LPa[1], log10(S5G_T_LPa[5]), col=rgb(.6,0,0), pch=16, lwd=2);

## TPC + TSR(Lmat+Linf+Lv) ====
par(mfg=c(1,6), mar = c(.5, .25, .5, 1.5))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(9E-8), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=seq(-7,-4), labels=c("", expression(10^{-6}), expression(10^{-5}), expression(10^{-4})), las=2, cex.axis=1.1)
# TPC + TSR(Lmat+Linf)
S5J_T_CR_II[c(1,77),c(3:9)] <- 1E-10; #Artificially replacing the initial null values of consumer by a very low number to draw a continuous line from the invasion threshold
lines(S5J_T_CR_II[,1], log10(S5J_T_CR_II[,5]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5J_T_PCR_IIb[,1], log10(S5J_T_PCR_IIb[,5]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5J_T_PCR_IIa[,1], log10(S5J_T_PCR_IIa[,5]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5J_T_BPE_a[1], log10(S5J_T_BPE_a[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5J_T_BPE_b[1], log10(S5J_T_BPE_b[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5J_T_BPE_c[1], log10(S5J_T_BPE_c[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5J_T_BPE_d[1], log10(S5J_T_BPE_d[5]), col=rgb(0,0,0), pch=16, lwd=2);

# TPC + TSR(Lv+Lmat+Linf)
S5M_T_CR_II[c(1,77),c(3:9)] <- 1E-10; #Artificially replacing the initial null values of consumer by a very low number to draw a continuous line from the invasion threshold
lines(S5M_T_CR_II[,1], log10(S5M_T_CR_II[,5]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5M_T_PCR_IIb[,1], log10(S5M_T_PCR_IIb[,5]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5M_T_PCR_IIa[,1], log10(S5M_T_PCR_IIa[,5]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5M_T_BPE_a[1], log10(S5M_T_BPE_a[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5M_T_BPE_b[1], log10(S5M_T_BPE_b[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5M_T_BPE_c[1], log10(S5M_T_BPE_c[5]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5M_T_BPE_d[1], log10(S5M_T_BPE_d[5]), col=rgb(0,0,0), pch=16, lwd=2);


### Row 2: Biomass Vuln Juv ====
## TSR(Lmat) alone, and +TSR(Lv) ====
par(mfg=c(2,1), mar = c(.5, 1, .5, .25))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-7), log10(9E-5)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=seq(-7,-4), labels=c("", expression(10^{-6}), expression(10^{-5}), expression(10^{-4})), las=2, cex.axis=1.1)
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TSR(Lmat)
lines(S5B_T_CR_II[,1], log10(S5B_T_CR_II[,7]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5B_T_PCR_IIa[1:44,1], log10(S5B_T_PCR_IIa[1:44,7]), type="l", col=rgb(.9,0,0), lwd=2, lty=2); # unstable branch
lines(S5B_T_PCR_IIa[44:dim(S5B_T_PCR_IIa)[1],1], log10(S5B_T_PCR_IIa[44:dim(S5B_T_PCR_IIa)[1],7]), type="l", col=rgb(.9,0,0), lwd=2, lty=1); # stable branch
points(S5B_T_LPa[1], log10(S5B_T_LPa[7]), col=rgb(.9,0,0), pch=16, lwd=2);
# TSR(Lv+Lmat)
lines(S5E_T_CR_II[,1], log10(S5E_T_CR_II[,7]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5E_T_PCR_IIa[,1], log10(S5E_T_PCR_IIa[,7]), type="l", col=rgb(.6,0,0), lwd=2, lty=4); # unstable branch
lines(S5E_T_PCR_IIb[,1], log10(S5E_T_PCR_IIb[,7]), type="l", col=rgb(.6,0,0), lwd=2, lty=3); # stable branch

## TPC +TSR(Lmat + Lv) ====
par(mfg=c(2,2), mar = c(.5, .25, .5, .75))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-7), log10(9E-5)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TPC + TSR(Lmat)
lines(S5H_T_CR_II[,1], log10(S5H_T_CR_II[,7]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5H_T_PCR_IIb[,1], log10(S5H_T_PCR_IIb[,7]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5H_T_PCR_IIa[,1], log10(S5H_T_PCR_IIa[,7]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5H_T_BPE_a[1], log10(S5H_T_BPE_a[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5H_T_BPE_b[1], log10(S5H_T_BPE_b[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5H_T_BPE_c[1], log10(S5H_T_BPE_c[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5H_T_BPE_d[1], log10(S5H_T_BPE_d[7]), col=rgb(0,0,0), pch=16, lwd=2);
# TPC + TSR(Lmat + Lv)
lines(S5K_T_CR_II[,1], log10(S5K_T_CR_II[,7]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5K_T_PCR_IIb[,1], log10(S5K_T_PCR_IIb[,7]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5K_T_PCR_IIa[,1], log10(S5K_T_PCR_IIa[,7]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5K_T_BPE_a[1], log10(S5K_T_BPE_a[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5K_T_BPE_b[1], log10(S5K_T_BPE_b[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5K_T_BPE_c[1], log10(S5K_T_BPE_c[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5K_T_BPE_d[1], log10(S5K_T_BPE_d[7]), col=rgb(0,0,0), pch=16, lwd=2);

## TSR(Linf) alone, and +TSR(Lv) ====
par(mfg=c(2,3), mar = c(.5, .75, .5, .25))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-7), log10(9E-5)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TSR(Linf)
lines(S5C_T_CR_II[,1], log10(S5C_T_CR_II[,7]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5C_T_PCR_IIa[1:48,1], log10(S5C_T_PCR_IIa[1:48,7]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5C_T_PCR_IIa[48:dim(S5C_T_PCR_IIa)[1],1], log10(S5C_T_PCR_IIa[48:dim(S5C_T_PCR_IIa)[1],7]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5C_T_LPa[1], log10(S5C_T_LPa[7]), col=rgb(.9,0,0), pch=16, lwd=2);
# TSR(Lv+Linf)
lines(S5F_T_CR_II[,1], log10(S5F_T_CR_II[,7]), type="l", col=rgb(0,0,0), lwd=2, lty=3); 
lines(S5F_T_PCR_II[1:264,1], log10(S5F_T_PCR_II[1:264,7]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5F_T_PCR_II[264:dim(S5F_T_PCR_II)[1],1], log10(S5F_T_PCR_II[264:dim(S5F_T_PCR_II)[1],7]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5F_T_LP[1], log10(S5F_T_LP[7]), col=rgb(.6,0,0), pch=16, lwd=2);

## TPC + TSR(Linf+Lv) ====
par(mfg=c(2,4), mar = c(.5, .25, .5, .75))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-7), log10(9E-5)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TPC + TSR(Linf)
lines(S5I_T_CR_II[,1], log10(S5I_T_CR_II[,7]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5I_T_PCR_IIb[,1], log10(S5I_T_PCR_IIb[,7]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5I_T_PCR_IIa[,1], log10(S5I_T_PCR_IIa[,7]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5I_T_BPE_a[1], log10(S5I_T_BPE_a[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5I_T_BPE_b[1], log10(S5I_T_BPE_b[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5I_T_BPE_c[1], log10(S5I_T_BPE_c[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5I_T_BPE_d[1], log10(S5I_T_BPE_d[7]), col=rgb(0,0,0), pch=16, lwd=2);
# TSR(Lv) + TSR+TPC(Linf)
lines(S5L_T_CR_II[,1], log10(S5L_T_CR_II[,7]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5L_T_PCR_IIb[,1], log10(S5L_T_PCR_IIb[,7]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5L_T_PCR_IIa[,1], log10(S5L_T_PCR_IIa[,7]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5L_T_BPE_a[1], log10(S5L_T_BPE_a[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5L_T_BPE_b[1], log10(S5L_T_BPE_b[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5L_T_BPE_c[1], log10(S5L_T_BPE_c[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5L_T_BPE_d[1], log10(S5L_T_BPE_d[7]), col=rgb(0,0,0), pch=16, lwd=1);

## TSR(Lmat+Linf) alone, and +TSR(Lv) ====
par(mfg=c(2,5), mar = c(.5, .75, .5, .25))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-7), log10(9E-5)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TSR(Lmat+Linf)
lines(S5D_T_CR_II[,1], log10(S5D_T_CR_II[,7]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5D_T_PCR_IIa[1:32,1], log10(S5D_T_PCR_IIa[1:32,7]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);#unstable branch
lines(S5D_T_PCR_IIa[32:45,1], log10(S5D_T_PCR_IIa[32:45,7]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);#stable branch
points(S5D_T_LPa[1], log10(S5D_T_LPa[7]), col=rgb(.9,0,0), pch=16, lwd=2);

# TSR(Lv+Lmat+Linf)
lines(S5G_T_CR_II[,1], log10(S5G_T_CR_II[,7]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5G_T_PCR_IIa[1:51,1], log10(S5G_T_PCR_IIa[1:51,7]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);#unstable branch
lines(S5G_T_PCR_IIa[51:dim(S5G_T_PCR_IIa)[1],1], log10(S5G_T_PCR_IIa[51:dim(S5G_T_PCR_IIa)[1],7]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);#stable branch
points(S5G_T_LPa[1], log10(S5G_T_LPa[7]), col=rgb(.6,0,0), pch=16, lwd=2);

## TPC + TSR(Lmat+Linf+Lv) ====
par(mfg=c(2,6), mar = c(.5, .25, .5, 1.5))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-7), log10(9E-5)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=seq(-7,-4), labels=c("", expression(10^{-6}), expression(10^{-5}), expression(10^{-4})), las=2, cex.axis=1.1)
# TPC + TSR(Lmat+Linf)
lines(S5J_T_CR_II[,1], log10(S5J_T_CR_II[,7]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5J_T_PCR_IIb[,1], log10(S5J_T_PCR_IIb[,7]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5J_T_PCR_IIa[,1], log10(S5J_T_PCR_IIa[,7]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5J_T_BPE_a[1], log10(S5J_T_BPE_a[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5J_T_BPE_b[1], log10(S5J_T_BPE_b[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5J_T_BPE_c[1], log10(S5J_T_BPE_c[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5J_T_BPE_d[1], log10(S5J_T_BPE_d[7]), col=rgb(0,0,0), pch=16, lwd=2);
# TPC + TSR(Lv+Lmat+Linf)
lines(S5M_T_CR_II[,1], log10(S5M_T_CR_II[,7]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5M_T_PCR_IIb[,1], log10(S5M_T_PCR_IIb[,7]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5M_T_PCR_IIa[,1], log10(S5M_T_PCR_IIa[,7]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5M_T_BPE_a[1], log10(S5M_T_BPE_a[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5M_T_BPE_b[1], log10(S5M_T_BPE_b[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5M_T_BPE_c[1], log10(S5M_T_BPE_c[7]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5M_T_BPE_d[1], log10(S5M_T_BPE_d[7]), col=rgb(0,0,0), pch=16, lwd=2);
### Row 3: Biomass Juv ====
## TSR(Lmat) alone, and +TSR(Lv) ====
par(mfg=c(3,1), mar = c(.5, 1, .5, .25))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-6), log10(8E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=seq(-7,-4), labels=c("", expression(10^{-6}), expression(10^{-5}), expression(10^{-4})), las=2, cex.axis=1.1)
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TSR(Lmat)
lines(S5B_T_CR_II[,1], log10(S5B_T_CR_II[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5B_T_PCR_IIa[1:44,1], log10(S5B_T_PCR_IIa[1:44,8]), type="l", col=rgb(.9,0,0), lwd=2, lty=2); # unstable branch
lines(S5B_T_PCR_IIa[44:dim(S5B_T_PCR_IIa)[1],1], log10(S5B_T_PCR_IIa[44:dim(S5B_T_PCR_IIa)[1],8]), type="l", col=rgb(.9,0,0), lwd=2, lty=1); # stable branch
points(S5B_T_LPa[1], log10(S5B_T_LPa[8]), col=rgb(.9,0,0), pch=16, lwd=2);
# TSR(Lv+Lmat)
lines(S5E_T_CR_II[,1], log10(S5E_T_CR_II[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5E_T_PCR_IIa[,1], log10(S5E_T_PCR_IIa[,8]), type="l", col=rgb(.6,0,0), lwd=2, lty=4); # unstable branch
lines(S5E_T_PCR_IIb[,1], log10(S5E_T_PCR_IIb[,8]), type="l", col=rgb(.6,0,0), lwd=2, lty=3); # stable branch

## TPC +TSR(Lmat + Lv) ====
par(mfg=c(3,2), mar = c(.5, .25, .5, .75))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-6), log10(8E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TPC + TSR(Lmat)
lines(S5H_T_CR_II[,1], log10(S5H_T_CR_II[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5H_T_PCR_IIb[,1], log10(S5H_T_PCR_IIb[,8]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5H_T_PCR_IIa[,1], log10(S5H_T_PCR_IIa[,8]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5H_T_BPE_a[1], log10(S5H_T_BPE_a[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5H_T_BPE_b[1], log10(S5H_T_BPE_b[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5H_T_BPE_c[1], log10(S5H_T_BPE_c[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5H_T_BPE_d[1], log10(S5H_T_BPE_d[8]), col=rgb(0,0,0), pch=16, lwd=2);

# TPC + TSR(Lmat + Lv)
lines(S5K_T_CR_II[,1], log10(S5K_T_CR_II[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5K_T_PCR_IIb[,1], log10(S5K_T_PCR_IIb[,8]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5K_T_PCR_IIa[,1], log10(S5K_T_PCR_IIa[,8]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5K_T_BPE_a[1], log10(S5K_T_BPE_a[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5K_T_BPE_b[1], log10(S5K_T_BPE_b[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5K_T_BPE_c[1], log10(S5K_T_BPE_c[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5K_T_BPE_d[1], log10(S5K_T_BPE_d[8]), col=rgb(0,0,0), pch=16, lwd=2);

## TSR(Linf) alone, and +TSR(Lv) ====
par(mfg=c(3,3), mar = c(.5, .75, .5, .25))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-6), log10(8E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TSR(Linf)
lines(S5C_T_CR_II[,1], log10(S5C_T_CR_II[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5C_T_PCR_IIa[1:48,1], log10(S5C_T_PCR_IIa[1:48,8]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5C_T_PCR_IIa[48:dim(S5C_T_PCR_IIa)[1],1], log10(S5C_T_PCR_IIa[48:dim(S5C_T_PCR_IIa)[1],8]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5C_T_LPa[1], log10(S5C_T_LPa[8]), col=rgb(.9,0,0), pch=16, lwd=2);
# TSR(Lv+Linf)
lines(S5F_T_CR_II[,1], log10(S5F_T_CR_II[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=3); 
lines(S5F_T_PCR_II[1:264,1], log10(S5F_T_PCR_II[1:264,8]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5F_T_PCR_II[264:dim(S5F_T_PCR_II)[1],1], log10(S5F_T_PCR_II[264:dim(S5F_T_PCR_II)[1],8]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5F_T_LP[1], log10(S5F_T_LP[8]), col=rgb(.6,0,0), pch=16, lwd=2);

## TPC + TSR(Linf+Lv) ====
par(mfg=c(3,4), mar = c(.5, .25, .5, .75))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-6), log10(8E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TPC + TSR(Linf)
lines(S5I_T_CR_II[,1], log10(S5I_T_CR_II[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5I_T_PCR_IIb[,1], log10(S5I_T_PCR_IIb[,8]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5I_T_PCR_IIa[,1], log10(S5I_T_PCR_IIa[,8]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5I_T_BPE_a[1], log10(S5I_T_BPE_a[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5I_T_BPE_b[1], log10(S5I_T_BPE_b[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5I_T_BPE_c[1], log10(S5I_T_BPE_c[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5I_T_BPE_d[1], log10(S5I_T_BPE_d[8]), col=rgb(0,0,0), pch=16, lwd=2);
# TSR(Lv) + TSR+TPC(Linf)
lines(S5L_T_CR_II[,1], log10(S5L_T_CR_II[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5L_T_PCR_IIb[,1], log10(S5L_T_PCR_IIb[,8]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5L_T_PCR_IIa[,1], log10(S5L_T_PCR_IIa[,8]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5L_T_BPE_a[1], log10(S5L_T_BPE_a[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5L_T_BPE_b[1], log10(S5L_T_BPE_b[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5L_T_BPE_c[1], log10(S5L_T_BPE_c[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5L_T_BPE_d[1], log10(S5L_T_BPE_d[8]), col=rgb(0,0,0), pch=16, lwd=1);

## TSR(Lmat+Linf) alone, and +TSR(Lv) ====
par(mfg=c(3,5), mar = c(.5, .75, .5, .25))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-6), log10(8E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TSR(Lmat+Linf)
lines(S5D_T_CR_II[,1], log10(S5D_T_CR_II[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5D_T_PCR_IIa[1:32,1], log10(S5D_T_PCR_IIa[1:32,8]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);#unstable branch
lines(S5D_T_PCR_IIa[32:45,1], log10(S5D_T_PCR_IIa[32:45,8]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);#stable branch
points(S5D_T_LPa[1], log10(S5D_T_LPa[8]), col=rgb(.9,0,0), pch=16, lwd=2);

# TSR(Lv+Lmat+Linf)
lines(S5G_T_CR_II[,1], log10(S5G_T_CR_II[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5G_T_PCR_IIa[1:51,1], log10(S5G_T_PCR_IIa[1:51,8]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);#unstable branch
lines(S5G_T_PCR_IIa[51:dim(S5G_T_PCR_IIa)[1],1], log10(S5G_T_PCR_IIa[51:dim(S5G_T_PCR_IIa)[1],8]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);#stable branch
points(S5G_T_LPa[1], log10(S5G_T_LPa[8]), col=rgb(.6,0,0), pch=16, lwd=2);

## TPC + TSR(Lmat+Linf+Lv) ====
par(mfg=c(3,6), mar = c(.5, .25, .5, 1.5))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(3E-6), log10(8E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=seq(-7,-4), labels=c("", expression(10^{-6}), expression(10^{-5}), expression(10^{-4})), las=2, cex.axis=1.1)
# TPC + TSR(Lmat+Linf)
lines(S5J_T_CR_II[,1], log10(S5J_T_CR_II[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5J_T_PCR_IIb[,1], log10(S5J_T_PCR_IIb[,8]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5J_T_PCR_IIa[,1], log10(S5J_T_PCR_IIa[,8]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5J_T_BPE_a[1], log10(S5J_T_BPE_a[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5J_T_BPE_b[1], log10(S5J_T_BPE_b[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5J_T_BPE_c[1], log10(S5J_T_BPE_c[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5J_T_BPE_d[1], log10(S5J_T_BPE_d[8]), col=rgb(0,0,0), pch=16, lwd=2);
# TPC + TSR(Lv+Lmat+Linf)
lines(S5M_T_CR_II[,1], log10(S5M_T_CR_II[,8]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5M_T_PCR_IIb[,1], log10(S5M_T_PCR_IIb[,8]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5M_T_PCR_IIa[,1], log10(S5M_T_PCR_IIa[,8]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5M_T_BPE_a[1], log10(S5M_T_BPE_a[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5M_T_BPE_b[1], log10(S5M_T_BPE_b[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5M_T_BPE_c[1], log10(S5M_T_BPE_c[8]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5M_T_BPE_d[1], log10(S5M_T_BPE_d[8]), col=rgb(0,0,0), pch=16, lwd=2);

### Row 4: Biomass Adults ====
## TSR(Lmat) alone, and +TSR(Lv) ====
par(mfg=c(4,1), mar = c(1, 1, .5, .25))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-7), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=1);
axis(1, at=c(5,15,25), label=T, las=1, cex.axis=1.1);
axis(2, at=seq(-7,-4), labels=c("", expression(10^{-6}), expression(10^{-5}), expression(10^{-4})), las=2, cex.axis=1.1)
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TSR(Lmat)
lines(S5B_T_CR_II[,1], log10(S5B_T_CR_II[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5B_T_PCR_IIa[1:44,1], log10(S5B_T_PCR_IIa[1:44,9]), type="l", col=rgb(.9,0,0), lwd=2, lty=2); # unstable branch
lines(S5B_T_PCR_IIa[44:dim(S5B_T_PCR_IIa)[1],1], log10(S5B_T_PCR_IIa[44:dim(S5B_T_PCR_IIa)[1],9]), type="l", col=rgb(.9,0,0), lwd=2, lty=1); # stable branch
points(S5B_T_LPa[1], log10(S5B_T_LPa[9]), col=rgb(.9,0,0), pch=16, lwd=2);
# TSR(Lv+Lmat)
lines(S5E_T_CR_II[,1], log10(S5E_T_CR_II[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5E_T_PCR_IIa[,1], log10(S5E_T_PCR_IIa[,9]), type="l", col=rgb(.6,0,0), lwd=2, lty=4); # unstable branch
lines(S5E_T_PCR_IIb[,1], log10(S5E_T_PCR_IIb[,9]), type="l", col=rgb(.6,0,0), lwd=2, lty=3); # stable branch

## TPC +TSR(Lmat + Lv) ====
par(mfg=c(4,2), mar = c(1, .25, .5, .75))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-7), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=1);
axis(1, at=c(5,15,25), label=T, las=1, cex.axis=1.1);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TPC + TSR(Lmat)
lines(S5H_T_CR_II[,1], log10(S5H_T_CR_II[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5H_T_PCR_IIb[,1], log10(S5H_T_PCR_IIb[,9]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5H_T_PCR_IIa[,1], log10(S5H_T_PCR_IIa[,9]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5H_T_BPE_a[1], log10(S5H_T_BPE_a[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5H_T_BPE_b[1], log10(S5H_T_BPE_b[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5H_T_BPE_c[1], log10(S5H_T_BPE_c[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5H_T_BPE_d[1], log10(S5H_T_BPE_d[9]), col=rgb(0,0,0), pch=16, lwd=2);

# TPC + TSR(Lmat + Lv)
lines(S5K_T_CR_II[,1], log10(S5K_T_CR_II[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5K_T_PCR_IIb[,1], log10(S5K_T_PCR_IIb[,9]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5K_T_PCR_IIa[,1], log10(S5K_T_PCR_IIa[,9]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5K_T_BPE_a[1], log10(S5K_T_BPE_a[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5K_T_BPE_b[1], log10(S5K_T_BPE_b[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5K_T_BPE_c[1], log10(S5K_T_BPE_c[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5K_T_BPE_d[1], log10(S5K_T_BPE_d[9]), col=rgb(0,0,0), pch=16, lwd=2);

## TSR(Linf) alone, and +TSR(Lv) ====
par(mfg=c(4,3), mar = c(1, .75, .5, .25))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-7), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=1);
axis(1, at=c(5,15,25), label=T, las=1, cex.axis=1.1);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TSR(Linf)
lines(S5C_T_CR_II[,1], log10(S5C_T_CR_II[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5C_T_PCR_IIa[1:48,1], log10(S5C_T_PCR_IIa[1:48,9]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5C_T_PCR_IIa[48:dim(S5C_T_PCR_IIa)[1],1], log10(S5C_T_PCR_IIa[48:dim(S5C_T_PCR_IIa)[1],9]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5C_T_LPa[1], log10(S5C_T_LPa[9]), col=rgb(.9,0,0), pch=16, lwd=2);
# TSR(Lv+Linf)
lines(S5F_T_CR_II[,1], log10(S5F_T_CR_II[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=3); 
lines(S5F_T_PCR_II[1:264,1], log10(S5F_T_PCR_II[1:264,9]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5F_T_PCR_II[264:dim(S5F_T_PCR_II)[1],1], log10(S5F_T_PCR_II[264:dim(S5F_T_PCR_II)[1],9]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5F_T_LP[1], log10(S5F_T_LP[9]), col=rgb(.6,0,0), pch=16, lwd=2);

## TPC + TSR(Linf+Lv) ====
par(mfg=c(4,4), mar = c(1, .25, .5, .75))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-7), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=1);
axis(1, at=c(5,15,25), label=T, las=1, cex.axis=1.1);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TPC + TSR(Linf)
lines(S5I_T_CR_II[,1], log10(S5I_T_CR_II[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5I_T_PCR_IIb[,1], log10(S5I_T_PCR_IIb[,9]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5I_T_PCR_IIa[,1], log10(S5I_T_PCR_IIa[,9]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5I_T_BPE_a[1], log10(S5I_T_BPE_a[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5I_T_BPE_b[1], log10(S5I_T_BPE_b[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5I_T_BPE_c[1], log10(S5I_T_BPE_c[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5I_T_BPE_d[1], log10(S5I_T_BPE_d[9]), col=rgb(0,0,0), pch=16, lwd=2);
# TSR(Lv) + TSR+TPC(Linf)
lines(S5L_T_CR_II[,1], log10(S5L_T_CR_II[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5L_T_PCR_IIb[,1], log10(S5L_T_PCR_IIb[,9]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5L_T_PCR_IIa[,1], log10(S5L_T_PCR_IIa[,9]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5L_T_BPE_a[1], log10(S5L_T_BPE_a[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5L_T_BPE_b[1], log10(S5L_T_BPE_b[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5L_T_BPE_c[1], log10(S5L_T_BPE_c[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5L_T_BPE_d[1], log10(S5L_T_BPE_d[9]), col=rgb(0,0,0), pch=16, lwd=1);

## TSR(Lmat+Linf) alone, and +TSR(Lv) ====
par(mfg=c(4,5), mar = c(1, .75, .5, .25))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-7), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(1, at=c(5,15,25), label=T, las=1, cex.axis=1.1);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
# TSR(Lmat+Linf)
lines(S5D_T_CR_II[,1], log10(S5D_T_CR_II[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5D_T_PCR_IIa[1:32,1], log10(S5D_T_PCR_IIa[1:32,9]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);#unstable branch
lines(S5D_T_PCR_IIa[32:45,1], log10(S5D_T_PCR_IIa[32:45,9]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);#stable branch
points(S5D_T_LPa[1], log10(S5D_T_LPa[9]), col=rgb(.9,0,0), pch=16, lwd=2);

# TSR(Lv+Lmat+Linf)
lines(S5G_T_CR_II[,1], log10(S5G_T_CR_II[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5G_T_PCR_IIa[1:51,1], log10(S5G_T_PCR_IIa[1:51,9]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);#unstable branch
lines(S5G_T_PCR_IIa[51:dim(S5G_T_PCR_IIa)[1],1], log10(S5G_T_PCR_IIa[51:dim(S5G_T_PCR_IIa)[1],9]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);#stable branch
points(S5G_T_LPa[1], log10(S5G_T_LPa[9]), col=rgb(.6,0,0), pch=16, lwd=2);

## TPC + TSR(Lmat+Linf+Lv) ====
par(mfg=c(4,6), mar = c(1, .25, .5, 1.5))
plot(1, 1, type="l", xlim=c(1, 30), xlab="", xaxt="n", xaxs="i", ylim=c(log10(2E-7), log10(5E-4)), ylab="", yaxt="n", yaxs="i");
axis(1, at=seq(0,30,by=5), label=F, las=2);
axis(1, at=c(5,15,25), label=T, las=1, cex.axis=1.1);
axis(2, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=c(log10(SEQ7), log10(SEQ6), log10(SEQ5), log10(SEQ4)), label=F, las=2);
axis(4, at=seq(-7,-4), labels=c("", expression(10^{-6}), expression(10^{-5}), expression(10^{-4})), las=2, cex.axis=1.1)
# TPC + TSR(Lmat+Linf)
lines(S5J_T_CR_II[,1], log10(S5J_T_CR_II[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=1);
lines(S5J_T_PCR_IIb[,1], log10(S5J_T_PCR_IIb[,9]), type="l", col=rgb(.9,0,0), lwd=2, lty=2);# unstable branch
lines(S5J_T_PCR_IIa[,1], log10(S5J_T_PCR_IIa[,9]), type="l", col=rgb(.9,0,0), lwd=2, lty=1);# stable branch
points(S5J_T_BPE_a[1], log10(S5J_T_BPE_a[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5J_T_BPE_b[1], log10(S5J_T_BPE_b[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5J_T_BPE_c[1], log10(S5J_T_BPE_c[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5J_T_BPE_d[1], log10(S5J_T_BPE_d[9]), col=rgb(0,0,0), pch=16, lwd=2);
# TPC + TSR(Lv+Lmat+Linf)
lines(S5M_T_CR_II[,1], log10(S5M_T_CR_II[,9]), type="l", col=rgb(0,0,0), lwd=2, lty=3);
lines(S5M_T_PCR_IIb[,1], log10(S5M_T_PCR_IIb[,9]), type="l", col=rgb(.6,0,0), lwd=2, lty=4);# unstable branch
lines(S5M_T_PCR_IIa[,1], log10(S5M_T_PCR_IIa[,9]), type="l", col=rgb(.6,0,0), lwd=2, lty=3);# stable branch
points(S5M_T_BPE_a[1], log10(S5M_T_BPE_a[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5M_T_BPE_b[1], log10(S5M_T_BPE_b[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5M_T_BPE_c[1], log10(S5M_T_BPE_c[9]), col=rgb(0,0,0), pch=16, lwd=2);
points(S5M_T_BPE_d[1], log10(S5M_T_BPE_d[9]), col=rgb(0,0,0), pch=16, lwd=2);
