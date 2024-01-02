library(bellreg)           # Used to generate random "Bell" data
library(lamW)              # Used for the Lambert function to get MLE

set.seed(123456)
nrep<- 50000
nboot<- 999
n<- 500
theta<- 0.75
sumhata<- msehata<- 0

#START THE SIMULATION LOOP
##########################

for (i in 1:nrep) {

sumboot<- 0
y<- rbell(n,theta)
ybar<-mean(y)
mle<- lambertW0(ybar)    # Data are positive, so only the principal branch of the Lambert function exits.

# NOW START THE BOOTSTRAP LOOP
# ============================

for(ii in 1:nboot) {

# Non-parametric bootstrap
yy<-sample(y,n,replace=T)

# OR

# Parametric bootstrap
#yy<- rbell(n,mle)

yybar<- mean(yy)
theta_boot<- lambertW0(yybar)
sumboot<- sumboot+theta_boot

}

# END OF THE BOOTSTRAP LOOP
# CALCULATE THE BOOTSTRAP BIAS
bootbias<- (sumboot/nboot)-mle

# COMPUTE THE BIAS-ADJUSTED ESTIMATORS:
# ===================================
theta_adj<- mle-bootbias
sumhata<- sumhata+theta_adj
msehata<- msehata+(theta_adj-theta)^2

}
# END OF THE SIMULATION LOOP


# CALCULATE THE BIASES OF THE BOOTSTRAP BIAS-ADJUSTED ESTIMATORS

biashata<- (sumhata/nrep)-theta
# CALCULATE THE % BIASES:
pbiashata<- 100*biashata/theta

# CALCULATE THE % MSE's & % MSE's OF THE BOOTSTRAP BIAS-ADJUSTED ESTIMATORS

msehata<- msehata/nrep
pmsehata=100*msehata/(theta^2)

cat("Theta =", theta,"; Sample Size =", n,"; Replications =", nrep, "; Boot-samples = ", nboot, "\n")
#============================================================================================================
#
cat("Simulated % Bias of Adjusted Theta-hat =", pbiashata, "\n")
#
cat("Simulated % MSE of Adjusted Theta-hat =", pmsehata, "\n")


#############################################

