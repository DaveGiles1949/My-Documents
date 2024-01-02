# NOTE: MLE = MOM

library(bellreg)           # Used to generate random "Bell" data
library(lamW)              # Used for the Lambert function to get MLE

set.seed(123456)
nrep<- 50000
n<- 75
theta<- 3.0

theta_hat<- vector()            # Cox-Snell
theta_tilde<- vector()          # MLE
theta_cup<- vector()            # Firth
theta_star<- vector()           # Godwin-Giles

fun_firth<- function(x,dat,n) {
n*(dat/x-exp(x))+(x+2)/(2*(x+1))
}                                     # Used to get Firth estimator

fun_gg<- function(x,n,mle) {
x-x*(x+2)/(2*n*exp(x)*(1+x)^2)-mle
}                                     # Used to get Godwin-Giles estimator

for (i in 1:nrep) {
y<- rbell(n,theta)
ybar<-mean(y)
theta_tilde[i]<- lambertW0(ybar)    # Data are positive, so only the principal branch of the Lambert function exits.
est_bias<- -theta_tilde[i]*(theta_tilde[i]+2)/(2*n*exp(theta_tilde[i])*(1+theta_tilde[i])^2)
theta_hat[i]<- theta_tilde[i]-est_bias
theta_cup[i]<-uniroot(fun_firth, dat=ybar,n=n, lower=0.0001, upper=theta+1)$root     # Upper limit set acccording to true theta value
theta_star[i]<- uniroot(fun_gg, n=n, mle=theta_tilde[i], lower=0.0001, upper=theta+1)$root
}

summary(theta_tilde)
summary(theta_hat)
summary(theta_cup)
summary(theta_star)

perc_bias_tilde<- 100*(mean(theta_tilde)-theta)/theta
perc_bias_hat<- 100*(mean(theta_hat)-theta)/theta
perc_bias_cup<- 100*(mean(theta_cup)-theta)/theta
perc_bias_star<- 100*(mean(theta_star)-theta)/theta

perc_mse_tilde<- 100*((mean(theta_tilde)-theta)^2 + var(theta_tilde))/theta^2
perc_mse_hat<- 100*((mean(theta_hat)-theta)^2 + var(theta_hat))/theta^2
perc_mse_cup<- 100*((mean(theta_cup)-theta)^2 + var(theta_cup))/theta^2
perc_mse_star<- 100*((mean(theta_star)-theta)^2 + var(theta_star))/theta^2

c(perc_bias_tilde,perc_bias_hat,perc_bias_cup, perc_bias_star)
c(perc_mse_tilde,perc_mse_hat,perc_mse_cup,perc_mse_star)
c(theta, n, nrep)

par(mfrow=c(4,1))
hist(theta_tilde)
hist(theta_hat)
hist(theta_cup)
hist(theta_star)

# END OF THE BASIC SIMULATION EXPERIMENT
#############################################
n<- 10
# plot the bias function as a function of theta:
dev.off()
t<- seq(1, 600)/100
bias_func<- -t*(t+2)/(2*n*exp(t)*(1+t)^2)
plot(t,bias_func, type="l", xlab=expression(paste(theta)), ylab="Bias", main=c(expression(paste("Figure 2: First-Order Bias of")~ tilde(theta))))
mtext("(n = 10)", side=3)

# Where is the bias function's turning point?
f<- function(x) {
 db<-    x^3+3*x^2+2*x-2  
}
uniroot(f,lower=0.1, upper=0.8)   # The turning point is at theta = 0.5213

#############################################

