require(FunWithNumbers)
require(sm)
require(mgcv)
require(sads)
require(stats4)
require(actuar)   # for the ZTP distribution

N<- c()           # sample size
n<-  9            # number of categories (probabilities)
n1<- n-1 
Nstep<- 2         # increment for N
Nvals<- 20

nrep<- 10000      # number of M.C. replications  
set.seed (321)
power_10<- matrix(nrow=1, ncol=Nvals)
power_5<-  matrix(nrow=1, ncol=Nvals)
power_1<-  matrix(nrow=1, ncol=Nvals)
p<- matrix(nrow=1,ncol=n)
unsq<- matrix(nrow=1,ncol=nrep)    # store the test statistic values here  
m<- matrix(nrow=n1, ncol=n1)
one<- c(rep(1,n1))
zero<- c(rep(0,n1))

s<- 2.0       # Change this value as required
lambda<-  1  # for the ZTP case
q<- c(0.0873	,	0.1133	,	0.1771		)	
#       10%      5%    1%           # Change these tabulated quantiles to match "s"

# START OF THE MONTE CARLO LOOP
# -----------------------------

for(jj in 1:Nvals) {
N[jj]<- 5+jj*Nstep      # smallest value of N is 10
for(ii in 1:nrep) {
H<- 0
#xb<- rzipf(N[jj],n,s)                       # The null hypothesis is TRUE
# Replace last line with one of the following for the null being FALSE
#xb<- rztpois(N[jj], lambda)                # Generate zero-truncated Poisson variates
xb <- benprob(N[jj])                       # Generate Benford first-digit variates

for (i in 1:n) {
H<- H+1/i^s
}
for (i in 1:n) {
p[i]<- 1/(i^s*H)
}

breaks<- c(0:n)
xbb <- binning(xb, breaks=breaks) 
freq<- xbb$table.freq                     # these are the empirical frequencies
term1<- (freq/N[jj] - p)                  # compare empirical relative frequencies & theoretical probs. 
S<- cumsum(term1)   
Ssq<- S^2
unsq[ii]<- (N[jj]/n)*(sum(Ssq)-Ssq[n]-(sum(S)-S[n])^2/n)

}     #END OF THE MONTE CARLO LOOP
#  ---------------------------

# we now have the "nrep" values for the test statistic.
# check to see what fractions of them exceed the various critical values

power_10[jj]<- mean(unsq>q[1])
power_5[jj]<- mean(unsq>q[2])
power_1[jj]<- mean(unsq>q[3])

}   # end of "N" loop

power_10
power_5
power_1


plot(N , power_10,col="black", lty=1,type="l",lwd=2,ylab="Power", ylim=c(0.0,1.0),main="Figure 4(d): Powers of Test \n(Benford Alternative; s = 2.0)" )  

 lines(N,power_5, col="red",lty=2,lwd=2)            
 lines(N, power_1, col="blue", lty=3,lwd=2)        
 
legend(33,0.9,legend = c(expression(paste(alpha, " = ", "0.10")),
                         expression(paste(alpha, " = ", "0.05")),
                         expression(paste(alpha, " = ", "0.01"))),
                         box.col="white", col=c("black","red", "blue"),
                         lty=c(1,2,3), lwd=c(2,2,2),ncol=1)
