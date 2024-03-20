## Compute a p-value, for a computed value of the test statistic:

require(mgcv)
U<- 0.0583      #  Computed value of test statistic
s<- 1.6547      #  (Estimated) value of "s"
n<- 5           #  Value of "n"
n1<- n-1
H<- 0
p<- c()

m<- matrix(nrow=n1, ncol=n1)
one<- c(rep(1,n1))
zero<- c(rep(0,n1))

for (i in 1:n) {
H<- H+1/i^s
}

for (i in 1:n) {
p[i]<- 1/(i^s*H)
}

for (i in 1:n1) {              
for (j in 1:n1) {
sum<- 0
for (k in 1:n1) {
sum<- sum+p[k]*(n-max(k,j))*min(k,j)
}
max1<- max(i,j)
min1<- min(i,j)
m[i,j]=(p[i]/(n^2))*((n-max1)*min1-sum)
}
}

lb<- eigen(m, only.values=TRUE)$values   ## weights

pval<- psum.chisq(U,lb,one,zero,tol=1e-06)
pval		# Note - this p-value will only be approximate if "s" is estimated
		# Thse approximation will be good, however if this (asymptotic) test is used with a sufficiently
		# large sample size and the MLE estimator has been used (given its consistency)	
