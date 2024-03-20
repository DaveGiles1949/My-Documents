require(mgcv)
require(formattable)      # Used to display output

q90<- c()
q95<- c()
q99<- c()
n<- c()
mini<- 0.     # mini and maxi provide a range that's anticipated to 
maxi<- .5     # "bracket" the values of the quantiles being computed

s<- 1.0

for (ii in seq(2,20,1)) {

n<- c(n,ii)
N<- ii
N1<- N-1
H<- 0
p<- c()

num_eig<- N1               # number of largest eigenvalues to be used (<= n1) 

m<- matrix(nrow=N1, ncol=N1)
one<- c(rep(1,num_eig))
zero<- c(rep(0,num_eig))

for (i in 1:N) {
H<- H+1/i^s
}

for (i in 1:N) {
p[i]<- 1/(i^s*H)
}

for (i in 1:N1) {              
for (j in 1:N1) {
sum<- 0
for (k in 1:N1) {
sum<- sum+p[k]*(N-max(k,j))*min(k,j)
}
max1<- max(i,j)
min1<- min(i,j)
m[i,j]=(p[i]/(N^2))*((N-max1)*min1-sum)
}
}

#lb <- eigs(m, num_eig, opts = list(retvec = FALSE))$values 

lb<- eigen(m, only.values=TRUE)$values   ## weights

f<- function(x,lb,df,nc) 1-psum.chisq(x,lb,df,nc,tol=1e-06) - 0.90
q90<- c(q90,uniroot(f, c(mini, maxi), tol = 0.00001,lb=lb,df=one,nc=zero )$root)
f<- function(x,lb,df,nc) 1-psum.chisq(x,lb,df,nc,tol=1e-06) - 0.95
q95<- c(q95,uniroot(f, c(mini, maxi), tol = 0.00001,lb=lb,df=one,nc=zero )$root)
f<- function(x,lb,df,nc) 1-psum.chisq(x,lb,df,nc,tol=1e-06) - 0.99
q99<- c(q99,uniroot(f, c(mini, maxi), tol = 0.00001,lb=lb,df=one,nc=zero )$root)

}


# Print the quantile results in a clean format

s
df1<-data.frame(n,q90,q95,q99)
df1$n<-formattable(df1$n,format="f",digits=0)
df1$q90<-formattable(df1$q90,format="f",digits=4)
df1$q95<-formattable(df1$q95,format="f",digits=4)
df1$q99<-formattable(df1$q99,format="f",digits=4)

print(df1, row.names = FALSE) 