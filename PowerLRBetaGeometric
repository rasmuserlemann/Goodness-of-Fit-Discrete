library(stats)
library(DiscreteWeibull)

loglike = function(theta, x){
  if (theta[1] <= 0 | theta[2] <= 0){
    return(10**(10))
  }
  if (theta[1] > 100 | theta[2] > 100){
    return(10**(10))
  }
  final = 0
  alpha = theta[1]
  beta = theta[2]
  for (data in x){
    final = final + log(beta(alpha+1, data+beta)/beta(alpha,beta))
  }
  return(-final)
}

likeG = function(x, p){
  final = 0
  for (data in x){
    final = final+log((p**data)*(1-p))
  }
  return(final)
}
likeBG = function(x, theta){
  final = 0
  alpha = theta[1]
  beta = theta[2]
  for (data in x){
    final = final+log(beta(alpha+1, data+beta)/beta(alpha,beta))
  }
  return(-final)
}

LR = function(x){
  mleBG = nlm(loglike, c(1, 1), x=x)$estimate
  mleG = length(x)/(sum(x)+length(x))
  return(-2*(likeG(x,mleG)-likeBG(x, mleBG)))
}

alg1 = function(n, t){
  N = 0
  V = 0
  I = t+n-1
  k = c()
  while (N<n-1){
    p = runif(1)
    if (p<(n-1-N)/(t+n-1-N-V)){
      k = c(I,k)
      N = N+1
      I = I-1
    }
    else{
      V = V+1
      I = I-1
    }
  }
  return(k)
}

randombars = function(n,t){
  s1 = alg1(n, t)
  final = c(s1[1]-1)
  for (ind in 2:length(s1)){
    val = s1[ind]-ind-(s1[ind-1]-(ind-1))
    final = append(final, val)
  }
  final = append(final, t-(s1[length(s1)]-(n-1)))
  return (final)
}

vvv = 800
sig = 0.05
sig2 = 0.1
count2 = 0
count3 = 0
for (ooo in 1:vvv){
  
#Alternative
n = 25
d = rpois(n, 2)
t = sum(d)

mmm = 500
CO = LR(d)
C = c()

for (pp in 1:mmm){
  s = randombars(n,t)
  C[pp] = LR(s)
}

condpvalue = length(C[C > CO])/length(C)
if (condpvalue < sig){
  count2 = count2 + 1
}
if (condpvalue < sig2){
  count3 = count3 + 1
}
}

cat("significance level", sig, "power is", 100*count2/vvv)
cat("significance level", sig2, "power is", 100*count3/vvv)
