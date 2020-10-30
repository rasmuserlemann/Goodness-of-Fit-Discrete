# NECESSARY R-PACKAGES

library(VGAM)
library(DiscreteWeibull)

###################################################################

#Drawing conditional simulations

simcount = function(n,t)
{
  barsvec=randombars(n,t)
  tell=c()
  mm=max(barsvec)
  for (k in 0:mm)
  {
    tell[k+1]=0
    for (i in 1:n)
    {
      tell[k+1] = tell[k+1] + 1*(barsvec[i]==k)
    }
  }
  return(tell)
}

# Draws from the bars and stars!
#
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

ct=function(d)
{
  ct=c()
  m = max(d)
  for (i in 0:m) ct[i+1]=length(d[d==i])
  return(ct) 
}


# Function that calculates r-th moment 
# of the data 'count', wher k is the maximum value of x

mom <- function(k,count,r)
{
  vec=c(0:k)
  momr = sum(vec^r*count)/sum(count)
  return(momr)
}


# TEST STATISTICS

# DATA SETS ARE USUALLY GIVEN BY VECTORS CALLED count

# Calculation of n, t and p-hat are done inside the functions.

# FUNCTIONS TO GO FROM count TO d AND VICE VERSA 

# Function that takes 'count' to 'd'

dt=function(count)
{
  m=length(count)-1
  y=rep(0,count[1])
  if (m > 0) 
  {
    for (i in 1:m) y=append(y,rep(i,count[i+1]))
  }
  return(y) 
}


# Function that takes 'count' to 'd'
dt=function(count)
{
  m=length(count)-1
  y=rep(0,count[1])
  return(y) 
}

# Function that takes 'd' to 'count' 
ct = function(dat){
  out = c()
  M = max(dat)
  for (v in 0:M)
    out = c(out, sum(dat == v))
  out
}

# FUNCTION FOR THE VECTOR Z FROM THE ARTICLE (NOT USED LATER
# SINCE IT IS INCORPORATED IN OTHER FUNCTIONS

Z <- function(count)
{
  m=length(count)
  n=sum(count) # number of obs
  vec=c(0:(m-1))
  t=sum(vec*count) #sum of xi
  phat = n/(n+t)
  O=cumsum(count)
  E=c()
  for (k in 1:m) E[k]=n*(1-(1-phat)^k)
  Z=O-E
  return(Z)
}

# THIS IS THE CRAMER VON MISES STATISTIC WHICH SUMS FROM 0 TO THE 
# UPPER NUMBER MU GIVEN IN THE PAPER 

CM <- function(count)
{
  m=length(count)
  n=sum(count) # number of obs
  vec=c(0:(m-1))
  t=sum(vec*count) #sum of xi
  phat = n/(n+t)
  O=cumsum(count)
  E=c()
  ph=c()
  for (k in 1:m) E[k]=n*(1-(1-phat)^k)
  for (k in 1:m) ph[k]=phat*(1-phat)^(k-1)
  Z=O-E
  mm = m
  while (ph[mm]*(1-phat) >= 0.001/n)
  {
    ph[mm+1] = ph[mm]*(1-phat)
    Z[mm+1] = Z[mm]-n*ph[mm+1]
    mm = mm+1
  }
  CM=sum(Z^2*ph)/n
  return(CM)
}


# THIS IS THE ANDERSON-DARLING STATISTIC WHICH SUMS FROM 0 TO THE 
# UPPER NUMBER MU GIVEN IN THE PAPER 

AD <- function(count)
{
  m=length(count)
  AD=0
  if (m > 1)
  {
    n=sum(count) # number of obs
    vec=c(0:(m-1))
    t=sum(vec*count) #sum of xi
    phat = n/(n+t)
    O=cumsum(count)
    E=c()
    ph=c()
    for (k in 1:m) E[k]=n*(1-(1-phat)^k)
    for (k in 1:m) ph[k]=phat*(1-phat)^(k-1)
    H=cumsum(ph)
    Z=O-E
    mm = m
    while (ph[mm]*(1-phat) >= 0.001/n)
    {
      ph[mm+1] = ph[mm]*(1-phat)
      Z[mm+1] = Z[mm]-n*ph[mm+1]
      mm = mm+1
    }
    H=cumsum(ph)
    AD=sum(Z^2*ph/(H*(1-H)))/n
  }
  return(AD)
}


# THIS IS THE KOLMOGOROV-SMIRNOV STATISTIC

KS <- function(count)
{
  m=length(count)
  n=sum(count) # number of obs
  vec=c(0:(m-1))
  t=sum(vec*count) #sum of xi
  phat = n/(n+t)
  O=cumsum(count)
  E=c()
  ph=c()
  for (k in 1:m) E[k]=n*(1-(1-phat)^k)
  for (k in 1:m) ph[k]=phat*(1-phat)^k
  Z=O-E
  maxz=0
  for (r in 1:m)
  {
    if (abs(Z[r]) > maxz) maxz = abs(Z[r])
  }
  return(maxz)
}


# THIS IS THE Rao statistic (FULL LIKELIHOOD APPROACH)

CR = function(count)
{
  d=dt(count)
  n=length(d)
  t = sum(d)
  phat=n/(n+t)
  sum1=0
  sum2 = 0
  for (i in 1:n)
  {
    if (d[i] != 0){
    sum1=sum1+d[i]*log(d[i])
    sum2 = sum2 + (d[i]+1)*log(d[i]+1)
    }
  }
  sc=sum1-sum2
  return(sc)
}

# THIS IS THE SCORE TEST FOR THE BETA-GEOMETRIC

SB <- function(count)
{
  k=length(count)-1
  m1=mom(k,count,1)
  m2=mom(k,count,2)
  m=length(count)
  score = (m2-m1-2*m1^2)
  return(score)
}


# THIS IS THE PREFERRED SCORE TEST FOR THE BETA-GEOMETRIC.
# IT IS CALLED SBabs, BUT CORRESPONDS TO THE SB0 IN THE PAPER
# (SEE LAST LINE)

SBabs <- function(count)
{
  k=length(count)-1
  m1=mom(k,count,1)
  m2=mom(k,count,2)
  m=length(count)
  score = (m2-m1-2*m1^2)
  return(max(score,0))
}


# THIS IS THE MOMENT ESTIMATOR FOR THETA IN THE BETA-GEOMETRIC
# RECALL THAT THETA=0 IF THE DISGTRIBUTION IS GEOMETRIC

theta <- function(count)
{
  k=length(count)-1
  theta=0
  if (k > 0)
  {
    m1=mom(k,count,1)
    m2=mom(k,count,2)
    theta=(m2-m1-2*m1^2)/(2*m2-m1^2+m1*m2)
  }
  return(theta)
}


# THIS IS The Weibull Score statistic - 
# NOTE THE ABSOLUTE VALUE IN LAST LINE

SW = function(count)
{
  d=dt(count)
  n=length(d)
  t = sum(d)
  phat=n/(n+t)
  sum1=0
  sum2 = 0
  for (i in 1:n)
  {
    if (d[i] != 0) sum1=sum1+d[i]*log(d[i])
    sum2 = sum2 + (1-phat)*(d[i]+1)*log(d[i]+1)
  }
  sc=sum1-sum2
  return(abs(sc))
}

# BELOW ARE THE ONE-SIDED TEST STATISTICS, 
# SWL IS FOR ALTERNATIVE BETA < 1, SWU IS FOR BETA > 1

SWL = function(count)
{
  d=dt(count)
  n=length(d)
  t = sum(d)
  phat=n/(n+t)
  sum1=0
  sum2 = 0
  for (i in 1:n)
  {
    if (d[i] != 0) sum1=sum1+d[i]*log(d[i])
    sum2 = sum2 + (1-phat)*(d[i]+1)*log(d[i]+1)
  }
  sc=sum1-sum2
  return(sc)
}

SWU = function(count)
{
  d=dt(count)
  n=length(d)
  t = sum(d)
  phat=n/(n+t)
  sum1=0
  sum2 = 0
  for (i in 1:n)
  {
    if (d[i] != 0) sum1=sum1+d[i]*log(d[i])
    sum2 = sum2 + (1-phat)*(d[i]+1)*log(d[i]+1)
  }
  sc=sum1-sum2
  return(-sc)
}


## THIS IS FOR P-VALUE CALCULATION FOR A GIVEN DATA VECTOR count
## NOTE THE NUMBER OF SIMULATIONS, 'ant' THAT HAS TO BE SET (E.G. 10 000)
## THE RESULT IS A VECTOR GIVING P-VALUES FOR THE 10 TEST STATISTICS 
## IN ORDER

pvalue = function(count,ant)
{
  m=length(count)
  n=sum(count) # number of obs
  vec=c(0:(m-1))
  t=sum(vec*count) #sum of xi
  W1=CM(count)
  W2=AD(count)
  W3=KS(count)
  W4=CR(count)
  W5=SB(count)
  W6=SBabs(count)
  W7=theta(count)
  W8=SW(count)
  W9=SWL(count)
  W10=SWU(count)
  tell1=0
  tell2=0
  tell3=0
  tell4=0
  tell5=0
  tell6=0
  tell7=0
  tell8=0
  tell9=0
  tell10=0
  for (i in 1:ant)
  {
    csim = simcount(n,t)
    if (CM(csim) >= W1) tell1=tell1+1
    if (AD(csim) >= W2) tell2=tell2+1
    if (KS(csim) >= W3) tell3=tell3+1
    if (CR(csim) >= W4) tell4=tell4+1
    if (SB(csim) >= W5) tell5=tell5+1
    if (SBabs(csim) >= W6) tell6=tell6+1
    if (theta(csim) >= W7) tell7=tell7+1
    if (SW(csim) >= W8) tell8=tell8+1
    if (SWL(csim) >= W9) tell9=tell9+1
    if (SWU(csim) >= W10) tell10=tell10+1
  }
  return(c(tell1,tell2,tell3,tell4,tell5,tell6,tell7,tell8,tell9,tell10)/ant)
}


## POWER SIMULATION FOR TESTS VERSUS BETA-GEOMETRIC.
## HERE a=alpha, b=beta FROM ARTICLE. level=significance level,
## number = number of samples from alternative distribution,
## antsim = number of simulations to get p-value for each data set

powerB = function(n,a,b,level,number,antsim)
{
  reject1=0
  reject2=0
  reject3=0
  reject4=0
  reject5=0
  reject6=0
  reject7=0
  reject8=0
  reject9=0
  reject10=0
  for (r in 1:number)
  {
    dat=rbetageom(n,a,b)
    n=length(dat)
    t=sum(dat)
    count=ct(dat)
    wobs1=CM(count)
    wobs2=AD(count)
    wobs3=KS(count)
    wobs4=CR(count)
    wobs5=SB(count)
    wobs6=SBabs(count)
    wobs7=theta(count)
    wobs8=SW(count)
    wobs9=SWL(count)
    wobs10=SWU(count)
    tell1=0
    tell2=0
    tell3=0
    tell4=0
    tell5=0
    tell6=0
    tell7=0
    tell8=0
    tell9=0
    tell10=0
    for (i in 1:antsim)
    {
      simd = simcount(n,t)
      wnull1= CM(simd)
      wnull2= AD(simd)
      wnull3= KS(simd)
      wnull4= CR(simd)
      wnull5= SB(simd)
      wnull6= SBabs(simd)
      wnull7= theta(simd)
      wnull8= SW(simd)
      wnull9= SWL(simd)
      wnull10= SWU(simd)
      if (wobs1 <= wnull1) tell1=tell1+1
      if (wobs2 <= wnull2) tell2=tell2+1
      if (wobs3 <= wnull3) tell3=tell3+1
      if (wobs4 <= wnull4) tell4=tell4+1
      if (wobs5 <= wnull5) tell5=tell5+1
      if (wobs6 <= wnull6) tell6=tell6+1
      if (wobs7 <= wnull7) tell7=tell7+1
      if (wobs8 <= wnull8) tell8=tell8+1
      if (wobs9 <= wnull9) tell9=tell9+1
      if (wobs10 <= wnull10) tell10=tell10+1
    }
    pval1 = tell1/antsim
    pval2 = tell2/antsim
    pval3 = tell3/antsim
    pval4 = tell4/antsim
    pval5 = tell5/antsim
    pval6 = tell6/antsim
    pval7 = tell7/antsim
    pval8 = tell8/antsim
    pval9 = tell9/antsim
    pval10 = tell10/antsim
    if (pval1 <= level) reject1=reject1+1
    if (pval2 <= level) reject2=reject2+1
    if (pval3 <= level) reject3=reject3+1
    if (pval4 <= level) reject4=reject4+1
    if (pval5 <= level) reject5=reject5+1
    if (pval6 <= level) reject6=reject6+1
    if (pval7 <= level) reject7=reject7+1
    if (pval8 <= level) reject8=reject8+1
    if (pval9 <= level) reject9=reject9+1
    if (pval10 <= level) reject10=reject10+1
  }
  power1=reject1/number
  power2=reject2/number
  power3=reject3/number
  power4=reject4/number
  power5=reject5/number
  power6=reject6/number
  power7=reject7/number
  power8=reject8/number
  power9=reject9/number
  power10=reject10/number
  return(c(power1,power2,power3,power4,power5,power6,power7,power8,power9,power10))
}



## POWER SIMULATION FOR TESTS VERSUS DISCRETE WEIBULL TYPE I.
## HERE q and beta ARE AS IN ARTICLE. level=significance level,
## number = number of samples from alternative distribution,
## antsim = number of simulations to get p-value for each data set

powerW = function(n,q,beta,level,number,antsim)
{
  reject1=0
  reject2=0
  reject3=0
  reject4=0
  reject5=0
  reject6=0
  reject7=0
  reject8=0
  reject9=0
  reject10=0
  for (r in 1:number)
  {
    dat=rdweibull(n,q,beta,zero=TRUE)
    n=length(dat)
    t=sum(dat)
    count=ct(dat)
    wobs1=CM(count)
    wobs2=AD(count)
    wobs3=KS(count)
    wobs4=CR(count)
    wobs5=SB(count)
    wobs6=SBabs(count)
    wobs7=theta(count)
    wobs8=SW(count)
    wobs9=SWL(count)
    wobs10=SWU(count)
    tell1=0
    tell2=0
    tell3=0
    tell4=0
    tell5=0
    tell6=0
    tell7=0
    tell8=0
    tell9=0
    tell10=0
    for (i in 1:antsim)
    {
      simd = simcount(n,t)
      wnull1= CM(simd)
      wnull2= AD(simd)
      wnull3= KS(simd)
      wnull4= CR(simd)
      wnull5= SB(simd)
      wnull6= SBabs(simd)
      wnull7= theta(simd)
      wnull8= SW(simd)
      wnull9= SWL(simd)
      wnull10= SWU(simd)
      if (wobs1 <= wnull1) tell1=tell1+1
      if (wobs2 <= wnull2) tell2=tell2+1
      if (wobs3 <= wnull3) tell3=tell3+1
      if (wobs4 <= wnull4) tell4=tell4+1
      if (wobs5 <= wnull5) tell5=tell5+1
      if (wobs6 <= wnull6) tell6=tell6+1
      if (wobs7 <= wnull7) tell7=tell7+1
      if (wobs8 <= wnull8) tell8=tell8+1
      if (wobs9 <= wnull9) tell9=tell9+1
      if (wobs10 <= wnull10) tell10=tell10+1
    }
    pval1 = tell1/antsim
    pval2 = tell2/antsim
    pval3 = tell3/antsim
    pval4 = tell4/antsim
    pval5 = tell5/antsim
    pval6 = tell6/antsim
    pval7 = tell7/antsim
    pval8 = tell8/antsim
    pval9 = tell9/antsim
    pval10 = tell10/antsim
    if (pval1 <= level) reject1=reject1+1
    if (pval2 <= level) reject2=reject2+1
    if (pval3 <= level) reject3=reject3+1
    if (pval4 <= level) reject4=reject4+1
    if (pval5 <= level) reject5=reject5+1
    if (pval6 <= level) reject6=reject6+1
    if (pval7 <= level) reject7=reject7+1
    if (pval8 <= level) reject8=reject8+1
    if (pval9 <= level) reject9=reject9+1
    if (pval10 <= level) reject10=reject10+1
  }
  power1=reject1/number
  power2=reject2/number
  power3=reject3/number
  power4=reject4/number
  power5=reject5/number
  power6=reject6/number
  power7=reject7/number
  power8=reject8/number
  power9=reject9/number
  power10=reject10/number
  return(c(power1,power2,power3,power4,power5,power6,power7,power8,power9,power10))
}


## POWER SIMULATION FOR TESTS UNDER NULL HYPOTHESIS
## HERE p IS THE PROBABILITY PARAMETER, level=significance level,
## number = number of initial samples,
## antsim = number of simulations to get p-value for each data set

powerG = function(n,p,level,number,antsim)
{
  reject1=0
  reject2=0
  reject3=0
  reject4=0
  reject5=0
  reject6=0
  reject7=0
  reject8=0
  reject9=0
  reject10=0
  for (r in 1:number)
  {
    dat=rgeom(n,p)
    n=length(dat)
    t=sum(dat)
    count=ct(dat)
    wobs1=CM(count)
    wobs2=AD(count)
    wobs3=KS(count)
    wobs4=CR(count)
    wobs5=SB(count)
    wobs6=SBabs(count)
    wobs7=theta(count)
    wobs8=SW(count)
    wobs9=SWL(count)
    wobs10=SWU(count)
    tell1=0
    tell2=0
    tell3=0
    tell4=0
    tell5=0
    tell6=0
    tell7=0
    tell8=0
    tell9=0
    tell10=0
    for (i in 1:antsim)
    {
      simd = simcount(n,t)
      wnull1= CM(simd)
      wnull2= AD(simd)
      wnull3= KS(simd)
      wnull4= CR(simd)
      wnull5= SB(simd)
      wnull6= SBabs(simd)
      wnull7= theta(simd)
      wnull8= SW(simd)
      wnull9= SWL(simd)
      wnull10= SWU(simd)
      if (wobs1 <= wnull1) tell1=tell1+1
      if (wobs2 <= wnull2) tell2=tell2+1
      if (wobs3 <= wnull3) tell3=tell3+1
      if (wobs4 <= wnull4) tell4=tell4+1
      if (wobs5 <= wnull5) tell5=tell5+1
      if (wobs6 <= wnull6) tell6=tell6+1
      if (wobs7 <= wnull7) tell7=tell7+1
      if (wobs8 <= wnull8) tell8=tell8+1
      if (wobs9 <= wnull9) tell9=tell9+1
      if (wobs10 <= wnull10) tell10=tell10+1
    }
    pval1 = tell1/antsim
    pval2 = tell2/antsim
    pval3 = tell3/antsim
    pval4 = tell4/antsim
    pval5 = tell5/antsim
    pval6 = tell6/antsim
    pval7 = tell7/antsim
    pval8 = tell8/antsim
    pval9 = tell9/antsim
    pval10 = tell10/antsim
    if (pval1 <= level) reject1=reject1+1
    if (pval2 <= level) reject2=reject2+1
    if (pval3 <= level) reject3=reject3+1
    if (pval4 <= level) reject4=reject4+1
    if (pval5 <= level) reject5=reject5+1
    if (pval6 <= level) reject6=reject6+1
    if (pval7 <= level) reject7=reject7+1
    if (pval8 <= level) reject8=reject8+1
    if (pval9 <= level) reject9=reject9+1
    if (pval10 <= level) reject10=reject10+1
  }
  power1=reject1/number
  power2=reject2/number
  power3=reject3/number
  power4=reject4/number
  power5=reject5/number
  power6=reject6/number
  power7=reject7/number
  power8=reject8/number
  power9=reject9/number
  power10=reject10/number
  return(c(power1,power2,power3,power4,power5,power6,power7,power8,power9,power10))
}

# EXAMPLE TO SHOW THAT THE SIZE (TYPE I ERROR PROBABILITIES CAN BE MUCH 
# LESS THAN =0.05 IF n IS SMALL, AND APPROX EQUAL TO 0.05 IF N IS LARGE.
# NOTE THAT THE RESULTS 0.057 > 0.05 ARE DUE TO MONTE CARLO ERROR. THEORETICALLY
# THEY ARE < 0.05 (BUT CLOSE IF n IS LARGE)

#Example
powerG(25,0.5,0.1,100,100)