#From Bolker book (pg 192):
#
#		2(negLL_restricted - negLL_full) ~ X^2

# Reduced base model
nLL0 = 89.1812
nLL0_p=56

# Model with extra params
nLL1 =	82.0095
nLL1_p = 58	

D=2*(nLL0-nLL1)
degF=nLL1_p-nLL0_p

prob=dchisq(D, df= degF)
prob

1-pchisq(D,2)

##Here's another way to compare two models from R with pos LL's:

logLik(model.OLS)
'log Lik.' -273.5342 (df=3)
logLik(model1)
'log Lik.' -115.7798 (df=4)
LR1 <- 2*(logLik(model1)-logLik(model.OLS))[1]
1-pchisq(LR1,1)
