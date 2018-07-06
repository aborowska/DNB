
setwd('C:/Users/aba228/Dropbox/New Projects/Integer/Codes/DNBSkellamCode/src')

set.seed(123)
library("stochvol")
#  data("exrates")
#  ret <- logret(exrates$USD, demean = TRUE)
ret <- read.csv('OutPut/__DNB_SimulatedData.csv')
vol_true <- read.csv('__DNB_SimulatedLogInt.csv')
par(mfrow = c(3, 1))
plot(ret[,1], type = "l" )
plot(ret[,2], type = "l" )
plot(vol_true)
type = "l" )


ret <- ret[,1] 
ret_short <- ret[1:1000] 
# sim <- svsim(500, mu = -9, phi = 0.99, sigma = 0.1)
# par(mfrow = c(1, 1))
# plot(sim)

res <- svsample(ret_short, priormu = c(-10, 1), priorphi = c(20, 1.1), priorsigma = 0.1)
summary(res, showlatent = FALSE)

volplot(res, forecast = 0)