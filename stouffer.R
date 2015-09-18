load("ch8.RData")

require(stats4)

dist <- distSel
diag(dist) <- Inf

rankMat=matrix(0, ncol=dim(distSel)[2], nrow=dim(distSel)[1])

for(i in 1:952){
  for(j in 1:952){
    rankMat[i,j]=PuSel[j]/sum(PuSel[dist[i,j]>dist[i,]])
  }}

plot(rankMat[1,], dist[1,], log='xy')

diag(rankMat)=0

rankMat[!is.finite(rankMat)]=0

par <- c(theta=-15, tau1=0, tau2=-1)

ySel <- t(MuSel)/PuSel
fadeout <- ySel==0
S <- t(SuSel)
beta <- t(exp(betaSel))

relevant <- which(fadeout[,-546])
obs <- fadeout[,-1][relevant]

likfn <- function (theta, tau1, tau2) {
  theta <- exp(theta)
  tau1<-exp(tau1)
  tau2<-exp(tau2)
  iotaSel <- theta*(PuSel^tau1)*(rankMat^tau2%*%ySel)
  lambda <- beta*S*iotaSel^alpha
  prob <- exp(-lambda[relevant])
  #print(sum(!is.finite(prob)))
  -sum(log(ifelse(obs,prob,1-prob)))
}

test2<-mle(likfn, start=list(theta=-1, tau1=-3, tau2=-0),
           control=list(trace=4), method='Nelder-Mead')
ctest2<-confint(test2)
ptest2<-profile(test2)

# Maximum likelihood estimation

# Call:
# mle(minuslogl = likfn, start = list(theta = -1, tau1 = -3, tau2 = -0), 
#     method = "Nelder-Mead", control = list(trace = 4))

# Coefficients:
#          Estimate  Std. Error
# theta -6.43423775 0.048529532
# tau1  -0.02826678 0.005209889
# tau2  -0.48567701 0.010815780

# -2 log L: 241369.7 

# ctest2
#             2.5 %      97.5 %
# theta -6.52823647 -6.33790100
# tau1  -0.03878094 -0.01834258
# tau2  -0.50796697 -0.46544077
# > exp(ctest2)2
# Error: unexpected numeric constant in "exp(ctest2)2"
# > exp(ctest2)
#             2.5 %      97.5 %
# theta 0.001461581 0.001768009
# tau1  0.961961416 0.981824625
# tau2  0.601717648 0.627858302

# ptest2
# An object of class "profile.mle"
# Slot "profile":
# $theta
#             z par.vals.theta par.vals.tau1 par.vals.tau2
# 1  -2.5967115    -6.55924154   -0.01750674   -0.49391818
# 2  -2.0831918    -6.53424078   -0.01966610   -0.49235265
# 3  -1.5694884    -6.50924003   -0.02197279   -0.49125771
# 4  -1.0552887    -6.48423927   -0.02396540   -0.48902428
# 5  -0.5393694    -6.45923851   -0.02626909   -0.48792731
# 6   0.0000000    -6.43423775   -0.02826678   -0.48567701
# 7   0.4859472    -6.40923699   -0.03061391   -0.48463738
# 8   1.0027096    -6.38423624   -0.03273088   -0.48283539
# 9   1.5196009    -6.35923548   -0.03492434   -0.48130784
# 10  2.0357328    -6.33423472   -0.03728752   -0.48026782
# 11  2.5522818    -6.30923396   -0.03943161   -0.47844298
# 12  3.0693440    -6.28423320   -0.04173636   -0.47718894

# $tau1
#             z par.vals.theta par.vals.tau1 par.vals.tau2
# 1  -3.0179350    -6.31402788   -0.04437053   -0.49667129
# 2  -2.5107652    -6.33365753   -0.04168657   -0.49484979
# 3  -2.0021611    -6.35334517   -0.03900261   -0.49293026
# 4  -1.4918167    -6.37487048   -0.03631865   -0.49159026
# 5  -0.9798143    -6.39478644   -0.03363470   -0.48970809
# 6  -0.4658035    -6.41391666   -0.03095074   -0.48753675
# 7   0.0000000    -6.43423775   -0.02826678   -0.48567701
# 8   0.5600955    -6.45487308   -0.02558283   -0.48405156
# 9   1.0779393    -6.47630106   -0.02289887   -0.48253871
# 10  1.5971668    -6.49655127   -0.02021491   -0.48072934
# 11  2.1170122    -6.51674698   -0.01753096   -0.47870588
# 12  2.6391898    -6.53813157   -0.01484700   -0.47728593

# $tau2
#             z par.vals.theta par.vals.tau1 par.vals.tau2
# 1  -2.9116923    -6.47517373   -0.03350713   -0.51910853
# 2  -2.4395026    -6.46722942   -0.03273689   -0.51353661
# 3  -1.9598592    -6.46074260   -0.03182402   -0.50796469
# 4  -1.4727723    -6.45396342   -0.03094010   -0.50239277
# 5  -0.9773638    -6.44712455   -0.03007367   -0.49682085
# 6  -0.4776925    -6.44115557   -0.02915374   -0.49124893
# 7   0.0000000    -6.43423775   -0.02826678   -0.48567701
# 8   0.5570162    -6.42807101   -0.02729137   -0.48010509
# 9   1.0807845    -6.41981669   -0.02660503   -0.47453316
# 10  1.6167464    -6.41111401   -0.02590194   -0.46896124
# 11  2.1608197    -6.40504918   -0.02497111   -0.46338932
# 12  2.7139213    -6.39801314   -0.02409489   -0.45781740


# Slot "summary":
# Maximum likelihood estimation

# Call:
# mle(minuslogl = likfn, start = list(theta = -1, tau1 = -3, tau2 = -0), 
#     method = "Nelder-Mead", control = list(trace = 4))

# Coefficients:
#          Estimate  Std. Error
# theta -6.43423775 0.048529532
# tau1  -0.02826678 0.005209889
# tau2  -0.48567701 0.010815780

# -2 log L: 241369.7 
