rm(list=ls(all=TRUE))
#===============================================================================#
RRcoefToVarComponent <- function(times,numCoeff, CoeffMatrix, LegendrePolynomials) {
  varCovr=matrix(NA, ncol=times*numOfTraits, nrow=times*numOfTraits)
  k1 = 1; k2 = times
  z1 = 1; z2 = numCoeff
  for (i in 1:(numOfTraits)){
    k3 = 1; k4 = times
    z3 = 1; z4 = numCoeff
    for (j in 1:(numOfTraits)){
      varCovr[(k1:k2),(k3:k4)] = LegendrePolynomials %*% CoeffMatrix[(z1:z2),(z3:z4)] %*% t(LegendrePolynomials)
      k3 = k3+times;    k4 = k4+times
      z3 = z3+numCoeff; z4 = z4+numCoeff
    }
    k1 = k1+times;    k2 = k2+times
    z1 = z1+numCoeff; z2 = z2+numCoeff
  }
  return(varCovr)
}

#===============================================================================#
Phi = as.matrix(read.table("Legendre_polynomials.txt"))
#===============================================================================#

# 1. One trait
RRcoeff = as.matrix(read.table("3_RR_coefficients_for_1_traits.txt", header = F))
#===============================================================================#
numCoeff = ncol(Phi)  # Number of coefficients (i.e., order of fit)
timePoints = nrow(Phi) # Number of time points (i.e, number of times of repeated records)
#===============================================================================#
numOfTraits =nrow(RRcoeff)/numCoeff # Number of traits
#===============================================================================#
varCovr = RRcoefToVarComponent(times=timePoints,numCoeff=numCoeff,CoeffMatrix=RRcoeff, LegendrePolynomials=Phi)
varCovr
correlation = cov2cor(varCovr)
correlation


# 2. Two traits
RRcoeff = as.matrix(read.table("3_RR_coefficients_for_2_traits.txt", header = F))
#===============================================================================#
numCoeff = ncol(Phi)
timePoints = nrow(Phi)
#===============================================================================#
numOfTraits =nrow(RRcoeff)/numCoeff
#===============================================================================#
varCovr = RRcoefToVarComponent(times=timePoints,numCoeff=numCoeff,CoeffMatrix=RRcoeff, LegendrePolynomials=Phi)
varCovr
correlation = cov2cor(varCovr)
correlation




