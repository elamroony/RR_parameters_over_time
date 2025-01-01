rm(list=ls(all=TRUE))
#===============================================================================#
Phi = as.matrix(read.table("Legendre_polynomials.txt"))
H = as.matrix(read.table("3_RR_coefficients_for_1_traits.txt", header =F)) 
#===============================================================================#

varCovr <- Phi %*% H %*% t(Phi)
varCovr
#The additive genetic variance:
diag(varCovr)

correlation = cov2cor(varCovr)
correlation


