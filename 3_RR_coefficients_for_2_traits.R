rm(list=ls(all=TRUE))
#===============================================================================#
Phi = as.matrix(read.table("Legendre_polynomials.txt"))
H = as.matrix(read.table("3_RR_coefficients_for_2_traits.txt", header =F)) 
#===============================================================================#

H1 = H[1:3,1:3] # Block 1
H1
H2 = H[1:3,4:6] # Block 2, same as transpose of H[4:6,1:3]
H2
H3 = H[4:6,4:6] # Block 3
H3

HH1 = Phi %*%H1 %*% t(Phi)
HH1
HH2 = Phi %*%H2 %*% t(Phi)
HH2
HH3 = Phi %*%H3 %*% t(Phi)
HH3


temp1 = cbind(HH1,HH2)
temp1
temp2 = cbind(t(HH2),HH3) # Make sure to transpose Block 2
temp2
varCovr = rbind(temp1,temp2)
varCovr
correlation = cov2cor(varCovr)
correlation








