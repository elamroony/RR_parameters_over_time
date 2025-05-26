import numpy as np
import pandas as pd

# Define the function
def RRcoefToVarComponent(times, numCoeff, CoeffMatrix, LegendrePolynomials, numOfTraits):
    varCovr = np.full((times * numOfTraits, times * numOfTraits), np.nan)
    k1, k2 = 0, times
    z1, z2 = 0, numCoeff
    for i in range(numOfTraits):
        k3, k4 = 0, times
        z3, z4 = 0, numCoeff
        for j in range(numOfTraits):
            varCovr[k1:k2, k3:k4] = LegendrePolynomials @ CoeffMatrix[z1:z2, z3:z4] @ LegendrePolynomials.T
            k3 += times
            k4 += times
            z3 += numCoeff
            z4 += numCoeff
        k1 += times
        k2 += times
        z1 += numCoeff
        z2 += numCoeff
    return varCovr

# Function to convert var/covar to correlations
def cov2cor(cov):
    stddev = np.sqrt(np.diag(cov))
    outer_stddev = np.outer(stddev, stddev)
    correlation = cov / outer_stddev
    correlation[cov == 0] = 0
    return correlation

# ================= Examples ==========================================#
# ------------ Single trait ---------------------
# Read Legendre_polynomials
Phi = pd.read_csv("Legendre_polynomials.txt", sep='\\s+', header=None).values

# Read RR coefficients:
RRcoeff = pd.read_csv("3_RR_coefficients_for_1_traits.txt", sep='\\s+', header=None).values

# Define parameters:
timePoints = Phi.shape[0]
numCoeff = Phi.shape[1]
numOfTraits = RRcoeff.shape[0] // numCoeff  

# Compute variance-covariance matrix
varCovr = RRcoefToVarComponent(timePoints, numCoeff, RRcoeff, Phi, numOfTraits)
print("Variance-Covariance Matrix:\n", varCovr)

# Compute correlation matrix
correlation = cov2cor(varCovr)
print("\nCorrelation Matrix:\n", correlation)


#======================================================#
# ------------ Two traits ---------------------
# Read Legendre_polynomials
Phi = pd.read_csv("Legendre_polynomials.txt", sep='\\s+', header=None).values

# Read RR coefficients:
RRcoeff = pd.read_csv("3_RR_coefficients_for_2_traits.txt", sep='\\s+', header=None).values

# Define parameters:
timePoints = Phi.shape[0]
numCoeff = Phi.shape[1]
numOfTraits = RRcoeff.shape[0] // numCoeff  

# Compute variance-covariance matrix
varCovr = RRcoefToVarComponent(timePoints, numCoeff, RRcoeff, Phi, numOfTraits)
print("Variance-Covariance Matrix:\n", varCovr)

# Compute correlation matrix
correlation = cov2cor(varCovr)
print("\nCorrelation Matrix:\n", correlation)
# =====================================================================#