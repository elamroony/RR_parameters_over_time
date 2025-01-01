import pandas as pd
import os
import numpy as np
#==========================================================================#
# Read Legendre Polynomials
LegendrePolynomials = pd.DataFrame([
[0.707107,   -1.224745,    1.581139],
[0.707107,   -0.612372,   -0.197642],
[0.707107,    0.000000,   -0.790569],
[0.707107,    0.612372,   -0.197642],
[0.707107,    1.224745,    1.581139]])
#==========================================================================#
# EBVs at each time point:
coeff=pd.DataFrame([6.33, 4.08,-1.83])
coeff
EBVs = LegendrePolynomials @ coeff
EBVs
#==========================================================================#

