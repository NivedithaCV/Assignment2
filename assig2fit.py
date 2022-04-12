import numpy as np
A=np.loadtxt("assign2fit.txt")

import math
import pandas as pd
import utilities as li
import numpy as np
import matplotlib.pyplot as plt
File_data = np.loadtxt("msfit.txt")
# A= li.readfile(File_data)
#t, N, dN = np.loadtxt("msfit.txt", unpack=True)

A=[[1,	106	,10],
[15,	80,	9],
[30	,98,	10],
[45,	75,	9],
[60	,74,	8],
[75,	73,	8],
[90	,49,	7],
[105,	38,	6],
[120,	37,	6],
[135,	22,	5]]

A=np.asarray(File_data)
x, y,dN = A[:, 0], A[:, 1] ,A[:,2]
X=x
Y= np.log(y)
dY= dN/y
B, A, dBn, dA,cov = li.LineFitWt(X, Y, dY)
redchisqr = li.redchisq(X, Y, dY, A, B)
#error
