
import utilities as ul
import numpy as np
r1,it=ul.rand(1000,65,0.0,1021,30)
r2,it_=ul.rand(3,65,9,15373,30)
def f2(x):
    k=4*(1-x**2)
    return k
v=ul.M_C(f2,10000)
volume=2*v
print("Volume by Monte Carlo merthod:", volume)
