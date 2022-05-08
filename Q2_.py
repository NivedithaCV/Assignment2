import utilities as ul
import numpy as np
r1,it=ul.rand(1000,65,0.0,1021,30)
r2,it_=ul.rand(3,65,10,16381,30)

def pi_calculation(n):
    first,it = ul.rand(3,3007,654946,18,n)
    second,it__= ul.rand(5,5003,956755,21,n)
    k_i=0
    for i in range (n):
        p=[]
        p.append(first[i])
        p.append(second[i])
        k=first[i]-0.5
        l=second[i]-0.5
        norm=((k*k)+(l*l))**(0.5)
        if norm<0.5:
            k_i+=1
    pi=4*((k_i)/(len(first)))
    return pi
pi=pi_calculation(1000000)
print("pi by points :",pi)

def fun(x):
    return np.sqrt(1-(x**2))
Pi=4*(ul.M_C(fun,10000))
print("pi by Monte carlo method is :",Pi)
