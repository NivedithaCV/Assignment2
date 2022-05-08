 # import numpy as np
#
#
# def read_txt(path,sym):
# # ti specify the path use '/' eg: c:/sample/reaadme.txt
# # mode a: appending text, r: reading txt, w: writing text
#     data= np.loadtxt(path, dtype= 'str' ,delimiter=sym )
#
#  area
import random

import numpy as np
import math
def rand(seed,a,c,m,n):

    rfun = np.zeros(n)
    rfun[0] = seed
    itr=[]
    itr.append(0)
    for i in range(1, n):
        rfun[i] = ((a*rfun[i-1]) + c) % m
        itr.append(i)
    return rfun/m,itr
def Matrix_read(Z,Y):
    M=[]
    N=[]
    A=[[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]]
    B=[[0],[0],[0],[0],[0],[0]]
    fhand=open(Z)
    ghand=open(Y)
    for line in fhand:
        line=line.rstrip()
        li=list(line.split(","))
        c=len(li)
        M.append(li)
    for k in ghand:
        k=k.rstrip()
        lis=list(k.split(","))
        d=len(lis)
        N.append(lis)
    r=len(M)

    for i in range(r):
        B[i][0]=int(N[i][0])
        for j in range(r):
            A[i][j]=int(M[i][j])
    A=np.asarray(A)
    B=np.asarray(B)
    return(A,B)

def readfile(f1):
    a = f1.read()
    a1 = [item.split('  ') for item in a.split('\n')]
    A = []
    for i in range(len(a1)):
        row = []
        for j in range(len(a1[0])):
            row.append(0)
        A.append(row)
    for i in range(len(a1)):
        for k in range(len(a1[0])):
            A[i][k] = float(a1[i][k])
    return A

def read1d(p):
    b = p.read()
    b1 = (b.split(' '))
    C=[]
    for i in range(len(b1)):
        C.append(0)
    for i in range(4):
        C[i] = float(b1[i])
    return C

def matmul(M,A):
    B = []
    for i in range(len(M)):
        row = []
        for j in range(len(A[0])):
            row.append(0)
        B.append(row)

    for x in range(len(M)):
        for y in range(len(A[0])):
            for z in range(len(M[0])):
                B[x][y] += M[x][z] * A[z][y]
    return B

#forward substitution - LY= B. Find for Y. For linear Eqn
# def fwrdsub(A , B):
#     global Y
#     Y = []
#     for k in range(len(A)):
#         Y.append(float(0))
#     for i in range(0, len(A)):
#         val = 0
#         for j in range(0,i):
#             val+=A[i][j]*Y[j]
#         Y[i] += (B[i] - val)
#     return Y

#partial pivoting
def partial_pivot(A, B, k):
    if np.abs(A[k][k]) < 1e-10:
        n = len(B)
        for i in range(k + 1, n):
            if abs(A[i][k]) > abs(A[k][k]) and abs(A[k][k]) == 0:
                A[k], A[i] = A[i], A[k]
                B[k], B[i] = B[i], B[k]
    return A, B

# def fwrdsub2D(L):
#     Y = [[0 for x in range(len(A))]
#          for y in range(len(A))]
#     I = np.identity(len(A))
#     n = len(L)
#     for i in range(n):
#         Y[i][i] = 1/L[i][i]
#     for i in range(1,n):
#         for j in range(i):
#             sum = 0
#             for k in range(j,i):
#                 sum += L[i][k]*Y[k][j]
#             Y[i][j] = -(sum)/L[i][i]
#     return Y
#
#backward substitution- UX = Y find for X.
# def bkwdsub(A , B):
#     global X
#     X = []
#     for k in range(len(A)):
#         X.append(float(0))
#     for i in reversed(range(len(A))):
#         val = 0
#         for j in reversed(range(0, len(A))):
#             if j > i:
#                 val += A[i][j]*X[j]
#         X[i] += (1/A[i][i])*(B[i] - val)
#     return X
#
# def LU_decomp (A, C):
#     #A=LU
#     for j in range(len(A)):
#         partial_pivot(A, C, j)
#         for i in range(len(A)):
#         #elements of upper triangle
#
#             if i <= j:
#                 sum = 0
#                 for k in range(i):
#                         sum += A[i][k] * A[k][j]
#                 A[i][j] = A[i][j] - sum
#             #elements of lower triangle
#
#             if i > j:
#                 sum = 0
#                 for k in range(j):
#                         sum += A[i][k] * A[k][j]
#                 A[i][j] = (1/A[j][j])*(A[i][j]-sum)
#     return A
#
#
#Gauss_Jordan elemination
def GaussJordan(A, B):
    n = len(B)
    #Aug = np.hstack([A,B.T])
    for k in range(n):
        partial_pivot(A, B, k)
        print(A,B)
        pivot = A[k][k]
        for i in range(k, n):
            A[k][i] = A[k][i]/pivot
        B[k] = B[k] / pivot

        for i in range(n):
            if abs(A[i][k]) < 1e-10 or i == k:
                continue
            else:
                term = A[i][k]
                for j in range(k, n):
                    A[i][j] = A[i][j] - term * A[k][j]
                B[i] = B[i] - term * B[k]
    return A,B

def partial_pivot(A, B, k):
    if np.abs(A[k][k]) < 1e-10:
        n = len(B)
        for i in range(k + 1, n):
            if abs(A[i][k]) > abs(A[k][k]) and abs(A[k][k]) == 0:
                A[k], A[i] = A[i], A[k]
                B[k], B[i] = B[i], B[k]
    return A, B


def cholesky(A):
    n = len(A)
    L = [[0 for x in range(n)]
         for y in range(n)]

    for i in range(len(A)):
        sum = 0
        for j in range(i):
            sum += L[i][j]**2
        L[i][i] = math.sqrt(A[i][i] - sum)

        for k in range(i+1,n):
            sum = 0
            for j in range(i):
                sum += L[k][j]*L[i][j]
            L[k][i] = (1/L[i][i])*(A[k][i] - sum)
    return L

# A = [[4, 12, -16],
#      [12, 37, -43],
#      [-16, -43, 98]]
# Y = cholesky(A)
# # print(Y)
# S = fwrdsub2D(Y)
# # print("Result is: ",S)
def gauss_sidel(A,B):
    N = len(A)
    x = [0 for x in range(len(B))]
    for i in range(1, N+1):
        x_new = [0 for x in range(len(B))]
        for i in range(A.shape[0]):
            s1 = np.dot(A[i, :i], x_new[:i])
            s2 = np.dot(A[i, i + 1:], x[i + 1:])
            x_new[i] = (B[i] - s1 - s2) / A[i, i]
        if np.linalg.norm(np.subtract(x,x_new)) < 1e-5:
            break
        x = x_new
    return x


def gauss_seidel(A, b, tolerance, max_iterations=10000):

    x = np.zeros_like(b, dtype=np.double)
    er = np.array([])
    kk = np.array([])
    # Iterate
    for k in range(max_iterations):
        x_old = x.copy()

        # Loop over rows
        for i in range(A.shape[0]):
            x[i] = (b[i] - np.dot(A[i, :i], x[:i]) - np.dot(A[i, (i+1):], x_old[(i+1):])) / A[i, i]
        er = np.append(er, np.linalg.norm(x - x_old, ord=np.inf) / np.linalg.norm(x, ord=np.inf))
        kk = np.append(kk, k)
        # Stop condition
        if np.linalg.norm(x - x_old, ord=np.inf) / np.linalg.norm(x, ord=np.inf) < tolerance:
            break

    return x, er, kk
def guass_s_inverse(A,tol):
    I = np.zeros((len(A), len(A)))  # Pre-allocate matrix
    for i in range(len(A)):
        b=[[0] for i in range(len(A))]
        b[i][0]=1
        column = gauss_seidel(A, b, tol)
        for j in range(1, len(b)):
            I[:, i] = column[0]


def conjugate_gradient(A,B,x):
    # residual
    r = B - np.matmul(A,x)
    p = r
    rs_old = np.matmul(np.transpose(r),r)
    for i in range(len(B)):
        Ap = np.matmul(A,p)
        alpha = rs_old/(np.matmul(np.transpose(p),Ap))
        x = x + alpha*p
        r = r - alpha*Ap
        rs_new = np.matmul(np.transpose(r),r)
        if np.sqrt(rs_new) < 1e-10:
            break
        else:
            p = r + (rs_new/rs_old)* p
            rs_old = rs_new
    return x
# A = np.array([[4, 1],[1, 3]])
# B = np.array([1, 2])
# x = np.array([2, 1])
# print(conjugate_gradient(A,B,x))
# Jacobi Method - All eigenvalues. For numpy arrays
def Jacobi(A):
    # largest off-diag element
    n = len(A)
    def maxind(A):
        Amax = 0
        for i in range(n-1):
            for j in range(i+1,n):
                if abs(A[i,j])>=Amax:
                    Amax = abs(A[i,j])
                    k = i
                    l = j
        return Amax, k,l

    # to make A[k,l] = 0 by rotation and define rotation matrix
    def rotate(A,p,k,l):
        A_diff = A[l,l]-A[k,k]
        if abs(A[k,l])< abs(A_diff)*1e-30:
            t = A[k,l]/A_diff
        else:
            phi = A_diff/(2*A[k,l])
            t = 1/(abs(phi)+np.sqrt(phi**2+1))
            if phi<0:
                t = -t
        c = 1/np.sqrt(t**2+1)
        s = t*c
        tau = s/(1+c)

        term = A[k,l]
        A[k,l] = 0
        A[k,k] = A[k,k] - t*term
        A[l,l] = A[l,l] + t*term
        for i in range(k):
            term = A[i,k]
            A[i,k] = term - s*(A[i,l] + tau*term )
            A[i,l] += s*(term- tau*A[i,l])
        for i in range(k+1,l):
            term = A[k,i]
            A[k,i] = term - s*(A[i,l] + tau*A[k,i])
            A[i,l] += s*(term - tau*A[i,l])
        for i in range(l+1, n):
            term = A[k,i]
            A[k,i] = term - s*(A[l,i] + tau*term)
            A[l,i] += s*(term - tau*A[l,i])
        for i in range(n):
            term = p[i,k]
            p[i,k] = term - s*(p[i,l] - tau*p[i,k])
            p[i,l] += s*(term - tau*p[i,l])

    p = np.identity(n)
    for i in range(4*n**2):
        Amax, k,l = maxind(A)
        if Amax < 1e-9:
            return np.diagonal(A)
        rotate(A,p,k,l)
    print("This method did not converge")

A = np.array([[1,np.sqrt(2),2],[np.sqrt(2),3,np.sqrt(2)],[2,np.sqrt(2),1]])
# print("Jacobi is:",Jacobi(A))




def jacobi_inv(a, b, x):
    # here x is the initial guess of x
    # here k is the iteration limit
    m = len(a)
    n = len(b)
    X = np.array([0.0 for i in range(5)])
    print(X)
    # X = np.zeros_like(b)
    k = 15
    for z in range(k):
        for i in range(n):
            sum = 0.0
            for j in range(n):
                if i != j:
                    sum += ((a[i][j]) * x[j])
            X[i] = (1 / (a[i][i])) * (b[i] - sum)

            if X[i] == X[i - 1]:
                break
        if np.allclose(x, X, atol=1e-10, rtol=0.):
            break
        x = X
    return X



def power_eig(A,x,n):
    e_value_=0
    for i in range(n):
        x=np.dot(A,x)

        x= x/ np.linalg.norm(x)
        e_value = np.dot(x.T, np.dot(A, x.T))/np.dot(x.T,x)
        error= np.abs((e_value-e_value_)/e_value)
        i=i+1
        e_value_= e_value
    return(x,e_value)

def power_method_for_non_dominant(A,x,y,n,m):
    #m is the number of non-dominant eigen vectors needed
    eigval=[]
    eigvec=[]
    for i in range(m):
        eigvalue, eigvect = power_eig(A, x, n)
        A=A-(eigvalue*(np.matmul(eigvect, (np.transpose(eigvect)))))
        eigval.append(eigvalue)
        eigvec.append(eigvect)
    return eigval,eigvec

def jackknife(x, func):
    """Jackknife estimate of the estimator func"""
    n = len(x)
    idx = np.arange(n)
    return np.sum(func(x[idx!=i]) for i in range(n))/float(n)
def jackknife_var(x, func):
    """Jackknife estiamte of the variance of the estimator func."""
    n = len(x)
    idx = np.arange(n)
    j_est = jackknife(x, func)
    return (n-1)/(n + 0.0) * np.sum((func(x[idx!=i]) - j_est)**2.0
                                    for i in range(n))
def LineFitWt(x, y, sig):
    """
    Fit to straight line.
    Inputs: x and y arrays and uncertainty array (unc) for y data.
    Ouputs: slope and y-intercept of best fit to data.
    """
    sig2 = sig**2
    x2=x**2
    S = (1./sig2).sum()
    S_x = (x/sig2).sum()
    S_y = (y/sig2).sum()
    S_xx= (x2/sig2).sum()
    S_xy=((x*y)/sig2).sum()
    Del = (S*S_xx)-(S_x**2)
    slope = ((S_xx*S_y)-(S_x*S_xy))/Del
    int = ((S*S_xy)-(S_x*S_y))/Del
    sig2_slope = S_xx/Del
    sig2_int = S/Del
    Cov = -S_x/Del
    #r2= S_xy/(S_xx*S_xy)
    return slope, int, np.sqrt(sig2_slope), np.sqrt(sig2_int),Cov#,r2

def redchisq(x, y, dy, slope, yint):
    chisq = (((y-yint-slope*x)/dy)**2).sum()
    return chisq/float(x.size-2)
import cmath
def compute_dft_complex(input):
	n = len(input)
	output = []
	for k in range(n):  # For each output element
		s = complex(0)
		for t in range(n):  # For each input element
			angle = 2j * cmath.pi * t * k / n
			s += input[t] * cmath.exp(-angle)
		output.append(s)
	return output


#####

# A=1,0,1,2         B=6
#   0,1,-2,0         -3
#   1,2,-1,0         -2
#   2,1,3,-2         0


def L_Udec(A,c_A):
    for j in range(c_A):
        for i in range(len(A)):

            #diagonal
            if i==j:
                sum=0
                for u in range(i):
                    sum=sum+A[i][u]*A[u][i]
                A[i][i]=A[i][i]-sum

                #elements of upper triangle
            if i<j:
                sum=0
                for k in range(i):
                    sum=sum+A[i][k]*A[k][j]
                A[i][j]=A[i][j]-sum

                #elements of lower triangle
            if i>j:
                sum=0
                for z in range(j):
                    sum=sum+A[i][z]*A[z][j]
                A[i][j]=(A[i][j]-sum)/A[j][j]
    return(A)

def forw_backw(A,B,col):
    for k in range(col):
        r=len(A)

        #forward substitution
        Y=[[0] for y in range(r)]
        for i in range(r):
            sum=0
            for k in range(i):
                sum=sum+A[i][k]*Y[k][0]
            Y[i][0]=B[i][0]-sum
        print("matrix Y",Y)

        #backward substitution
        X=[[0] for w in range(r)]
        for l in range(r-1,-1,-1):
            sum=0
            for m in range(l+1,r):
                sum=sum+A[l][m]*X[m][0]
            X[l][0]=(Y[l][0]-sum)/A[l][l]
    print("solution matrix is",X)


def jk_resampling(data):
    n = data.shape[0]
    for i in range(n):
        resamples[i] = np.delete(data, i)

    return resamples

def jackknife_stats(data, confidence_level=0.95):
    N=data.shape[0]
    samples=jk_resampling(data)
    mean_sample = samples.mean(axis=1)
    fk_mean= np.mean(mean_sample)
    for i in range(len(mean_sample)):
        s=s+(mean_sample[i]-fk_mean)**2
    var2= (N-1*s)/N
    error2 =var2/(N-1)
    return fk_mean, var2 ,error2


#Takes an nxn matrix and an 1xn vector
#and solves for x by iterating until the
#given tolerance (tol) is met. tol_type
#is a single character to indicate absolute
#or relative convergence ('a' or 'r').
#Prints the solution and required iterations
#to meet the tolerance.

def go(A, b, tol_type, tol):

    shape = np.shape(A)
    m = shape[0]
    n = shape[1]

    if m != n:
        exit(1)

    if m != np.shape(b)[0]:
        exit(1)


    if tol_type == 'a':
        numpy_solution = np.linalg.solve(A, b);

    old_x = np.zeros(np.shape(b))
    new_x = np.zeros(np.shape(b))
    prev_x = np.zeros(np.shape(b))
#We'll hold the difference (x(i) - x) here.
    diff = np.zeros(np.shape(b))

    error = tol + 1;

    num_iterations = 0
    while error > tol:
        num_iterations += 1
        for i in range(0, m):
            prev_x[i] = new_x[i]
            sum = b[i]
            for j in range(0, m):
                if i != j:
                    sum = sum - A[i][j]*old_x[j][0]
            sum = sum / A[i][i]
            new_x[i][0] = sum
        old_x = new_x
        if tol_type == 'a':
            diff = np.subtract(new_x, numpy_solution)
            error = np.linalg.norm(diff) / np.linalg.norm(new_x)
        if tol_type == 'r':
            diff = np.subtract(new_x, prev_x)
            error = np.linalg.norm(diff) / np.linalg.norm(new_x)
    print(new_x)

    return
def M_C(f, n):
    x,itr= rand(3,687,2,186845,n)
    x=np.array(x)
    sum = 0.0
    for i in range(n):
        sum += f(x[i])
    sum/=n
    return sum
