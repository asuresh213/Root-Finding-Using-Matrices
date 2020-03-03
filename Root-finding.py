import numpy as np
import matplotlib.pyplot as plt
import random
from scipy import interpolate

# -------------------------------------------------------------------------------
# Finding zeroes for randomly generated polynomial
n = 15  # degree of polynomial
coeffs = [1]
f = 0
x = np.linspace(-10, 10, num=1000)
A = np.zeros((n, n))
# assigning random integer coefficients.
for i in range(n):
    coeffs.append(random.randrange(-10, 10))
print(coeffs)

# creating upper hessenberg matrix
for i in range(n):
    for j in range(n):
        if(j == n-1):
            A[i][j] = -coeffs[n-i]
        if(i == j+1):
            A[i][j] = 1

print(np.linalg.eig(A)[0])

# plotting the polynomial
for i in range(len(coeffs)):
    f += coeffs[i]*(x**(n-i))

plt.plot(x, f)
plt.show()

# ------------------------------------------------------------------------------

# Using taylor approximation of non-polynomial function to identify zeroes
# -------------------------------------------------------------------------------
coeffs = []
tempsum = 0

# Defining a factorial function


def fact(n):
    temp = 1
    for i in range(1, n+1):
        temp *= i
    return temp

# Getting the taylor polynomial


def p(x, nmax):
    tempsum = 0
    for n in range(nmax+1):
        nfact = fact(2*n)
        tempsum += (-1)**n * ((x**(2*n))/nfact)
    return tempsum

# Getting coefficients of the taylor polynomial (This couldve been avoided)


def get_coeffs(nmax):
    coeffs = []
    for n in range(nmax):
        if(n % 2 == 0):
            if(n % 4 == 0):
                coeffs.append(1/fact(n))
            else:
                coeffs.append(-1/fact(n))
        else:
            coeffs.append(0)
    return np.array(coeffs)


x = np.arange(-4*np.pi, 4*np.pi, 0.01)  # initalizing x values
y = np.cos(x)
t = p(x, 10)  # obtaining taylor polynomial for first 6 terms
degree = 150 # and obviously, with degree 12
coeffs = fact(degree)*get_coeffs(degree)  # making the polynomial monic
dim = len(coeffs)
C = np.zeros((dim, dim))
# creating the upper hessenberg matrix
for i in range(dim):
    for j in range(dim):
        if(j == dim-1):
            C[i][j] = -coeffs[i]
        if(i == j+1):
            C[i][j] = 1

ev = []
# solving for real eigenvalues
for i in np.linalg.eig(C)[0]:
    if(np.abs(i.imag) == 0):
        ev.append(i.real/np.pi)
        print("Zero", i.real, "Eigenvalue:", i.real/np.pi)

print(sorted(ev))
per_error = []
for i in ev:
    temp = int(i)+0.5
    tempcalc = 100*(np.abs(temp - i)/np.abs(temp))
    per_error.append(tempcalc)
print(np.average(per_error))

# plotting
plt.plot(x, y, label="cos(x)")
plt.plot(x, t, label="Approximation upto 12 terms")
plt.axis([-4*np.pi, 4*np.pi, -2.5, 2.5])
plt.legend()
plt.show()
# -------------------------------------------------------------------------------

# Rootfinding for the nth energy state of the Quantum Harmonic Oscillator
# -------------------------------------------------------------------------------

# Funtion to define all hermite polynomials


def Hermite(n, x):
    if(n == 0):
        return [1, [1]]
    if(n == 1):
        return [2*x, [1, 0]]
    if(n == 2):
        return [4*(x**2) - 2, [1, 0, (-1/2)]]
    if(n == 3):
        return [8*(x**3) - 12*x, [1, 0, (-12/8), 0]]
    if(n == 4):
        return [16*(x**4)-48*(x**2)+12, [1, 0, -48/16, 0, 12/16]]
    if(n == 5):
        return [32*(x**5) - 160*(x**3) + 120*x, [1, 0, -160/32, 0, 120/32, 0]]
    if(n == 6):
        return [64*(x**6) - 480*(x**4) + 720*(x**2) - 120, [1, 0, -480/64, 0, 720/64, 0 - 120/64]]
    if(n == 7):
        return [128*(x**7)-1344*(x**5) + 3360*(x**3) - 1680*x, [1, 0, -1344/128, 0, 3360/128, 0, -1680/128]]
    if(n == 8):
        return [256*(x**8) - 3584*(x**6) + 13440*(x**4) - 13440*(x**2) + 1680, [1, 0, -3584/256, 0, 13440/256, 0, -13440/256, 0, 1680]]
    if(n == 9):
        return [512*(x**9)-9216*(x**7)+48384*(x**5) - 80640*(x**3) + 30240*x, [1, 0, -9216/512, 0, 48384/512, 0, 80640/512, 0, 30240/512]]
    if(n == 10):
        return [1024*(y**10)-23040*(y**8)+161280*(y**6)-403200*(y**4)+302400*(y**2)-30240, [1, 0, -23040/1024, 0, 161280/1024, 0, -403200/1024, 0, 302400/1024, 0, -30240]]


# psi_n = 0 \implies H_n = 0. Solving for H_n = 0
coeffs = []
n = 10
y = []
x = np.linspace(-10, 10, num=1000000)
B = np.zeros((n, n))
hbar = 1  # because theoretical physics
m = 1  # normalized mass
omega = np.sqrt(1/m)  # k = 1 in the oscillator

# set up for psi_n
alpha = (m*omega)/hbar
y = np.sqrt(alpha)*x

# defining psi_n
psi_n = ((alpha/np.pi)**(1/4))*(1/np.sqrt((2**n)*np.math.factorial(n))) * \
    (Hermite(n, y)[0])*(np.e**(-0.5*y**2))

# obtaining the coefficients for the Hermite polynomials
coeffs = Hermite(n, y)[1]
# creating upper hessenberg matrix
for i in range(n):
    for j in range(n):
        if(j == n-1):
            B[i][j] = -coeffs[n-i]
        if(i == j+1):
            B[i][j] = 1
# solving for eigenvalues.
print(np.linalg.eig(B)[0])

# proof of concept! Very tedious to actually perform.
# interpolate through n data points to obtain n cubic functions.
# In each domain see if the resp. cubic function changes sign
# if it does, do rootfinding with the restricted domain.
# Repeat process for all cubic functions to obtain all roots.
# Can be optimized with further work.

#fitted = 0
#spl1 = interpolate.CubicSpline(y, psi_n)
#for i in range(4):
    #fitted += spl1.c[i, 500]*(y**(3-i))

# plotting
plt.plot(y, psi_n, label="original")
#plt.plot(y, fitted, label="proof of concept")
plt.axis([-10, 10, -1, 1])
plt.legend()
plt.show()
