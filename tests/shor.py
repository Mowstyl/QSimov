# Functions related with Shor algorithm

import numpy as np
from math import sqrt
from cmath import exp, pi

# Function f takes N and a as parameters and returns a function f_a,N(x) that returns a^x Mod N
def f(N, a):
	return lambda x: pow(a, x, N)

# Now using map we can apply a function f(N, a) to a list of values [x1, ..., xn].
# i. e.: list(map(f(15, 2), [0, 1, 2, 3]))

# We are trying to find the period of a function f(N, a). We are going to define a function that tries to find the period.
# We call period to the smallest r > 0 with a^r Mod N = 1.
# We won't use map because we don't need to store all the results in memory, we only need to check them one by one.
def findPeriod(N, a):
	f_aN = f(N, a)
	for i in range(1, N):
		if f_aN(i) == 1:
			return i
	return None

# Function that returns an array of nth roots of the unity
def nroots(n, rc = 14): # Rounds to 14 decimal places by default
	c = 2j * pi / n
	gen = (exp(k * c) for k in range(n))
	return [round(i.real, rc) + round(i.imag, rc) * 1j for i in gen]

# Function that creates a vandermonde matrix for the given values
def vandermonde(*args):
	return np.array([[pow(arg, i) for i in range(len(args))] for arg in args])


# Function that creates a vandermonde matrix of the given size with the roots of the unity
def vanderomega(size, rc = 14):
	return vandermonde(*nroots(size, rc))

def DFT(size, rc = 14):
	return 1/sqrt(size) * vanderomega(size, rc)
