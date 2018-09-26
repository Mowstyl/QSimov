# -*- coding: utf-8 -*-

import cmath as cm
import numpy as np
from structures.qregistry import *
from structures.qgate import *
from structures.qcircuit import *

# np.zeros((h,w), dtype=complex) Inicializa una matriz de numeros complejos con alto h y ancho w
# La suma de matrices se realiza con +. A + B
# La multiplicacion por un escalar se hace con *. n * A
# Para multiplicar las matrices A y B se usa np.dot(A,B)
# El producto Kronecker de A y B esta definido con np.kron(A,B)


def _hMat(n): # Devuelve una matriz que al ser multiplicada por 1/sqrt(2^n) resulta en la puerta Hadamard para n bits
	H = np.ones((2,2), dtype=complex)
	H[1,1] = -1
	if n > 1:
		H = np.kron(H, _hMat(n - 1))
	return H

def H(n): # Devuelve una puerta Hadamard para n QuBits
	H = QGate("H")
	H.addLine(_hMat(n))
	H.setMult(1 / np.sqrt(2**n))
	return H

def PauliX(): # Also known as NOT
	px = QGate("NOT")
	m = np.array([0,1,1,0], dtype=complex)
	m.shape = (2,2)
	px.addLine(m)
	return px

def PauliY():
	py = QGate("Y")
	m = np.array([0,-1j,1j,0], dtype=complex)
	m.shape = (2,2)
	py.addLine(m)
	return py

def PauliZ():
	pz = QGate("Z")
	m = np.array([1,0,0,-1], dtype=complex)
	m.shape = (2,2)
	pz.addLine(m)
	return pz

def SqrtNOT(): # Square root of NOT gate, usually seen in its controlled form C-√NOT. Sometimes called C-√X gate.
	v = QGate("√NOT")
	m = np.array([1, -1j, -1j, 1], dtype=complex)
	m.shape = (2,2)
	v.addLine(m)
	v.setMult((1 + 1j)/2)
	return v

def ControlledU(gate): # Returns a controlled version of the given gate
	g = gate
	name = "U"
	if type(gate) == QGate:
		g = gate.m
		name = gate.name
	gdim = g.shape[0]
	m = np.eye(gdim*2, dtype=complex)
	m[gdim:,gdim:] = g
	cu = QGate("C-" + name)
	cu.addLine(m)
	return cu

def CNOT(): # Returns a CNOT gate for two QuBits, also called Feynman gate
	#return ControlledU(PauliX())
	cn = QGate("C-NOT")
	m = np.zeros((4,4), dtype=complex)
	m[0,0] = 1
	m[1,1] = 1
	m[2,3] = 1
	m[3,2] = 1
	cn.addLine(m)
	return cn

def NOTC(): # Returns a CNOT gate for two QuBits, first QuBit objective and second one control
	#return SWAP() @ CNOT() @ SWAP()
	nc = QGate("NOT-C")
	m = np.zeros((4,4), dtype=complex)
	m[0,0] = 1
	m[3,1] = 1
	m[2,2] = 1
	m[1,3] = 1
	nc.addLine(m)
	return nc

def SWAP(): # SWAP gate for 2 qubits
	sw = QGate("SWAP")
	#m = np.zeros((4,4), dtype=complex)
	#m[0,0] = 1
	#m[1,2] = 1
	#m[2,1] = 1
	#m[3,3] = 1
	#sw.addLine(m)
	sw.addLine(CNOT())
	sw.addLine(NOTC())
	sw.addLine(CNOT())
	return sw

def SqrtSWAP(): # Square root of SWAP gate for 2 qubits
	sw = QGate("√SWAP")
	m = np.zeros((4,4), dtype=complex)
	m[0,0] = 1
	m[1,1] = 0.5 * (1+1j)
	m[1,2] = 0.5 * (1-1j)
	m[2,1] = 0.5 * (1-1j)
	m[2,2] = 0.5 * (1+1j)
	m[3,3] = 1
	sw.addLine(m)
	return sw

def Toffoli(): # Returns a CCNOT gate for three QuBits. A, B, C -> P = A, Q = B, R = AB XOR C.
	''' # This does the same as the line below. Circuit with the implementation of Toffoli gate using SWAP, CNOT, Controlled-SNot and Controlled-SNot+
	# Gates needed (without control SWAPs): 5
	gate = np.kron(I(1), ControlledU(V()))
	gate = np.dot(gate, np.kron(SWAP(), I(1)))
	gate = np.dot(gate, np.kron(I(1), ControlledU(V())))
	gate = np.dot(gate, np.kron(CNOT(), I(1)))
	gate = np.dot(gate, np.kron(I(1), ControlledU(Dagger(V()))))
	gate = np.dot(gate, np.kron(CNOT(), I(1)))
	gate = np.dot(gate, np.kron(SWAP(), I(1)))
	return gate
	'''
	return ControlledU(CNOT())

def Fredkin(): # Returns a CSWAP gate for three QuBits
	return ControlledU(SWAP())

def Deutsch(angle): # Returns Deutsh gate with specified angle. D(pi/2) = Toffoli
	d = np.eye(8, dtype=complex)
	can = np.cos(angle)
	san = np.sin(angle)
	d[6,6] = can * 1j
	d[6,7] = san
	d[7,6] = san
	d[7,7] = can * 1j
	g = QGate("D-" + str(angle))
	g.addLine(d)
	return g

def getSC(number): # Gets the number of significative ciphers of a given number
	return len(str(number).replace('.', ''))

def setSC(number, sc): # Returns the specified number with the specified significative ciphers
	res = 0
	num = str(number).split('.')
	i = len(num[0])
	d = 0
	if i >= sc:
		diff = i - sc
		res = int(num[0][0:sc]+"0"*diff)
	elif len(num) == 2:
		d = len(num[1])
		tsc = min(sc - i, d)
		diff = 0
		if sc - i > d:
			diff = sc - i - d
		res = float(num[0] + '.' + num[1][0:tsc]+"0"*diff)
		if d > tsc and num[1][tsc] >= '5':
			res += 10**-tsc
	return res

def toComp(angle, sc=None): # Returns a complex number with module 1 and the specified phase.
	while angle >= 2*np.pi:
		angle -= 2*np.pi
	if sc == None:
		sc = getSC(angle)
	res = np.around(np.cos(angle), decimals=sc-1) + np.around(np.sin(angle), decimals=sc-1)*1j
	return res

def PhaseShift(angle): # Phase shift (R) gate, rotates qubit with specified angle (in radians)
	ps = np.array([1, 0, 0, toComp(angle, 16)], dtype=complex)
	ps.shape = (2,2)
	g = QGate("R(" + str(angle) + ")")
	return ps

def Peres(): # A, B, C -> P = A, Q = A XOR B, R = AB XOR C. Peres gate.
	''' # Implementation of Peres gate with smaller gates.
	# Gates needed (without control SWAPs): 4
	p = QGate("Peres")
	p.addLine(SWAP(), I(1))
	p.addLine(I(1), ControlledU(SqrtNOT()))
	p.addLine(SWAP(), I(1))
	p.addLine(I(1), ControlledU(SqrtNOT()))
	p.addLine(CNOT(), I(1))
	p.addLine(I(1), ControlledU(Dagger(SqrtNOT())))
	return p
	'''
	p = QGate("Peres")
	p.addLine(Toffoli())
	p.addLine(CNOT(), I(1))
	return p

def R(): # A, B, C -> P = A XOR B, Q = A, R = AB XOR ¬C. R gate.
	# Optimized implementation with smaller gates
	# Gates needed (without control SWAPs): 6
	r = QGate("R")
	r.addLine(SWAP(), PauliX())
	r.addLine(I(1), ControlledU(SqrtNOT()))
	r.addLine(SWAP(), I(1))
	r.addLine(I(1), ControlledU(SqrtNOT()))
	r.addLine(CNOT(), I(1))
	r.addLine(I(1), ControlledU(Dagger(SqrtNOT())))
	r.addLine(SWAP(), I(1))
	return r

def TR(): # A, B, C -> P = A, Q = A XOR B, R = A¬B XOR C. TR gate.
	# Implementation of TR gate with smaller gates.
	# Gates needed (without control SWAPs): 6
	tr = QGate("TR")
	tr.addLine(I(1), PauliX(), I(1))
	tr.addLine(SWAP(), I(1))
	tr.addLine(I(1), ControlledU(SqrtNOT()))
	tr.addLine(SWAP(), I(1))
	tr.addLine(I(1), ControlledU(SqrtNOT()))
	tr.addLine(CNOT(), I(1))
	tr.addLine(I(1), ControlledU(Dagger(SqrtNOT())))
	tr.addLine(I(1), PauliX(), I(1))
	return tr

def URG(): # A, B, C -> P = (A+B) XOR C, Q = B, R = AB XOR C.
	# Implementation of URG gate with smaller gates.
	# Gates needed (without control SWAPs): 8
	urg = QGate("URG")
	urg.addLine(I(1), ControlledU(SqrtNOT()))
	urg.addLine(SWAP(), I(1))
	urg.addLine(I(1), ControlledU(SqrtNOT()))
	urg.addLine(CNOT(), I(1))
	urg.addLine(I(1), ControlledU(Dagger(SqrtNOT())))
	urg.addLine(CNOT(), I(1))
	urg.addLine(I(1), SWAP())
	urg.addLine(CNOT(), I(1))
	urg.addLine(I(1), CNOT())
	urg.addLine(CNOT(), I(1))
	urg.addLine(I(1), SWAP())
	urg.addLine(SWAP(), I(1))
	return urg

def BJN(): # A, B, C -> P = A, Q = B, R = (A+B) XOR C. BJN gate.
	# Implementation of TR gate with smaller gates.
	# Gates needed (without control SWAPs): 5
	bjn = QGate("BJN")
	bjn.addLine(SWAP(), I(1))
	bjn.addLine(I(1), ControlledU(SqrtNOT()))
	bjn.addLine(SWAP(), I(1))
	bjn.addLine(I(1), ControlledU(SqrtNOT()))
	bjn.addLine(CNOT(), I(1))
	bjn.addLine(I(1), ControlledU(SqrtNOT()))
	bjn.addLine(CNOT(), I(1))
	return bjn

def BlochCoords(qbit):
	alpha = qbit[0][0]
	pcorr = cm.rect(1, -cm.phase(alpha))
	alpha *= pcorr
	alpha = alpha.real
	beta = qbit[0][1] * pcorr
	theta = np.arccos(alpha) * 2
	s = np.sin(theta/2)
	if (s != 0):
		phi = np.log(beta/s)/1j
		phi = phi.real
	else:
		phi = 0.
	return (theta, phi)

def getTruthTable(gate, ancilla=None, garbage=0, iterations=1): # Prints the truth table of the given gate.
	# You can set the ancilla bits to not include them in the table, with the list of values they must have.
	# For example, if you have two 0s and one 1 as ancilla bits, ancilla[0,0,1]. It always takes the last bits as the ancilla ones!
	# The garbage=n removes the last n bits from the truth table, considering them garbage.
	# For example, if you have 6 outputs and the last 4 outputs are garbage, only the value of the first two would be printed.
	# Always removes the last n bits!
	num = int(np.log2(gate.shape[0]))
	mesd = {}
	for iteration in range(iterations):
		for i in range(0, gate.shape[0]):
			nbin = [int(x) for x in bin(i)[2:]]
			qinit = [0 for j in range(num - len(nbin))]
			qinit += nbin
			if ancilla == None or qinit[-len(ancilla):] == ancilla:
				qr = QRegistry(qinit)
				qr.ApplyGate(gate)
				mes = qr.Measure([1 for j in range(num-garbage)])
				if ancilla != None:
					ini = qinit[:-len(ancilla)]
				else:
					ini = qinit
				if str(ini) not in mesd:
					mesd[str(ini)] = np.zeros(num)
				mesd[str(ini)] = [x + y for x, y in zip(mesd[str(ini)], mes)]
	for k in mesd:
		print(k + " -> " + str(["P(1)=" + str(v/iterations) if v/iterations != 1.0 and v/iterations != 0.0 else int(v/iterations) for v in mesd[k]]))

def QEq(q1, q2):
	return np.array_equal(q1,q2) and str(q1) == str(q2)

def HalfSubstractor(): # A, B, 0 -> P = A-B, Q = Borrow, R = B = Garbage
	hs = QGate("Half Substractor")
	hs.addLine(SWAP(), I(1))
	hs.addLine(TR())
	hs.addLine(SWAP(), I(1))
	hs.addLine(I(1), SWAP())
	return hs

def Substractor(): # A, B, Bin, 0, 0, 0 -> P = A-B, Q = Borrow, R = B1 = Garbage, S = B1B2 = Garbage, T = Bin = Garbage, U = B = Garbage
	# Can be used as a comparator. Q will be 0 if A>=B, 1 otherwise.
	fs = QGate("Substractor")
	fs.addLine(I(2), SWAP(), I(2))
	fs.addLine(HalfSubstractor(), I(3))
	fs.addLine(I(2), SWAP(), I(2))
	fs.addLine(I(1), SWAP(), SWAP(), I(1))
	fs.addLine(I(2), SWAP(), SWAP())
	fs.addLine(HalfSubstractor(), I(3))
	fs.addLine(I(2), SWAP(), I(2))
	fs.addLine(I(3), SWAP(), I(1))
	fs.addLine(I(1), URG(), I(2))
	return fs

# Function that returns the 2^nth root of the unity
def nroot(n, rc = 14): # Rounds to 14 decimal places by default
	r = cm.exp(2j * cm.pi / pow(2, n))
	return round(r.real, rc) + round(r.imag, rc) * 1j

def RUnity(m, rc = 14):
	ru = np.eye(2, dtype=complex)
	ru[1,1] = nroot(m, rc)
	g = QGate("RU" + str(m))
	g.addLine(ru)
	return g

def QFT(size, rc = 14):
	'''
	size = 4
	uft = np.kron(Hadamard(1), I(size - 1))
	uft = np.dot(uft, np.kron(np.dot(SWAP(), ControlledU(RUnity(2, rc))), I(size - 2)))
	uft = np.dot(uft, np.kron(Hadamard(1), np.kron(np.dot(SWAP(), ControlledU(RUnity(3, rc))), I(size - 3))))
	uft = np.dot(uft, np.kron(np.dot(SWAP(), ControlledU(RUnity(2, rc))), np.dot(SWAP(), ControlledU(RUnity(4, rc)))))
	uft = np.dot(uft, np.kron(Hadamard(1), np.kron(np.dot(SWAP(), ControlledU(RUnity(3, rc))), I(size - 3))))
	uft = np.dot(uft, np.kron(np.dot(SWAP(), ControlledU(RUnity(2, rc))), I(size - 2)))
	uft = np.dot(uft, np.kron(Hadamard(1), I(size - 1)))
	uft = np.dot(uft, np.kron(SWAP(), I(size - 2)))
	uft = np.dot(uft, np.kron(I(size - 3), np.kron(SWAP(), I(size - 3))))
	uft = np.dot(uft, np.kron(SWAP(), SWAP()))
	uft = np.dot(uft, np.kron(I(size - 3), np.kron(SWAP(), I(size - 3))))
	uft = np.dot(uft, np.kron(SWAP(), I(size - 2)))

	return uft
	'''
	from tests.shor import DFT
	return DFT(pow(2, size))
