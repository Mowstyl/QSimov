import cmath as cm
import numpy as np
import random as rnd
from structures.qgate import _getMatrix

class QRegistry:
	def __init__(self, nqbits, **kwargs):
		# nqbits -> number of QuBits in the registry.
		# Seed for the Pseudo Random Number Generation can be specified with seed = <seed> as an argument.
		self.state = np.zeros(2**nqbits, dtype=complex)
		self.state[0] = 1
		self.state.shape = (1, 2**nqbits)

	def measure(self, msk, remove = False): # List of numbers with the QuBits that should be measured. 0 means not measuring that qubit, 1 otherwise. remove = True if you want to remove a QuBit from the registry after measuring
		if (type(msk) != list or len(msk) != int(np.log2(self.state.size)) or \
			not all(type(num) == int and (num == 0 or num == 1) for num in msk)):
			raise ValueError('Not valid mask')
		mask = []
		for i in range(len(msk)):
			if msk[i] == 1:
				mask.append(i)
		tq = int(np.log2(self.state.size))
		if (not all(num < tq and num > -1 for num in mask)):
			raise ValueError('Out of range')
		mes = []
		for qbit in mask:
			r = rnd.random()
			p = 0
			max = 2**(tq - (qbit + 1))
			cnt = 0
			rdy = True
			for i in range(0, self.state.size):
				if (cnt == max):
					rdy = not rdy
					cnt = 0
				if (rdy):
					p += cm.polar(self.state[0,i])[0]**2
				cnt += 1
			if (r < p):
				me = 0
			else:
				me = 1
			mes.append(me)
			self.collapse((tq - (qbit + 1)), me, remove)
		return mes

	def applyGate(self, *gates): # Applies a quantum gate to the registry.
		gate = _getMatrix(gates[0])
		for g in list(gates)[1:]:
			gate = np.kron(gate, _getMatrix(g))
		self.state = np.transpose(np.dot(gate, ket(self.state)))

	def collapse(self, qbit, mes, remove): # Collapses a qubit from the registry. qbit is the id of the qubit, numerated as q0..qn in the registry. mes is the value obtained when measuring it. remove indicates whether it should be removed from the registry.
		max = 2**qbit
		cnt = 0
		rdy = mes == 1
		mfd = []
		for i in range(0, self.state.size):
			if (cnt == max):
				rdy = not rdy
				cnt = 0
			if (rdy):
				self.state[0, i] = 0
				mfd.append(i)
			cnt += 1
		if (remove):
			for qbit in mfd[::-1]:
				self.state = np.delete(self.state, qbit, 1)
		normalize(self.state)
	def densityMatrix(self):
		return np.dot(ket(self.state), bra(self.state))
	def vnEntropy(self, **kwargs):
		base = kwargs.get('base', "e")
		#dm = self.densityMatrix()
		#evalues, m = np.linalg.eig(dm)
		entropy = 0
		#for e in evalues:
		#	if e != 0:
		#		entropy += e * np.log(e)
		for amp in self.state[0]:
			p = cm.polar(amp)[0]**2
			if p > 0:
				if base == "e":
					entropy += p * np.log(p)
				elif type(base) == int or type(base) == float:
					entropy += p * np.log(p)/np.log(base)
		return -entropy

def prob(q, x): # Devuelve la probabilidad de obtener x al medir el qbit q
	p = 0
	if (x < q.size):
		p = cm.polar(q[0,x])[0]**2
	return p

def bra(v): # Devuelve el vector pasado como parametro en forma de fila conjugado. <v|
	b = v[:]
	s = v.shape
	if s[0] != 1:
		b = np.matrix.getH(b)
	else:
		b = np.conjugate(b)
	return b

def ket(v): # Devuelve el vector pasado como parametro en forma de columna. |v>
	k = v[:]
	s = v.shape
	if s[1] != 1:
		k = np.transpose(k)
	return k

def superposition(x, y): # Devuelve el estado compuesto por los dos QuBits.
	z = np.kron(x, y)
	normalize(z)
	return z

def normalize(state): # Funcion que asegura que se cumpla la propiedad que dice que |a|^2 + |b|^2 = 1 para cualquier QuBit. Si no se cumple, modifica el QuBit para que la cumpla si se puede.
	sqs = 0
	for i in range(0, state.size):
		sqs += cm.polar(state[0, i])[0]**2
	sqs = np.sqrt(sqs)
	if (sqs == 0):
		raise ValueError('Impossible QuBit')
	if (sqs != 1):
		for bs in state:
			bs /= sqs

def QBit(a,b): # Devuelve un QuBit con a y b. q = a|0> + b|1>, con a y b complejos
	q = np.array([a,b], dtype=complex)
	q.shape = (1,2)
	normalize(q)
	return q

def QZero(): # Devuelve un QuBit en el estado 0
	q = np.array([complex(1,0),complex(0,0)])
	q.shape = (1,2)
	return q

def QOne(): # Devuelve un QuBit en el estado 1
	q = np.array([complex(0,0),complex(1,0)])
	q.shape = (1,2)
	return q
