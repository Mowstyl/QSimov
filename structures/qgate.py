import numpy as np
import structures.funmatrix as fm
import ctypes as ct

__qsimov__ = ct.CDLL("libqsimov.dll")
__cIdentity__ = __qsimov__.Identity
__cIdentity__.argtypes = [ct.c_int]
__cIdentity__.restype = ct.c_void_p

class QGate(object):
	def __init__(self, name="UNNAMED"):
		self.m = 1
		self.simple = True
		self.lines = []
		self.name = name

	def __getitem__(self, key):
		return fm.getItem(self.m, key)

	def __setitem__(self, key, value):
		self.m[key] = value

	def __delitem__(self, key):
		del self.m[key]

	def __repr__(self):
		return self.name

	def __str__(self):
		return self.name

	def __lt__(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		return self.m.__lt__(m)

	def __le_(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		return self.m.__le__(m)

	def __eq__(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		return self.m.__eq__(m)

	def __ne_(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		return self.m.__ne__(m)

	def __gt__(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		return self.m.__gt__(m)

	def __ge_(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		return self.m.__ge__(m)

	def __add__(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		sol = QGate()
		sol.addLine(fm.madd(self.m, m))
		return sol

	def __sub__(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		sol = QGate()
		sol.addLine(fm.msub(self.m, m))
		return sol

	def __mod__(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		sol = QGate()
		sol.addLine(self.m.__mod__(m))
		return sol

	def __mul__(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		sol = QGate()
		sol.addLine(fm.ewmul(self.m, m))
		return sol

	def __rmul__(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		sol = QGate()
		sol.addLine(fm.ewmul(self.m, m))
		return sol

	def __imul__(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		sol = QGate()
		sol.addLine(fm.ewmul(m, self.m))
		return sol

	def __matmul__(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		sol = QGate()
		sol.addLine(fm.matmul(self.m, m))
		return sol

	def __pow__(self, other):
		m = other
		if type(other) == QGate:
			m = other.m
		sol = QGate()
		sol.addLine(fm.kron(self.m, m))
		return sol

	def addLine(self, *args):
		self.lines.append(list(args))
		if self.simple and (len(list(args)) > 1 or len(self.lines) > 1):
			self.simple = False
		aux = args[0]
		for gate in args[1:]:
			g = gate
			if type(gate) == QGate:
				g = gate.m
			aux = fm.kron(aux, g)
		if (self.m != 1):
			self.m = fm.matmul(aux, self.m)
		else:
			self.m = aux

	def setName(self, name):
		self.name = name

def I(n): # Returns Identity Matrix for the specified number of QuBits
	return ct.c_void_p(__cIdentity__(ct.c_int(n)))

def _getMatrix(gate):
	m = gate
	if type(gate) == QGate:
		m = gate.m
	return m

def unitaryMatrix(mat, decimals=10):
	mustbei = np.around(np.dot(_getMatrix(mat), _getMatrix(dagger(mat))), decimals=decimals)
	return (mustbei == I(int(np.log2(mustbei.shape[0])))).all()

def normalizeGate(mat):
	det = np.linalg.det(mat)
	if det != 0:
		return mat/det
	else:
		return None

def transpose(gate): # Returns the Transpose of the given matrix
	if type(gate) == QGate:
		t = QGate(gate.name + "T")
		if type(gate.m) == fm.FunctionalMatrix:
			t.addLine(gate.m.transpose())
		else:
			t.addLine(np.matrix.transpose(gate.m))
	elif type(gate) == ct.c_void_p:
		t = QGate("UT")
		t.addLine(fm.transpose(gate))
	else:
		t = QGate("UT")
		t.addLine(np.matrix.transpose(gate))
	return t

def dagger(gate): # Returns the Hermitian Conjugate or Conjugate Transpose of the given matrix
	if type(gate) == QGate:
		t = QGate(gate.name + "†")
		if gate.simple:
			t.addLine(dagger(gate.m))
		else:
			lines = gate.lines[::-1]
			for line in lines:
				t.addLine(*[dagger(g).m for g in line])
			t.setMult(gate.mult)
	else:
		t = QGate("U†")
		if type(gate) == ct.c_void_p:
			t.addLine(fm.dagger(gate))
		else:
			t.addLine(np.matrix.getH(gate))
	return t

def invert(gate): # Returns the inverse of the given matrix
	if type(gate) == QGate:
		t = QGate(gate.name + "-¹")
		t.addLine(np.linalg.inv(gate.m))
	else:
		t = QGate("U-¹")
		t.addLine(np.linalg.inv(gate))
	return t
