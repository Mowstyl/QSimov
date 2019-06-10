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
        if type(aux) == QGate:
                aux = aux.m
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

    def transpose(self): # Returns the Transpose of the given matrix
        self.name = self.name + "T"
        fm.transpose(self.m)
        return self

    def dagger(self): # Returns the Hermitian Conjugate or Conjugate Transpose of the given matrix
        self.name = self.name + "†"
        fm.dagger(self.m)
        return self

    def invert(self):
        return self.dagger()

def I(n=1): # Returns Identity Matrix for the specified number of QuBits
    ig = QGate("I(" + str(n) + ")")
    ig.addLine(ct.c_void_p(__cIdentity__(ct.c_int(n))))
    return ig

def _getMatrix(gate):
    m = gate
    if type(gate) == QGate:
        m = gate.m
    return m
