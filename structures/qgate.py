from structures.funmatrix import Funmatrix
from structures.measure import Measure
from functools import reduce
import numpy as np
import ctypes as ct

__qsimov__ = ct.CDLL("libqsimov.dll")
__cIdentity__ = __qsimov__.Identity
__cIdentity__.argtypes = [ct.c_int]
__cIdentity__.restype = ct.c_void_p

class QGate(object):
    def __init__(self, name="UNNAMED", size=None):
        self.simple = True
        self.lines = []
        self.name = name
        self.size = size

    def __getitem__(self, key):
        return self.getMatrix()[key]

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def addLine(self, *args):
        if (self.simple and self.lines == [] and len(args) == 1):
            size = int(np.log2(args[0].getMatrix().shape()[0]))
            if self.size == None:
                self.size = size
            if self.size == size:
                if type(args[0]) == QGate:
                    self.lines = args[0].lines
                else:
                    self.lines = [[args[0]]]
            else:
                raise ValueError("This gate requires a " + str(self.size) + " QuBit matrix. Received " + str(size) + " QuBit matrix.")
        else:
            size = sum(map(lambda gate: int(np.log2(gate.getMatrix().shape()[0])), args))
            if self.size == None:
                self.size = size
            if self.size == size:
                self.simple = False
                self.lines += joinGates(args)
            else:
                raise ValueError("This gate requires a " + str(self.size) + " QuBit matrix. Received " + str(size) + " QuBit matrix.")

    def setName(self, name):
        self.name = name

    def dagger(self): # Returns the Hermitian Conjugate or Conjugate Transpose of the given matrix
        invgate = QGate(self.name + "-1")
        if self.simple:
            invgate.lines = [[self.lines[0][0].invert()]]
        else:
            invgate.lines = [[gate.invert() for gate in line] for line in self.lines[::-1]]
        return invgate

    def invert(self):
        return self.dagger()

    def getMatrix(self):
        if self.simple:
            return self.lines[0][0]
        else:
            return reduce(lambda gate1, gate2: gate1 @ gate2, map(lambda line: reduce(lambda gate1, gate2: gate1.getMatrix() ** gate2.getMatrix(), line), self.lines))

def I(n=1): # Returns Identity Matrix for the specified number of QuBits
    ig = QGate("I(" + str(n) + ")")
    ig.addLine(Funmatrix(ct.c_void_p(__cIdentity__(ct.c_int(n))), ig.name))
    return ig

def joinGates(gates):
    maxgatelen = max(map(getLines, gates))
    newlines = []
    for i in range(maxgatelen):
        newline = []
        for gate in gates:
            if type(gate) != Measure:
                if i < len(gate.lines):
                    newline += gate.lines[i]
                else:
                    newline += [I(gate.size)]
            else:
                if i < 1:
                    newline += [gate]
                else:
                    newline += [I(getGateSize(gate))]
        newlines += [newline]
    return newlines

def getGateSize(gate):
    if type(gate) == str:
        gate = qapi.getGate(gate)
    elif type(gate) == Measure:
        return len(gate.mask)
    return int(np.log2(gate.getMatrix().shape()[0]))

def getLines(gate):
    if type(gate) == Measure:
        return 1
    return len(gate.lines)
