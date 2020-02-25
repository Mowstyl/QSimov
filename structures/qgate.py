import numpy as np
import ctypes as ct
import platform as plat
import connectors.parser as prs
from structures.funmatrix import Funmatrix
from structures.measure import Measure
from functools import reduce
from os.path import dirname, abspath, sep

# DLL Load
if plat.system() == "Windows":
    extension = ".dll"
else:
    extension = ".so"
__rootfolder__ = dirname(dirname(abspath(__file__)))
__libfolder__ = __rootfolder__ + sep + "lib"
__qsimovpath__ = __libfolder__ + sep + "libqsimov" + extension
__qsimov__ = ct.CDLL(__qsimovpath__)

c_double_p = ct.POINTER(ct.c_double)

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

__cGetQGate__ = __qsimov__.getQGate
__cGetQGate__.argtypes = [ct.c_char_p, ct.c_double, ct.c_double, ct.c_double, ct.c_int]
__cGetQGate__.restype = ct.c_void_p

__cGetQGateQubits__ = __qsimov__.getQGateQubits
__cGetQGateQubits__.argtypes = [ct.c_void_p]
__cGetQGateQubits__.restype = ct.c_int

__cGetQGateSize__ = __qsimov__.getQGateSize
__cGetQGateSize__.argtypes = [ct.c_void_p]
__cGetQGateSize__.restype = ct.c_uint

__cGet2dGate__ = __qsimov__.get2dGate
__cGet2dGate__.argtypes = [ct.c_void_p]
__cGet2dGate__.restype = c_double_p
def getGate(gatename):
    name, arg1, arg2, arg3, invert = prs.getGateData(gatename)
    qgate = ct.c_void_p(__cGetQGate__(ct.c_char_p(name.encode()), ct.c_double(arg1), ct.c_double(arg2), ct.c_double(arg3), ct.c_int(int(invert))))
    if (qgate or False) == qgate: # NOTNULLPTR and True returns True, NOTNULLPTR or False returns NOTNULLPTR
        nqubits = int(__cGetQGateQubits__(qgate))
        size = int(__cGetQGateSize__(qgate))
        plainmatrix2d = __cGet2dGate__(qgate)[:size*size*2]
        rematrix2d = plainmatrix2d[:size*size]
        immatrix2d = plainmatrix2d[size*size:size*size*2]
        matrix = np.array([complex(rematrix2d[i], immatrix2d[i]) for i in range(size*size)])
        matrix = matrix.reshape(size, size)
        # print("Gate for " + str(nqubits) + " qubits with size " + str(size))
    else:
        matrix = None
        print("Error while getting the specified gate!")

    return matrix

def getGateData(gatename):
    name, arg1, arg2, arg3, invert = prs.getGateData(gatename)
    if "SWAP" in name or name == "XX" or name == "YY" or name == "ZZ":
        shape = (2, 4)
    else:
        qgate = ct.c_void_p(__cGetQGate__(ct.c_char_p(name.encode()), ct.c_double(arg1), ct.c_double(arg2), ct.c_double(arg3), ct.c_int(int(invert))))
        if (qgate or False) == qgate: # NOTNULLPTR and True returns True, NOTNULLPTR or False returns NOTNULLPTR
            nqubits = int(__cGetQGateQubits__(qgate))
            size = int(__cGetQGateSize__(qgate))
            shape = (nqubits, size)
        else:
            shape = None
            print("Error while getting the specified gate!")

    return shape

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
                    newline += [None for i in range(gate.size)]
            else:
                if i < 1:
                    newline += [gate]
                else:
                    newline += [None for i in range(getGateSize(gate))]
        newlines += [newline]
    return newlines

def getGateSize(gate):
    if type(gate) == str:
        return getGateData(gate)
    elif type(gate) == Measure:
        return len(gate.mask)
    return gate.size

def getLines(gate):
    if type(gate) == Measure:
        return 1
    return len(gate.lines)
