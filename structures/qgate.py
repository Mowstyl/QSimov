import numpy as np
import ctypes as ct
import platform as plat
import connectors.parser as prs
from collections.abc import Iterable
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
        self.empty = True
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
        args = [_rebuildGateName(gate) isinstance(type(gate), str) or (isinstance(gate, Iterable) and isinstance(type(gate[0]), str)) else gate for gate in args]
        size = sum([getGateSize(gate) for gate in args])
        if size > 0:
            if (self.empty):
                self.size = size
                self.empty = False
            if (self.size != size):
                raise ValueError("This gate requires a " + str(self.size) + " QuBit gate. Received " + str(size) + " QuBit gate.")
            parties = set()
            for i in range(len(args)):
                myparty = _getQuBitArg(args[i][0])
                if myparty is None:
                    myparty = set([i])
                if args[i][0] is not None:
                    if len(myparty.intersection(args[i][1])) == 0:
                        myparty = myparty.union(args[i][1])
                    else:
                        raise ValueError("You can't apply a gate to a qubit and use it as a control: " + str(myparty.intersection(args[i][1])))
                if args[i][0] is not None:
                    if len(myparty.intersection(args[i][2])) == 0:
                        myparty = myparty.union(args[i][2])
                    else:
                        raise ValueError("You can't apply a gate to a qubit and use it as a control, or use it as control and anticontrol at the same time: " + str(myparty.intersection(args[i][2])))
                if len(parties.intersection(myparty)) == 0:
                    parties = parties.union(myparty)
                else:
                    raise ValueError("You can't apply two or more gates to the same qubit in the same line: " + str(parties.intersection(myparty)))

            self.lines += [args]
        else:
            print ("No gates. Just Monika.")

    def setName(self, name):
        self.name = name

    def dagger(self): # Returns the Hermitian Conjugate or Conjugate Transpose of the given matrix
        invgate = QGate(self.name + "-1")
        for line in self.lines[::-1]:
            invgate.addLine(*[[_invertStrGate(gate[0]), gate[1], gate[2]] if type(gate) == list else gate.dagger() for gate in line])
        return invgate

    def invert(self):
        return self.dagger()

    def _applyGate(self, registry, qubit, controls, anticontrols):
        for line in self.lines:
            for i in range(len(line)):
                if line[i][1] is not None:
                    if controls is None:
                        controls = [qubit + c for c in line[i][1]]
                    else:
                        controls += [qubit + c for c in line[i][1]]
                if line[i][2] is not None:
                    if anticontrols is None:
                        anticontrols = [qubit + ac for ac in line[i][2]]
                    else:
                        anticontrols += [qubit + ac for ac in line[i][2]]
                if type(line[i][0]) == str:
                    registry.applyGate(_addQubitOffset(line[i][0], qubit), qubit=i+qubit, control=controls, anticontrol=anticontrols)
                else:
                    line[i][0]._applyGate(registry, qubit+i, controls, anticontrols)

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

def _rebuildGateName(gate):
    cons = []
    acons = []
    gatename = gate
    if type(gate) == list:
        gatename = gate[0]
        if (len(gate) > 1):
            if gate[1] is None or len(gate[1]) == 0:
                cons = None
            else:
                cons = set(gate[1])
        if (len(gate) > 2):
            if gate[2] is None or len(gate[2]) == 0:
                acons = None
            else:
                acons = set(gate[2])
    name, arg1, arg2, arg3, invert = prs.getGateData(gatename)
    if arg1 is not None:
        name += "(" + str(arg1)
    if arg2 is not None:
        name += str(arg2)
    if arg3 is not None:
        name += str(arg3) + ")"
    if invert:
        name += "-1"
    return [name, cons, acons]

def _getQuBitArg(gatename):
    args = None
    qubitargs = ["XX", "YY", "ZZ"]
    name, arg1, arg2, arg3, invert = prs.getGateData(gatename)
    if name in qubitargs or "SWAP" in name:
        if name in qubitargs:
            args = set([arg2, arg3])
        else:
            args = set([arg1, arg2])
    return args

def _addQubitOffset(gatename, offset):
    qubitargs = ["XX", "YY", "ZZ"]
    name, arg1, arg2, arg3, invert = prs.getGateData(gatename)
    if name in qubitargs or "SWAP" in name:
        arg2 += offset
        if name in qubitargs:
            arg3 += offset
        else:
            arg1 += offset
        if arg1 is not None:
            name += "(" + str(arg1)
        if arg2 is not None:
            name += str(arg2)
        if arg3 is not None:
            name += str(arg3) + ")"
        if invert:
            name += "-1"
        gatename = name
    return gatename


def _invertStrGate(gatename):
    selfinvert = ["X", "Y", "Z", "H", "SWAP"]
    name, arg1, arg2, arg3, invert = prs.getGateData(gatename)
    if name not in selfinvert:
        if invert:
            gatename = gatename[:-2]
        else:
            gatename += "-1"
    return gatename

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
    size = 0
    if type(gate) == str:
        size = getGateData(gate)[0]
    if type(gate) == list:
        cacs = 0
        if len(gate) > 1:
            if gate[1] is not None:
                cacs += len(gate[1])
        if len(gate) > 2:
            if gate[2] is not None:
                cacs += len(gate[2])
        size = getGateData(gate)[0] + cacs
    elif type(gate) == Measure:
        size = len(gate.mask)
    elif type(gate) == QGate:
        size = gate.size
    else:
        raise ValueError(str(gate) + " is not a gate!")
    return size

def getLines(gate):
    if type(gate) == Measure:
        return 1
    return len(gate.lines)
