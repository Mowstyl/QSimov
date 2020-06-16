import numpy as np
import ctypes as ct
import ctypes.util
import platform as plat
import os
from os.path import sep
from qsimov.structures.funmatrix import Funmatrix
from qsimov.structures.qgate import QGate
from qsimov.connectors.parser import getGateData
from collections.abc import Iterable

# DLL Load
if plat.system() == "Windows":
    __libc__ = ct.cdll.msvcrt
    extension = ".dll"
else:
    __libc__ = ct.cdll.LoadLibrary(ctypes.util.find_library("c"))
    extension = ".so"
__rootfolder__ = os.getcwd() + sep + "qsimov"
__libfolder__ = __rootfolder__ + sep + "lib"
__funmatpath__ = __libfolder__ + sep + "libfunmat" + extension
__qsimovpath__ = __libfolder__ + sep + "libqsimov" + extension
__funmat__ = ct.CDLL(__funmatpath__)
__qsimov__ = ct.CDLL(__qsimovpath__)

# C Pointer types
c_double_p = ct.POINTER(ct.c_double)
c_int_p = ct.POINTER(ct.c_int)

free = __libc__.free
free.argtypes = [ct.c_void_p]
free.restype = ct.c_void_p

__new_QRegistry__ = __qsimov__.new_QRegistry
__new_QRegistry__.argtypes = [ct.c_uint]
__new_QRegistry__.restype = ct.c_void_p

__cToString__ = __qsimov__.QR_toString
__cToString__.argtypes = [ct.c_void_p]
__cToString__.restype = ct.c_char_p

'''
__cApplyGateQubit__ = __qsimov__.applyGateQubit
__cApplyGateQubit__.argtypes = [ct.c_void_p, ct.c_char_p, ct.c_double, ct.c_double, ct.c_double, ct.c_int, ct.c_int]
__cApplyGateQubit__.restype = ct.c_int
'''

__cApplyGate__ = __qsimov__.applyGate
__cApplyGate__.argtypes = [ct.c_void_p, ct.c_char_p, ct.c_double, ct.c_double, ct.c_double, ct.c_int, ct.c_int, c_int_p, ct.c_int, c_int_p, ct.c_int]
__cApplyGate__.restype = ct.c_int

__cBlochCoords__ = __qsimov__.blochCoords
__cBlochCoords__.argtypes = [ct.c_void_p]
__cBlochCoords__.restype = c_double_p

__cHopfCoords__ = __qsimov__.hopfCoords
__cHopfCoords__.argtypes = [ct.c_void_p]
__cHopfCoords__.restype = c_double_p

__cMeasure__ = __qsimov__.measure
__cMeasure__.argtypes = [ct.c_void_p, c_int_p, ct.c_int, c_int_p, ct.c_int]
__cMeasure__.restype = ct.c_int

__cGetState__ = __qsimov__.getState
__cGetState__.argtypes = [ct.c_void_p]
__cGetState__.restype = c_double_p

__cSetState__ = __qsimov__.setState
__cSetState__.argtypes = [ct.c_void_p, c_double_p, ct.c_int]
__cSetState__.restype = ct.c_int

__cGetNQubits__ = __qsimov__.getNQubits
__cGetNQubits__.argtypes = [ct.c_void_p]
__cGetNQubits__.restype = ct.c_int

__cGetSize__ = __qsimov__.getSize
__cGetSize__.argtypes = [ct.c_void_p]
__cGetSize__.restype = ct.c_int

__cDensityMat__ = __qsimov__.DensityMat
__cDensityMat__.argtypes = [ct.c_void_p]
__cDensityMat__.restype = ct.c_void_p

__cFreeState__ = __qsimov__.freeState
__cFreeState__.argtypes = [ct.c_void_p]
__cFreeState__.restype = ct.c_int

__partialTrace__ = __funmat__.partial_trace
__partialTrace__.argtypes = [ct.c_void_p, ct.c_int]
__partialTrace__.restype = ct.c_void_p


class QRegistry:
    def __init__(self, nqbits):
        # nqbits -> number of QuBits in the registry.
        # Seed for the Pseudo Random Number Generation can be specified with seed = <seed> as an argument.
        self.reg = __new_QRegistry__(nqbits)

    def __del__(self):
        __cFreeState__(self.reg)
        self.reg = None

    def getSize(self):
        return int(__cGetSize__(self.reg))

    def getNQubits(self):
        return int(__cGetNQubits__(self.reg))

    def toString(self):
        return __cToString__(self.reg).decode("ascii")

    def measure(self, msk, remove=False):  # List of numbers with the QuBits that should be measured. 0 means not measuring that qubit, 1 otherwise. remove = True if you want to remove a QuBit from the registry after measuring
        nqubits = self.getNQubits()
        if (not isinstance(msk, Iterable) or len(msk) != nqubits or
                not all(type(num) == int and (num == 0 or num == 1) for num in msk)):
            raise ValueError('Not valid mask')
        mask = []
        for i in range(nqubits):
            if msk[-1-i] == 1:
                mask.append(i)
        if (not all(num < nqubits and num > -1 for num in mask)):
            raise ValueError('Out of range')

        nq = len(mask)
        if nq > 0:
            int_array = ct.c_int * nq
            result = int_array(*[0 for i in mask])
            rem = 0
            if (remove):
                rem = 1
            if int(__cMeasure__(self.reg, int_array(*mask), ct.c_int(nq), result, ct.c_int(rem))) == 0:
                print("Error measuring!")
            else:
                return list(result)[::-1]
        else:
            return []

    def applyGate(self, gate, qubit=0, control=None, anticontrol=None):
        if np.issubdtype(type(qubit), np.integer) and qubit < self.getNQubits() and qubit >= 0:
            if not isinstance(control, Iterable) and type(control) != int and not (control is None):
                raise ValueError("Control must be an int, a list of ints or None!")
            elif not isinstance(anticontrol, Iterable) and type(anticontrol) != int and not (anticontrol is None):
                raise ValueError("Anticontrol must be an int, a list of ints or None!")
            elif isinstance(control, Iterable) and isinstance(anticontrol, Iterable) and len(set(control) & set(anticontrol)) > 0:
                raise ValueError("A control can't be an anticontrol!")
            else:
                allOk = True

                clen = 0
                if np.issubdtype(type(control), np.integer):
                    int_array = ct.c_int * 1
                    control = int_array(control)
                    clen = 1
                elif control is None or len(control) == 0:
                    control = c_int_p()
                    clen = 0
                else:
                    control = set(control)
                    clen = len(control)
                    int_array = ct.c_int * clen
                    control = int_array(*control)
                    allOk = all(0 <= id < self.getNQubits() for id in control)

                aclen = 0
                if np.issubdtype(type(anticontrol), np.integer):
                    int_array = ct.c_int * 1
                    anticontrol = int_array(control)
                    aclen = 1
                elif anticontrol is None or len(anticontrol) == 0:
                    anticontrol = c_int_p()
                    aclen = 0
                else:
                    anticontrol = set(anticontrol)
                    aclen = len(anticontrol)
                    int_array = ct.c_int * aclen
                    anticontrol = int_array(*anticontrol)
                    allOk = allOk and all(0 <= id < self.getNQubits() for id in anticontrol)

                if allOk:
                    if type(gate) == str:
                        name, arg1, arg2, arg3, invert = getGateData(gate)
                        if arg1 is None:
                            arg1 = 0
                        if arg2 is None:
                            arg2 = 0
                        if arg3 is None:
                            arg3 = 0
                        if (name.lower() == "swap"):
                            __SWAPQubits__(self.reg, arg1, arg2, False, False, invert, control, clen, anticontrol, aclen)
                        elif (name.lower() == "iswap"):
                            __SWAPQubits__(self.reg, arg1, arg2, True, False, invert, control, clen, anticontrol, aclen)
                        elif (name.lower() == "sqrtswap"):
                            __SWAPQubits__(self.reg, arg1, arg2, False, True, invert, control, clen, anticontrol, aclen)
                        elif (name.lower() == "xx"):
                            __IsingQubits__(self.reg, 0, arg1, arg2, arg3, invert, control, clen, anticontrol, aclen)
                        elif (name.lower() == "yy"):
                            __IsingQubits__(self.reg, 1, arg1, arg2, arg3, invert, control, clen, anticontrol, aclen)
                        elif (name.lower() == "zz"):
                            __IsingQubits__(self.reg, 2, arg1, arg2, arg3, invert, control, clen, anticontrol, aclen)
                        else:
                            if int(__cApplyGate__(self.reg, ct.c_char_p(name.encode()), ct.c_double(arg1), ct.c_double(arg2), ct.c_double(arg3), ct.c_int(int(invert)), ct.c_int(qubit), control, ct.c_int(clen), anticontrol, ct.c_int(aclen))) == 0:
                                print("Error applying gate to specified QuBit!")
                    elif type(gate) == QGate:
                        gate._applyGate(self, qubit, control[:clen], anticontrol[:aclen])
                else:
                    print("The ids must be between 0 and " + str(self.getNQubits))
        else:
            if not np.issubdtype(type(qubit), np.integer):
                print("Qubit must be of integer type!")
            elif qubit >= self.getNQubits() or qubit < 0:
                print("The specified qubit doesn't exist!")

    def getState(self):
        rawState = __cGetState__(self.reg)
        size = self.getSize()
        state = np.array([complex(rawState[i], rawState[i + size]) for i in range(size)])
        free(rawState)

        return state

    def setState(self, newState, nqubits):
        size = 2 << (nqubits - 1)  # 2^n_qubits
        statetype = ct.c_double * (size * 2)
        cstate = statetype(*[newState[i].real if i < size else newState[i - size].imag for i in range(size * 2)])
        result = __cSetState__(self.reg, cstate, ct.c_int(nqubits)) == 1

        return result

    def densityMatrix(self):
        return Funmatrix(ct.c_void_p(__cDensityMat__(self.reg)), "Rho")

    def reducedDensityMatrix(self, elem):
        nq = self.getNQubits()
        if 0 <= elem < nq:
            # elem = nq - elem - 1
            return Funmatrix(ct.c_void_p(__partialTrace__(ct.c_void_p(__cDensityMat__(self.reg)), ct.c_int(elem))), "Tr_" + str(elem) + "(Rho)")
        else:
            print("The specified QuBit doesn't exist in this registry!")

    def reducedTrace(self, elem):
        rho_a = np.array(self.reducedDensityMatrix(elem)[:])
        rt = (rho_a @ rho_a).trace().real
        return rt if rt <= 1.0 else 1.0

    def vnEntropy(self, **kwargs):
        base = kwargs.get('base', 2)
        # dm = self.densityMatrix()
        # evalues, m = np.linalg.eig(dm)
        entropy = 0
        # for e in evalues:
        #     if e != 0:
        #         entropy += e * np.log(e)
        size = self.getSize()
        for i in range(size):
            p = self.prob(i)
            if p > 0:
                if base == "e":
                    entropy += p * np.log(p)
                elif type(base) == int or type(base) == float:
                    entropy += p * np.log(p)/np.log(base)
        return -entropy

    def hopfCoords(self):
        if self.getNQubits() == 1:
            return __cHopfCoords__(self.reg)[:3]
        else:
            # print("You can only use 1 qubit registries!")
            return None

    def blochCoords(self):
        if self.getNQubits() == 1:
            return __cBlochCoords__(self.reg)[:2]
        else:
            # print("You can only use 1 qubit registries!")
            return None

    def bra(self):  # Devuelve el vector de estado en forma de fila conjugado. <v|
        k = np.array(self.getState())
        k.shape = (1, k.shape[0])
        return np.conjugate(k)

    def ket(self):  # Devuelve el vector de estado en forma de columna. |v>
        k = np.array(self.getState())
        k.shape = (k.shape[0], 1)
        return k

    def prob(self, x):  # Devuelve la probabilidad de obtener x al medir el registro
        p = 0
        if (x < self.getSize()):
            p = self.densityMatrix()[x, x].real
        return p


def _calculateState(tnum, cnum, size, fun):
    indexes = [i for i in range(cnum, size, tnum)]
    return [(index, fun(index, 0)) for index in indexes]


__cJoinStates__ = __qsimov__.joinStates
__cJoinStates__.argtypes = [ct.c_void_p, ct.c_void_p]
__cJoinStates__.restype = ct.c_void_p


def superposition(a, b):  # Devuelve el estado compuesto por los dos QuBits.
    r = QRegistry(1)
    __cFreeState__(r.reg)
    r.reg = ct.c_void_p(__cJoinStates__(a.reg, b.reg))
    return r


__cSWAPQubits__ = __qsimov__.SWAPQubits
__cSWAPQubits__.argtypes = [ct.c_void_p, ct.c_int, ct.c_int, ct.c_int, ct.c_int, ct.c_int, c_int_p, ct.c_int, c_int_p, ct.c_int]
__cSWAPQubits__.restype = ct.c_int


def __SWAPQubits__(reg, qubit1, qubit2, imaginary, sqrt, invert, control, clen, anticontrol, aclen):  # Applies a quantum gate to the registry
    if int(__cSWAPQubits__(reg, ct.c_int(qubit1), ct.c_int(qubit2), ct.c_int(int(imaginary)), ct.c_int(int(sqrt)), ct.c_int(int(invert)), control, ct.c_int(clen), anticontrol, ct.c_int(aclen))) == 0:
        print("Error swapping specified QuBits!")


__cIsingQubits__ = __qsimov__.IsingQubits
__cIsingQubits__.argtypes = [ct.c_void_p, ct.c_int, ct.c_double, ct.c_int, ct.c_int, ct.c_int, c_int_p, ct.c_int, c_int_p, ct.c_int]
__cIsingQubits__.restype = ct.c_int


def __IsingQubits__(reg, type, angle, qubit1, qubit2, invert, control, clen, anticontrol, aclen):  # Applies a quantum gate to the registry
    if int(__cIsingQubits__(reg, ct.c_int(type), ct.c_double(angle), ct.c_int(qubit1), ct.c_int(qubit2), ct.c_int(int(invert)), control, ct.c_int(clen), anticontrol, ct.c_int(aclen))) == 0:
        print("Error coupling specified QuBits!")
