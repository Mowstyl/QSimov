"""Module that provides the lowest level data structure: QRegistry.

Data Structures:
    QRegistry: Quantum Registry, base of all quantum related operations

Functions:
    superposition: join two registries into one by calculating tensor product.
"""
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
__rootfolder__ = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
__libfolder__ = __rootfolder__ + sep + "lib"
__funmatpath__ = __libfolder__ + sep + "libfunmat" + extension
__qsimovpath__ = __libfolder__ + sep + "libqsimov" + extension
if hasattr(os, "add_dll_directory"):
    __funmatpath__ = "libfunmat" + extension
    __qsimovpath__ = "libqsimov" + extension
__funmat__ = ct.CDLL(__funmatpath__)
__qsimov__ = ct.CDLL(__qsimovpath__)

# C Pointer types
c_double_p = ct.POINTER(ct.c_double)
c_int_p = ct.POINTER(ct.c_int)

free = __libc__.free
free.argtypes = [ct.c_void_p]
free.restype = ct.c_void_p

__new_QRegistry__ = __qsimov__.new_QRegistry  # C QRegistry Constructor
__new_QRegistry__.argtypes = [ct.c_uint]      # Number of qubits
__new_QRegistry__.restype = ct.c_void_p       # Pointer to the C QRegistry
"""C QRegistry Constructor.
Positional arguments:
    unsigned integer -> number of qubits
Return:
    Pointer to the C QRegistry
"""


__cToString__ = __qsimov__.QR_toString
__cToString__.argtypes = [ct.c_void_p]
__cToString__.restype = ct.c_char_p
"""C QRegistry toString function.
Positional arguments:
    pointer -> C QRegistry
Return:
    String with the QRegistry representation
"""


__cApplyGate__ = __qsimov__.applyGate
__cApplyGate__.argtypes = [ct.c_void_p, ct.c_char_p, ct.c_double, ct.c_double,
                           ct.c_double, ct.c_int, ct.c_int, c_int_p, ct.c_int,
                           c_int_p, ct.c_int]
__cApplyGate__.restype = ct.c_int
"""C QRegistry applyGate function.
Positional arguments:
    pointer -> C QRegistry
    string -> Name of the gate
    double -> First parameter of the gate (optional)
    double -> Second parameter of the gate (optional)
    double -> Third parameter of the gate (optional)
    int -> Boolean. Whether or not to invert the gate. 1 -> Invert
    int -> id of the less significant qubit this gate will be applied to
    int array -> array containing the ids of the control qubits (optional)
    int -> size of the control ids array
    int array -> array containing the ids of the anticontrol qubits (optional)
    int -> size of the anticontrol ids array
Return:
    int: 0 -> Failure. 1 -> Success
"""


__cSWAPQubits__ = __qsimov__.SWAPQubits
__cSWAPQubits__.argtypes = [ct.c_void_p, ct.c_int, ct.c_int, ct.c_int,
                            ct.c_int, ct.c_int, c_int_p, ct.c_int,
                            c_int_p, ct.c_int]
__cSWAPQubits__.restype = ct.c_int
"""C QRegistry SWAPQubits function.
Positional arguments:
    pointer -> C QRegistry
    int -> id of the first qubit to swap
    int -> id of the second qubit to swap
    int -> Boolean. Whether or not this is the ISWAP gate.
        0 -> SWAP
        1 -> ISWAP
    int -> Boolean. Whether or not this is the sqrt of the gate. 1 -> sqrt
    int -> Boolean. Whether or not to invert the gate. 1 -> Invert
    int array -> array containing the ids of the control qubits (optional)
    int -> size of the control ids array
    int array -> array containing the ids of the anticontrol qubits (optional)
    int -> size of the anticontrol ids array
Return:
    int: 0 -> Failure. 1 -> Success
"""


__cIsingQubits__ = __qsimov__.IsingQubits
__cIsingQubits__.argtypes = [ct.c_void_p, ct.c_int, ct.c_double, ct.c_int,
                             ct.c_int, ct.c_int, c_int_p, ct.c_int,
                             c_int_p, ct.c_int]
__cIsingQubits__.restype = ct.c_int
"""C QRegistry IsingQubits function.
Positional arguments:
    pointer -> C QRegistry
    int -> type of the Ising gate
        0 -> XX
        1 -> YY
        2 -> ZZ
    double -> angle parameter of the Ising gate
    int -> id of the first qubit
    int -> id of the second qubit
    int -> Boolean. Whether or not to invert the gate. 1 -> Invert
    int array -> array containing the ids of the control qubits (optional)
    int -> size of the control ids array
    int array -> array containing the ids of the anticontrol qubits (optional)
    int -> size of the anticontrol ids array
Return:
    int: 0 -> Failure. 1 -> Success
"""


__cBlochCoords__ = __qsimov__.blochCoords
__cBlochCoords__.argtypes = [ct.c_void_p]
__cBlochCoords__.restype = c_double_p
"""C QRegistry blochCoords function.
Positional arguments:
    pointer -> C QRegistry
Return:
    double[2]: the polar coordinates of the qubit in the Bloch sphere
"""


__cMeasure__ = __qsimov__.measure
__cMeasure__.argtypes = [ct.c_void_p, c_int_p, ct.c_int, c_int_p, ct.c_int]
__cMeasure__.restype = ct.c_int
"""C QRegistry measure function.
Positional arguments:
    pointer -> C QRegistry
    int array -> ids of the qubits we want to measure
    int -> size of the qubit ids array
    int array -> array in which the result will be stored
        Same size as previous array
        Sorted from most significant to least significant qubit
Return:
    int: 0 -> Failure. 1 -> Success
"""

__cGetState__ = __qsimov__.getState
__cGetState__.argtypes = [ct.c_void_p]
__cGetState__.restype = c_double_p
"""C QRegistry getState function.
Positional arguments:
    pointer -> C QRegistry
Return:
    double array: array of doubles containing the state
        The first half of the array contains the real part of the amplitudes
        The second half contains the imaginary part of the amplitudes
"""

__cSetState__ = __qsimov__.setState
__cSetState__.argtypes = [ct.c_void_p, c_double_p, ct.c_int]
__cSetState__.restype = ct.c_int
"""C QRegistry setState function.
Positional arguments:
    pointer -> C QRegistry
    double array: array of doubles containing the state
        The first half of the array contains the real part of the amplitudes
        The second half contains the imaginary part of the amplitudes
    int: number of qubits the new state has
Return:
    int: 0 -> Failure. 1 -> Success
"""

__cGetNQubits__ = __qsimov__.getNQubits
__cGetNQubits__.argtypes = [ct.c_void_p]
__cGetNQubits__.restype = ct.c_int
"""C QRegistry getNQubits function.
Positional arguments:
    pointer -> C QRegistry
Return:
    int: number of qubits in the registry
"""

__cGetSize__ = __qsimov__.getSize
__cGetSize__.argtypes = [ct.c_void_p]
__cGetSize__.restype = ct.c_int
"""C QRegistry getSize function.
Positional arguments:
    pointer -> C QRegistry
Return:
    int: number of amplitudes in the state vector of the registry
"""

__cDensityMat__ = __qsimov__.DensityMat
__cDensityMat__.argtypes = [ct.c_void_p]
__cDensityMat__.restype = ct.c_void_p
"""C QRegistry DensityMat function.
Positional arguments:
    pointer -> C QRegistry
Return:
    pointer to a functional matrix representing the density matrix
"""

__cFreeState__ = __qsimov__.freeState
__cFreeState__.argtypes = [ct.c_void_p]
__cFreeState__.restype = ct.c_int
"""C QRegistry freeState function.
Positional arguments:
    pointer -> C QRegistry
Return:
    int: 0 -> Failure. 1 -> Success
"""

__partialTrace__ = __funmat__.partial_trace
__partialTrace__.argtypes = [ct.c_void_p, ct.c_int]
__partialTrace__.restype = ct.c_void_p
"""C QRegistry partial_trace function.
Positional arguments:
    pointer -> C QRegistry
    int -> id of the qubit to trace out
Return:
    pointer to a functional matrix representing the reduced density matrix
"""


class QRegistry:
    """Quantum Registry, base of all quantum related operations."""

    def __init__(self, nqbits):
        """Initialize QRegistry to state 0.

        nqbits -> number of QuBits in the registry.
        """
        self.reg = __new_QRegistry__(nqbits)

    def __del__(self):
        """Release memory held by the QRegistry."""
        __cFreeState__(self.reg)
        self.reg = None

    def getSize(self):
        """Return the number of elements in the state vector."""
        return int(__cGetSize__(self.reg))

    def getNQubits(self):
        """Return the number of qubits in this registry."""
        return int(__cGetNQubits__(self.reg))

    def toString(self):
        """Return string representation of this registry."""
        return __cToString__(self.reg).decode("ascii")

    def measure(self, msk, remove=False):
        """Measure specified qubits of this registry and collapse.

        Positional arguments:
            msk -> List of numbers with the QuBits that should be measured
                0 means not measuring that qubit, 1 otherwise
        Keyworded arguments:
            remove = True if you want to remove a QuBit from the registry
                after measuring

        Return:
            List with the value obtained after each measure
        """
        nqubits = self.getNQubits()
        if (not isinstance(msk, Iterable)
                or len(msk) != nqubits
                or not all(type(num) == int
                           and (num == 0 or num == 1)
                           for num in msk)):
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
            if int(__cMeasure__(self.reg, int_array(*mask), ct.c_int(nq),
                                result, ct.c_int(rem))) == 0:
                print("Error measuring!")
            else:
                return list(result)[::-1]
        else:
            return []

    def applyGate(self, gate, qubit=0, control=None, anticontrol=None,
                  optimize=True):
        """Apply specified gate to specified qubit with specified controls."""
        if (np.issubdtype(type(qubit), np.integer)
                and qubit < self.getNQubits() and qubit >= 0):
            if (not isinstance(control, Iterable)
                    and type(control) != int
                    and not (control is None)):
                raise ValueError("Control must be an int, " +
                                 "a list of ints or None!")
            elif (not isinstance(anticontrol, Iterable)
                  and type(anticontrol) != int
                  and not (anticontrol is None)):
                raise ValueError("Anticontrol must be an int, " +
                                 "a list of ints or None!")
            elif (isinstance(control, Iterable)
                  and isinstance(anticontrol, Iterable)
                  and len(set(control) & set(anticontrol)) > 0):
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
                    allOk = allOk and all(0 <= id < self.getNQubits()
                                          for id in anticontrol)

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
                            __SWAPQubits__(self.reg, arg1, arg2, False, False,
                                           invert, control, clen,
                                           anticontrol, aclen)
                        elif (name.lower() == "iswap"):
                            __SWAPQubits__(self.reg, arg1, arg2, True, False,
                                           invert, control, clen,
                                           anticontrol, aclen)
                        elif (name.lower() == "sqrtswap"):
                            __SWAPQubits__(self.reg, arg1, arg2, False, True,
                                           invert, control, clen,
                                           anticontrol, aclen)
                        elif (name.lower() == "xx"):
                            __IsingQubits__(self.reg, 0, arg1, arg2, arg3,
                                            invert, control, clen,
                                            anticontrol, aclen)
                        elif (name.lower() == "yy"):
                            __IsingQubits__(self.reg, 1, arg1, arg2, arg3,
                                            invert, control, clen,
                                            anticontrol, aclen)
                        elif (name.lower() == "zz"):
                            __IsingQubits__(self.reg, 2, arg1, arg2, arg3,
                                            invert, control, clen,
                                            anticontrol, aclen)
                        else:
                            if int(__cApplyGate__(self.reg,
                                                  ct.c_char_p(name.encode()),
                                                  ct.c_double(arg1),
                                                  ct.c_double(arg2),
                                                  ct.c_double(arg3),
                                                  ct.c_int(int(invert)),
                                                  ct.c_int(qubit),
                                                  control,
                                                  ct.c_int(clen),
                                                  anticontrol,
                                                  ct.c_int(aclen))) == 0:
                                print("Error applying gate to specified QuBit")
                    elif type(gate) == QGate:
                        gate._applyGate(self, qubit, control[:clen],
                                        anticontrol[:aclen], optimize=optimize)
                else:
                    print("The ids must be between 0 and " +
                          str(self.getNQubits))
        else:
            if not np.issubdtype(type(qubit), np.integer):
                print("Qubit must be of integer type!")
            elif qubit >= self.getNQubits() or qubit < 0:
                print("The specified qubit doesn't exist!")

    def getState(self):
        """Return the state vector of the registry."""
        rawState = __cGetState__(self.reg)
        size = self.getSize()
        state = np.array([complex(rawState[i], rawState[i + size])
                          for i in range(size)])
        free(rawState)

        return state

    def setState(self, newState, nqubits):
        """Set the state vector of the registry."""
        size = 2 << (nqubits - 1)  # 2^n_qubits
        statetype = ct.c_double * (size * 2)
        cstate = statetype(*[newState[i].real if i < size
                             else newState[i - size].imag
                             for i in range(size * 2)])
        result = __cSetState__(self.reg, cstate, ct.c_int(nqubits)) == 1

        return result

    def densityMatrix(self):
        """Return functional matrix of the density matrix."""
        return Funmatrix(ct.c_void_p(__cDensityMat__(self.reg)), "Rho")

    def reducedDensityMatrix(self, elem):
        """Return functional matrix of the reduced density matrix."""
        nq = self.getNQubits()
        if 0 <= elem < nq:
            # elem = nq - elem - 1
            c_density_mat = ct.c_void_p(__cDensityMat__(self.reg))
            c_reduced_mat = ct.c_void_p(__partialTrace__(c_density_mat,
                                                         ct.c_int(elem)))
            return Funmatrix(c_reduced_mat, "Tr_" + str(elem) + "(Rho)")
        else:
            print("The specified QuBit doesn't exist in this registry!")

    def reducedTrace(self, elem):
        """Trace of the square reduced density matrix."""
        rho_a = np.array(self.reducedDensityMatrix(elem)[:])
        rt = (rho_a @ rho_a).trace().real
        return rt if rt <= 1.0 else 1.0

    def vnEntropy(self, **kwargs):
        """Calculate Von Newmann Entropy.

        Keyworded arguments:
            base: the base of the logarithm. Default = 2
                The string "e" is a valid value
        """
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

    def blochCoords(self):
        """Get the polar coordinates of ONE qubit in the bloch sphere."""
        if self.getNQubits() == 1:
            return __cBlochCoords__(self.reg)[:2]
        else:
            # print("You can only use 1 qubit registries!")
            return None

    def bra(self):
        """Get the conjugated row form state vector (bra <v|)."""
        k = np.array(self.getState())
        k.shape = (1, k.shape[0])
        return np.conjugate(k)

    def ket(self):
        """Get the column form state vector (ket |v>)."""
        k = np.array(self.getState())
        k.shape = (k.shape[0], 1)
        return k

    def prob(self, x):
        """Get the odds of measuring a certain state x (not qubit)."""
        p = 0
        if (x < self.getSize()):
            p = self.densityMatrix()[x, x].real
        return p


def _calculateState(tnum, cnum, size, fun):
    """TODO: Find out the purpose of this function."""
    indexes = [i for i in range(cnum, size, tnum)]
    return [(index, fun(index, 0)) for index in indexes]


__cJoinStates__ = __qsimov__.joinStates
__cJoinStates__.argtypes = [ct.c_void_p, ct.c_void_p]
__cJoinStates__.restype = ct.c_void_p


def superposition(a, b):  # Devuelve el estado compuesto por los dos QuBits.
    """Join two registries into one by calculating tensor product."""
    r = QRegistry(1)
    __cFreeState__(r.reg)
    r.reg = ct.c_void_p(__cJoinStates__(a.reg, b.reg))
    return r


def __SWAPQubits__(reg, qubit1, qubit2, imaginary, sqrt, invert,
                   control, clen, anticontrol, aclen):
    """Apply a type of SWAP gate to the registry."""
    if int(__cSWAPQubits__(reg, ct.c_int(qubit1), ct.c_int(qubit2),
                           ct.c_int(int(imaginary)), ct.c_int(int(sqrt)),
                           ct.c_int(int(invert)), control, ct.c_int(clen),
                           anticontrol, ct.c_int(aclen))) == 0:
        print("Error swapping specified QuBits!")


def __IsingQubits__(reg, type, angle, qubit1, qubit2, invert,
                    control, clen, anticontrol, aclen):
    """Apply a type of Ising gate to the registry."""
    if int(__cIsingQubits__(reg, ct.c_int(type), ct.c_double(angle),
                            ct.c_int(qubit1), ct.c_int(qubit2),
                            ct.c_int(int(invert)), control, ct.c_int(clen),
                            anticontrol, ct.c_int(aclen))) == 0:
        print("Error coupling specified QuBits!")
