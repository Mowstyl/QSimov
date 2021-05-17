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
from qsimov.connectors.parser import get_gate_data
from collections.abc import Iterable

# DLL Load
if plat.system() == "Windows":
    _libc = ct.cdll.msvcrt
    extension = ".dll"
else:
    _libc = ct.cdll.LoadLibrary(ctypes.util.find_library("c"))
    extension = ".so"
_root_folder = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
_lib_folder = _root_folder + sep + "lib"
_funmat_path = _lib_folder + sep + "libfunmat" + extension
_qsimov_path = _lib_folder + sep + "libqsimov" + extension
if hasattr(os, "add_dll_directory"):
    _funmat_path = "libfunmat" + extension
    _qsimov_path = "libqsimov" + extension
_funmat = ct.CDLL(_funmat_path)
_qsimov = ct.CDLL(_qsimov_path)

# C Pointer types
c_double_p = ct.POINTER(ct.c_double)
c_int_p = ct.POINTER(ct.c_int)

_cFree = _libc.free
_cFree.argtypes = [ct.c_void_p]
_cFree.restype = ct.c_void_p
"""C free function.
Positional arguments:
    pointer to release
Return:
    nothing
"""

_new_QRegistry = _qsimov.new_QRegistry  # C QRegistry Constructor
_new_QRegistry.argtypes = [ct.c_uint]      # Number of qubits
_new_QRegistry.restype = ct.c_void_p       # Pointer to the C QRegistry
"""C QRegistry Constructor.
Positional arguments:
    unsigned integer -> number of qubits
Return:
    Pointer to the C QRegistry
"""

'''
_cGetRegMemory = _qsimov.getRegMemory
_cGetRegMemory.argtypes = [ct.c_void_p]
_cGetRegMemory.restype = ct.c_int
'''
"""C sizeof function.
Positional arguments:
    pointer to reg
Return:
    size of state in bytes
"""


_cToString = _qsimov.QR_toString
_cToString.argtypes = [ct.c_void_p]
_cToString.restype = ct.c_char_p
"""C QRegistry toString function.
Positional arguments:
    pointer -> C QRegistry
Return:
    String with the QRegistry representation
"""


_cApplyGate = _qsimov.applyGate
_cApplyGate.argtypes = [ct.c_void_p, ct.c_char_p, ct.c_double, ct.c_double,
                        ct.c_double, ct.c_int, ct.c_int, c_int_p, ct.c_int,
                        c_int_p, ct.c_int]
_cApplyGate.restype = ct.c_int
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


_cSWAPQubits = _qsimov.SWAPQubits
_cSWAPQubits.argtypes = [ct.c_void_p, ct.c_int, ct.c_int, ct.c_int, ct.c_int,
                         ct.c_int, c_int_p, ct.c_int, c_int_p, ct.c_int]
_cSWAPQubits.restype = ct.c_int
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


_cIsingQubits = _qsimov.IsingQubits
_cIsingQubits.argtypes = [ct.c_void_p, ct.c_int, ct.c_double, ct.c_int,
                          ct.c_int, ct.c_int, c_int_p, ct.c_int,
                          c_int_p, ct.c_int]
_cIsingQubits.restype = ct.c_int
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


_cBlochCoords = _qsimov.blochCoords
_cBlochCoords.argtypes = [ct.c_void_p]
_cBlochCoords.restype = c_double_p
"""C QRegistry blochCoords function.
Positional arguments:
    pointer -> C QRegistry
Return:
    double[2]: the polar coordinates of the qubit in the Bloch sphere
"""


_cMeasure = _qsimov.measure
_cMeasure.argtypes = [ct.c_void_p, c_int_p, ct.c_int, c_int_p, ct.c_int]
_cMeasure.restype = ct.c_int
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

_cGetState = _qsimov.getState
_cGetState.argtypes = [ct.c_void_p]
_cGetState.restype = c_double_p
"""C QRegistry getState function.
Positional arguments:
    pointer -> C QRegistry
Return:
    double array: array of doubles containing the state
        The first half of the array contains the real part of the amplitudes
        The second half contains the imaginary part of the amplitudes
"""

_cSetState = _qsimov.setState
_cSetState.argtypes = [ct.c_void_p, c_double_p, ct.c_int]
_cSetState.restype = ct.c_int
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

_cGetNQubits = _qsimov.getNQubits
_cGetNQubits.argtypes = [ct.c_void_p]
_cGetNQubits.restype = ct.c_int
"""C QRegistry getNQubits function.
Positional arguments:
    pointer -> C QRegistry
Return:
    int: number of qubits in the registry
"""

_cGetSize = _qsimov.getSize
_cGetSize.argtypes = [ct.c_void_p]
_cGetSize.restype = ct.c_int
"""C QRegistry getSize function.
Positional arguments:
    pointer -> C QRegistry
Return:
    int: number of amplitudes in the state vector of the registry
"""

_cDensityMat = _qsimov.DensityMat
_cDensityMat.argtypes = [ct.c_void_p]
_cDensityMat.restype = ct.c_void_p
"""C QRegistry DensityMat function.
Positional arguments:
    pointer -> C QRegistry
Return:
    pointer to a functional matrix representing the density matrix
"""

_cFreeState = _qsimov.freeState
_cFreeState.argtypes = [ct.c_void_p]
_cFreeState.restype = ct.c_int
"""C QRegistry freeState function.
Positional arguments:
    pointer -> C QRegistry
Return:
    int: 0 -> Failure. 1 -> Success
"""

_partialTrace = _funmat.partial_trace
_partialTrace.argtypes = [ct.c_void_p, ct.c_int]
_partialTrace.restype = ct.c_void_p
"""C QRegistry partial_trace function.
Positional arguments:
    pointer -> C QRegistry
    int -> id of the qubit to trace out
Return:
    pointer to a functional matrix representing the reduced density matrix
"""


_cJoinStates = _qsimov.joinStates
_cJoinStates.argtypes = [ct.c_void_p, ct.c_void_p]
_cJoinStates.restype = ct.c_void_p
"""C QRegistry joinStates function.
Positional arguments:
    pointer -> C QRegistry A
    pointer -> C QRegistry B
Return:
    pointer to a C Qregistry A tensor product B
"""


class QRegistry:
    """Quantum Registry, base of all quantum related operations."""

    def __init__(self, nqbits):
        """Initialize QRegistry to state 0.

        nqbits -> number of QuBits in the registry.
        """
        self.reg = ct.c_void_p(_new_QRegistry(nqbits))

    def __del__(self):
        """Release memory held by the QRegistry."""
        _cFreeState(self.reg)
        self.reg = None

    '''
    def get_memory(self):
        """Return memory allocated by C structure in bytes."""
        return int(_cGetRegMemory(self.reg))
    '''

    def toString(self):
        """Use state_string method instead. DEPRECATED."""
        print("Method QRegistry.toString is deprecated.",
              "Please use state_string if you seek the same functionality")
        return self.state_string()

    def state_string(self):
        """Return string representation of the state of this registry."""
        print("state_string might be removed in future versions")
        return _cToString(self.reg).decode("ascii")

    def getSize(self):
        """Use get_state_size method instead. DEPRECATED."""
        print("Method QRegistry.getSize is deprecated.",
              "Please use get_state_size if you seek the same functionality")
        return self.get_state_size()

    def get_state_size(self):
        """Return the number of elements in the state vector."""
        return int(_cGetSize(self.reg))

    def getNQubits(self):
        """Use get_num_qubits method instead. DEPRECATED."""
        print("Method QRegistry.getNQubits is deprecated.",
              "Please use get_num_qubits if you seek the same functionality")
        return self.get_num_qubits()

    def get_num_qubits(self):
        """Return the number of qubits in this registry."""
        return int(_cGetNQubits(self.reg))

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
        nqubits = self.get_num_qubits()
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
            if int(_cMeasure(self.reg, int_array(*mask), ct.c_int(nq),
                             result, ct.c_int(rem))) == 0:
                print("Error measuring!")
            else:
                return list(result)[::-1]
        else:
            return []

    def applyGate(self, *args, **kwargs):
        """Use apply_gate method instead. DEPRECATED."""
        print("Method QRegistry.applyGate is deprecated.",
              "Please use apply_gate if you seek the same functionality")
        return self.apply_gate(*args, **kwargs)

    def apply_gate(self, gate, qubit=0, control=None, anticontrol=None,
                   optimize=True):
        """Apply specified gate to specified qubit with specified controls.

        Positional arguments:
            gate: string with the name of the gate to apply, or a QGate
        Keyworded arguments:
            qubit: id of the least significant qubit the gate will target
            control: id or list of ids of the qubit that will act as controls
            anticontrol: id or list of ids of the qubit that will act as
                         anticontrols
            optimize: only for QGates. Whether to use optimized lines or
                      user defined lines
        """
        if (np.issubdtype(type(qubit), np.integer)
                and qubit < self.get_num_qubits() and qubit >= 0):
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
                    allOk = all(0 <= id < self.get_num_qubits()
                                for id in control)

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
                    allOk = allOk and all(0 <= id < self.get_num_qubits()
                                          for id in anticontrol)

                if allOk:
                    if type(gate) == str:
                        name, arg1, arg2, arg3, invert = get_gate_data(gate)
                        if arg1 is None:
                            arg1 = 0
                        if arg2 is None:
                            arg2 = 0
                        if arg3 is None:
                            arg3 = 0
                        if (name.lower() == "swap"):
                            _swap_qubits(self.reg, arg1, arg2, False, False,
                                         invert, control, clen,
                                         anticontrol, aclen)
                        elif (name.lower() == "iswap"):
                            _swap_qubits(self.reg, arg1, arg2, True, False,
                                         invert, control, clen,
                                         anticontrol, aclen)
                        elif (name.lower() == "sqrtswap"):
                            _swap_qubits(self.reg, arg1, arg2, False, True,
                                         invert, control, clen,
                                         anticontrol, aclen)
                        elif (name.lower() == "xx"):
                            _ising_qubits(self.reg, 0, arg1, arg2, arg3,
                                          invert, control, clen,
                                          anticontrol, aclen)
                        elif (name.lower() == "yy"):
                            _ising_qubits(self.reg, 1, arg1, arg2, arg3,
                                          invert, control, clen,
                                          anticontrol, aclen)
                        elif (name.lower() == "zz"):
                            _ising_qubits(self.reg, 2, arg1, arg2, arg3,
                                          invert, control, clen,
                                          anticontrol, aclen)
                        else:
                            if int(_cApplyGate(self.reg,
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
                    else:
                        gate._apply_gate(self, qubit,
                                         control[:clen], anticontrol[:aclen],
                                         optimize=optimize)
                else:
                    print("The ids must be between 0 and " +
                          str(self.get_num_qubits()))
        else:
            if not np.issubdtype(type(qubit), np.integer):
                print("Qubit must be of integer type!")
            elif qubit >= self.get_num_qubits() or qubit < 0:
                print("The specified qubit doesn't exist!")

    def getState(self):
        """Use get_state method instead. DEPRECATED."""
        print("Method QRegistry.getState is deprecated.",
              "Please use get_state if you seek the same functionality")
        return self.get_state()

    def get_state(self):
        """Return the state vector of the registry."""
        rawState = _cGetState(self.reg)
        size = self.get_state_size()
        state = np.array([complex(rawState[i], rawState[i + size])
                          for i in range(size)])
        _cFree(rawState)

        return state

    def setState(self, newState, nqubits):
        """Use _set_state method instead. DEPRECATED."""
        print("Method QRegistry.setState is deprecated.",
              "Please use _set_state if you seek the same functionality")
        return self._set_state(newState, nqubits)

    def _set_state(self, newState, nqubits):
        """Set the state vector of the registry."""
        size = 2 << (nqubits - 1)  # 2^n_qubits
        statetype = ct.c_double * (size * 2)
        cstate = statetype(*[newState[i].real if i < size
                             else newState[i - size].imag
                             for i in range(size * 2)])
        result = _cSetState(self.reg, cstate, ct.c_int(nqubits)) == 1

        return result

    def densityMatrix(self):
        """Use density_matrix method instead. DEPRECATED."""
        print("Method QRegistry.densityMatrix is deprecated.",
              "Please use density_matrix if you seek the same functionality")
        return self.density_matrix()

    def density_matrix(self):
        """Return functional matrix of the density matrix."""
        return Funmatrix(ct.c_void_p(_cDensityMat(self.reg)), "Rho")

    def reducedDensityMatrix(self, elem):
        """Use reduced_density_matrix method instead. DEPRECATED."""
        print("Method QRegistry.reducedDensityMatrix is deprecated.",
              "Please use reduced_density_matrix if you seek the",
              "same functionality")
        return self.reduced_density_matrix(elem)

    def reduced_density_matrix(self, elem):
        """Return functional matrix of the reduced density matrix."""
        nq = self.get_num_qubits()
        if 0 <= elem < nq:
            # elem = nq - elem - 1
            c_density_mat = ct.c_void_p(_cDensityMat(self.reg))
            c_reduced_mat = ct.c_void_p(_partialTrace(c_density_mat,
                                                      ct.c_int(elem)))
            return Funmatrix(c_reduced_mat, "Tr_" + str(elem) + "(Rho)")
        else:
            print("The specified QuBit doesn't exist in this registry!")

    def reducedTrace(self, elem):
        """Use reduced_trace method instead. DEPRECATED."""
        print("Method QRegistry.reducedTrace is deprecated.",
              "Please use reduced_trace if you seek the same functionality")
        return self.reduced_trace(elem)

    def reduced_trace(self, elem):
        """Trace of the square reduced density matrix."""
        rho_a = np.array(self.reduced_density_matrix(elem)[:])
        rt = (rho_a @ rho_a).trace().real
        return rt if rt <= 1.0 else 1.0

    def vnEntropy(self, **kwargs):
        """Use get_entropy method instead. DEPRECATED."""
        print("Method QRegistry.vnEntropy is deprecated.",
              "Please use get_entropy if you seek the same functionality")
        return self.get_entropy(**kwargs)

    def get_entropy(self, **kwargs):
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
        size = self.get_state_size()
        for i in range(size):
            p = self.prob(i)
            if p > 0:
                if base == "e":
                    entropy += p * np.log(p)
                elif type(base) == int or type(base) == float:
                    entropy += p * np.log(p)/np.log(base)
        return -entropy

    def blochCoords(self):
        """Use get_bloch_coords method instead. DEPRECATED."""
        print("Method QRegistry.blochCoords is deprecated.",
              "Please use get_bloch_coords if you seek the same functionality")
        return self.get_bloch_coords()

    def get_bloch_coords(self):
        """Get the polar coordinates of ONE qubit in the bloch sphere."""
        if self.get_num_qubits() == 1:
            return _cBlochCoords(self.reg)[:2]
        else:
            print("You can only get bloch coordinates for 1 qubit registries")
            return None

    def bra(self):
        """Get the conjugated row form state vector (bra <v|)."""
        k = np.array(self.get_state())
        k.shape = (1, k.shape[0])
        return np.conjugate(k)

    def ket(self):
        """Get the column form state vector (ket |v>)."""
        k = np.array(self.get_state())
        k.shape = (k.shape[0], 1)
        return k

    def prob(self, x):
        """Get the odds of measuring a certain state x (not qubit)."""
        p = 0
        if (x < self.get_state_size()):
            p = self.density_matrix()[x, x].real
        return p


def superposition(a, b):
    """Join two registries into one by calculating tensor product."""
    r = QRegistry(1)
    _cFreeState(r.reg)
    r.reg = ct.c_void_p(_cJoinStates(a.reg, b.reg))
    return r


def _swap_qubits(reg, qubit1, qubit2, imaginary, sqrt, invert,
                 control, clen, anticontrol, aclen):
    """Apply a type of SWAP gate to the registry."""
    if int(_cSWAPQubits(reg, ct.c_int(qubit1), ct.c_int(qubit2),
                        ct.c_int(int(imaginary)), ct.c_int(int(sqrt)),
                        ct.c_int(int(invert)), control, ct.c_int(clen),
                        anticontrol, ct.c_int(aclen))) == 0:
        print("Error swapping specified QuBits!")


def _ising_qubits(reg, type, angle, qubit1, qubit2, invert,
                  control, clen, anticontrol, aclen):
    """Apply a type of Ising gate to the registry."""
    if int(_cIsingQubits(reg, ct.c_int(type), ct.c_double(angle),
                         ct.c_int(qubit1), ct.c_int(qubit2),
                         ct.c_int(int(invert)), control, ct.c_int(clen),
                         anticontrol, ct.c_int(aclen))) == 0:
        print("Error coupling specified QuBits!")
