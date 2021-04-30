"""Module with functional matrices related stuff."""
import ctypes as ct
import qsimov.connectors.parser as prs
import platform as plat
import os
import numpy as np
from os.path import sep

# DLL Load
if plat.system() == "Windows":
    extension = ".dll"
else:
    extension = ".so"
_root_folder = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
_lib_folder = _root_folder + sep + "lib"
_funmat_path = _lib_folder + sep + "libfunmat" + extension
if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(_root_folder)
    os.add_dll_directory(_lib_folder)
    _funmat_path = "libfunmat" + extension
_funmat = ct.CDLL(_funmat_path)

c_double_p = ct.POINTER(ct.c_double)

_cGetItem = _funmat.getitemaux
_cGetItem.argtypes = [ct.c_void_p, ct.c_uint, ct.c_uint,
                      c_double_p, c_double_p]
_cGetItem.restype = ct.c_int
"""C FunMatrix getitemaux function.
Positional arguments:
    pointer: FunMatrix struct pointer
    unsigned int: index of the row to get
    unsigned int: index of the column to get
    double pointer: real part of the result (output)
    double pointer: imaginary part of the result (output)
Return:
    int: 0 -> Failure. 1 -> Success
"""

_cGetRows = _funmat.rows
_cGetRows.argtypes = [ct.c_void_p]
_cGetRows.restype = ct.c_int
"""C FunMatrix rows function.
Positional arguments:
    pointer: FunMatrix struct pointer
Return:
    int: number of rows of the FunMatrix
"""

_cGetCols = _funmat.columns
_cGetCols.argtypes = [ct.c_void_p]
_cGetCols.restype = ct.c_int
"""C FunMatrix columns function.
Positional arguments:
    pointer: FunMatrix struct pointer
Return:
    int: number of columns of the FunMatrix
"""

_cMadd = _funmat.madd
_cMadd.argtypes = [ct.c_void_p, ct.c_void_p]
_cMadd.restype = ct.c_void_p
"""C FunMatrix madd function.
Positional arguments:
    pointer: FunMatrix struct pointer A
    pointer: FunMatrix struct pointer B
Return:
    pointer: FunMatrix struct pointer A + B
"""

_cMsub = _funmat.msub
_cMsub.argtypes = [ct.c_void_p, ct.c_void_p]
_cMsub.restype = ct.c_void_p
"""C FunMatrix msub function.
Positional arguments:
    pointer: FunMatrix struct pointer A
    pointer: FunMatrix struct pointer B
Return:
    pointer: FunMatrix struct pointer A - B
"""

_cMprod = _funmat.mprodaux
_cMprod.argtypes = [ct.c_double, ct.c_double, ct.c_void_p]
_cMprod.restype = ct.c_void_p
"""C FunMatrix mprodaux function.
Positional arguments:
    double: real part of the scalar value s
    double: imaginary part of the scalar value s
    pointer: FunMatrix struct pointer A
Return:
    pointer: FunMatrix struct pointer s * A
"""

_cEwmul = _funmat.ewmul
_cEwmul.argtypes = [ct.c_void_p, ct.c_void_p]
_cEwmul.restype = ct.c_void_p
"""C FunMatrix ewmul function.
Positional arguments:
    pointer: FunMatrix struct pointer A
    pointer: FunMatrix struct pointer B
Return:
    pointer: FunMatrix struct pointer A * B (entity wise multiplication)
"""

_cMdiv = _funmat.mdivaux
_cMdiv.argtypes = [ct.c_double, ct.c_double, ct.c_void_p]
_cMdiv.restype = ct.c_void_p
"""C FunMatrix mdivaux function.
Positional arguments:
    double: real part of the scalar value s
    double: imaginary part of the scalar value s
    pointer: FunMatrix struct pointer A
Return:
    pointer: FunMatrix struct pointer A/s
"""

_cMatmul = _funmat.matmul
_cMatmul.argtypes = [ct.c_void_p, ct.c_void_p]
_cMatmul.restype = ct.c_void_p
"""C FunMatrix matmul function.
Positional arguments:
    pointer: FunMatrix struct pointer A
    pointer: FunMatrix struct pointer B
Return:
    pointer: FunMatrix struct pointer A . B (matrix multiplication)
"""

_cKron = _funmat.kron
_cKron.argtypes = [ct.c_void_p, ct.c_void_p]
_cKron.restype = ct.c_void_p
"""C FunMatrix kron function.
Positional arguments:
    pointer: FunMatrix struct pointer A
    pointer: FunMatrix struct pointer B
Return:
    pointer: FunMatrix struct pointer A x B (tensor product)
"""

_cTranspose = _funmat.transpose
_cTranspose.argtypes = [ct.c_void_p]
_cTranspose.restype = ct.c_void_p
"""C FunMatrix transpose function.
Positional arguments:
    pointer: FunMatrix struct pointer A
Return:
    pointer: FunMatrix struct pointer A transposed
"""

_cDagger = _funmat.dagger
_cDagger.argtypes = [ct.c_void_p]
_cDagger.restype = ct.c_void_p
"""C FunMatrix dagger function.
Positional arguments:
    pointer: FunMatrix struct pointer A
Return:
    pointer: FunMatrix struct pointer A transposed and conjugated
"""

_cGetMemory = _funmat.getMemory
_cGetMemory.argtypes = [ct.c_void_p]
_cGetMemory.restype = ct.c_int
"""C FunMatrix getMemory function.
Positional arguments:
    pointer: FunMatrix struct pointer A
Return:
    int: size allocated by the FunMatrix structure in bytes
"""


class Funmatrix(object):
    """Functional Matrices related stuff in python."""

    def __init__(self, fmatrix, name="Lazy Matrix"):
        """Functional Matrix constructor.

        Positional arguments:
            fmatrix: pointer to C FunMatrix structure
        Keyworded arguments:
            name: string with the name of the matrix
        """
        self.m = fmatrix
        self.name = name
        self._shape = (int(_cGetRows(fmatrix)), int(_cGetCols(fmatrix)))

    def __repr__(self):
        """Return the string representation of the FunMatrix."""
        return self.name

    def __str__(self):
        """Return the string representation of the FunMatrix."""
        return self.name

    def _get_single_item(self, i, j):
        """Return the item in (i, j) position of the matrix."""
        res = None
        re = ct.c_double(0.0)
        im = ct.c_double(0.0)
        try:
            aux = _cGetItem(self.m, ct.c_uint(i), ct.c_uint(j),
                            c_double_p(re), c_double_p(im))
            if (aux == 0):
                print("Error getting the specified item!")
            else:
                res = complex(re.value, im.value)
        except Exception as e:
            print(e)
        return res

    def __getitem__(self, key):
        """Return the items of the matrix specified in key."""
        rows, cols = self.shape

        if (type(key) == int):
            if (key >= rows):
                raise IndexError("index " + str(key) + " is out of bounds " +
                                 "for axis 0 with shape " + str(rows))
            if (key < 0):
                key = rows + key
            return np.array([self._get_single_item(key, j)
                             for j in range(cols)])
        elif (type(key) == slice):
            start = key.start if key.start is not None else 0
            if (start < 0):
                start = rows + start
            stop = key.stop if key.stop is not None else rows
            if (stop < 0):
                stop = rows + stop
            step = key.step if key.step is not None else 1
            return np.array([self.__getitem__(i)
                             for i in range(start, stop, step)])
        elif (type(key) == tuple):
            if (type(key[1]) is None):
                return self.__getitem__(key[0])
            elif (type(key[0]) is None):
                if (type(key[1]) == int):
                    if (key[1] >= cols):
                        raise IndexError("index " + str(key[1]) + " is out " +
                                         "of bounds for axis 1 with shape " +
                                         str(cols))
                    if (key[1] < 0):
                        key = (key[0], cols + key[1])
                    item_list = [np.around(self._get_single_item(i, key[1]),
                                           decimals=self.dec)
                                 for i in range(rows)]
                    return np.array(item_list)
                elif (type(key[1]) == slice):
                    start = key[1].start if key[1].start is not None else 0
                    if (start < 0):
                        start = cols + start
                    stop = key[1].stop if key[1].stop is not None else cols
                    if (stop < 0):
                        stop = cols + stop
                    step = key[1].step if key[1].step is not None else 1
                    return np.array([self.__getitem__(None, i)
                                     for i in range(start, stop, step)])
            elif (type(key[0]) == int and type(key[1]) == int):
                if (key[0] >= rows):
                    raise IndexError("index " + str(key[0]) + " is out of " +
                                     "bounds for axis 0 with shape " +
                                     str(rows))
                if (key[1] >= cols):
                    raise IndexError("index " + str(key[1]) + " is out of " +
                                     "bounds for axis 1 with shape " +
                                     str(cols))
                if (key[0] < 0):
                    key = (rows + key[0], key[1])
                if (key[1] < 0):
                    key = (key[0], cols + key[1])
                return self._get_single_item(key[0], key[1])
            elif (type(key[0]) == slice and type(key[1]) == int):
                if (key[1] >= cols):
                    raise IndexError("index " + str(key[1]) + " is out of " +
                                     "bounds for axis 1 with shape " +
                                     str(cols))
                start = key[0].start if key[0].start is not None else 0
                if (start < 0):
                    start = rows + start
                stop = key[0].stop if key[0].stop is not None else rows
                if (stop < 0):
                    stop = rows + stop
                step = key[0].step if key[0].step is not None else 1
                if (key[1] < 0):
                    key = (key[0], cols + key[1])
                return np.array([self._get_single_item(i, key[1])
                                 for i in range(start, stop, step)])
            elif (type(key[0]) == int and type(key[1]) == slice):
                if (key[0] >= rows):
                    raise IndexError("index " + str(key[0]) + " is out of " +
                                     "bounds for axis 0 with shape " +
                                     str(rows))
                start = key[1].start if key[1].start is not None else 0
                if (start < 0):
                    start = cols + start
                stop = key[1].stop if key[1].stop is not None else cols
                if (stop < 0):
                    stop = cols + stop
                step = key[1].step if key[1].step is not None else 1
                if (key[0] < 0):
                    key = (rows + key[0], key[1])
                return np.array([self._get_single_item(key[0], i)
                                 for i in range(start, stop, step)])
            elif (type(key[0]) == slice and type(key[1]) == slice):
                start0 = key[0].start if key[0].start is not None else 0
                if (start0 < 0):
                    start0 = rows + start0
                stop0 = key[0].stop if key[0].stop is not None else rows
                if (stop0 < 0):
                    stop0 = rows + stop0
                step0 = key[0].step if key[0].step is not None else 1
                start1 = key[1].start if key[1].start is not None else 0
                if (start1 < 0):
                    start1 = cols + start1
                stop1 = key[1].stop if key[1].stop is not None else cols
                if (stop1 < 0):
                    stop1 = cols + stop1
                step1 = key[1].step if key[1].step is not None else 1
                return np.array([[self._get_single_item(i, j)
                                  for i in range(start0, stop0, step0)]
                                 for j in range(start1, stop1, step1)])

    def __add__(self, other):
        """Return the result of adding self + other matrices."""
        if type(other) != Funmatrix:
            raise TypeError("unsupported operand type(s) for +: '" +
                            type(self).__name__ + "' and '" +
                            type(self).__name__ + "'")
        res = None
        try:
            res = Funmatrix(ct.c_void_p(_cMadd(self.m, other.m)))
        except Exception as e:
            print(e)
        return res

    def __sub__(self, other):
        """Return the result of substracting self - other matrices."""
        if type(other) != Funmatrix:
            raise TypeError("unsupported operand type(s) for -: '" +
                            type(self).__name__ + "' and '" +
                            type(self).__name__ + "'")
        res = None
        try:
            res = Funmatrix(ct.c_void_p(_cMsub(self.m, other.m)))
        except Exception as e:
            print(e)
        return res

    def __mul__(self, other):
        """Calculate entity wise self * other or scalar product."""
        if type(other) == Funmatrix:
            res = None
            try:
                res = Funmatrix(ct.c_void_p(_cEwmul(self.m, other.m)))
            except Exception as e:
                print(e)
            return res
        elif (type(other) == int
                or type(other) == float
                or type(other) == complex):
            res = None
            re = 0
            im = 0
            if type(other) == int or type(other) == float:
                re = other
            else:
                re = other.real
                im = other.imag
            try:
                res = Funmatrix(ct.c_void_p(_cMprod(ct.c_double(re),
                                                    ct.c_double(im),
                                                    self.m)))
            except Exception as e:
                print(e)
            return res
        else:
            raise TypeError("unsupported operand type(s) for *: '" +
                            type(self).__name__ + "' and '" +
                            type(self).__name__ + "'")

    def __rmul__(self, other):
        """Calculate entity wise other * self or scalar product."""
        if type(other) == Funmatrix:
            res = None
            try:
                res = Funmatrix(ct.c_void_p(_cEwmul(other.m, self.m)))
            except Exception as e:
                print(e)
            return res
        elif (type(other) == int
                or type(other) == float
                or type(other) == complex):
            res = None
            re = 0
            im = 0
            if type(other) == int or type(other) == float:
                re = other
            else:
                re = other.real
                im = other.imag
            try:
                res = Funmatrix(ct.c_void_p(_cMprod(ct.c_double(re),
                                                    ct.c_double(im),
                                                    self.m)))
            except Exception as e:
                print(e)
            return res
        else:
            raise TypeError("unsupported operand type(s) for *: '" +
                            type(self).__name__ + "' and '" +
                            type(self).__name__ + "'")

    def __truediv__(self, other):
        """Calculate scalar division self / scalar."""
        if (type(other) == int
                or type(other) == float
                or type(other) == complex):
            res = None
            re = 0
            im = 0
            if type(other) == int or type(other) == float:
                re = other
            else:
                re = other.real
                im = other.imag
            res = None
            try:
                res = Funmatrix(ct.c_void_p(_cMdiv(ct.c_double(re),
                                                   ct.c_double(im),
                                                   self.m)))
            except Exception as e:
                print(e)
            return res
        else:
            raise TypeError("unsupported operand type(s) for /: '" +
                            type(self).__name__ + "' and '" +
                            type(self).__name__ + "'")

    def __matmul__(self, other):
        """Calculate matrix multiplication self . other."""
        if type(other) != Funmatrix:
            raise TypeError("unsupported operand type(s) for @: '" +
                            type(self).__name__ + "' and '" +
                            type(self).__name__ + "'")
        res = None
        try:
            res = Funmatrix(ct.c_void_p(_cMatmul(self.m, other.m)))
        except Exception as e:
            print(e)
        return res

    def __pow__(self, other):
        """Calculate matrix tensor product self x other."""
        if type(other) != Funmatrix:
            raise TypeError("unsupported operand type(s) for **: '" +
                            type(self).__name__ + "' and '" +
                            type(self).__name__ + "'")
        res = None
        try:
            res = Funmatrix(ct.c_void_p(_cKron(self.m, other.m)))
        except Exception as e:
            print(e)
        return res

    def transpose(self):
        """Return the transpose of this matrix."""
        res = None
        try:
            res = Funmatrix(ct.c_void_p(_cTranspose(self.m)))
        except Exception as e:
            print(e)
        return res

    def dagger(self):
        """Return the conjugate (Hermitian) transpose of this matrix."""
        res = None
        if prs.getGroups(self.name)[1] in ["I", "H", "X", "NOT",
                                           "Y", "Z", "SWAP"]:
            res = self
        else:
            try:
                res = Funmatrix(ct.c_void_p(_cDagger(self.m)),
                                self.name + "-1")
            except Exception as e:
                print(e)
        return res

    def get_matrix(self):
        """Return self."""
        return self

    def get_memory(self):
        """Return memory allocated by C structure in bytes."""
        return int(_cGetMemory(self.m))

    @property
    def shape(self):
        """Return tuple with (number of rows, number of columns)."""
        return self._shape
