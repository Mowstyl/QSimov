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
__libfolder__ = os.getcwd() + sep + "qsimov" + sep + "lib"
__funmatpath__ = __libfolder__ + sep + "libfunmat" + extension
if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(os.getcwd() + sep + "qsimov")
    os.add_dll_directory(__libfolder__)
    __funmatpath__ = "libfunmat" + extension
__funmat__ = ct.CDLL(__funmatpath__)

c_double_p = ct.POINTER(ct.c_double)

__cGetItem__ = __funmat__.getitemaux
__cGetItem__.argtypes = [ct.c_void_p, ct.c_uint, ct.c_uint, c_double_p, c_double_p]
__cGetItem__.restype = ct.c_int

__cGetRows__ = __funmat__.rows
__cGetRows__.argtypes = [ct.c_void_p]
__cGetRows__.restype = ct.c_int

__cGetCols__ = __funmat__.columns
__cGetCols__.argtypes = [ct.c_void_p]
__cGetCols__.restype = ct.c_int

__cMadd__ = __funmat__.madd
__cMadd__.argtypes = [ct.c_void_p, ct.c_void_p]
__cMadd__.restype = ct.c_void_p

__cMsub__ = __funmat__.msub
__cMsub__.argtypes = [ct.c_void_p, ct.c_void_p]
__cMsub__.restype = ct.c_void_p

__cMprod__ = __funmat__.mprodaux
__cMprod__.argtypes = [ct.c_double, ct.c_double, ct.c_void_p]
__cMprod__.restype = ct.c_void_p

__cEwmul__ = __funmat__.ewmul
__cEwmul__.argtypes = [ct.c_void_p, ct.c_void_p]
__cEwmul__.restype = ct.c_void_p

__cMdiv__ = __funmat__.mdivaux
__cMdiv__.argtypes = [ct.c_double, ct.c_double, ct.c_void_p]
__cMdiv__.restype = ct.c_void_p

__cMatmul__ = __funmat__.matmul
__cMatmul__.argtypes = [ct.c_void_p, ct.c_void_p]
__cMatmul__.restype = ct.c_void_p

__cKron__ = __funmat__.kron
__cKron__.argtypes = [ct.c_void_p, ct.c_void_p]
__cKron__.restype = ct.c_void_p

__cTranspose__ = __funmat__.transpose
__cTranspose__.argtypes = [ct.c_void_p]
__cTranspose__.restype = ct.c_void_p

__cDagger__ = __funmat__.dagger
__cDagger__.argtypes = [ct.c_void_p]
__cDagger__.restype = ct.c_void_p

__cGetMemory__ = __funmat__.getMemory
__cGetMemory__.argtypes = [ct.c_void_p]
__cGetMemory__.restype = ct.c_int


class Funmatrix(object):
    def __init__(self, fmatrix, name="Lazy Matrix"):
        self.m = fmatrix
        self.name = name
        self._shape = (int(__cGetRows__(fmatrix)), int(__cGetCols__(fmatrix)))

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def __getSingleItem__(self, i, j):
        res = None
        re = ct.c_double(0.0)
        im = ct.c_double(0.0)
        try:
            aux = __cGetItem__(self.m, ct.c_uint(i), ct.c_uint(j), c_double_p(re), c_double_p(im))
            if (aux == 0):
                print("Error getting the specified item!")
            else:
                res = complex(re.value, im.value)
        except Exception as e:
            print(e)
        return res

    def __getitem__(self, key):
        rows, cols = self.shape

        if (type(key) == int):
            if (key >= rows):
                raise IndexError("index " + str(key) + " is out of bounds for axis 0 with shape " + str(rows))
            if (key < 0):
                key = rows + key
            return np.array([self.__getSingleItem__(key, j) for j in range(cols)])
        elif (type(key) == slice):
            start = key.start if key.start is not None else 0
            if (start < 0):
                start = rows + start
            stop = key.stop if key.stop is not None else rows
            if (stop < 0):
                stop = rows + stop
            step = key.step if key.step is not None else 1
            return np.array([self.__getitem__(i) for i in range(start, stop, step)])
        elif (type(key) == tuple):
            if (type(key[1]) is None):
                return self.__getitem__(key[0])
            elif (type(key[0]) is None):
                if (type(key[1]) == int):
                    if (key[1] >= cols):
                        raise IndexError("index " + str(key[1]) + " is out of bounds for axis 1 with shape " + str(cols))
                    if (key[1] < 0):
                        key = (key[0], cols + key[1])
                    return np.array([np.around(self.__getSingleItem__(i, key[1]), decimals=self.dec) for i in range(rows)])
                elif (type(key[1]) == slice):
                    start = key[1].start if key[1].start is not None else 0
                    if (start < 0):
                        start = cols + start
                    stop = key[1].stop if key[1].stop is not None else cols
                    if (stop < 0):
                        stop = cols + stop
                    step = key[1].step if key[1].step is not None else 1
                    return np.array([self.__getitem__(None, i) for i in range(start, stop, step)])
            elif (type(key[0]) == int and type(key[1]) == int):
                if (key[0] >= rows):
                    raise IndexError("index " + str(key[0]) + " is out of bounds for axis 0 with shape " + str(rows))
                if (key[1] >= cols):
                    raise IndexError("index " + str(key[1]) + " is out of bounds for axis 1 with shape " + str(cols))
                if (key[0] < 0):
                    key = (rows + key[0], key[1])
                if (key[1] < 0):
                    key = (key[0], cols + key[1])
                return self.__getSingleItem__(key[0], key[1])
            elif (type(key[0]) == slice and type(key[1]) == int):
                if (key[1] >= cols):
                    raise IndexError("index " + str(key[1]) + " is out of bounds for axis 1 with shape " + str(cols))
                start = key[0].start if key[0].start is not None else 0
                if (start < 0):
                    start = rows + start
                stop = key[0].stop if key[0].stop is not None else rows
                if (stop < 0):
                    stop = rows + stop
                step = key[0].step if key[0].step is not None else 1
                if (key[1] < 0):
                    key = (key[0], cols + key[1])
                return np.array([self.__getSingleItem__(i, key[1]) for i in range(start, stop, step)])
            elif (type(key[0]) == int and type(key[1]) == slice):
                if (key[0] >= rows):
                    raise IndexError("index " + str(key[0]) + " is out of bounds for axis 0 with shape " + str(rows))
                start = key[1].start if key[1].start is not None else 0
                if (start < 0):
                    start = cols + start
                stop = key[1].stop if key[1].stop is not None else cols
                if (stop < 0):
                    stop = cols + stop
                step = key[1].step if key[1].step is not None else 1
                if (key[0] < 0):
                    key = (rows + key[0], key[1])
                return np.array([self.__getSingleItem__(key[0], i) for i in range(start, stop, step)])
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
                return np.array([[self.__getSingleItem__(i, j) for i in range(start0, stop0, step0)] for j in range(start1, stop1, step1)])

    def __add__(self, other):
        if type(other) != Funmatrix:
            raise TypeError("unsupported operand type(s) for +: '" + type(self).__name__ + "' and '" + type(self).__name__ + "'")
        res = None
        try:
            res = Funmatrix(ct.c_void_p(__cMadd__(self.m, other.m)))
        except Exception as e:
            print(e)
        return res

    def __sub__(self, other):
        if type(other) != Funmatrix:
            raise TypeError("unsupported operand type(s) for -: '" + type(self).__name__ + "' and '" + type(self).__name__ + "'")
        res = None
        try:
            res = Funmatrix(ct.c_void_p(__cMsub__(self.m, other.m)))
        except Exception as e:
            print(e)
        return res

    def __mul__(self, other):
        if type(other) == Funmatrix:
            res = None
            try:
                res = Funmatrix(ct.c_void_p(__cEwmul__(self.m, other.m)))
            except Exception as e:
                print(e)
            return res
        elif type(other) == int or type(other) == float or type(other) == complex:
            res = None
            re = 0
            im = 0
            if type(other) == int or type(other) == float:
                re = other
            else:
                re = other.real
                im = other.imag
            try:
                res = Funmatrix(ct.c_void_p(__cMprod__(ct.c_double(re), ct.c_double(im), self.m)))
            except Exception as e:
                print(e)
            return res
        else:
            raise TypeError("unsupported operand type(s) for *: '" + type(self).__name__ + "' and '" + type(self).__name__ + "'")

    def __rmul__(self, other):
        if type(other) == Funmatrix:
            res = None
            try:
                res = Funmatrix(ct.c_void_p(__cEwmul__(other.m, self.m)))
            except Exception as e:
                print(e)
            return res
        elif type(other) == int or type(other) == float or type(other) == complex:
            res = None
            re = 0
            im = 0
            if type(other) == int or type(other) == float:
                re = other
            else:
                re = other.real
                im = other.imag
            try:
                res = Funmatrix(ct.c_void_p(__cMprod__(ct.c_double(re), ct.c_double(im), self.m)))
            except Exception as e:
                print(e)
            return res
        else:
            raise TypeError("unsupported operand type(s) for *: '" + type(self).__name__ + "' and '" + type(self).__name__ + "'")

    def __truediv__(self, other):
        if type(other) == int or type(other) == float or type(other) == complex:
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
                res = Funmatrix(ct.c_void_p(__cMdiv__(ct.c_double(re), ct.c_double(im), self.m)))
            except Exception as e:
                print(e)
            return res
        else:
            raise TypeError("unsupported operand type(s) for /: '" + type(self).__name__ + "' and '" + type(self).__name__ + "'")

    def __matmul__(self, other):
        if type(other) != Funmatrix:
            raise TypeError("unsupported operand type(s) for @: '" + type(self).__name__ + "' and '" + type(self).__name__ + "'")
        res = None
        try:
            res = Funmatrix(ct.c_void_p(__cMatmul__(self.m, other.m)))
        except Exception as e:
            print(e)
        return res

    def __pow__(self, other):
        if type(other) != Funmatrix:
            raise TypeError("unsupported operand type(s) for **: '" + type(self).__name__ + "' and '" + type(self).__name__ + "'")
        res = None
        try:
            res = Funmatrix(ct.c_void_p(__cKron__(self.m, other.m)))
        except Exception as e:
            print(e)
        return res

    def transpose(self):
        res = None
        try:
            res = Funmatrix(ct.c_void_p(__cTranspose__(self.m)))
        except Exception as e:
            print(e)
        return res

    def invert(self):
        res = None
        if prs.getGroups(self.name)[1] in ["I", "H", "X", "NOT", "Y", "Z", "SWAP"]:
            res = self
        else:
            try:
                res = Funmatrix(ct.c_void_p(__cDagger__(self.m)), self.name + "-1")
            except Exception as e:
                print(e)
        return res

    def getMatrix(self):
        return self

    def getsizeof(self):
        return int(__cGetMemory__(self.m))

    @property
    def shape(self):
        return self._shape
