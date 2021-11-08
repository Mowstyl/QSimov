"""Module with functional matrices related stuff."""
import doki
import numpy as np


class Funmatrix(object):
    """Functional Matrices related stuff in python."""

    def __init__(self, fmatrix, name="Lazy Matrix", verbose=False):
        """Functional Matrix constructor.

        Positional arguments:
            fmatrix: pointer to C FunMatrix structure
        Keyworded arguments:
            name: string with the name of the matrix
        """
        self.m = fmatrix
        self.name = name
        self._shape = doki.funmatrix_shape(fmatrix, verbose)

    def __repr__(self):
        """Return the string representation of the FunMatrix."""
        return self.name

    def __str__(self):
        """Return the string representation of the FunMatrix."""
        return self.name

    def _get_single_item(self, i, j, verbose=False):
        """Return the item in (i, j) position of the matrix."""
        return doki.funmatrix_get(self.m, i, j, verbose)

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
                            type(other).__name__ + "'")
        if self.shape != other.shape:
            raise ValueError("Both matrices must have the same shape")
        return Funmatrix(doki.funmatrix_add(self.m, other.m, False))

    def __sub__(self, other):
        """Return the result of substracting self - other matrices."""
        if type(other) != Funmatrix:
            raise TypeError("unsupported operand type(s) for -: '" +
                            type(self).__name__ + "' and '" +
                            type(other).__name__ + "'")
        if self.shape != other.shape:
            raise ValueError("Both matrices must have the same shape")
        return Funmatrix(doki.funmatrix_sub(self.m, other.m, False))

    def __mul__(self, other):
        """Calculate entity wise self * other or scalar product."""
        if type(other) == Funmatrix:
            if self.shape != other.shape:
                raise ValueError("Both matrices must have the same shape")
            return Funmatrix(doki.funmatrix_ewmul(self.m, other.m, False))
        elif np.isscalar(other):
            scalar = complex(other)
            return Funmatrix(doki.funmatrix_scalar_mul(self.m, scalar, False))
        else:
            raise TypeError("unsupported operand type(s) for *: '" +
                            type(self).__name__ + "' and '" +
                            type(other).__name__ + "'")

    def __rmul__(self, other):
        """Calculate entity wise other * self or scalar product."""
        if type(other) == Funmatrix:
            if self.shape != other.shape:
                raise ValueError("Both matrices must have the same shape")
            return Funmatrix(doki.funmatrix_ewmul(other.m, self.m, False))
        elif np.isscalar(other):
            scalar = complex(other)
            return Funmatrix(doki.funmatrix_scalar_mul(self.m, scalar, False))
        else:
            raise TypeError("unsupported operand type(s) for *: '" +
                            type(other).__name__ + "' and '" +
                            type(self).__name__ + "'")

    def __truediv__(self, other):
        """Calculate scalar division self / scalar."""
        if np.isscalar(other):
            scalar = complex(other)
            return Funmatrix(doki.funmatrix_scalar_div(self.m, scalar, False))
        else:
            raise TypeError("unsupported operand type(s) for /: '" +
                            type(self).__name__ + "' and '" +
                            type(other).__name__ + "'")

    def __matmul__(self, other):
        """Calculate matrix multiplication self . other."""
        if type(other) != Funmatrix:
            raise TypeError("unsupported operand type(s) for @: '" +
                            type(self).__name__ + "' and '" +
                            type(other).__name__ + "'")
        if self.shape[1] != other.shape[0]:
            raise ValueError("Matrix product is only defined when the number" +
                             " of columns of the first operand equals the" +
                             " number of rows of the second operand")
        return Funmatrix(doki.funmatrix_matmul(self.m, other.m, False))

    def __pow__(self, other):
        """Calculate matrix tensor product self x other or usual power."""
        if type(other) == Funmatrix:
            return Funmatrix(doki.funmatrix_kron(self.m, other.m, False))
        elif np.issubdtype(other, np.integer):
            exponent = int(other)
            if exponent < 0:
                raise ValueError("Negative exponents are not supported")
            res = self
            if exponent == 0:
                res = 1
            for i in range(exponent - 1):
                res = res @ self
            return res
        else:
            raise TypeError("unsupported operand type(s) for **: '" +
                            type(self).__name__ + "' and '" +
                            type(self).__name__ + "'")

    def transpose(self):
        """Return the transpose of this matrix."""
        return Funmatrix(doki.funmatrix_transpose(self.m, False))

    def dagger(self):
        """Return the conjugate (Hermitian) transpose of this matrix."""
        return Funmatrix(doki.funmatrix_dagger(self.m, False))

    def partial_trace(self, elem):
        """Return the partial trace of this matrix, after tracing out elem."""
        if self.shape[0] != self.shape[1]:
            raise ValueError("Partial trace undefined for non square matrices")
        num_elems = np.log2(self.shape[0])
        if not np.allclose(num_elems % 1, 0):
            raise NotImplementedError("Partial trace unsupported for " +
                                      "shapes that are not a power of 2")
        num_elems = int(num_elems)
        if not np.allclose(elem % 1, 0):
            raise ValueError("elem must be an integer")
        elem = int(elem)
        if elem < 0 or elem >= num_elems:
            raise ValueError(f"elem must be in [0, {num_elems - 1}] for" +
                             " this matrix")
        return Funmatrix(doki.funmatrix_partialtrace(self.m, elem, False))

    def trace(self):
        """Return the trace of this matrix."""
        if self.shape[0] != self.shape[1]:
            raise ValueError("Trace undefined for non square matrices")
        return doki.funmatrix_trace(self.m, False)

    def get_matrix(self):
        """Return self."""
        return self

    @property
    def shape(self):
        """Return tuple with (number of rows, number of columns)."""
        return self._shape
