import ctypes as ct

c_double_p = ct.POINTER(ct.c_double)

__funmat__ = ct.CDLL("libfunmat.dll")
__cGetItem__ = __funmat__.getitemaux
__cGetItem__.argtypes = [ct.c_void_p, ct.c_uint, ct.c_uint, c_double_p, c_double_p]
__cGetItem__.restype = ct.c_int

def __getItem__(mat, i, j):
    res = None
    re = ct.c_double(0.0)
    im = ct.c_double(0.0)
    try:
        aux = __cGetItem__(mat, ct.c_uint(i), ct.c_uint(j), c_double_p(re), c_double_p(im))
        if (aux == 0):
            print ("Error getting the specified item!")
        else:
            res = complex(re.value, im.value)
    except Exception as e:
        print (e)
    return res

__cGetRows__ = __funmat__.rows
__cGetRows__.argtypes = [ct.c_void_p]
__cGetRows__.restype = ct.c_int

__cGetCols__ = __funmat__.columns
__cGetCols__.argtypes = [ct.c_void_p]
__cGetCols__.restype = ct.c_int

def getItem(mat, key):
    rows = int(__cGetRows__(mat))
    cols = int(__cGetCols__(mat))

    if (type(key) == int):
        if (key >= rows):
            raise IndexError("index " + str(key) + " is out of bounds for axis 0 with shape " + str(rows))
        if (key < 0):
            key = rows + key
        return [__getItem__(mat, key, j) for j in range(cols)]
    elif (type(key) == slice):
        start = key.start if key.start != None else 0
        if (start < 0):
            start = rows + start
        stop = key.stop if key.stop != None else rows
        if (stop < 0):
            stop = rows + stop
        step = key.step if key.step != None else 1
        return [getItem(mat, i) for i in range(start, stop, step)]
    elif (type(key) == tuple):
        if (type(key[1]) == None):
            return getItem(mat, key[0])
        elif (type(key[0]) == None):
            if (type(key[1]) == int):
                if (key[1] >= cols):
                    raise IndexError("index " + str(key[1]) + " is out of bounds for axis 1 with shape " + str(cols))
                if (key[1] < 0):
                    key = (key[0], cols + key[1])
                return [np.around(__getItem__(mat, i, key[1]), decimals=self.dec) for i in range(rows)]
            elif (type(key[1]) == slice):
                start = key[1].start if key[1].start != None else 0
                if (start < 0):
                    start = cols + start
                stop = key[1].stop if key[1].stop != None else cols
                if (stop < 0):
                    stop = cols + stop
                step = key[1].step if key[1].step != None else 1
                return [getItem(mat, None, i) for i in range(start, stop, step)]
        elif (type(key[0]) == int and type(key[1]) == int):
            if (key[0] >= rows):
                raise IndexError("index " + str(key[0]) + " is out of bounds for axis 0 with shape " + str(rows))
            if (key[1] >= cols):
                raise IndexError("index " + str(key[1]) + " is out of bounds for axis 1 with shape " + str(cols))
            if (key[0] < 0):
                key = (rows + key[0], key[1])
            if (key[1] < 0):
                key = (key[0], cols + key[1])
            return __getItem__(mat, key[0], key[1])
        elif (type(key[0]) == slice and type(key[1]) == int):
            if (key[1] >= cols):
                raise IndexError("index " + str(key[1]) + " is out of bounds for axis 1 with shape " + str(cols))
            start = key[0].start if key[0].start != None else 0
            if (start < 0):
                start = rows + start
            stop = key[0].stop if key[0].stop != None else rows
            if (stop < 0):
                stop = rows + stop
            step = key[0].step if key[0].step != None else 1
            if (key[1] < 0):
                key = (key[0], cols + key[1])
            return [__getItem__(mat, i, key[1]) for i in range(start, stop, step)]
        elif (type(key[0]) == int and type(key[1]) == slice):
            if (key[0] >= rows):
                raise IndexError("index " + str(key[0]) + " is out of bounds for axis 0 with shape " + str(rows))
            start = key[1].start if key[1].start != None else 0
            if (start < 0):
                start = cols + start
            stop = key[1].stop if key[1].stop != None else cols
            if (stop < 0):
                stop = cols + stop
            step = key[1].step if key[1].step != None else 1
            if (key[0] < 0):
                key = (rows + key[0], key[1])
            return [__getItem__(mat, key[0], i) for i in range(start, stop, step)]
        elif (type(key[0]) == slice and type(key[1]) == slice):
            start0 = key[0].start if key[0].start != None else 0
            if (start0 < 0):
                start0 = rows + start0
            stop0 = key[0].stop if key[0].stop != None else rows
            if (stop0 < 0):
                stop0 = rows + stop0
            step0 = key[0].step if key[0].step != None else 1
            start1 = key[1].start if key[1].start != None else 0
            if (start1 < 0):
                start1 = cols + start1
            stop1 = key[1].stop if key[1].stop != None else cols
            if (stop1 < 0):
                stop1 = cols + stop1
            step1 = key[1].step if key[1].step != None else 1
            return [[__getItem__(mat, i, j) for i in range(start0, stop0, step0)] for j in range(start1, stop1, step1)]

__cMadd__ = __funmat__.madd
__cMadd__.argtypes = [ct.c_void_p, ct.c_void_p]
__cMadd__.restype = ct.c_void_p

def madd(a, b):
    res = None
    try:
        res = __cMadd__(a, b)
    except Exception as e:
        print (e)
    return res

__cMsub__ = __funmat__.msub
__cMsub__.argtypes = [ct.c_void_p, ct.c_void_p]
__cMsub__.restype = ct.c_void_p

def msub(a, b):
    res = None
    try:
        res = __cMsub__(a, b)
    except Exception as e:
        print (e)
    return res

__cMprod__ = __funmat__.mprodaux
__cMprod__.argtypes = [ct.c_double, ct.c_double, ct.c_void_p]
__cMprod__.restype = ct.c_void_p

def mprod(r, a):
    res = None
    try:
        res = __cMprod__(ct.c_double(r.real), ct.c_double(r.imag), a)
    except Exception as e:
        print (e)
    return res

__cMdiv__ = __funmat__.mdivaux
__cMdiv__.argtypes = [ct.c_double, ct.c_double, ct.c_void_p]
__cMdiv__.restype = ct.c_void_p

def mdiv(r, a):
    res = None
    try:
        res = __cMdiv__(ct.c_double(r.real), ct.c_double(r.imag), a)
    except Exception as e:
        print (e)
    return res

__cMatmul__ = __funmat__.matmul
__cMatmul__.argtypes = [ct.c_void_p, ct.c_void_p]
__cMatmul__.restype = ct.c_void_p

def matmul(a, b):
    res = None
    try:
        res = __cMatmul__(a, b)
    except Exception as e:
        print (e)
    return res

__cEwmul__ = __funmat__.ewmul
__cEwmul__.argtypes = [ct.c_void_p, ct.c_void_p]
__cEwmul__.restype = ct.c_void_p

def ewmul(a, b):
    res = None
    try:
        res = __cEwmul__(a, b)
    except Exception as e:
        print (e)
    return res

__cKron__ = __funmat__.kron
__cKron__.argtypes = [ct.c_void_p, ct.c_void_p]
__cKron__.restype = ct.c_void_p

def kron(a, b):
    res = None
    try:
        res = __cKron__(a, b)
    except Exception as e:
        print (e)
    return res

__cTranspose__ = __funmat__.transpose
__cTranspose__.argtypes = [ct.c_void_p]

def transpose(m):
    try:
        __cTranspose__(m)
    except Exception as e:
        print (e)

__cDagger__ = __funmat__.dagger
__cDagger__.argtypes = [ct.c_void_p]

def dagger(m):
    try:
        __cDagger__(m)
    except Exception as e:
        print (e)
