from structures.qgate import QGate, I
from structures.qregistry import QRegistry, superposition
from structures.measure import Measure
from structures.funmatrix import Funmatrix
import ctypes as ct
import connectors.parser as prs
import structures.funmatrix as fm
import gc

__qsimov__ = ct.CDLL("libqsimov.dll")
c_double_p = ct.POINTER(ct.c_double)

__cH__ = __qsimov__.H
__cH__.argtypes = [ct.c_int]
__cH__.restype = ct.c_void_p
def H(n=1): # Devuelve una puerta Hadamard para n QuBits
    H = QGate("H(" + str(n) + ")")
    H.addLine(Funmatrix(ct.c_void_p(__cH__(ct.c_int(n))), H.name))
    return H

__cX__ = __qsimov__.X
__cX__.restype = ct.c_void_p
def PauliX(): # Also known as NOT
    px = QGate("X")
    px.addLine(Funmatrix(ct.c_void_p(__cX__()), px.name))
    return px

__cY__ = __qsimov__.Y
__cY__.restype = ct.c_void_p
def PauliY():
    py = QGate("Y")
    py.addLine(Funmatrix(ct.c_void_p(__cY__()), py.name))
    return py

__cZ__ = __qsimov__.Z
__cZ__.restype = ct.c_void_p
def PauliZ():
    pz = QGate("Z")
    pz.addLine(Funmatrix(ct.c_void_p(__cZ__()), pz.name))
    return pz

__cRX__ = __qsimov__.RX
__cRX__.argtypes = [ct.c_double]
__cRX__.restype = ct.c_void_p
def Rx(theta):
    rx = QGate("Rx(" + str(theta) + ")")
    rx.addLine(Funmatrix(ct.c_void_p(__cRX__(ct.c_double(theta))), rx.name))
    return rx

__cRY__ = __qsimov__.RY
__cRY__.argtypes = [ct.c_double]
__cRY__.restype = ct.c_void_p
def Ry(theta):
    ry = QGate("Ry(" + str(theta) + ")")
    ry.addLine(Funmatrix(ct.c_void_p(__cRY__(ct.c_double(theta))), ry.name))
    return ry

__cRZ__ = __qsimov__.RZ
__cRZ__.argtypes = [ct.c_double]
__cRZ__.restype = ct.c_void_p
def Rz(theta):
    rz = QGate("Rz(" + str(theta) + ")")
    rz.addLine(Funmatrix(ct.c_void_p(__cRZ__(ct.c_double(theta))), rz.name))
    return rz

__cSqrtX__ = __qsimov__.SqrtX
__cSqrtX__.restype = ct.c_void_p
def SqrtNOT(): # Square root of NOT gate, usually seen in its controlled form C-√NOT. Sometimes called √X or even V gate.
    v = QGate("SqrtX")
    v.addLine(Funmatrix(ct.c_void_p(__cSqrtX__()), "√NOT"))
    return v

__cCU__ = __qsimov__.CU
__cCU__.argtypes = [ct.c_void_p]
__cCU__.restype = ct.c_void_p
def CU(gate): # Returns a controlled version of the given gate [I 0; 0 U]
    cu = QGate("C-" + gate.name)
    cu.addLine(Funmatrix(ct.c_void_p(__cCU__(gate.getMatrix().m)), cu.name))
    return cu

__cCNOT__ = __qsimov__.CNOT
__cCNOT__.argtypes = [ct.c_int, ct.c_int]
__cCNOT__.restype = ct.c_void_p
# control is the id of the control qubit, objective is the id of the objective qubit
def CNOT(control=0, objective=1): # Returns a CNOT gate for two QuBits, also called Feynman gate
    cn = QGate("C-NOT(" + str(control) + ", " + str(objective) + ")")
    cn.addLine(Funmatrix(ct.c_void_p(__cCNOT__(ct.c_int(abs(objective - control)), ct.c_int(control > objective))), cn.name))
    return cn

__cSWAP__ = __qsimov__.SWAP
__cSWAP__.restype = ct.c_void_p
def SWAP(): # SWAP gate for 2 qubits
    sw = QGate("SWAP")
    sw.addLine(Funmatrix(ct.c_void_p(__cSWAP__()), sw.name))
    return sw

__cISWAP__ = __qsimov__.ISWAP
__cISWAP__.restype = ct.c_void_p
def ISWAP(): # ISWAP gate for 2 qubits
    isw = QGate("ISWAP")
    isw.addLine(Funmatrix(ct.c_void_p(__cISWAP__()), isw.name))
    return isw

def SqrtSWAP(): # Square root of SWAP gate for 2 qubits
    m = np.zeros((4,4), dtype=complex)
    m[0,0] = 1
    m[1,1] = 0.5 * (1+1j)
    m[1,2] = 0.5 * (1-1j)
    m[2,1] = 0.5 * (1-1j)
    m[2,2] = 0.5 * (1+1j)
    m[3,3] = 1
    return customGate("SqrtSWAP", m)

def Toffoli(): # Returns a CCNOT gate for three QuBits. A, B, C -> P = A, Q = B, R = AB XOR C.
    ''' # This does the same as the line below. Circuit with the implementation of Toffoli gate using SWAP, CNOT, Controlled-SNot and Controlled-SNot+
    # Gates needed (without control SWAPs): 5
    gate = np.kron(I(1), ControlledU(V()))
    gate = np.dot(gate, np.kron(SWAP(), I(1)))
    gate = np.dot(gate, np.kron(I(1), ControlledU(V())))
    gate = np.dot(gate, np.kron(CNOT(), I(1)))
    gate = np.dot(gate, np.kron(I(1), ControlledU(Dagger(V()))))
    gate = np.dot(gate, np.kron(CNOT(), I(1)))
    gate = np.dot(gate, np.kron(SWAP(), I(1)))
    return gate
    '''
    return CU(CNOT())

def Fredkin(): # Returns a CSWAP gate for three QuBits
    return CU(SWAP())

def Deutsch(angle): # Returns Deutsch gate with specified angle. D(pi/2) = Toffoli
    d = np.eye(8, dtype=complex)
    can = np.cos(angle)
    san = np.sin(angle)
    d[6,6] = can * 1j
    d[6,7] = san
    d[7,6] = san
    d[7,7] = can * 1j
    return customGate("D-" + str(angle), d)

__cU__ = __qsimov__.U
__cU__.argtypes = [ct.c_double, ct.c_double, ct.c_double]
__cU__.restype = ct.c_void_p
def u3(theta, phi, lamb): # U gate
    g = QGate("U(" + str(theta) + ", " + str(phi) + ", " + str(lamb) + ")")
    g.addLine(Funmatrix(ct.c_void_p(__cU__(ct.c_double(theta), ct.c_double(phi), ct.c_double(lamb))), g.name))
    return g

__cU2__ = __qsimov__.U2
__cU2__.argtypes = [ct.c_double, ct.c_double]
__cU2__.restype = ct.c_void_p
def u2(phi, lamb): # Equivalent to U(pi/2, phi, lambda)
    g = QGate("U(pi/2, " + str(phi) + ", " + str(lamb) + ")")
    g.addLine(Funmatrix(ct.c_void_p(__cU2__(ct.c_double(phi), ct.c_double(lamb))), g.name))
    return g

__cU1__ = __qsimov__.U1
__cU1__.argtypes = [ct.c_double]
__cU1__.restype = ct.c_void_p
def u1(angle): # Phase shift (R) gate, rotates qubit with specified angle (in radians). Equivalent to U(0, 0, lambda)
    g = QGate("U(0, 0, " + str(angle) + ")")
    g.addLine(Funmatrix(ct.c_void_p(__cU1__(ct.c_double(angle))), g.name))
    return g

__cPyCustomGate__ = __qsimov__.PyCustomGate
__cPyCustomGate__.argtypes = [c_double_p, c_double_p, ct.c_uint, ct.c_uint]
__cPyCustomGate__.restype = ct.c_void_p
def customGate(name, mat):
    g = QGate(name)
    size = mat.size
    nrows = mat.shape[0]
    mat = mat.reshape(size)
    array_type = ct.c_double * size
    remat = array_type(*[e.real for e in mat])
    immat = array_type(*[e.imag for e in mat])
    g.addLine(Funmatrix(ct.c_void_p(__cPyCustomGate__(remat, immat, ct.c_uint(nrows), ct.c_uint(size))), g.name))
    return g

__cIAA__ = __qsimov__.IAA
__cIAA__.argtypes = [ct.c_int]
__cIAA__.restype = ct.c_void_p
def IAA(n): # Inversion about the average
    g = QGate("Inversion about the average")
    g.addLine(Funmatrix(ct.c_void_p(__cIAA__(ct.c_int(n))), g.name))
    return g

def Peres(): # A, B, C -> P = A, Q = A XOR B, R = AB XOR C. Peres gate.
    ''' # Implementation of Peres gate with smaller gates.
    # Gates needed (without control SWAPs): 4
    p = QGate("Peres")
    p.addLine(SWAP(), I(1))
    p.addLine(I(1), ControlledU(SqrtNOT()))
    p.addLine(SWAP(), I(1))
    p.addLine(I(1), ControlledU(SqrtNOT()))
    p.addLine(CNOT(), I(1))
    p.addLine(I(1), ControlledU(Dagger(SqrtNOT())))
    return p
    '''
    p = QGate("Peres")
    p.addLine(Toffoli())
    p.addLine(CNOT(), I(1))
    return p

def R(): # A, B, C -> P = A XOR B, Q = A, R = AB XOR ¬C. R gate.
    # Optimized implementation with smaller gates
    # Gates needed (without control SWAPs): 6
    r = QGate("R")
    r.addLine(SWAP(), PauliX())
    r.addLine(I(1), ControlledU(SqrtNOT()))
    r.addLine(SWAP(), I(1))
    r.addLine(I(1), ControlledU(SqrtNOT()))
    r.addLine(CNOT(), I(1))
    r.addLine(I(1), ControlledU(Dagger(SqrtNOT())))
    r.addLine(SWAP(), I(1))
    return r

def TR(): # A, B, C -> P = A, Q = A XOR B, R = A¬B XOR C. TR gate.
    # Implementation of TR gate with smaller gates.
    # Gates needed (without control SWAPs): 6
    tr = QGate("TR")
    tr.addLine(I(1), PauliX(), I(1))
    tr.addLine(SWAP(), I(1))
    tr.addLine(I(1), ControlledU(SqrtNOT()))
    tr.addLine(SWAP(), I(1))
    tr.addLine(I(1), ControlledU(SqrtNOT()))
    tr.addLine(CNOT(), I(1))
    tr.addLine(I(1), ControlledU(Dagger(SqrtNOT())))
    tr.addLine(I(1), PauliX(), I(1))
    return tr

def URG(): # A, B, C -> P = (A+B) XOR C, Q = B, R = AB XOR C.
    # Implementation of URG gate with smaller gates.
    # Gates needed (without control SWAPs): 8
    urg = QGate("URG")
    urg.addLine(I(1), ControlledU(SqrtNOT()))
    urg.addLine(SWAP(), I(1))
    urg.addLine(I(1), ControlledU(SqrtNOT()))
    urg.addLine(CNOT(), I(1))
    urg.addLine(I(1), ControlledU(Dagger(SqrtNOT())))
    urg.addLine(CNOT(), I(1))
    urg.addLine(I(1), SWAP())
    urg.addLine(CNOT(), I(1))
    urg.addLine(I(1), CNOT())
    urg.addLine(CNOT(), I(1))
    urg.addLine(I(1), SWAP())
    urg.addLine(SWAP(), I(1))
    return urg

def BJN(): # A, B, C -> P = A, Q = B, R = (A+B) XOR C. BJN gate.
    # Implementation of TR gate with smaller gates.
    # Gates needed (without control SWAPs): 5
    bjn = QGate("BJN")
    bjn.addLine(SWAP(), I(1))
    bjn.addLine(I(1), ControlledU(SqrtNOT()))
    bjn.addLine(SWAP(), I(1))
    bjn.addLine(I(1), ControlledU(SqrtNOT()))
    bjn.addLine(CNOT(), I(1))
    bjn.addLine(I(1), ControlledU(SqrtNOT()))
    bjn.addLine(CNOT(), I(1))
    return bjn

def HalfSubstractor(): # A, B, 0 -> P = A-B, Q = Borrow, R = B = Garbage
    hs = QGate("Half Substractor")
    hs.addLine(SWAP(), I(1))
    hs.addLine(TR())
    hs.addLine(SWAP(), I(1))
    hs.addLine(I(1), SWAP())
    return hs

def Substractor(): # A, B, Bin, 0, 0, 0 -> P = A-B, Q = Borrow, R = B1 = Garbage, S = B1B2 = Garbage, T = Bin = Garbage, U = B = Garbage
    # Can be used as a comparator. Q will be 0 if A>=B, 1 otherwise.
    fs = QGate("Substractor")
    fs.addLine(I(2), SWAP(), I(2))
    fs.addLine(HalfSubstractor(), I(3))
    fs.addLine(I(2), SWAP(), I(2))
    fs.addLine(I(1), SWAP(), SWAP(), I(1))
    fs.addLine(I(2), SWAP(), SWAP())
    fs.addLine(HalfSubstractor(), I(3))
    fs.addLine(I(2), SWAP(), I(2))
    fs.addLine(I(3), SWAP(), I(1))
    fs.addLine(I(1), URG(), I(2))
    return fs

# Function that returns the 2^nth root of the unity
def nroot(n, rc = 14): # Rounds to 14 decimal places by default
    r = cm.exp(2j * cm.pi / pow(2, n))
    return round(r.real, rc) + round(r.imag, rc) * 1j

def RUnity(m, rc = 14):
    ru = np.eye(2, dtype=complex)
    ru[1,1] = nroot(m, rc)
    g = QGate("RU" + str(m))
    g.addLine(ru)
    return g

__cQFT__ = __qsimov__.QFT
__cQFT__.argtypes = [ct.c_int]
__cQFT__.restype = ct.c_void_p
def QFT(size):
    '''
    size = 4
    uft = np.kron(Hadamard(1), I(size - 1))
    uft = np.dot(uft, np.kron(np.dot(SWAP(), ControlledU(RUnity(2, rc))), I(size - 2)))
    uft = np.dot(uft, np.kron(Hadamard(1), np.kron(np.dot(SWAP(), ControlledU(RUnity(3, rc))), I(size - 3))))
    uft = np.dot(uft, np.kron(np.dot(SWAP(), ControlledU(RUnity(2, rc))), np.dot(SWAP(), ControlledU(RUnity(4, rc)))))
    uft = np.dot(uft, np.kron(Hadamard(1), np.kron(np.dot(SWAP(), ControlledU(RUnity(3, rc))), I(size - 3))))
    uft = np.dot(uft, np.kron(np.dot(SWAP(), ControlledU(RUnity(2, rc))), I(size - 2)))
    uft = np.dot(uft, np.kron(Hadamard(1), I(size - 1)))
    uft = np.dot(uft, np.kron(SWAP(), I(size - 2)))
    uft = np.dot(uft, np.kron(I(size - 3), np.kron(SWAP(), I(size - 3))))
    uft = np.dot(uft, np.kron(SWAP(), SWAP()))
    uft = np.dot(uft, np.kron(I(size - 3), np.kron(SWAP(), I(size - 3))))
    uft = np.dot(uft, np.kron(SWAP(), I(size - 2)))

    return uft
    '''
    qft = QGate("QFT(" + str(size) + ")")
    qft.addLine(Funmatrix(ct.c_void_p(__cQFT__(size)), qft.name))
    return qft

__gateDict__ = {}
__gateDict__["x"] = (PauliX, 0, 0)
__gateDict__["not"] = (PauliX, 0, 0)
__gateDict__["sqrtnot"] = (SqrtNOT, 0, 0)
__gateDict__["sqrtx"] = (SqrtNOT, 0, 0)
__gateDict__["v"] = (SqrtNOT, 0, 0)
__gateDict__["y"] = (PauliY, 0, 0)
__gateDict__["z"] = (PauliZ, 0, 0)
__gateDict__["rx"] = (Rx, 1, 1)
__gateDict__["ry"] = (Ry, 1, 1)
__gateDict__["rz"] = (Rz, 1, 1)
__gateDict__["h"] = (H, 0, 1)
__gateDict__["i"] = (I, 0, 1)
__gateDict__["swap"] = (SWAP, 0, 0)
__gateDict__["fredkin"] = (Fredkin, 0, 0)
__gateDict__["iswap"] = (ISWAP, 0, 0)
__gateDict__["sqrtswap"] = (SqrtSWAP, 0, 0)
__gateDict__["u"] = (u3, 1, 3)
__gateDict__["u3"] = (u3, 3, 3)
__gateDict__["u2"] = (u2, 2, 2)
__gateDict__["u1"] = (u1, 1, 1)
__gateDict__["iaa"] = (IAA, 1, 1)
__gateDict__["qft"] = (QFT, 1, 1)
__gateDict__["d"] = (Deutsch, 1, 1)
__gateDict__["toffoli"] = (Toffoli, 0, 0)

def getGate(gate):
    if type(gate) == str:
        ncontrols, gatename, invert, nargs, args = prs.getGroups(gate) # TODO: invert
        gatename = gatename.lower()
        if gatename in __gateDict__:
            gatemet, minargs, maxargs = __gateDict__[gatename]
            if gatename == "not":
                if ncontrols == 1:
                    gatemet, minargs, maxargs = CNOT, 0, 2
                    ncontrols = 0
                elif ncontrols >= 2:
                    gatemet, minargs, maxargs = Toffoli, 0, 0
                    ncontrols -= 2
            if gatename == "swap" and ncontrols >= 1:
                gatemet, minargs, maxargs = Fredkin, 0, 0
                ncontrols -= 1
            if gatename == "u":
                if nargs == 3:
                    minargs = 3
                elif nargs == 2:
                    gatemet = u2
                    minargs, maxargs = 2, 2
                elif nargs == 1:
                    gatemet = u1
                    minargs, maxargs = 1, 1

            if minargs <= nargs <= maxargs: # Adoro Python
                if nargs == 0:
                    gate = gatemet()
                else:
                    gate = gatemet(*args)
                for i in range(ncontrols):
                    gate = CU(gate)
                return gate
            else:
                raise ValueError(gatename + " gate number of args must be between " + str(minargs) + " and " + str(maxargs))
        else:
            raise ValueError(gatename + " can't be used with QSimovAPI")
    else:
        return gate

def _executeOnce(qregistry, lines, ancilla=None): # You can pass a QRegistry or an array to build a new QRegistry. When the second option is used, the ancilliary qubits will be added to the specified list.
    g = []
    r = qregistry
    firstGate = False
    if type(qregistry) == int:
        r = QRegistry(qregistry)
    elif type(qregistry) == list:
        r = QRegistry(len(qregistry))
        if any(qregistry):
            r.applyGate(*[I(1) if i == 0 else PauliX() for i in qregistry])
    elif ancilla is not None:
        raise ValueError("Can not use ancilla with precreated registry!")

    if ancilla is not None and len(ancilla) > 0:
        a = QRegistry(len(ancilla))
        if any(ancilla):
            a.applyGate(*[I(1) if i == 0 else PauliX() for i in ancilla])
        raux = superposition(r, a)
        del r
        del a
        r = raux
    try:
        mres = []
        for line in lines:
            g = line[0]
            if type(g) != Measure:
                g = getGate(g)
                if type(g) == QGate:
                    g = g.getMatrix()
                for rawgate in line[1:]:
                    gate = getGate(rawgate)
                    if type(gate) == QGate:
                        gate = gate.getMatrix()
                    g = g ** gate
                r.applyGate(g)
                del g
            else:
                r = g.check(r)
                mres += r[1]
                r = r[0]
            gc.collect()
    finally:
        gc.collect()
    return (r, mres)

def execute(qregistry, iterations=1, lines=[], ancilla=None):
    sol = [_executeOnce(qregistry, lines, ancilla) for i in range(iterations)]
    if iterations == 1:
        sol = sol[0]
    return sol
