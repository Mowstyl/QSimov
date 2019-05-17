from structures.qgate import QGate, _getMatrix, I
from structures.qregistry import *
import structures.funmatrix as fm
import gc
import ctypes as ct

__qsimov__ = ct.CDLL("libqsimov.dll")
__cX__ = __qsimov__.X
__cX__.restype = ct.c_void_p
__pX__ = ct.c_void_p(__cX__())

class Measure(object):
    def __init__(self, mask, conds=[], remove=False):
        # Mask is a list of 0 and 1, where 0 means not to measure a certain QuBit and 1 means that the QuBit has to be measured
        self.mask = mask
        self.conds = conds
        self.remove = remove

    def __repr__(self):
        return ["Measure" if i == 1 else "I" for i in self.mask]

    def __str__(self):
        return self.__repr__()

    def _mesToList(self, mresults):
        lin = 0
        mres = []
        for m in self.mask:
            tap = None
            if m == 1:
                tap = mresults[lin]
                lin += 1
            mres.append(tap)
        return mres

    def check(self, qregistry):
        res = qregistry.measure(self.mask, remove=self.remove)
        res = self._mesToList(res)
        # print ("Measure result: " + str(res))
        r = (qregistry, [res])
        for cond in self.conds:
            aux = cond.evaluate(r[0], res)
            if type(aux) == tuple:
                r = (aux[0], r[1] + aux[1])
            else:
                r = (aux, r[1])
        return r


class Condition(object):
    def __init__(self, cond, ifcase, elcase, typeif, typeel):
        # cond is an array of what we expect to have measured in each QuBit. None if we don't care about a certain value. Example: [0, 1, None, None, 1].
        # ifcase and elcase can be Conditions or QCircuits to be applied to the registry. They can also be functions that take the registry and the result as a parameter.
        # typeif and typeel are integers from -1 to 2.
        # -1 means that nothing has to be done in that case.
        # 0 means that ifcase/elcase is another Condition.
        # 1 means that ifcase/elcase is a QGate to be applied.
        # 2 means that ifcase/elcase is a QCircuit.
        self.cond = cond
        self.ifcase = ifcase
        self.elcase = elcase
        self.typeif = typeif
        self.typeel = typeel

    def evaluate(self, qregistry, mresults):
        case = self.elcase
        t = self.typeel
        if _specialCompare(self.cond, mresults):
            case = self.ifcase
            t = self.typeif
        if t == 0: # Condition
            r = case.evaluate(qregistry, mresults)
        elif t == 1: # QGate
            r = qregistry
            r.applyGate(case)
        elif t == 2: # QCircuit
            r = case._executeOnce(qregistry)
        else: # Do nothing
            r = qregistry
        return r

class QCircuit(object):
    def __init__(self, name="UNNAMED", ancilla=[]): # You can choose whether to save the circuit and apply gates separately on each computation (faster circuit creation) or to precompute the matrixes (faster execution)
        self.name = name
        self.matrix = [1]
        self.measure = []
        self.lines = []
        self.plan = [0]
        self.ancilla = ancilla

    def addLine(self, *args):
        try:
            self.lines.append(list(args))
        finally:
            gc.collect()

    def _executeOnce(self, qregistry): # You can pass a QRegistry or an array to build a new QRegistry. When the second option is used, the ancilliary qubits will be added to the specified list.
        g = []
        r = qregistry
        firstGate = False
        if type(qregistry) == int:
            r = QRegistry(qregistry)
        elif type(qregistry) == list:
            r = QRegistry(len(qregistry))
            if any(qregistry):
                r.applyGate(*[I(1) if i == 0 else __pX__ for i in qregistry])
        elif self.ancilla is not None:
            raise ValueError("Can not use ancilla with precreated registry!")

        if self.ancilla is not None and len(self.ancilla) > 0:
            a = QRegistry(len(self.ancilla))
            if any(self.ancilla):
                a.applyGate(*[I(1) if i == 0 else __pX__ for i in self.ancilla])
            raux = superposition(r, a)
            #freeRegistry(r)
            #freeRegistry(a)
            r = raux

        try:
            print(r.getState())
            mres = []
            for line in self.lines:
                g = line[0]
                if type(g) != Measure:
                    g = _getMatrix(g)
                    for gate in line[1:]:
                        g = fm.kron(g, _getMatrix(gate))
                    r.applyGate(g)
                    print(r.getState())
                    del g
                else:
                    r = g.check(r)
                    mres += r[1]
                    r = r[0]
                gc.collect()
        finally:
            print ("Done!")
            gc.collect()
        return (r, mres)

    def execute(self, qregistry, iterations=1, qmachine=None, args=None):
        if type(qregistry) == QRegistry and iterations > 1:
            raise ValueError("Can not do more than one iteration with a precreated registry!")
        elif (qmachine == None or qmachine == "local"):
            sol = [self._executeOnce(qregistry) for i in range(iterations)]
            if iterations == 1:
                sol = sol[0]
            return sol
        else:
            raise ValueError("Unsupported qmachine!")

def _specialCompare(a, b):
    same = len(a) == len(b)
    if (same):
        for i in range(len(a)):
            if a[i] is not None and a[i] != b[i]:
                same = False
                break
    return same
