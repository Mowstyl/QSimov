from structures.qgate import QGate
from structures.qregistry import QRegistry, superposition
from structures.measure import Measure
from structures.funmatrix import Funmatrix
import ctypes as ct
import connectors.parser as prs
import structures.funmatrix as fm
import gc

def _executeOnce(qsystem, lines, ancilla=None): # You can pass a QRegistry or an array to build a new QRegistry. When the second option is used, the ancilliary qubits will be added to the specified list.
    g = []
    r = qregistry
    firstGate = False
    if type(qregistry) == int:
        r = QSystem(qsystem)
    elif type(qregistry) == list:
        r = QSystem(len(qsystem))
        if any(qsystem):
            r.applyGate(*[None if i == 0 else "X" for i in qregistry])
    elif ancilla is not None:
        raise ValueError("Can not use ancilla with precreated registry!")

    if ancilla is not None and len(ancilla) > 0:
        a = QSystem(len(ancilla))
        if any(ancilla):
            a.applyGate(*[None if i == 0 else "X" for i in ancilla])
        raux = superposition(r, a)
        del r
        del a
        r = raux
    try:
        mres = []
        for line in lines:
            g = line[0]
            if type(g) != Measure:
                #g = getGate(g)   # TODO
                if type(g) == QGate:
                    g = g.getMatrix()
                for rawgate in line[1:]:
                    #gate = getGate(rawgate)
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
