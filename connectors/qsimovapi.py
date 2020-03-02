from structures.qgate import QGate
from structures.qregistry import QRegistry, superposition
from structures.qsystem import QSystem, joinSystems
from structures.measure import Measure
from structures.funmatrix import Funmatrix
from collections.abc import Iterable
import ctypes as ct
import connectors.parser as prs
import structures.funmatrix as fm
import gc
import numpy as np

def _executeOnce(qsystem, lines, ancilla=None, useSystem=True): # You can pass a QRegistry or an array to build a new QRegistry. When the second option is used, the ancilliary qubits will be added to the specified list.
    g = []
    r = qsystem
    firstGate = False
    QDataStruct = QSystem if useSystem else QRegistry
    if np.issubdtype(type(qsystem), np.integer):
        r = QDataStruct(qsystem)
    elif isinstance(qsystem, Iterable):
        r = QDataStruct(len(qsystem))
        for i in range(len(qsystem)):
            if qsystem[i] != 0:
                r.applyGate("X", qubit=i)
    elif ancilla is not None:
        raise ValueError("Can not use ancilla with precreated registry!")

    if ancilla is not None and len(ancilla) > 0:
        joinFunct = joinSystems
        if not useSystem:
            joinFunct = superposition
        a = QDataStruct(len(ancilla))
        for i in range(len(ancilla)):
            if ancilla[i] != 0:
                a.applyGate("X", qubit=i)
        raux = joinFunct(r, a)
        if not useSystem:
            del r
            del a
        r = raux
    try:
        mres = []
        for line in lines:
            g = line[0]
            if type(g) != Measure:
                for i in range(len(line)):
                    controls = line[i][1]
                    anticontrols = line[i][2]
                    if type(line[i][0]) == str:
                        r.applyGate(line[i][0], qubit=i, control=controls, anticontrol=anticontrols)
                    else:
                        line[i][0]._applyGate(r, i, controls, anticontrols)
            else:
                r = g.check(r)
                mres += r[1]
                r = r[0]
            gc.collect()
    finally:
        gc.collect()
    return (r, mres)

def execute(qregistry, iterations=1, lines=[], ancilla=None, useSystem=True):
    sol = [_executeOnce(qregistry, lines, ancilla, useSystem) for i in range(iterations)]
    if iterations == 1:
        sol = sol[0]
    return sol
