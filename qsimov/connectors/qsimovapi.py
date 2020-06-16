from qsimov.structures.qgate import QGate, getGateSize
from qsimov.structures.qregistry import QRegistry, superposition
from qsimov.structures.qsystem import QSystem, joinSystems
from qsimov.structures.measure import Measure
from collections.abc import Iterable
import gc
import numpy as np


def _executeOnce(qsystem, lines, ancilla=None, useSystem=True):  # You can pass a QRegistry or an array to build a new QRegistry. When the second option is used, the ancilliary qubits will be added to the specified list.
    # print(qsystem)
    g = []
    r = qsystem
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
            def joinFunct(r1, r2):
                return superposition(r2, r1)
        a = QDataStruct(len(ancilla))
        for i in range(len(ancilla)):
            if ancilla[i] != 0:
                a.applyGate("X", qubit=i)
        raux = joinFunct(r, a)
        del r
        del a
        r = raux
    try:
        mres = []
        for line in lines:
            g = line[0]
            if type(g) != Measure:
                currbit = 0
                for i in range(len(line)):
                    if line[i] is not None and line[i][0] is not None and (isinstance(line[i][0], QGate) or line[i][0].lower() != "i"):
                        controls = line[i][1]
                        anticontrols = line[i][2]
                        if type(line[i][0]) == str:
                            r.applyGate(line[i][0], qubit=currbit, control=controls, anticontrol=anticontrols)
                        else:
                            line[i][0]._applyGate(r, currbit, controls, anticontrols)
                        currbit += getGateSize(line[i][0])
                    else:
                        currbit += 1
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
