"""Module that executes a circuit using Doki (QSimov core).

This module provides functions that execute a given set of quantum gates
on a given quantum system, using Doki (the QSimov core written in C) to
perform the simulation.
"""

from qsimov.structures.qgate import QGate, get_gate_qubits
from qsimov.structures.qregistry import QRegistry, superposition
from qsimov.structures.qsystem import QSystem, join_systems
from qsimov.structures.measure import Measure
from collections.abc import Iterable
import gc
import numpy as np


def _execute_once(qsystem, lines, ancilla=None, useSystem=True, optimize=True):
    """Execute the gates in lines on a qsystem once.

    You can pass a QRegistry or an array to build a new QRegistry.
    When the second option is used,
    the ancilliary qubits will be added to the specified list.
    """
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
                r.apply_gate("X", qubit=i)
    elif ancilla is not None:
        raise ValueError("Can not use ancilla with precreated registry!")

    if ancilla is not None and len(ancilla) > 0:
        join_funct = join_systems
        if not useSystem:
            def join_funct(r1, r2):
                return superposition(r2, r1)
        a = QDataStruct(len(ancilla))
        for i in range(len(ancilla)):
            if ancilla[i] != 0:
                a.apply_gate("X", qubit=i)
        raux = join_funct(r, a)
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
                    if line[i] is not None and line[i][0] is not None and \
                            (isinstance(line[i][0], QGate) or
                             line[i][0].lower() != "i"):
                        controls = line[i][1]
                        anticontrols = line[i][2]
                        if type(line[i][0]) == str:
                            r.apply_gate(line[i][0],
                                         qubit=currbit,
                                         control=controls,
                                         anticontrol=anticontrols)
                        else:
                            line[i][0]._apply_gate(r, currbit,
                                                   controls,
                                                   anticontrols,
                                                   optimize=optimize)
                        currbit += get_gate_qubits(line[i][0])
                    else:
                        currbit += 1
            else:
                r = g.check(r, optimize=optimize)
                mres += r[1]
                r = r[0]
            gc.collect()
    finally:
        gc.collect()
    return (r, mres)


def execute(qregistry, iterations=1, lines=[], ancilla=None,
            useSystem=True, optimize=True):
    """Execute the gates in lines on a qsystem.

    Gets repeated the specified number of iterations.
    Returns the result of each iteration.
    """
    sol = [_execute_once(qregistry, lines, ancilla, useSystem, optimize)
           for i in range(iterations)]
    if iterations == 1:
        sol = sol[0]
    return sol
