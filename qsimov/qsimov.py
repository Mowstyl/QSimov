# -*- coding: utf-8 -*-
"""Base module for the Quantum Computing Simulator named QSimov.

Data Structures:
    Funmatrix: Functional Matrices related stuff in python
    QRegistry: Quantum Registry, base of all quantum related operations
    QSystem: Quantum System, preferred over QRegistry (can save a lot of space)
    QGate: Quantum Gate, built from elemental gates or other QGates
    QCircuit: Quantum Circuit, built from gates and measurements
    Measure: Data structure that represents a measurement in a circuit
    Condition: A condition to be evaluated after a measurement

Functions:
    setRandomSeed: sets the seed for measurements.
        Automatically called when this module is imported
    getGate: returns data about an elemental gate (not QGate)
    QEq: UNSTABLE. Compares two items
    getTruthTable: UNSTABLE. Returns the truth table of a gate

Data structures and functions marked as UNSTABLE might not work
Avoid using them as they might be removed in future versions
"""

import numpy as np
import ctypes as ct
import ctypes.util
import time as t
import platform as plat
from qsimov.structures.funmatrix import Funmatrix
from qsimov.structures.qregistry import QRegistry, superposition
from qsimov.structures.qsystem import QSystem, join_systems
from qsimov.structures.qgate import QGate, getGate, get_gate
from qsimov.structures.qcircuit import QCircuit
from qsimov.structures.measure import Measure
from qsimov.structures.condition import Condition

# np.zeros((h,w), dtype=complex) Inicializa una matriz de numeros complejos
# con alto h y ancho w
# La suma de matrices se realiza con +. A + B
# La multiplicacion por un escalar se hace con *. n * A
# Para multiplicar las matrices A y B se usa np.dot(A,B)
# El producto Kronecker de A y B esta definido con np.kron(A,B)

# Lib C functions
if plat.system() == "Windows":
    __libc__ = ct.cdll.msvcrt
else:
    __libc__ = ct.cdll.LoadLibrary(ctypes.util.find_library("c"))
__cSrand__ = __libc__.srand
__cSrand__.argtypes = [ct.c_uint]

_seed = None


def get_seed():
    """Return the current seed."""
    global _seed
    return _seed


def setRandomSeed(seed, debug=False):
    """Use set_seed method instead. DEPRECATED."""
    print("Method QSimov.setRandomSeed is deprecated.",
          "Please use set_seed if you seek the same functionality")
    return set_seed(seed, debug=debug)


def set_seed(seed, debug=False):
    """Set the seed used in measurements.

    seed: The seed to use
    debug: Whether to print the seed or not (default: False)
    """
    if seed is None:
        seed = ct.c_uint(int(t.time()))
    else:
        seed = ct.c_uint(int(seed))
    if debug:
        print("Seed: " + str(seed.value))
    __cSrand__(seed)
    global _seed
    _seed = seed.value


def getTruthTable(gate, ancilla=None, garbage=0, iterations=1):
    """Print the truth table of the given gate.

    You can set the ancilla bits to not include them in the table,
    with the list of values they must have.
    For example, if you have two 0s and one 1 as ancilla bits, ancilla[0,0,1].
    It always takes the last bits as the ancilla ones!
    The garbage=n removes the last n bits from the truth table,
    considering them garbage.
    For example, if you have 6 outputs and the last 4 outputs are garbage,
    only the value of the first two would be printed.
    Always removes the last n bits!
    """
    if (type(gate) == QGate):
        gate = gate.m
    num = int(np.log2(gate.getRows()))
    mesd = {}
    for iteration in range(iterations):
        for i in range(0, gate.getRows()):
            nbin = [int(x) for x in bin(i)[2:]]
            qinit = [0 for j in range(num - len(nbin))]
            qinit += nbin
            if ancilla is None or qinit[-len(ancilla):] == ancilla:
                qr = QSystem(len(qinit))
                qr.apply_gate(*[None if j == 0 else "X" for j in qinit])
                qr.apply_gate(gate)
                mes = qr.measure([1 for j in range(num-garbage)])
                if ancilla is not None:
                    ini = qinit[:-len(ancilla)]
                else:
                    ini = qinit
                if str(ini) not in mesd:
                    mesd[str(ini)] = np.zeros(num)
                mesd[str(ini)] = [x + y for x, y in zip(mesd[str(ini)], mes)]
    for k in mesd:
        for i in range(len(mesd[k])):
            mesd[k][i] /= iterations
            if (mesd[k][i] == 1.0 or mesd[k][i] == 0.0):
                mesd[k][i] = int(mesd[k][i])
        print(k + " -> " + str(["P(1)=" + str(v)
                                if type(v) != int and type(v) != int
                                else v
                                for v in mesd[k]]))
    return mesd


def QEq(q1, q2):
    """Return if q1 and q2 are equal and represented the same."""
    return np.array_equal(q1, q2) and str(q1) == str(q2)


set_seed(None, debug=False)
