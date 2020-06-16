# -*- coding: utf-8 -*-

import numpy as np
import ctypes as ct
import ctypes.util
import time as t
import platform as plat
from qsimov.structures.funmatrix import Funmatrix
from qsimov.structures.qregistry import QRegistry
from qsimov.structures.qsystem import QSystem
from qsimov.structures.qgate import QGate, getGate
from qsimov.structures.qcircuit import QCircuit
from qsimov.structures.measure import Measure
from qsimov.structures.condition import Condition

# np.zeros((h,w), dtype=complex) Inicializa una matriz de numeros complejos con alto h y ancho w
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


def setRandomSeed(seed, debug=False):
    if seed is None:
        seed = ct.c_uint(int(t.time()))
    else:
        seed = ct.c_uint(int(seed))
    if debug:
        print("Seed: " + str(seed.value))
    __cSrand__(seed)


def getTruthTable(gate, ancilla=None, garbage=0, iterations=1):  # Prints the truth table of the given gate.
    # You can set the ancilla bits to not include them in the table, with the list of values they must have.
    # For example, if you have two 0s and one 1 as ancilla bits, ancilla[0,0,1]. It always takes the last bits as the ancilla ones!
    # The garbage=n removes the last n bits from the truth table, considering them garbage.
    # For example, if you have 6 outputs and the last 4 outputs are garbage, only the value of the first two would be printed.
    # Always removes the last n bits!
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
                qr.applyGate(*[None if j == 0 else "X" for j in qinit])
                qr.applyGate(gate)
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
        print(k + " -> " + str(["P(1)=" + str(v) if type(v) != int and type(v) != int else v for v in mesd[k]]))
    return mesd


def QEq(q1, q2):
    return np.array_equal(q1, q2) and str(q1) == str(q2)


setRandomSeed(None)
