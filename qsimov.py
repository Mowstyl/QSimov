# -*- coding: utf-8 -*-

import cmath as cm
import numpy as np
import ctypes as ct
import time as t
import re
import structures.funmatrix as fm
from structures.qregistry import *
from structures.qgate import *
from structures.qcircuit import *
from connectors.qsimovapi import *


# np.zeros((h,w), dtype=complex) Inicializa una matriz de numeros complejos con alto h y ancho w
# La suma de matrices se realiza con +. A + B
# La multiplicacion por un escalar se hace con *. n * A
# Para multiplicar las matrices A y B se usa np.dot(A,B)
# El producto Kronecker de A y B esta definido con np.kron(A,B)

# Lib C functions
__libc__ = ct.cdll.msvcrt
__cSrand__ = __libc__.srand
__cSrand__.argtypes = [ct.c_uint]
def setRandomSeed(seed, debug=False):
    if (seed == None):
        seed = ct.c_uint(int(t.time()))
    else:
        seed = ct.c_uint(int(seed))
    if debug:
        print ("Seed: " + str(seed.value))
    __cSrand__(seed)

__rep__ = re.compile("^((?:C-)+)?([a-zA-Z0-9]+)(\((?:(?:(?:[a-zA-Z]+)|(?:[0-9]+(?:\.[0-9]+)?))\,\s*)*(?:(?:(?:[a-zA-Z]+)|(?:[0-9]+(?:\.[0-9]+)?)))\))?$")
def getGroups(str_gate):
    res = __rep__.match(str_gate)
    return parseGroups(res.groups()) if res is not None else None

def parseGroups(groups):
    errored = False
    g1 = groups[0].count("C") if groups[0] is not None else 0
    g2 = groups[1]
    if groups[2] is not None:
        aux = groups[2][1:-1].split(",")
        g3 = len(aux)
        g4 = []
        for attr in aux:
            attr = attr.strip()
            if "." in attr:
                attr = float(attr)
            elif attr[0] in "0123456789":
                attr = int(attr)
            elif attr.lower() == "pi":
                attr = np.pi
            elif attr.lower() == "tau":
                attr = 2 * np.pi
            elif attr.lower() == "e":
                attr = np.e
            else:
                print (attr)
                errored = True
                break
            g4.append(attr)
    else:
        g3 = 0
        g4 = None

    if not errored:
        return (g1, g2, g3, g4)
    else:
        return None

def getTruthTable(gate, ancilla=None, garbage=0, iterations=1): # Prints the truth table of the given gate.
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
            if ancilla == None or qinit[-len(ancilla):] == ancilla:
                qr = QRegistry(len(qinit))
                qr.applyGate(*[I(1) if j == 0 else PauliX() for j in qinit])
                qr.applyGate(gate)
                mes = qr.measure([1 for j in range(num-garbage)])
                if ancilla != None:
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
    return np.array_equal(q1,q2) and str(q1) == str(q2)

setRandomSeed(None)
