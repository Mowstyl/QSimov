# -*- coding: utf-8 -*-
'''
QSimov: A Quantum Computing Toolkit.
Copyright (C) 2017  Hernán Indíbil de la Cruz Calvo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

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
from qsimov.structures.qsystem import QSystem
from qsimov.structures.qgate import QGate

# np.zeros((h,w), dtype=complex) Inicializa una matriz de numeros complejos
# con alto h y ancho w
# La suma de matrices se realiza con +. A + B
# La multiplicacion por un escalar se hace con *. n * A
# Para multiplicar las matrices A y B se usa np.dot(A,B)
# El producto Kronecker de A y B esta definido con np.kron(A,B)


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
