#!/usr/bin/python
"""Deutsch-Jozsa with n qubit."""

from qsimov import Drewom, QCircuit, QGate
import sys
import time as t


def DJAlgCircuit(size, U_f):
    """Return Deutsch-Josza algorithm circuit.

    U_f is the oracle, having x1..xn and y as input qubits.
    x1..xn and y take only values either 0 or 1.
    Once applied it returns x1..xn and y', where y' = f(x1..xn) XOR y.
    Size is n + 1, where n is the number of input bits of f.
    """
    # The last qubit is defined as an ancilla qubit with value 1
    c = QCircuit(size, size-1, "Deutsch-Josza Algorithm", ancilla=[1])

    # We apply a Hadamard gate to all the qubits
    for i in range(size):
        c.add_operation("H", targets=i)

    # We apply the U_f oracle
    c.add_operation(U_f, targets=[i for i in range(size)])

    # We apply a Hadamard gate to all qubits but the last one
    for i in range(size - 1):
        c.add_operation("H", targets=i)

    # We measure x1..xn qubits.
    # If all of them are 0, f is constant, otherwise f is balanced.
    targets = [i for i in range(size - 1)]
    c.add_operation("MEASURE", targets=targets, outputs=targets)

    return c


def geq_zero(size):
    """DJ Oracle for f(x) = x >= 0. Balanced."""
    gate = QGate(size, 0, "U_(x>=0)")

    gate.add_operation("X", targets=size-1, anticontrols=size-2)

    return gate


def sample_main(nq, use_system=False, iterations=1):
    """Execute DJ with specified qubits."""
    # The number of qubits (x1..xn, y)
    gate = geq_zero(nq)  # The U_f oracle
    circuit = DJAlgCircuit(nq, gate)  # The Deutsch-Jozsa algorithm circuit
    # We specify useSystem = False to disable optimizations
    # other than Functional Matrices
    import numpy as np
    executor = Drewom(qmachine="doki",
                      extra={"num_threads": -1,
                             "random_generator": np.random.rand,
                             "use_system": use_system})
    init = t.time()
    mes = executor.execute(circuit, iterations=iterations)
    end = t.time()
    # print("RawMes:", mes)
    is_bal_iters = [any(iter[0][:-1]) for iter in mes]
    is_bal = all(is_bal_iters)
    if not is_bal and any(is_bal_iters):
        print("There is no consensus! BUG!")
    print("Is balanced?:", is_bal)
    print("Elapsed time:", end - init, "s")


if __name__ == "__main__":
    argv = sys.argv[1:]
    if 1 != len(argv) or int(argv[0]) < 2:  # We only have one argument
        print("Syntax: " + sys.argv[0] + " <number of qubits (min 2)>")
    else:
        sample_main(int(argv[0]), False, 1)
