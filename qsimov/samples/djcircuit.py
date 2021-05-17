#!/usr/bin/python
"""Deutsch-Jozsa with n qubit."""

from qsimov import QCircuit, QGate, Measure
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
    c = QCircuit("Deutsch-Josza Algorithm", ancilla=[1])

    # We apply a Hadamard gate to all the qubits
    c.add_line(*["H" for i in range(size)])
    c.add_line(U_f)  # Se aplica el oraculo

    # We apply a Hadamard gate to all qubits but the last one
    c.add_line(*("H" for i in range(size-1)), "I")

    # We measure x1..xn qubits.
    # If all of them are 0, f is constant, otherwise f is balanced.
    c.add_line(Measure([1 for i in range(size - 1)] + [0]))

    return c


def geq_zero(size):
    """DJ Oracle for f(x) = x >= 0. Balanced."""
    gate = QGate("U_(x>=0)")

    gate.add_line(*[None for i in range(size-1)], ["X", None, [size-2]])

    return gate


def sample_main():
    """Execute DJ with specified qubits."""
    argv = sys.argv[1:]
    if 1 != len(argv) or int(argv[0]) < 2:  # We only have one argument
        print("Syntax: " + sys.argv[0] + " <number of qubits (min 2)>")
        return
    # The number of qubits (x1..xn, y)
    nq = int(argv[0])
    gate = geq_zero(nq)  # The U_f oracle
    circuit = DJAlgCircuit(nq, gate)  # The Deutsch-Jozsa algorithm circuit
    # We specify useSystem = False to disable optimizations
    # other than Functional Matrices
    init = t.time()
    _, mes = circuit.execute([0 for i in range(nq - 1)],
                             args={"useSystem": False})
    end = t.time()
    mes = mes[0]
    print("Is balanced?:", any(mes[:-1]))
    print("Elapsed time:", end - init, "s")


if __name__ == "__main__":
    sample_main()
