#!/usr/bin/python
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

"""Test module for QSimov."""

import sys
import qsimov as qj
# import webbrowser as wb
import random as rnd
import numpy as np
from operator import add

from qsimov.samples.state_gate import first_column_gate


def int_state_test(nq, negative=False, verbose=False, QItem=qj.QRegistry):
    """Test application of one gate of one qubit."""
    state_size = 2**nq

    passed = 0
    total = state_size
    if verbose:
        if not negative:
            print(" Simple positive state tests:")
        else:
            print(" Simple negative state tests:")
    val = 1
    if negative:
        val = -1
    v = np.zeros(state_size)
    for i in range(state_size):
        v[i-1] = 0
        v[i] = val
        gate = first_column_gate(nq, v)
        system = QItem(nq)
        system.apply_gate(gate)
        allOk = np.allclose(v, system.get_state())
        if not allOk:
            if verbose:
                print(v)
                print(system.get_state())
                print(v == system.get_state())
                print("  Michael Bay visited your simulator...")
            del v
            del system
            break
        passed += 1
        del system
    if allOk:
        if verbose:
            print("  Noice")
        del v
    return (passed, total)


def state_tests(minqubits, maxqubits, seed=None, verbose=False):
    """Test high level structures: QGate and QCircuit."""
    if not (seed is None):
        qj.set_seed(seed)
        rnd.seed(seed)
        np.random.seed(seed)
    result = [(0, 0) for i in range(4)]  # We have 2 tests

    if verbose:
        print("Testing state gate")
    for nq in range(minqubits, maxqubits + 1):
        if verbose:
            print("Testing with " + str(nq) + " qubit circuits")
        # Deutsch-Josza algorithm with QRegistry tests
        result[0] = map(add, result[0],
                        int_state_test(nq, verbose=verbose))
        # Deutsch-Josza algorithm with QRegistry tests
        result[1] = map(add, result[1],
                        int_state_test(nq, verbose=verbose,
                                       QItem=qj.QSystem))
        # Deutsch-Josza algorithm with QRegistry tests
        result[2] = map(add, result[2],
                        int_state_test(nq, negative=True, verbose=verbose,
                                       QItem=qj.QSystem))
        # Deutsch-Josza algorithm with QRegistry tests
        result[3] = map(add, result[3],
                        int_state_test(nq, negative=True, verbose=verbose,
                                       QItem=qj.QSystem))

    for i in range(0, 4):
        result[i] = tuple(result[i])
    return result


def main():
    """Execute all tests."""
    argv = sys.argv[1:]
    if 2 <= len(argv) <= 4 and int(argv[0]) >= 2:
        results = []
        if len(argv) == 2:
            seed = rnd.randrange(2**32 - 1)
            print("Seed: " + str(seed))
            print("\tTesting states...")
            results += state_tests(int(argv[0]), int(argv[1]), seed=seed)
        elif len(argv) == 3:
            print("Seed: " + str(int(argv[2])))
            print("\tTesting states...")
            results += state_tests(int(argv[0]), int(argv[1]),
                                   seed=int(argv[2]))
        else:
            print("Seed: " + str(int(argv[2])))
            print("\tTesting states...")
            results += state_tests(int(argv[0]), int(argv[1]),
                                   seed=int(argv[2]),
                                   verbose=bool(argv[3]))
        passed = [int(result[0] == result[1]) for result in results]
        noice = sum(passed)
        total = len(results)
        print("Passed: " + str(noice) + "/" + str(total))
        if noice != total:
            for testid in range(total):
                if passed[testid] == 0:
                    if testid < 2:
                        print("Unsigned int test", str(testid), "failed!")
                    elif testid < 2 + 2:
                        print("Signed int test", str(testid - 2), "failed!")
                        # We currently ignore signed tests
                        noice += 1
                    elif testid < 15 + 37 + 39:
                        print("QSystem Test", str(testid - (15 + 37)),
                              "failed!")
                    else:
                        print("High level Test", str(testid - (15 + 37 + 39)),
                              "failed!")
        if noice == total:
            print("PEACE AND TRANQUILITY")
            # wb.open_new_tab("https://youtu.be/SHvhps47Lmc")
        else:
            print("SORROW")
            # wb.open_new_tab("https://youtu.be/4Js-XbNj6Tk?t=37")
        # We assert so the test fails if we failed something
        assert noice == total
    else:
        print("Syntax: " + sys.argv[0] + " <minimum number of qubits (min 2)>",
              "<maximum number of qubits> <seed (optional)>")


if __name__ == "__main__":
    main()
