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

import doki
import numpy as np
import qsimov as qj
import random as rnd
import sympy as sp
import sys
# import webbrowser as wb

# from operator import add
from qsimov.samples.djcircuit import DJAlgCircuit
from sympy.physics.quantum import TensorProduct


# th, ph, la = sp.symbols("θ φ λ", real=True)


def Bal(n, controlId=0):
    """Return Deutsch-Jozsa oracle for balanced function."""
    gate = qj.QGate(n, 0, "Balanced")
    gate.add_operation("X", targets=n-1, controls=controlId)
    return gate


def Const(n, twice=False):
    """Return Deutsch-Jozsa oracle for constant function."""
    gate = qj.QGate(n, 0, "Constant")
    gate.add_operation("X", targets=n-1)
    if twice:
        gate.add_operation("X", targets=n-1)
    return gate


def TeleportationCircuit(gate):
    """Return teleportation algorithm circuit.

    Positional arguments:
        gate: gate to apply to the qubit that is going to be sent
            (so we don't send a 0 or a 1)
    Return:
        QCircuit with teleportation algorithm
    """
    qc = qj.QCircuit(3, 2, "Teleportation", ancilla=(0, 0))
    qc.add_operation("H", targets=[1])
    qc.add_operation("X", targets=[2], controls=[1])
    # Aqui es donde trabajamos con el qubit Q que queremos enviar.
    # Se le aplica la puerta pasada como parámetro.
    qc.add_operation(gate, targets=0)
    # Una vez terminado todo lo que queremos hacerle al QuBit,
    # procedemos a preparar el envio
    # Se aplica una puerta C-NOT sobre Q (control) y B (objetivo).
    qc.add_operation("X", targets=1, controls=0)
    # Se aplica una puerta Hadamard sobre Q.
    qc.add_operation("H", targets=[0])

    qc.add_operation("MEASURE", targets=(0, 1), outputs=(0, 1))
    qg = qj.QGate(1, 1, "CX")
    qg.add_operation("X", targets=0, c_controls=0)
    qc.add_operation(qg, targets=2, c_targets=1)
    qc.add_operation("Z", targets=2, c_controls=0)

    return qc


def entangle_gate():
    """Return a QGate that creates a Bell pair."""
    e = qj.QGate(2, 0, "Entangle")

    e.add_operation("H", targets=0)
    e.add_operation("X", targets=1, controls=0)

    return e


def entangle_system(s, id, id2):
    """Entangle specified qubits of a system."""
    aux = s.apply_gate("H", targets=id)
    res = aux.apply_gate("X", targets=id2, controls=id)
    del aux
    return res


def inversion_tests(verbose=False):
    """Test gate inversion."""
    if verbose:
        print("  Testing gate inversion...")
    e = entangle_gate()  # TODO: Intensive inversion test

    ei = e.invert()
    if e.get_operations() != ei.get_operations()[::-1]:
        if verbose:
            print(e.get_operations())
            print(ei.get_operations())
            print(ei.get_operations()[::-1])
            print(e.get_operations() == ei.get_operations()[::-1])
            print([e.get_operations()[i] != ei.get_operations()[::-1][i]
                   for i in range(len(e.get_operations()))])
            print("    Michael Bay visited your simulator...")
        raise AssertionError("Failed invert test")

    ed = e.dagger()
    if e.get_operations() != ed.get_operations()[::-1]:
        if verbose:
            print(e.get_operations())
            print(ed.get_operations())
            print(ed.get_operations()[::-1])
            print(e.get_operations() != ed.get_operations()[::-1])
            print([e.get_operations()[i] != ed.get_operations()[::-1][i]
                   for i in range(len(e.get_operations()))])
            print("    Michael Bay visited your simulator...")
        raise AssertionError("Failed dagger test")

    r = qj.QRegistry(2)
    er = qj.QRegistry(2)
    s = qj.QSystem(2)
    es = qj.QSystem(2)

    aux = r.apply_gate(e)
    del r
    r = aux.apply_gate(ed)
    del aux
    aux = s.apply_gate(e)
    del s
    s = aux.apply_gate(ed)
    del aux
    if any(r.get_state() != er.get_state()) \
            or any(s.get_state() != es.get_state()):
        if verbose:
            print(r.get_state())
            print(er.get_state())
            print(s.get_state())
            print(es.get_state())
            print(r.get_state() == er.get_state())
            print(s.get_state() == es.get_state())
            print("    Michael Bay visited your simulator...")
        del r
        del er
        del s
        del es
        raise AssertionError("Error comparing gate+inversion result")
    if verbose:
        print("    Noice")


def entangle_test(QItem, id1, id2, verbose):
    e = entangle_gate()
    r = QItem(3)
    rg = QItem(3)
    aux = entangle_system(r, id1, id2)
    del r
    r = aux
    aux = rg.apply_gate(e, targets=[id1, id2])
    del rg
    rg = aux

    if any(r.get_state() != rg.get_state()):
        if verbose:
            print(id1)
            print(id2)
            print(e.get_operations())
            print(r.get_state())
            print(rg.get_state())
            print(r.get_state() == rg.get_state())
            print("    Michael Bay visited your simulator...")
        del r
        del rg
        raise AssertionError(f"Error entangling {id1} and {id2}")
    del r
    del rg


def entangle_tests(verbose=False, useSystem=False):
    """Test entanglement."""
    if verbose:
        print("  Testing entanglement gate...")

    if useSystem:
        QItem = qj.QSystem
    else:
        QItem = qj.QRegistry

    entangle_test(QItem, 0, 1, verbose)
    entangle_test(QItem, 0, 2, verbose)
    entangle_test(QItem, 1, 2, verbose)
    if verbose:
        print("    Noice")


def one_gate_tests(nq, verbose=False, QItem=qj.QRegistry):
    """Test application of one gate of one qubit."""
    if verbose:
        print(" One qubit gate tests:")
    size = 2**nq
    for id in range(nq):
        invert = bool(rnd.randint(0, 1))
        rands = np.random.rand(3) * 2 * np.pi - np.pi
        gate = "U(" + ",".join([str(angle)
                                for angle in rands]) + ")"
        if invert:
            gate += "-1"
        sgate = qj.SimpleGate(gate)
        if verbose:
            print("  Testing gate " + gate + " to qubit " + str(id) + "...")
        # print("    Gate: " + str(numpygate))
        b = QItem(nq, verbose=False)
        a = doki.registry_new(nq, False)
        a2 = doki.registry_apply(a, sgate.gate, [id], set(), set(),
                                 -1, False)
        del a
        if verbose:
            print("    Doki done")
        b2 = b.apply_gate(gate, targets=id)
        del b
        if verbose:
            print("    QSimov done")
        a2_state = np.array([doki.registry_get(a2, i, False, False)
                             for i in range(size)])
        b2_state = b2.get_state()
        if not np.allclose(a2_state, b2_state):
            if verbose:
                print("Expected:", a2_state)
                print("Received:", b2_state)
                print(a2_state == b2_state)
                print("Gate:", sgate.matrix)
                print("    Michael Bay visited your simulator...")
            del a2
            del b2
            raise AssertionError("Error comparing states after applying gate")
        if verbose:
            print("    Noice")
        del a2
        del b2


def TwoU_np(angle1_1, angle1_2, angle1_3,
            angle2_1, angle2_2, angle2_3):
    """Return numpy two qubit gate that may entangle."""
    gate1str = f"U({angle1_1},{angle1_2},{angle1_3})"
    gate2str = f"U({angle2_1},{angle2_2},{angle2_3})"
    gate1aux = qj.SimpleGate(gate1str)
    gate2aux = qj.SimpleGate(gate2str)
    gate1 = TensorProduct(gate1aux.matrix, sp.eye(2))
    gate2 = sp.eye(4)
    gate2[2, 2] = gate2aux.matrix[0, 0]
    gate2[2, 3] = gate2aux.matrix[0, 1]
    gate2[3, 2] = gate2aux.matrix[1, 0]
    gate2[3, 3] = gate2aux.matrix[1, 1]
    return gate2 @ gate1


def _add_two_U():
    """Add the TwoU gate to the list of available gates."""
    qj.add_gate("TwoU", TwoU_np, 6, 6,
                is_own_inverse=False, overwrite=True)


def two_gate_tests(nq, verbose=False, QItem=qj.QRegistry):
    """Test application of one gate of two qubits."""
    if verbose:
        print(" Two qubit gate tests:")
    if nq < 2:
        raise ValueError("Can't apply 2 qubit gates to 1 qubit structure")
    _add_two_U()
    size = 2**nq
    for id1 in range(nq):
        for id2 in range(nq):
            if id1 == id2:
                continue
            invert = bool(rnd.randint(0, 1))
            rands = np.random.rand(6) * 2 * np.pi - np.pi
            gate = "TwoU(" + ",".join([str(angle)
                                       for angle in rands]) + ")"
            if invert:
                gate += "-1"
            sgate = qj.SimpleGate(gate)
            if verbose:
                print(f"  Testing gate {gate} to qubits {id1} and {id2}...")
            b = QItem(nq, verbose=False)
            a = doki.registry_new(nq, False)
            a2 = doki.registry_apply(a, sgate.gate, [id1, id2], set(), set(),
                                     -1, verbose)
            del a
            if verbose:
                print("    Doki done")
            b2 = b.apply_gate(gate, targets=[id1, id2])
            del b
            if verbose:
                print("    QSimov done")
            a2_state = np.array([doki.registry_get(a2, i, False, False)
                                 for i in range(size)])
            b2_state = b2.get_state()
            if not np.allclose(a2_state, b2_state):
                if verbose:
                    print("Expected:", a2_state)
                    print("Received:", b2_state)
                    print(b2.qubitMap)
                    print(b2.regs[b2.qubitMap[0]][0].get_state())
                    print(b2.regs[b2.qubitMap[1]][0].get_state())
                    print(a2_state == b2_state)
                    print("    Michael Bay visited your simulator...")
                del a2
                del b2
                raise AssertionError("Error comparing states after " +
                                     "applying gate")
            if verbose:
                print("    Noice")
            del a2
            del b2


def controlled_gate_tests(nq, verbose=False, QItem=qj.QRegistry):
    """Test application of controlled gates."""
    if verbose:
        print(" Controlled gate tests:")
    size = 2**nq
    for id in range(nq):
        isControl = [rnd.randint(-1, 1) for i in range(nq - 1)]
        controls = {i if i < id else i + 1
                    for i in range(nq - 1) if isControl[i] == 1}
        anticontrols = {i if i < id else i + 1
                        for i in range(nq - 1) if isControl[i] == -1}
        invert = bool(rnd.randint(0, 1))
        rands = np.random.rand(3) * 2 * np.pi - np.pi
        gate = "U(" + ",".join([str(angle)
                                for angle in rands]) + ")"
        if invert:
            gate += "-1"
        sgate = qj.SimpleGate(gate)
        if verbose:
            print(f"  Testing gate {gate} to qubit {id} and {isControl} ...")
        # print("    Gate: " + str(numpygate))
        b = QItem(nq, verbose=False)
        a = doki.registry_new(nq, False)
        a2 = doki.registry_apply(a, sgate.gate, [id], controls, anticontrols,
                                 -1, False)
        del a
        if verbose:
            print("    Doki done")
        b2 = b.apply_gate(gate, targets=id,
                          controls=controls, anticontrols=anticontrols)
        del b
        if verbose:
            print("    QSimov done")
        a2_state = np.array([doki.registry_get(a2, i, False, False)
                             for i in range(size)])
        b2_state = b2.get_state()
        if not np.allclose(a2_state, b2_state):
            if verbose:
                print("Expected:", a2_state)
                print("Received:", b2_state)
                print(a2_state == b2_state)
                print("Target:", id)
                print("Cs:", controls)
                print("ACs:", anticontrols)
                print("    Michael Bay visited your simulator...")
            del a2
            del b2
            raise AssertionError("Error comparing states after applying gate")
        if verbose:
            print("    Noice")
        del a2
        del b2


def measure_registry_tests(nq, verbose=False):
    """Test measurement with QRegistry."""
    if verbose:
        print(" Measure QRegistry tests:")
    for id in range(nq):
        reg = qj.QRegistry(nq)
        reg2 = reg.apply_gate("X", targets=id)
        del reg
        reg = reg2
        aux1, mes = reg.measure({id})
        if nq > 1:
            aux2, mes2 = aux1.measure({i for i in range(nq) if i != id})
            if aux2 is None:
                raise AssertionError("registry is never None after measure")
        else:
            aux2 = None
            mes2 = None
        del reg
        if (not mes[id]
                or (mes2 is not None and any(mes2))
                or aux1 is None
                or aux1.num_qubits != nq-1
                or aux1.num_bits != 1
                or (aux2 is not None
                    and (aux2.num_qubits != 0 or aux2.num_bits != nq))):
            if verbose:
                print("M1:", mes)
                print("M2:", mes2)
                print("Check1:", not mes[id])
                print("Check2:", any(mes2))
                print("Check3:", aux1.num_qubits != nq-1)
                print("Check4:", aux2 is not None)
                print("    Michael Bay visited your simulator...")
            del aux1
            del aux2
            raise AssertionError("Error measuring states")
        del aux1
        del aux2
    if verbose:
        print("    Noice")


def compare_state(r, state, rdm0, rdm1, srt0=1, srt1=1, verbose=False):
    """Compare states, reduced density matrices and reduced traces.

    Positional arguments:
        r: QRegistry
        state: numpy array with the expected state of r
        rdm0: numpy array with the expected reduced density matrix after
              tracing out qubit 0.
        rdm1: numpy array with the expected reduced density matrix after
              tracing out qubit 1.
    Keyworded arguments:
        srt0: expected value for squared reduced trace after tracing out 0
        srt1: expected value for squared reduced trace after tracing out 1
        verbose: if messages with extra information should be printed
    """
    if not np.allclose(r.get_state(), state):
        if verbose:
            print("State error")
            print(r.get_state())
            print(state)
            print(r.get_state() == state)
            print("    Michael Bay visited your simulator...")
        return False

    dm = np.dot(state.reshape((state.size, 1)),  # Ket
                state.conj().reshape((1, state.size)))  # Bra
    qdm = r.density_matrix()
    np_qdm = qdm[:]
    if not np.allclose(np_qdm, dm):
        if verbose:
            print("Density matrix error")
            print(np_qdm)
            print(dm)
            print(np_qdm == dm)
            print("    Michael Bay visited your simulator...")
        return False

    qrdm0 = qdm.partial_trace(0)
    np_qrdm0 = qrdm0[:]
    if not np.allclose(np_qrdm0, rdm0):
        if verbose:
            print("RDM0")
            print(np_qrdm0)
            print(rdm0)
            print(np_qrdm0 == rdm0)
            print("    Michael Bay visited your simulator...")
        return False

    qrdm1 = qdm.partial_trace(1)
    np_qrdm1 = qrdm1[:]
    if not np.allclose(np_qrdm1, rdm1):
        if verbose:
            print("RDM1")
            print(np_qrdm1)
            print(rdm1)
            print(np_qrdm1 == rdm1)
            print("    Michael Bay visited your simulator...")
        return False

    qsrt0 = (qrdm0 @ qrdm0).trace()
    if not np.allclose(qsrt0, srt0):
        if verbose:
            print("SRT0")
            print(qrdm0[:])
            print(qsrt0)
            print(srt0)
            print(qsrt0 == srt0)
            print("    Michael Bay visited your simulator...")
        return False

    qsrt1 = (qrdm1 @ qrdm1).trace()
    if not np.allclose(qsrt1, srt1):
        if verbose:
            print("SRT1")
            print(qsrt1)
            print(srt1)
            print(qsrt1 == srt1)
            print("    Michael Bay visited your simulator...")
        return False
    return True


def tool_test(verbose=False):
    """Test QRegistry states, density matrix, reduced dm and reduced trace."""
    if verbose:
        print(" Tools for QRegistry:")
    reg = qj.QRegistry(2)
    state = np.array([1, 0, 0, 0])
    rdm0 = np.array([1, 0, 0, 0]).reshape((2, 2))
    rdm1 = rdm0[:]
    if not compare_state(reg, state, rdm0, rdm1, verbose=verbose):
        del reg
        raise AssertionError("Error on first step checking tools")
    del state
    del rdm0
    del rdm1

    reg2 = reg.apply_gate("H", targets=0)
    del reg
    reg = reg2
    state = np.array([1/np.sqrt(2), 1/np.sqrt(2), 0, 0])
    rdm0 = np.array([1, 0, 0, 0]).reshape((2, 2))
    rdm1 = np.array([0.5, 0.5, 0.5, 0.5]).reshape((2, 2))
    if not compare_state(reg, state, rdm0, rdm1, verbose=verbose):
        del reg
        raise AssertionError("Error on second step checking tools")
    del state
    del rdm0
    del rdm1

    reg2 = reg.apply_gate("X", targets=1, controls=0)
    del reg
    reg = reg2
    state = np.array([1/np.sqrt(2), 0, 0, 1/np.sqrt(2)])
    rdm0 = np.eye(2) * 0.5
    rdm1 = rdm0[:]
    if not compare_state(reg, state, rdm0, rdm1, srt0=0.5, srt1=0.5,
                         verbose=verbose):
        del reg
        raise AssertionError("Error on third step checking tools")

    if verbose:
        print("    Noice")
    del reg


def measure_system_minitest(verbose=False):
    s = qj.QSystem(2)
    s = s.apply_gate("X", targets=0)
    s, res = s.measure({0})
    if res != [1, None]:
        if verbose:
            print(res)
        raise AssertionError()
    s, res = s.measure({1})
    if res != [None, 0]:
        if verbose:
            print(res)
        raise AssertionError()
    s, res = s.measure({0, 1})
    if res != [1, 0]:
        if verbose:
            print(res)
        raise AssertionError()
    s, res = s.measure({1, 0})
    if res != [1, 0]:
        if verbose:
            print(res)
        raise AssertionError()
    s = s.apply_gate("X", targets=1, controls=0)
    s, res = s.measure({1, 0})
    if res != [1, 1]:
        if verbose:
            print(res)
        raise AssertionError()
    s, res = s.measure({0})
    if res != [1, None]:
        if verbose:
            print(res)
        raise AssertionError()
    s, res = s.measure({1})
    if res != [None, 1]:
        if verbose:
            print(res)
        raise AssertionError()
    s, res = s.measure({0, 1})
    if res != [1, 1]:
        if verbose:
            print(res)
        raise AssertionError()
    s = s.apply_gate("X", targets=0)
    s = s.apply_gate("X", targets=1)
    s = s.apply_gate("H", targets=0)
    s = s.apply_gate("X", targets=1, controls=0)
    s, res = s.measure({0, 1})
    if res[0] != res[1]:
        if verbose:
            print(res)
        raise AssertionError()


def measure_system_tests(nq, entangle=False, remove=False, verbose=False):
    """Test measurement with QSystem."""
    if verbose:
        print(f" Measure QSystem tests with entangle={entangle}:")
    for id in range(nq):
        reg = qj.QSystem(nq)
        if entangle:
            for control in range(1, nq, 2):
                reg2 = reg.apply_gate("X", targets=control-1, controls=control)
            if nq % 2 == 1:
                reg2 = reg.apply_gate("X", targets=nq-2, controls=nq-1)
            del reg
            reg = reg2
        reg2 = reg.apply_gate("X", targets=id)
        del reg
        reg = reg2
        aux1, mes = reg.measure({id})
        mes2 = None
        if nq > 1:
            aux2, mes2 = aux1.measure({i for i in range(nq) if i != id})
        del reg
        check_bad_measure = not mes[id]
        check_bad_measures = (nq > 1 and any(mes2[i]
                              for i in range(nq) if i != id))
        check_classic = aux1.get_classic(id) is not None
        check_not_classics = not all(aux1.get_classic(i) is None
                                     for i in range(nq) if i != id)
        check_all_classic = nq > 1 and not any(aux2.get_classic(i) is None
                                               for i in range(nq))
        if (check_bad_measure or check_bad_measures or check_classic
                or check_not_classics or check_all_classic):
            if verbose:
                print("M1:", mes)
                print("M2:", mes2)
                print("Check1:", check_bad_measure)
                print("Check2:", check_bad_measures)
                print("Check3:", check_classic)
                print("Check4:", check_not_classics)
                print("Check5:", check_all_classic)
                print("    Michael Bay visited your simulator...")
            del aux1
            if nq > 1:
                del aux2
            raise AssertionError("Error measuring states")
        del aux1
        if nq > 1:
            del aux2
    if verbose:
        print("    Noice")


def add_operation_tests(qdesign, verbose=False):
    """Test add_line method of the given qstruct object."""
    if verbose:
        print(" add_line tests with " + qdesign.__name__ + ":")
    qdes = qdesign(5, 0, "Test")
    cons = {1, 3}
    acons = {2, 4}
    qdes.add_operation("X", targets=0, controls=cons, anticontrols=acons)
    if len(qdes.get_operations()) != 1:
        raise AssertionError("Wrong operations list size: " +
                             f"{len(qdes.get_operations())}")
    op_data = qdes.get_operations()[0]
    if op_data["gate"]._str != "X" or op_data["targets"] != [0] \
            or op_data["controls"] != cons or op_data["anticontrols"] != acons:
        if verbose:
            print(op_data)
            print("    Michael Bay visited your simulator...")
        raise AssertionError("Wrong operation added")
    if verbose:
        print("    Noice")


def _deutsch_aux(executor, nq, gate):
    circuit = DJAlgCircuit(nq, gate)
    mess = executor.execute(circuit)
    mes = mess[0]

    reg2 = qj.QSystem(nq)  # Qubits (x1, ..., xn, y) initialized to zero
    aux = reg2.apply_gate("X", targets=nq-1)  # Qubit y set to one
    del reg2
    reg2 = aux
    # Apply Hadamard to all qubits
    for i in range(nq):
        aux = reg2.apply_gate("H", targets=i)
        del reg2
        reg2 = aux
    # Applied U (oracle)
    aux = reg2.apply_gate(gate, targets=[i for i in range(nq)])
    del reg2
    reg2 = aux
    # Applied Hadamard to (x1, ..., xn), nothing applied to y qubit
    for i in range(nq - 1):
        aux = reg2.apply_gate("H", targets=i)
        del reg2
        reg2 = aux
    # We measure (x1, ..., xn) qubits
    aux, mes2 = reg2.measure({i for i in range(nq - 1)})
    del reg2
    del aux
    # If any qubit (x1, ..., xn) is 1, balanced. Otherwise constant.

    return mes, mes2[:-1]


def deutschTests(nq, verbose=False, useSystem=False, optimize=False):
    """Test Deutsch-Jozsa algorithm for the specified number of qubits."""
    if verbose:
        print(" Deutsch circuit (" + (qj.QSystem.__name__ if useSystem
                                      else qj.QRegistry.__name__) + "):")
    executor = qj.Drewom(qmachine="doki",
                         extra={"num_threads": -1,
                                "random_generator": np.random.rand,
                                "use_system": useSystem})
    for id in range(nq - 1):
        gate = Bal(nq, id)
        mes, mes2 = _deutsch_aux(executor, nq, gate)
        if not mes == mes2 or not any(mes):
            if verbose:
                print(mes)
                print(mes2)
                print(mes == mes2)
                print("    Michael Bay visited your simulator...")
            raise AssertionError("Error checking DJ results")
    for id in range(2):
        gate = Const(nq, twice=(id == 1))
        mes, mes2 = _deutsch_aux(executor, nq, gate)
        if not mes == mes2 or any(mes):
            if verbose:
                print(mes)
                print(mes2)
                print(mes == mes2)
                print("    Michael Bay visited your simulator...")
            raise AssertionError("Error checking DJ results")
    if verbose:
        print("    Noice")


def teleportation_tests(verbose=False, useSystem=False, num_shots=16, optimize=False):
    """Execute teleportation algorithm related tests."""
    rands = np.random.rand(3) * 2 * np.pi - np.pi
    gate = "U(" + ",".join([str(angle)
                            for angle in rands]) + ")"
    initialValue = rnd.randrange(2)
    if verbose:
        print(" Teleportation circuit (" + (qj.QSystem.__name__ if useSystem
                                            else qj.QRegistry.__name__) + "):")
        print("  Gate: " + gate)
        print("  Initial value: " + str(initialValue))
    executor = qj.Drewom(qmachine="doki",
                         extra={"num_threads": -1,
                                "random_generator": np.random.rand,
                                "use_system": useSystem,
                                "return_struct": True})
    circuit = TeleportationCircuit(gate)
    mess = executor.execute(circuit, shots=num_shots)
    for i in range(num_shots):
        reg, mes = mess[i]

        reg2 = qj.QRegistry(1)
        aux = reg2.apply_gate(gate)
        del reg2
        reg2 = aux
        rstate = None
        if useSystem:
            rstate = reg.regs[reg.qubitMap[2]][0].get_state()
        else:
            rstate = reg.get_state()

        if not np.allclose(rstate, reg2.get_state()):
            if verbose:
                print("Ops:", circuit.get_operations())
                print(reg.get_state())
                print(reg2.get_state())
                print(reg.get_state() == reg2.get_state())
                print(mes)
                print("    Michael Bay visited your simulator...")
            del reg
            del reg2
            raise AssertionError("Error checking teleportation result!")
        else:
            if verbose:
                print("    Noice")
        del reg
        del reg2


def get_gate(gate_name):
    try:
        return [qj.SimpleGate(gate_name), ""]
    except ValueError as ex:
        s = ex.args[0].split("and")
        max_args = int(s[-1].split(".")[-2].strip())
        min_args = int(s[-2].strip().split(" ")[-1])
        for num_args in range(min_args, max_args + 1):
            args = []
            for i in range(num_args):
                args.append(str(rnd.random()))
            str_args = f"({','.join(args)})"
            return [qj.SimpleGate(gate_name + str_args), str_args]


def all_gate_tests(seed=None, verbose=False):
    """Execute all gate tests."""
    gate_names = qj.get_available_gates()
    gates = {}
    for gate_name in gate_names:
        gates[gate_name] = get_gate(gate_name)
    gate_aliases = qj.get_gate_aliases()
    for alias in gate_aliases:
        gate_names = gate_aliases[alias]
        if type(gate_names) == str:
            gate_names = (gate_names,)
        for gate_name in gate_names:
            gate, str_args = gates[gate_name]
            new_gate = qj.SimpleGate(alias + str_args)
            if not np.array_equal(gate.get_matrix(), new_gate.get_matrix()):
                if verbose:
                    print("Gate:", gate_name)
                    print("Alias:", alias)
                    print("Args:", str_args)
                    print("Expected:", gate.get_matrix())
                    print("Found:", new_gate.get_matrix())
                raise AssertionError("Error comparing gate with alias")


def data_structure_tests(minqubits, maxqubits, seed=None, verbose=False,
                         QItem=qj.QRegistry):
    """Execute all data structure tests."""
    if not (seed is None):
        rnd.seed(seed)
        np.random.seed(seed)
    for nq in range(minqubits, maxqubits + 1):
        if verbose:
            print("Testing with " + str(nq) + " qubit " + QItem.__name__)
        one_gate_tests(nq, verbose=verbose, QItem=QItem)
        if nq >= 2:
            two_gate_tests(nq, verbose=verbose, QItem=QItem)
            controlled_gate_tests(nq, verbose=verbose, QItem=QItem)
        if QItem == qj.QRegistry:
            measure_registry_tests(nq, verbose=verbose)
        else:
            measure_system_tests(nq, entangle=False, verbose=verbose)
            if nq >= 2:
                measure_system_tests(nq, entangle=True, verbose=verbose)
    if QItem == qj.QRegistry:
        # get_state, density_matrix,
        # reduced_density_matrix and reduced_trace tests
        tool_test(verbose=verbose)


def high_level_tests(minqubits, maxqubits, seed=None, verbose=False):
    """Test high level structures: QGate and QCircuit."""
    if not (seed is None):
        rnd.seed(seed)
        np.random.seed(seed)

    if verbose:
        print("Testing QGate inversion and application")
    # Dagger/Inversion QGate tests
    inversion_tests(verbose=verbose)
    # Entanglement QGate with QRegistry tests
    entangle_tests(verbose=verbose, useSystem=False)
    # Entanglement QGate with QSystem tests
    entangle_tests(verbose=verbose, useSystem=True)
    if maxqubits > 1:
        for nq in range(2, maxqubits + 1):
            if verbose:
                print("Testing Deutsch with " + str(nq) + " qubit circuits")
            # Deutsch-Josza algorithm with QRegistry tests
            deutschTests(nq, verbose=verbose, useSystem=False)
            # Deutsch-Josza algorithm with QSystem tests
            deutschTests(nq, verbose=verbose, useSystem=True)
    # Teleportation algorithm with QRegistry tests
    teleportation_tests(verbose=verbose, useSystem=False)
    # Teleportation algorithm with QSystem tests
    teleportation_tests(verbose=verbose, useSystem=True)
    # Control and anticontrol check for QGate
    add_operation_tests(qj.QGate, verbose=verbose)
    # Control and anticontrol check for QCircuit
    add_operation_tests(qj.QCircuit, verbose=verbose)


def main():
    """Execute all tests."""
    argv = sys.argv[1:]
    if 2 <= len(argv) <= 4:
        minqubits = int(argv[0])
        if minqubits < 1:
            print("minimum number of qubits must be at least 1")
        maxqubits = int(argv[1])
        if maxqubits < minqubits:
            print("minimum number of qubits cannot be greater than maximum")
        verbose = False
        seed = None
        verbose_location = 3
        if len(argv) == 4 or (len(argv) == 3 and argv[2].isnumeric()):
            seed = int(argv[2])
        elif len(argv) == 3:
            verbose_location = 2
        if len(argv) == verbose_location + 1:
            verbose = argv[verbose_location].lower() == "true"
            if not verbose and argv[verbose_location].lower() != "false":
                print('verbose must either be "true" or "false" (without ")')

        if seed is None:
            seed = rnd.randrange(2**32 - 1)
        print("Seed:", seed)
        print("\tTesting Gates...")
        all_gate_tests(seed=seed, verbose=verbose)
        print("\tTesting QRegistry...")
        data_structure_tests(minqubits, maxqubits, seed=seed, verbose=verbose,
                             QItem=qj.QRegistry)
        print("\tTesting QSystem...")
        data_structure_tests(minqubits, maxqubits, seed=seed, verbose=verbose,
                             QItem=qj.QSystem)
        print("\tTesting QGate and QCircuit...")
        high_level_tests(minqubits, maxqubits, seed=seed, verbose=verbose)
        print("PEACE AND TRANQUILITY")
    else:
        print("Syntax: " + sys.argv[0] + " <minimum number of qubits (min 1)>",
              "<maximum number of qubits> <seed (optional)> " +
              "<verbose (optional)")


if __name__ == "__main__":
    main()
