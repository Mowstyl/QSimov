#!/usr/bin/python
"""Test module for QSimov."""

import sys
import qsimov as qj
# import webbrowser as wb
import random as rnd
import numpy as np
from operator import add


def gate_to_id(numpygate, id, nq):
    """Return the first column of the matrix result of applying gate->qubit.

    Positional arguments:
        numpygate: numpy array with the gate matrix
        id: identifier of the qubit the gate will be applied to
        nq: number of qubits in the system
    Return:
        The tensor product between I_n, numpygate, I_m will be calculated.
        I_n and I_m are identity matrices for n and m qubits, the qubits not
        affected by numpygate.
        The first column of this matrix is then returned.
    """
    nqg = int(np.log2(numpygate.shape[0]))
    if 0 < id < nq - nqg:
        a = np.kron(np.eye(2**(nq - id - nqg)),
                    np.kron(numpygate, np.eye(2**id)))
    elif id == 0:
        a = np.kron(np.eye(2**(nq - nqg)), numpygate)
    else:
        a = np.kron(numpygate, np.eye(2**id))
    # print("    NpGate: " + str(a))
    return a[:, 0]


def gate_to_reg_id(numpygate, id, nq, reg):
    """Return result of applying a gate to a registry.

    Positional arguments:
        numpygate: numpy array with the gate matrix
        id: identifier of the qubit the gate will be applied to
        nq: number of qubits in the system
        reg: the system state vector
    Return:
        The product (I_n x numpygate x I_m) . |reg>
        x -> Tensor product
        . -> Matrix product
    """
    nqg = int(np.log2(numpygate.shape[0]))
    if 0 < id < nq - nqg:
        a = np.kron(np.eye(2**(nq - id - nqg)),
                    np.kron(numpygate, np.eye(2**id)))
    elif id == 0:
        a = np.kron(np.eye(2**(nq - nqg)), numpygate)
    else:
        a = np.kron(numpygate, np.eye(2**id))
    return np.dot(a, reg.reshape(reg.size, 1)).reshape(reg.size)


def H():
    """Return numpy array with Hadamard gate for 1 qubit."""
    gate = np.ones(4, dtype=complex).reshape(2, 2)
    aux = 1/np.sqrt(2)
    gate[1, 1] = -1
    return gate * aux


def X():
    """Return numpy array with Pauli X gate (aka NOT)."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    gate[0, 1] = 1
    gate[1, 0] = 1
    return gate


def Y():
    """Return numpy array with PauliY gate."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    gate[0, 1] = -1j
    gate[1, 0] = 1j
    return gate


def Z():
    """Return numpy array with PauliZ gate."""
    gate = np.eye(2, dtype=complex)
    gate[1, 1] = -1
    return gate


def SqrtX(invert):
    """Return numpy array with sqrt(NOT) gate."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    if not invert:
        gate[0, 0] = 0.5 + 0.5j
        gate[0, 1] = 0.5 - 0.5j
        gate[1, 0] = 0.5 - 0.5j
        gate[1, 1] = 0.5 + 0.5j
    else:
        gate[0, 0] = 0.5 - 0.5j
        gate[0, 1] = 0.5 + 0.5j
        gate[1, 0] = 0.5 + 0.5j
        gate[1, 1] = 0.5 - 0.5j
    return gate


def R(angle, invert):
    """Return numpy array with R gate."""
    gate = np.eye(2, dtype=complex)
    if not invert:
        gate[1, 1] = np.cos(angle) + np.sin(angle) * 1j
    else:
        gate[1, 1] = np.cos(angle) - np.sin(angle) * 1j
    return gate


def RX(angle, invert):
    """Return numpy array with rotation gate around X axis."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    cosan = np.cos(angle/2)
    sinan = np.sin(angle/2)
    if not invert:
        mult = -1j
    else:
        mult = 1j
    gate[0, 0] = cosan
    gate[0, 1] = mult * sinan
    gate[1, 0] = mult * sinan
    gate[1, 1] = cosan
    return gate


def RY(angle, invert):
    """Return numpy array with rotation gate around Y axis."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    cosan = np.cos(angle/2)
    sinan = np.sin(angle/2)
    gate[0, 0] = cosan
    gate[1, 1] = cosan
    if not invert:
        gate[0, 1] = -sinan
        gate[1, 0] = sinan
    else:
        gate[0, 1] = sinan
        gate[1, 0] = -sinan
    return gate


def RZ(angle, invert):
    """Return numpy array with rotation gate around Z axis."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    if not invert:
        gate[0, 0] = np.cos(-angle/2) + np.sin(-angle/2) * 1j
        gate[1, 1] = np.cos(angle/2) + np.sin(angle/2) * 1j
    else:
        gate[0, 0] = np.cos(-angle/2) - np.sin(-angle/2) * 1j
        gate[1, 1] = np.cos(angle/2) - np.sin(angle/2) * 1j
    return gate


def RUnity(n, invert):
    """Return numpy array with nth root of unity rotation gate."""
    return R(2*np.pi/(2**n), invert)


def HalfDeutsch(angle, invert):
    """Return numpy array with a portion of the Deutsch gate.

    This gate, when double controlled, is called Deutsch gate.
    """
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    cosan = np.cos(angle) * 1j
    sinan = np.sin(angle)
    if invert:
        cosan = -cosan
    gate[0, 0] = cosan
    gate[0, 1] = sinan
    gate[1, 0] = sinan
    gate[1, 1] = cosan
    return gate


def U(angle1, angle2, angle3, invert):
    """Return numpy array with U gate (IBM)."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    cosan = np.cos(angle1/2)
    sinan = np.sin(angle1/2)
    mult = 1
    if invert:
        mult = -1
    gate[0, 0] = cosan
    if not invert:
        gate[0, 1] = -sinan * np.cos(angle3) - sinan * np.sin(angle3) * 1j
        gate[1, 0] = sinan * np.cos(angle2) + sinan * np.sin(angle2) * 1j
    else:
        gate[0, 1] = sinan * np.cos(angle2) - sinan * np.sin(angle2) * 1j
        gate[1, 0] = -sinan * np.cos(angle3) + sinan * np.sin(angle3) * 1j
    gate[1, 1] = cosan * np.cos(angle2+angle3) \
        + mult * cosan * np.sin(angle2 + angle3) * 1j
    return gate


def U3(angle1, angle2, angle3, invert):
    """Return numpy array with U gate (IBM)."""
    return U(angle1, angle2, angle3, invert)


def U2(angle1, angle2, invert):
    """Return numpy array with U2 gate (IBM)."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    mult = 1
    if invert:
        mult = -1
    gate[0, 0] = 1
    if not invert:
        gate[0, 1] = -np.cos(angle2) - np.sin(angle2) * 1j
        gate[1, 0] = np.cos(angle1) + np.sin(angle1) * 1j
    else:
        gate[0, 1] = -np.cos(angle2) + np.sin(angle2) * 1j
        gate[1, 0] = np.cos(angle1) - np.sin(angle1) * 1j
    gate[1, 1] = np.cos(angle1+angle2) + mult * np.sin(angle1+angle2) * 1j
    return gate * (1/np.sqrt(2))


def U1(angle, invert):
    """Return numpy array with U1 gate (IBM)."""
    return R(angle, invert)


def SWAP():
    """Return numpy array with SWAP gate."""
    gate = np.zeros(4 * 4, dtype=complex)
    gate = gate.reshape(4, 4)

    gate[0][0] = 1
    gate[1][2] = 1
    gate[2][1] = 1
    gate[3][3] = 1

    return gate


def iSWAP(invert):
    """Return numpy array with iSWAP gate."""
    gate = np.zeros(4 * 4, dtype=complex)
    gate = gate.reshape(4, 4)

    gate[0][0] = 1
    if not invert:
        gate[1][2] = 1j
        gate[2][1] = 1j
    else:
        gate[1][2] = -1j
        gate[2][1] = -1j
    gate[3][3] = 1

    return gate


def sqrtSWAP(invert):
    """Return numpy array with sqrt(SWAP) gate."""
    gate = np.zeros(4 * 4, dtype=complex)
    gate = gate.reshape(4, 4)

    gate[0][0] = 1
    if not invert:
        gate[1][1] = 0.5 + 0.5j
        gate[1][2] = 0.5 - 0.5j
        gate[2][1] = 0.5 - 0.5j
        gate[2][2] = 0.5 + 0.5j
    else:
        gate[1][1] = 0.5 - 0.5j
        gate[1][2] = 0.5 + 0.5j
        gate[2][1] = 0.5 + 0.5j
        gate[2][2] = 0.5 - 0.5j
    gate[3][3] = 1

    return gate


def xx(angle, invert):
    """Return numpy array with Ising Coupling XX gate."""
    gate = np.eye(4, dtype=complex)
    if not invert:
        gate[0, 3] = np.sin(angle)-np.cos(angle)*1j
        gate[1, 2] = -1j
        gate[2, 1] = -1j
        gate[3, 0] = np.sin(-angle)-np.cos(-angle)*1j
    else:
        gate[0, 3] = np.sin(-angle)+np.cos(-angle)*1j
        gate[1, 2] = 1j
        gate[2, 1] = 1j
        gate[3, 0] = np.sin(angle)+np.cos(angle)*1j
    return gate*(1/np.sqrt(2))


def yy(angle, invert):
    """Return numpy array with Ising Coupling YY gate."""
    gate = np.eye(4, dtype=complex)
    gate = gate * np.cos(angle)
    ansin = np.sin(angle) * 1j
    if not invert:
        gate[0, 3] = ansin
        gate[1, 2] = -ansin
        gate[2, 1] = -ansin
        gate[3, 0] = ansin
    else:
        gate[0, 3] = -ansin
        gate[1, 2] = ansin
        gate[2, 1] = ansin
        gate[3, 0] = -ansin
    return gate


def zz(angle, invert):
    """Return numpy array with Ising Coupling ZZ gate."""
    gate = np.eye(4, dtype=complex)
    gate = gate * np.cos(angle)
    phi2 = angle/2
    if not invert:
        gate[0, 0] = np.cos(phi2) + np.sin(phi2) * 1j
        gate[1, 1] = np.cos(-phi2) + np.sin(-phi2) * 1j
        gate[2, 2] = gate[1, 1]
        gate[3, 3] = gate[0, 0]
    else:
        gate[0, 0] = np.cos(phi2) - np.sin(phi2) * 1j
        gate[1, 1] = np.cos(-phi2) - np.sin(-phi2) * 1j
        gate[2, 2] = gate[1, 1]
        gate[3, 3] = gate[0, 0]
    return gate


def swap_downstairs(id1, id2, nq, reg):
    """Swap qubit id1 with next qubit until reaches id2 (id1 < id2)."""
    swap = SWAP()
    for i in range(id1, id2 + 1):
        reg = gate_to_reg_id(swap, i, nq, reg)
    return reg


def swap_upstairs(id1, id2, nq, reg):
    """Swap qubit id1 with next qubit until reaches id2 (id1 > id2)."""
    swap = SWAP()
    for i in range(id2, id1 - 1, -1):
        reg = gate_to_reg_id(swap, i, nq, reg)
    return reg


def swap_downstairs_list(id1, id2, li):
    """Swap list element id1 with the next until reaches id2 (id1 > id2)."""
    for i in range(id1, id2 + 1):
        li[i], li[i+1] = li[i+1], li[i]
    return li


def swap_upstairs_list(id1, id2, li):
    """Swap list element id1 with the next until reaches id2 (id1 < id2)."""
    for i in range(id2, id1 - 1, -1):
        li[i], li[i+1] = li[i+1], li[i]
    return li


def sparseTwoGate(gate, id1, id2, nq, reg):
    """TODO: Find out or recall what the hell is this doing."""
    if id2 < id1:
        id1, id2 = id2, id1
    if id1 < 0 or id2 >= nq:
        reg = None
    else:
        if id2 - id1 > 1:
            reg = swap_downstairs(id1, id2 - 1, nq, reg)
        reg = gate_to_reg_id(gate, id2 - 1, nq, reg)
        if id2 - id1 > 1:
            reg = swap_upstairs(id1, id2 - 1, nq, reg)
    return reg


def CU(gate, ncontrols):
    """Return n-controlled version of given gate."""
    nqgate = int(np.log2(gate.shape[0]))
    cu = np.eye(2**(nqgate+ncontrols), dtype=complex)
    aux = cu.shape[0] - gate.shape[0]
    for i in range(gate.shape[0]):
        for j in range(gate.shape[1]):
            cu[aux + i, aux + j] = gate[i, j]

    return cu


def negateQubits(qubits, nq, reg):
    """Apply X gate to qubit ids specified."""
    for id in qubits:
        reg = gate_to_reg_id(X(), id, nq, reg)
    return reg


def applyCACU(gate, id, controls, anticontrols, nq, reg):
    """Apply gate with specified controls and anticontrols."""
    cset = set(controls)
    acset = set(anticontrols)
    cuac = list(cset.union(acset))
    if type(id) == list:
        extended_cuac = id + cuac
    else:
        extended_cuac = [id] + cuac
    qubitIds = [i for i in range(nq)]

    reg = negateQubits(acset, nq, reg)
    for i in range(len(extended_cuac)):
        if qubitIds[i] != extended_cuac[i]:
            indaux = qubitIds.index(extended_cuac[i])
            reg = swap_upstairs(i, indaux - 1, nq, reg)
            qubitIds = swap_upstairs_list(i, indaux - 1, qubitIds)
    reg = gate_to_reg_id(CU(gate, len(cuac)), 0, nq, reg)
    for i in range(nq):
        if qubitIds[i] != i:
            indaux = qubitIds.index(i)
            reg = swap_upstairs(i, indaux - 1, nq, reg)
            qubitIds = swap_upstairs_list(i, indaux - 1, qubitIds)
    reg = negateQubits(acset, nq, reg)
    return reg


def Bal(n, controlId=0):
    """Return Deutsch-Jozsa oracle for balanced function."""
    gate = qj.QGate("Balanced")
    gate.add_line(*[None for i in range(n-1)], ("X", controlId, None))
    return gate


def Const(n, twice=False):
    """Return Deutsch-Jozsa oracle for constant function."""
    gate = qj.QGate("Constant")
    gate.add_line(*(None for i in range(n-1)), "X")
    if twice:
        gate.add_line(*[None for i in range(n-1)], "X")
    return gate


def DJAlgCircuit(size, U_f):
    """Devuelve el circuito del algoritmo de Deutsch-Josza.

    U_f es el oraculo, que debe tener x1..xn e y como qubits.
    Tras aplicarlo el qubit y debe valer f(x1..xn) XOR y.
    El argumento size es n + 1, donde n es el numero de bits de entrada de f.
    """
    # El ultimo QuBit al ejecutar el algoritmo es de ancilla, con su valor a 1
    c = qj.QCircuit("Deutsch-Josza Algorithm", ancilla=[1])

    # Se aplica una compuerta hadamard a todos los qubits
    c.add_line(*["H" for i in range(size)])
    c.add_line(U_f)  # Se aplica el oraculo

    # Se aplica una puerta Hadamard a todos los qubits excepto al ultimo
    c.add_line(*("H" for i in range(size-1)), "I")

    # Se miden los qubit x, si es igual a 0 la funcion es constante.
    # En caso contrario no lo es.
    c.add_line(qj.Measure([1 for i in range(size - 1)] + [0]))

    return c


def TeleportationCircuit(gate, remove=True):
    """Return teleportation algorithm circuit.

    Positional arguments:
        gate: gate to apply to the qubit that is going to be sent
            (so we don't send a 0 or a 1)
    Keyworded arguments:
        remove: whether or not to remove the measured qubits from the system
    Return:
        QCircuit with teleportation algorithm
    """
    qc = qj.QCircuit("Teleportation", ancilla=(0, 0))
    qc.add_line(None, "H", None)
    qc.add_line(None, None, ("X", 1, None))
    # Aqui es donde trabajamos con el qubit Q que queremos enviar.
    # Se le aplica la puerta pasada como parámetro.
    qc.add_line(gate, None, None)
    # Una vez terminado todo lo que queremos hacerle al QuBit,
    # procedemos a preparar el envio
    # Se aplica una puerta C-NOT sobre Q (control) y B (objetivo).
    qc.add_line(None, ["X", 0, None], None)
    # Se aplica una puerta Hadamard sobre Q.
    qc.add_line("H", None, None)

    gate1 = "X"
    gate2 = "Z"
    if not remove:
        gate1 = {"gate": "X", "qubit": 2}
        gate2 = {"gate": "Z", "qubit": 2}
    c1 = qj.Condition([None, 1, None], gate1, None, 1, -1)
    c2 = qj.Condition((1, None, None), gate2, None, 1, -1)
    m = qj.Measure((1, 1, 0), conds=[c1, c2], remove=remove)
    qc.add_line(m)

    return qc


def entangle_gate():
    """Return a QGate that creates a Bell pair."""
    e = qj.QGate("Entangle")

    e.add_line("H", "I")
    e.add_line(None, ("X", 0))

    return e


def entangle_system(s, id, id2):
    """Entangle specified qubits of a system."""
    s.apply_gate("H", qubit=id)
    s.apply_gate("X", qubit=id2, control=id)


def inversion_tests(verbose=False):
    """Test gate inversion."""
    passed = 0
    total = 3
    if verbose:
        print("  Testing gate inversion...")

    e = entangle_gate()  # TODO: Intensive inversion test
    ei = e.invert()

    if e.lines != ei.lines[::-1]:
        if verbose:
            print(e.lines)
            print(ei.lines)
            print("    Michael Bay visited your simulator...")
        return (passed, total)

    passed += 1
    if verbose:
        print("    Noice")

    ed = e.dagger()

    if e.lines != ed.lines[::-1]:
        if verbose:
            print(e.lines)
            print(ed.lines)
            print("    Michael Bay visited your simulator...")
        return (passed, total)

    passed += 1
    if verbose:
        print("    Noice")

    r = qj.QRegistry(2)
    er = qj.QRegistry(2)
    s = qj.QSystem(2)
    es = qj.QSystem(2)

    r.apply_gate(e)
    r.apply_gate(ed)
    s.apply_gate(e)
    s.apply_gate(ed)
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
        return (passed, total)

    passed += 1
    if verbose:
        print("    Noice")

    return (passed, total)


def entangle_tests(verbose=False, useSystem=False):
    """Test entanglement."""
    passed = 0
    total = 2
    if verbose:
        print("  Testing entanglement gate...")

    if useSystem:
        QItem = qj.QSystem
    else:
        QItem = qj.QRegistry

    e = entangle_gate()

    r = QItem(3)
    rg = QItem(3)
    entangle_system(r, 0, 1)
    rg.apply_gate(e)

    if any(r.get_state() != rg.get_state()):
        if verbose:
            print(r.get_state())
            print(rg.get_state())
            print(r.get_state() == rg.get_state())
            print("    Michael Bay visited your simulator...")
        del r
        del rg
        return (passed, total)

    passed += 1
    del r
    del rg
    if verbose:
        print("    Noice")

    r = qj.QRegistry(3)
    rg1 = QItem(3)
    entangle_system(r, 1, 2)
    rg1.apply_gate(e, qubit=1)
    if useSystem:
        rg2 = QItem(3)
        rg2.apply_gate("I", e)

    if (any(r.get_state() != rg1.get_state())
            or (useSystem and any(r.get_state() != rg2.get_state()))):
        if verbose:
            print(r.get_state())
            print(rg1.get_state())
            print(r.get_state() == rg1.get_state())
            if useSystem:
                print(rg2.get_state())
                print(r.get_state() == rg2.get_state())
            print("    Michael Bay visited your simulator...")
        del r
        del rg1
        if useSystem:
            del rg2
        return (passed, total)

    passed += 1
    del r
    del rg1
    if useSystem:
        del rg2
    if verbose:
        print("    Noice")

    return (passed, total)


# TODO: Test for gate with no controls in controlled QGate


def gate_tests(gatename, verbose=False, hasInv=False, nArgs=0):
    """Test gate application."""
    passed = 0
    if not hasInv:
        total = 1
    elif nArgs == 0:
        total = 2
    else:
        total = 2 * 5 * nArgs
    if verbose:
        print("  Testing gate " + gatename + "...")

    if not hasInv:
        numpygate = qj.get_gate(gatename)
        numpygate2 = globals()[gatename]()
        allOk = np.allclose(numpygate, numpygate2)
        if not allOk:
            if verbose:
                print(numpygate)
                print(numpygate2)
                print(numpygate == numpygate2)
                print("    Michael Bay visited your simulator...")
            return (passed, total)
        passed += 1
        if verbose:
            print("    Noice")
    elif nArgs == 0:
        numpygate = qj.get_gate(gatename)
        numpygate2 = globals()[gatename](False)
        allOk = np.allclose(numpygate, numpygate2)
        if not allOk:
            if verbose:
                print(numpygate)
                print(numpygate2)
                print(numpygate == numpygate2)
                print("    Michael Bay visited your simulator...")
            return (passed, total)
        passed += 1
        if verbose:
            print("    Noice, testing inverse")
        numpygate = qj.get_gate(gatename + "-1")
        numpygate2 = globals()[gatename](True)
        allOk = np.allclose(numpygate, numpygate2)
        if not allOk:
            if verbose:
                print(numpygate)
                print(numpygate2)
                print(numpygate == numpygate2)
                print("    Michael Bay visited your simulator...")
            return (passed, total)
        passed += 1
        if verbose:
            print("    Noice")
    else:
        for i in range(5 * nArgs):
            if gatename == "RUnity":
                args = [rnd.randint(1, 10)]
            else:
                args = [rnd.random() for i in range(nArgs)]
            if verbose:
                print("    args: " + str(args))
            args.append(False)
            argstr = str(args[0])
            for arg in args[1:-1]:
                argstr += "," + str(arg)
            argstr = "(" + argstr + ")"
            numpygate = qj.get_gate(gatename + argstr)
            numpygate2 = globals()[gatename](*args)
            allOk = np.allclose(numpygate, numpygate2)
            if not allOk:
                if verbose:
                    print(numpygate)
                    print(numpygate2)
                    print(numpygate == numpygate2)
                    print("    Michael Bay visited your simulator...")
                return (passed, total)
            passed += 1
            if verbose:
                print("    Noice")
        for i in range(5 * nArgs):
            if gatename == "RUnity":
                args = [rnd.randint(1, 10)]
            else:
                args = [rnd.random() for i in range(nArgs)]
            args.append(True)
            argstr = str(args[0])
            for arg in args[1:-1]:
                argstr += "," + str(arg)
            argstr = "(" + argstr + ")"
            numpygate = qj.get_gate(gatename + argstr + "-1")
            numpygate2 = globals()[gatename](*args)
            allOk = np.allclose(numpygate, numpygate2)
            if not allOk:
                if verbose:
                    print(numpygate)
                    print(numpygate2)
                    print(numpygate == numpygate2)
                    print("    Michael Bay visited your simulator...")
                return (passed, total)
            passed += 1
            if verbose:
                print("    Noice")

    return (passed, total)


def one_gate_tests(nq, verbose=False, QItem=qj.QRegistry):
    """Test application of one gate of one qubit."""
    passed = 0
    total = nq
    if verbose:
        print(" One gate tests:")
    for id in range(nq):
        gate = "U(" + str(rnd.random()) + "," + str(rnd.random()) \
            + "," + str(rnd.random()) + ")"
        numpygate = qj.get_gate(gate)
        if verbose:
            print("  Testing gate " + gate + " to qubit " + str(id) + "...")
        # print("    Gate: " + str(numpygate))
        b = QItem(nq)
        a = gate_to_id(numpygate, id, nq)
        if verbose:
            print("    Numpy done")
        b.apply_gate(gate, id)
        if verbose:
            print("    QSimov done")
        allOk = np.allclose(a, b.get_state())
        if not allOk:
            if verbose:
                print(a)
                print(b.get_state())
                print(a == b.get_state())
                print("    Michael Bay visited your simulator...")
            del a
            del b
            break
        passed += 1
        if verbose:
            print("    Noice")
        del a
        del b
    return (passed, total)


def simple_swap_tests(nq, imaginary, sqrt, invert,
                      verbose=False, QItem=qj.QRegistry):
    """Test swapping qubits next to each other."""
    passed = 0
    total = nq - 1
    gate2 = "U(" + str(rnd.random()) + "," + str(rnd.random()) \
        + "," + str(rnd.random()) + ")"
    numpygate2 = qj.get_gate(gate2)
    a = gate_to_id(numpygate2, 0, nq)
    b = QItem(nq)
    b.apply_gate(gate2)
    if verbose:
        if not imaginary and not sqrt:
            print(" Simple SWAP tests:")
            print("  SWAP gate:")
        elif imaginary:
            if not invert:
                print("  ISWAP gate:")
            else:
                print("  Inverse ISWAP gate:")
        else:
            if not invert:
                print("  sqrtSWAP gate:")
            else:
                print("  Inverse sqrtSWAP gate:")
        print("   Applied gate " + gate2 + " to qubit 0")
    for id in range(nq - 1):
        gate = "SWAP"
        numpygate = SWAP()
        invstring = ""
        verbstring = ""
        if invert:
            invstring = "-1"
            verbstring = "Un"
        if imaginary:
            gate = "i" + gate
            numpygate = iSWAP(invert)
        elif sqrt:
            gate = "sqrt" + gate
            numpygate = sqrtSWAP(invert)
        if verbose:
            print("   " + verbstring + gate + "ping qubits", str(id), "and",
                  str(id + 1) + "...")
        a = gate_to_reg_id(numpygate, id, nq, a)
        if verbose:
            print("    Numpy done")
        b.apply_gate(gate + "(" + str(id) + "," + str(id + 1) + ")"
                     + invstring)
        if verbose:
            print("    QSimov done")
        try:
            allOk = np.allclose(a, b.get_state())
        except ValueError:
            allOk = False
        if not allOk:
            if verbose:
                print(a)
                print(b.get_state())
                print(a == b.get_state())
                print("    Michael Bay visited your simulator...")
            del a
            del b
            break
        passed += 1
        if verbose:
            print("    Noice")
    if allOk:
        del a
        del b
    return (passed, total)


def sparse_swap_tests(nq, imaginary, sqrt, invert,
                      verbose=False, QItem=qj.QRegistry):
    """Test swapping qubits at any distance."""
    passed = 0
    total = ((nq - 1) * (nq - 2))//2
    # for i in range(2, nq):
    #     total += (nq - i)
    gate2 = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," \
        + str(rnd.random()) + ")"
    numpygate2 = qj.get_gate(gate2)
    if verbose:
        if not imaginary and not sqrt:
            print(" Sparse SWAP tests:")
            print("  SWAP gate:")
        elif imaginary:
            if not invert:
                print("  ISWAP gate:")
            else:
                print("  Inverse ISWAP gate:")
        else:
            if not invert:
                print("  sqrtSWAP gate:")
            else:
                print("  Inverse sqrtSWAP gate:")
    for id1 in range(nq - 1):
        for id2 in range(id1 + 2, nq):
            a = gate_to_id(numpygate2, id1, nq)
            b = QItem(nq)
            b.apply_gate(gate2, id1)
            if verbose:
                print("   Applied gate " + gate2 + " to qubit 0")
                # print("    Gate: " + str(numpygate2))
            gate = "SWAP"
            numpygate = SWAP()
            invstring = ""
            verbstring = ""
            if invert:
                invstring = "-1"
                verbstring = "Un"
            if imaginary:
                gate = "i" + gate
                numpygate = iSWAP(invert)
            elif sqrt:
                gate = "sqrt" + gate
                numpygate = sqrtSWAP(invert)
            if verbose:
                print("   " + verbstring + gate + "ping qubits", str(id1),
                      "and", str(id2) + "...")
            a = sparseTwoGate(numpygate, id1, id2, nq, a)
            if verbose:
                print("    Numpy done")
            b.apply_gate(gate + "(" + str(id1) + "," + str(id2) + ")"
                         + invstring)
            if verbose:
                print("    QSimov done")
            allOk = np.allclose(a, b.get_state())
            if not allOk:
                if verbose:
                    print(a)
                    print(b.get_state())
                    print(a == b.get_state())
                    print("    Michael Bay visited your simulator...")
                del a
                del b
                break
            passed += 1
            if verbose:
                print("    Noice")
    if allOk:
        del a
        del b
    return (passed, total)


def simple_ising_tests(nq, type, invert, verbose=False, QItem=qj.QRegistry):
    """Test coupling qubits next to each other."""
    passed = 0
    total = nq - 1
    gate = "XX"
    gatefun = xx
    invstring = ""
    if not invert:
        if type == 0:
            if verbose:
                print(" Simple Ising coupling tests:")
                print("  Ising (XX) coupling gate:")
        elif type == 1:
            gate = "YY"
            gatefun = yy
            if verbose:
                print("  Ising (YY) coupling gate:")
        else:
            gate = "ZZ"
            gatefun = zz
            if verbose:
                print("  Ising (ZZ) coupling gate:")
    else:
        invstring = "-1"
        if type == 0:
            if verbose:
                print("  Inverse of Ising (XX) coupling gate:")
        elif type == 1:
            gate = "YY"
            gatefun = yy
            if verbose:
                print("  Inverse of Ising (YY) coupling gate:")
        else:
            gate = "ZZ"
            gatefun = zz
            if verbose:
                print("  Inverse of Ising (ZZ) coupling gate:")

    gate2 = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," \
        + str(rnd.random()) + ")"
    numpygate2 = qj.get_gate(gate2)
    a = gate_to_id(numpygate2, 0, nq)
    b = QItem(nq)
    b.apply_gate(gate2)
    if verbose:
        print("   Applied gate " + gate2 + " to qubit 0")
        # print("    Gate: " + str(numpygate2))
    for id in range(nq - 1):
        angle = rnd.random()
        numpygate = gatefun(angle, invert)
        if verbose:
            print("   Coupling qubits", str(id), "and", str(id + 1),
                  "with angle", str(angle) + "...")
            # print("    Gate: " + str(numpygate))
        a = gate_to_reg_id(numpygate, id, nq, a)
        if verbose:
            print("    Numpy done")
        b.apply_gate(gate + "(" + str(angle) + "," + str(id) + ","
                     + str(id + 1) + ")" + invstring)
        if verbose:
            print("    QSimov done")
        allOk = np.allclose(a, b.get_state())
        if not allOk:
            if verbose:
                print(a)
                print(b.get_state())
                print(a == b.get_state())
                print("    Michael Bay visited your simulator...")
            del a
            del b
            break
        passed += 1
        if verbose:
            print("    Noice")
    if allOk:
        del a
        del b
    return (passed, total)


def sparse_ising_tests(nq, type, invert, verbose=False, QItem=qj.QRegistry):
    """Test swapping qubits at any distance."""
    passed = 0
    total = ((nq - 1) * (nq - 2))//2
    # for i in range(2, nq):
    #     total += (nq - i)
    gate = "XX"
    gatefun = xx
    invstring = ""
    if not invert:
        if type == 0:
            if verbose:
                print(" Sparse Ising coupling tests:")
                print("  Ising (XX) coupling gate:")
        elif type == 1:
            gate = "YY"
            gatefun = yy
            if verbose:
                print("  Ising (YY) coupling gate:")
        else:
            gate = "ZZ"
            gatefun = zz
            if verbose:
                print("  Ising (ZZ) coupling gate:")
    else:
        invstring = "-1"
        if type == 0:
            if verbose:
                print("  Inverse of Ising (XX) coupling gate:")
        elif type == 1:
            gate = "YY"
            gatefun = yy
            if verbose:
                print("  Inverse of Ising (YY) coupling gate:")
        else:
            gate = "ZZ"
            gatefun = zz
            if verbose:
                print("  Inverse of Ising (ZZ) coupling gate:")
    for id1 in range(nq - 1):
        for id2 in range(id1 + 2, nq):
            gate2 = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," \
                    + str(rnd.random()) + ")"
            numpygate2 = qj.get_gate(gate2)
            angle = rnd.random()
            numpygate = gatefun(angle, invert)
            a = gate_to_id(numpygate2, id1, nq)
            b = QItem(nq)
            b.apply_gate(gate2, id1)
            if verbose:
                print("   Applied gate", gate2, "to qubit 0")
            if verbose:
                print("   Coupling qubits", str(id1), "and", str(id2),
                      "with angle", str(angle) + "...")
            a = sparseTwoGate(numpygate, id1, id2, nq, a)
            if verbose:
                print("    Numpy done")
            b.apply_gate(gate + "(" + str(angle) + "," + str(id1) + ","
                         + str(id2) + ")" + invstring)
            if verbose:
                print("    QSimov done")
            allOk = np.allclose(a, b.get_state())
            if not allOk:
                if verbose:
                    print(a)
                    print(b.get_state())
                    print(a == b.get_state())
                    print("    Michael Bay visited your simulator...")
                del a
                del b
                break
            passed += 1
            if verbose:
                print("    Noice")
    if allOk:
        del a
        del b
    return (passed, total)


def controlled_gate_tests(nq, verbose=False, QItem=qj.QRegistry):
    """Test application of controlled gates."""
    total = nq - 1
    passed = 0
    isControl = bool(rnd.randint(0, 1))
    qubitIds = np.random.permutation(nq)
    lastid = qubitIds[0]
    control = []
    anticontrol = []
    gate = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," \
           + str(rnd.random()) + ")"
    if verbose:
        print(" Controlled gate tests:")
        print("  Gate: " + gate + " to qubit " + str(lastid))
    numpygate = qj.get_gate(gate)
    a = gate_to_id(numpygate, lastid, nq)
    b = QItem(nq)
    b.apply_gate(gate, lastid)
    for id in qubitIds[1:]:
        if isControl:
            control.append(lastid)
        else:
            anticontrol.append(lastid)
        if verbose:
            print("   id: " + str(id))
            print("   controls: " + str(control))
            print("   anticontrols: " + str(anticontrol))
        a = applyCACU(numpygate, id, control, anticontrol, nq, a)
        b.apply_gate(gate, qubit=id, control=control, anticontrol=anticontrol)
        isControl = not isControl
        lastid = id
        allOk = np.allclose(a, b.get_state())
        if not allOk:
            if verbose:
                print(a)
                print(b.get_state())
                print(a == b.get_state())
                print("    Michael Bay visited your simulator...")
            break
        passed += 1
        if verbose:
            print("    Noice")
    del a
    del b
    return (passed, total)


def controlled_swap_tests(nq, imaginary, sqrt, invert,
                          verbose=False, QItem=qj.QRegistry):
    """Test swapping qubits with control qubits."""
    total = nq - 1
    passed = 0
    isControl = bool(rnd.randint(0, 1))
    qubitIds = np.random.permutation(nq)
    lastid = qubitIds[:2]
    control = []
    anticontrol = []
    gate = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," \
           + str(rnd.random()) + ")"
    if verbose:
        print(" Controlled gate tests:")
        print("  Gate: " + gate + " to qubit " + str(lastid))
    numpygate = qj.get_gate(gate)
    a = gate_to_id(numpygate, lastid[0], nq)
    b = QItem(nq)
    b.apply_gate(gate, lastid[0])

    gate = "SWAP"
    numpygate = SWAP()
    invstring = ""
    if invert:
        invstring = "-1"
    if imaginary:
        gate = "i" + gate
        numpygate = iSWAP(invert)
    elif sqrt:
        gate = "sqrt" + gate
        numpygate = sqrtSWAP(invert)

    a = sparseTwoGate(numpygate, lastid[0], lastid[1], nq, a)
    b.apply_gate(gate + "(" + str(lastid[0]) + "," + str(lastid[1]) + ")"
                 + invstring)
    allOk = np.allclose(a, b.get_state())
    if not allOk:
        if verbose:
            print(a)
            print(b.get_state())
            print(a == b.get_state())
            print("    Michael Bay visited your simulator...")
    passed += 1
    if verbose:
        print("    Noice")

    if allOk:
        for id in qubitIds[2:]:
            if isControl:
                control.append(lastid[0])
            else:
                anticontrol.append(lastid[0])
            lastid = [lastid[1], id]
            isControl = not isControl
            if verbose:
                print("   ids: " + str(lastid))
                print("   controls: " + str(control))
                print("   anticontrols: " + str(anticontrol))
            a = applyCACU(numpygate, lastid, control, anticontrol, nq, a)
            b.apply_gate(gate + "(" + str(lastid[0]) + "," + str(lastid[1])
                         + ")" + invstring,
                         control=control, anticontrol=anticontrol)
            allOk = np.allclose(a, b.get_state())
            if not allOk:
                if verbose:
                    print(a)
                    print(b.get_state())
                    print(a == b.get_state())
                    print("    Michael Bay visited your simulator...")
                break
            passed += 1
            if verbose:
                print("    Noice")
    del a
    del b
    return (passed, total)


def controlled_ising_tests(nq, type, invert,
                           verbose=False, QItem=qj.QRegistry):
    """Test coupling qubits with control qubits."""
    total = nq - 1
    passed = 0
    isControl = bool(rnd.randint(0, 1))
    qubitIds = np.random.permutation(nq)
    lastid = qubitIds[:2]
    control = []
    anticontrol = []
    gate = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," \
           + str(rnd.random()) + ")"
    if verbose:
        print(" Controlled gate tests:")
        print("  Gate: " + gate + " to qubit " + str(lastid))
    numpygate = qj.get_gate(gate)
    a = gate_to_id(numpygate, lastid[0], nq)
    b = QItem(nq)
    b.apply_gate(gate, lastid[0])

    gate = "XX"
    gatefun = xx
    invstring = ""
    if type == 1:
        gate = "YY"
        gatefun = yy
    elif type == 2:
        gate = "ZZ"
        gatefun = zz
    if invert:
        invstring = "-1"
    angle = rnd.random()
    a = sparseTwoGate(gatefun(angle, invert), lastid[0], lastid[1], nq, a)
    b.apply_gate(gate + "(" + str(angle) + "," + str(lastid[0]) + ","
                 + str(lastid[1]) + ")" + invstring)
    allOk = np.allclose(a, b.get_state())
    if not allOk:
        if verbose:
            print(a)
            print(b.get_state())
            print(a == b.get_state())
            print("    Michael Bay visited your simulator...")
    passed += 1
    if verbose:
        print("    Noice")

    if allOk:
        for id in qubitIds[2:]:
            if isControl:
                control.append(lastid[0])
            else:
                anticontrol.append(lastid[0])
            lastid = [lastid[1], id]
            isControl = not isControl
            if verbose:
                print("   ids: " + str(lastid))
                print("   controls: " + str(control))
                print("   anticontrols: " + str(anticontrol))
            a = applyCACU(gatefun(angle, invert), lastid,
                          control, anticontrol, nq, a)
            b.apply_gate(gate + "(" + str(angle) + "," + str(lastid[0]) + ","
                         + str(lastid[1]) + ")" + invstring,
                         control=control, anticontrol=anticontrol)
            allOk = np.allclose(a, b.get_state())
            if not allOk:
                if verbose:
                    print(a)
                    print(b.get_state())
                    print(a == b.get_state())
                    print("    Michael Bay visited your simulator...")
                break
            passed += 1
            if verbose:
                print("    Noice")
    del a
    del b
    return (passed, total)


def measure_registry_tests(nq, remove=False, verbose=False):
    """Test measurement with QRegistry."""
    passed = 0
    total = nq
    if verbose:
        print(" Measure QRegistry tests with remove=" + str(remove) + ":")
    for id in range(nq):
        reg = qj.QRegistry(nq)
        reg.apply_gate("X", qubit=id)
        eres = [0 if i != id else 1 for i in range(nq)]
        mes = reg.measure(eres, remove=remove)
        mes2 = reg.measure([1 for i in range(nq if not remove else nq-1)])
        if remove:
            eres = [0 for i in range(nq-1)]

        if (not mes[0] == 1
                or mes2 != eres
                or (remove and reg.get_state().size != 2**(nq-1))
                or (not remove and reg.get_state().size != 2**nq)):
            if verbose:
                print(eres)
                print(mes)
                print(mes2)
                print(not mes[0] == 1)
                print(mes2 != eres)
                print(reg.get_state())
                print(reg.get_state().size)
                if (remove):
                    print(2**(nq-1))
                else:
                    print(2**nq)
                print(not mes[0] == 1)
                print(mes2 != eres)
                print(remove and reg.get_state().size != 2**(nq-1))
                print(not remove and reg.get_state().size != 2**nq)
                print("    Michael Bay visited your simulator...")
            del reg
            break
        passed += 1
        if verbose:
            print("    Noice")
        del reg

    return (passed, total)


def compare_state(r, state, rdm0, rdm1, rt0=1, rt1=1, verbose=False):
    """Compare states, reduced density matrices and reduced traces.

    Positional arguments:
        r: QRegistry
        state: numpy array with the expected state of r
        rdm0: numpy array with the expected reduced density matrix after
              tracing out qubit 0.
        rdm1: numpy array with the expected reduced density matrix after
              tracing out qubit 1.
    Keyworded arguments:
        rt0: expected value for reduced trace after tracing out qubit 0
        rt1: expected value for reduced trace after tracing out qubit 1
        verbose: if messages with extra information should be printed
    """
    if not np.allclose(r.get_state(), state):
        if verbose:
            print(r.get_state())
            print(state)
            print(r.get_state() == state)
            print("    Michael Bay visited your simulator...")
        return False

    dm = state * state.reshape((4, 1))
    if not np.allclose(r.density_matrix()[:], dm):
        if verbose:
            print(r.density_matrix()[:])
            print(dm)
            print(r.density_matrix()[:] == dm)
            print("    Michael Bay visited your simulator...")
        return False

    if not np.allclose(r.reduced_density_matrix(0)[:], rdm0):
        if verbose:
            print("RDM0")
            print(r.reduced_density_matrix(0)[:])
            print(rdm0)
            print(r.reduced_density_matrix(0)[:] == rdm0)
            print("    Michael Bay visited your simulator...")
        return False

    if not np.allclose(r.reduced_density_matrix(1)[:], rdm1):
        if verbose:
            print("RDM1")
            print(r.reduced_density_matrix(1)[:])
            print(rdm1)
            print(r.reduced_density_matrix(1)[:] == rdm1)
            print("    Michael Bay visited your simulator...")
        return False

    if not np.allclose(r.reduced_trace(0), rt0):
        if verbose:
            print("RT0")
            print(r.reduced_trace(0))
            print(rt0)
            print(r.reduced_trace(0) == rt0)
            print("    Michael Bay visited your simulator...")
        return False

    if not np.allclose(r.reduced_trace(1), rt1):
        if verbose:
            print("RT1")
            print(r.reduced_trace(1))
            print(rt1)
            print(r.reduced_trace(1) == rt1)
            print("    Michael Bay visited your simulator...")
        return False
    return True


def tool_test(verbose=False):
    """Test QRegistry states, density matrix, reduced dm and reduced trace."""
    passed = 0
    total = 3
    if verbose:
        print(" Tools for QRegistry:")
    reg = qj.QRegistry(2)
    state = np.array([1, 0, 0, 0])
    rdm0 = np.array([1, 0, 0, 0]).reshape((2, 2))
    rdm1 = rdm0[:]
    if not compare_state(reg, state, rdm0, rdm1, verbose=verbose):
        del reg
        return (passed, total)
    passed += 1

    del state
    del rdm0
    del rdm1

    reg.apply_gate("H")
    state = np.array([1/np.sqrt(2), 1/np.sqrt(2), 0, 0])
    rdm0 = np.array([1, 0, 0, 0]).reshape((2, 2))
    rdm1 = np.array([0.5, 0.5, 0.5, 0.5]).reshape((2, 2))
    if not compare_state(reg, state, rdm0, rdm1, verbose=verbose):
        del reg
        return (passed, total)
    passed += 1

    del state
    del rdm0
    del rdm1

    reg.apply_gate("X", qubit=1, control=0)
    state = np.array([1/np.sqrt(2), 0, 0, 1/np.sqrt(2)])
    rdm0 = np.eye(2) * 0.5
    rdm1 = rdm0[:]
    if not compare_state(reg, state, rdm0, rdm1, rt0=0.5, rt1=0.5,
                         verbose=verbose):
        del reg
        return (passed, total)
    passed += 1

    if verbose:
        print("    Noice")
    del reg

    return (passed, total)


def measure_system_tests(nq, entangle=False, remove=False, verbose=False):
    """Test measurement with QSystem."""
    passed = 0
    total = nq
    if verbose:
        print(" Measure QSystem tests with remove=" + str(remove) + ":")
    for id in range(nq):
        reg = qj.QSystem(nq)
        if entangle:
            for control in range(1, nq, 2):
                reg.apply_gate("X", qubit=control-1, control=control)
            if nq % 2 == 1:
                reg.apply_gate("X", qubit=nq-2, control=nq-1)
        reg.apply_gate("X", qubit=id)
        eres = [0 if i != id else 1 for i in range(nq)]
        mes = reg.measure(eres)
        mes2 = reg.measure([1 for i in range(nq)])
        r2 = qj.QRegistry(nq)
        r2.apply_gate("X", qubit=id)
        r2.measure(eres)

        if (not mes[0] == 1
                or mes2 != eres
                or not all(reg.get_state() == r2.get_state())):
            if verbose:
                print(eres)
                print(mes)
                print(mes2)
                print(mes[0] == 1)
                print(mes2 != eres)
                print("    Michael Bay visited your simulator...")
            del reg
            break
        passed += 1
        if verbose:
            print("    Noice")
        del reg

    return (passed, total)


def add_line_tests(qstruct, verbose=False):
    """Test add_line method of the given qstruct object."""
    passed = 0
    total = 1
    if verbose:
        print(" add_line tests with " + qstruct.__name__ + ":")
    qstr = qstruct()
    cons = set([1, 3])
    acons = set([2, 4])
    qstr.add_line(["X", cons, acons])
    auxg, auxc, auxa = qstr.lines[0][0]

    if auxg != "X" or auxc != cons or auxa != acons:
        if verbose:
            print(auxg)
            print(auxc)
            print(auxa)
            print("    Michael Bay visited your simulator...")
    else:
        passed += 1

    if verbose:
        print("    Noice")

    return (passed, total)


def deutschTests(nq, verbose=False, useSystem=False, optimize=False):
    """Test Deutsch-Jozsa algorithm for the specified number of qubits."""
    passed = 0
    total = nq - 1 + 2
    allOk = True
    if verbose:
        print(" Deutsch circuit (" + (qj.QSystem.__name__ if useSystem
                                      else qj.QRegistry.__name__) + "):")
    for id in range(nq - 1):
        gate = Bal(nq, id)
        circuit = DJAlgCircuit(nq, gate)
        reg, mes = circuit.execute([0 for i in range(nq - 1)],
                                   args={"useSystem": useSystem},
                                   optimize=optimize)
        mes = mes[0]

        reg2 = qj.QSystem(nq)  # Qubits (x1, ..., xn, y) initialized to zero
        reg2.apply_gate("X", qubit=nq-1)  # Qubit y set to one
        reg2.apply_gate("H")  # Applied Hadamard to first qubit (TEST)
        # And then to the rest of the qubits
        reg2.apply_gate(None, *["H" for i in range(nq-1)])
        reg2.apply_gate(gate)  # Applied U (oracle)
        # Applied Hadamard to (x1, ..., xn), nothing applied to y qubit
        reg2.apply_gate(*["H" for i in range(nq-1)], "I")
        # We measure (x1, ..., xn) qubits
        mes2 = reg2.measure([1 for i in range(nq - 1)] + [0])
        # And we add a None to the measure result, representing y not being
        # measured. QCircuit does this automatically.
        mes2 += [None]
        # If any qubit (x1, ..., xn) is 1, balanced. Otherwise constant.

        if not all(reg.get_state() == reg2.get_state()) or not mes == mes2:
            if verbose:
                print(reg.get_state())
                print(reg2.get_state())
                print(reg.get_state() == reg2.get_state())
                print(mes)
                print(mes2)
                print(mes == mes2)
                allOk = False
                print("    Michael Bay visited your simulator...")
            del reg
            del reg2
            break
        passed += 1
        if verbose:
            print("    Noice")
        del reg
        del reg2
    if allOk:
        for id in range(2):
            gate = Const(nq, twice=id == 1)
            circuit = DJAlgCircuit(nq, gate)
            reg, mes = circuit.execute([0 for i in range(nq - 1)],
                                       args={"useSystem": useSystem},
                                       optimize=optimize)
            mes = mes[0]

            reg2 = qj.QSystem(nq)
            reg2.apply_gate("X", qubit=nq-1)
            reg2.apply_gate(*["H" for i in range(nq)])
            reg2.apply_gate(gate)
            reg2.apply_gate(*["H" for i in range(nq-1)], "I")
            mes2 = reg2.measure([1 for i in range(nq - 1)] + [0])
            mes2 += [None]

            if not all(reg.get_state() == reg2.get_state()) or not mes == mes2:
                if verbose:
                    print(reg.get_state())
                    print(reg2.get_state())
                    print(reg.get_state() == reg2.get_state())
                    print(mes)
                    print(mes2)
                    print(mes == mes2)
                    allOk = False
                    print("    Michael Bay visited your simulator...")
                del reg
                del reg2
                break
            passed += 1
            if verbose:
                print("    Noice")
            del reg
            del reg2

    return (passed, total)


def teleportation_tests(verbose=False, useSystem=False, remove=False,
                        optimize=False):
    """Execute teleportation algorithm related tests."""
    passed = 0
    total = 1
    gate = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," \
           + str(rnd.random()) + ")"
    initialValue = rnd.randrange(2)
    if verbose:
        print(" Teleportation circuit (" + (qj.QSystem.__name__ if useSystem
                                            else qj.QRegistry.__name__) + ")" +
              (" with remove" if remove else "") + ":")
        print("  Gate: " + gate)
        print("  Initial value: " + str(initialValue))
    gate = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," \
           + str(rnd.random()) + ")"
    initialValue = rnd.randrange(2)
    circuit = TeleportationCircuit(gate, remove=remove)
    reg, mes = circuit.execute([initialValue], args={"useSystem": useSystem},
                               optimize=optimize)
    mes = mes[0]

    if remove:
        reg2 = qj.QRegistry(1)
        if initialValue == 1:
            reg2.apply_gate("X")
        reg2.apply_gate(gate)
    else:
        reg2 = qj.QRegistry(3)
        if initialValue == 1:
            reg2.apply_gate("X", qubit=2)
        reg2.apply_gate(gate, qubit=2)
        if mes[0] == 1:
            reg2.apply_gate("X", qubit=0)
        if mes[1] == 1:
            reg2.apply_gate("X", qubit=1)

    if not np.allclose(reg.get_state(), reg2.get_state()):
        if verbose:
            print(reg.get_state())
            print(reg2.get_state())
            print(reg.get_state() == reg2.get_state())
            print(mes)
            print("    Michael Bay visited your simulator...")
    else:
        passed += 1
        if verbose:
            print("    Noice")
    del reg
    del reg2

    return (passed, total)


def all_gate_tests(seed=None, verbose=False):
    """Execute all gate tests."""
    if not (seed is None):
        qj.set_seed(seed)
        rnd.seed(seed)
        np.random.seed(seed)
    result = [(0, 0) for i in range(15)]  # We have 15 tests

    # H gate tests
    result[0] = gate_tests("H", verbose=verbose, hasInv=False, nArgs=0)
    # X gate tests
    result[1] = gate_tests("X", verbose=verbose, hasInv=False, nArgs=0)
    # Y gate tests
    result[2] = gate_tests("Y", verbose=verbose, hasInv=False, nArgs=0)
    # Z gate tests
    result[3] = gate_tests("Z", verbose=verbose, hasInv=False, nArgs=0)
    # SqrtX gate tests
    result[4] = gate_tests("SqrtX", verbose=verbose, hasInv=True, nArgs=0)
    # RX gate tests
    result[5] = gate_tests("RX", verbose=verbose, hasInv=True, nArgs=1)
    # RY gate tests
    result[6] = gate_tests("RY", verbose=verbose, hasInv=True, nArgs=1)
    # RZ gate tests
    result[7] = gate_tests("RZ", verbose=verbose, hasInv=True, nArgs=1)
    # Phase shift gate tests
    result[8] = gate_tests("R", verbose=verbose, hasInv=True, nArgs=1)
    # Roots of unity gate tests
    result[9] = gate_tests("RUnity", verbose=verbose, hasInv=True, nArgs=1)
    # Partial Deutsch gate tests
    result[10] = gate_tests("HalfDeutsch", verbose=verbose,
                            hasInv=True, nArgs=1)
    # U gate tests
    result[11] = gate_tests("U", verbose=verbose, hasInv=True, nArgs=3)
    # U3 gate tests
    result[12] = gate_tests("U3", verbose=verbose, hasInv=True, nArgs=3)
    # U2 gate tests
    result[13] = gate_tests("U2", verbose=verbose, hasInv=True, nArgs=2)
    # U1 gate tests
    result[14] = gate_tests("U1", verbose=verbose, hasInv=True, nArgs=1)

    return result


def data_structure_tests(minqubits, maxqubits, seed=None, verbose=False,
                         QItem=qj.QRegistry):
    """Execute all data structure tests."""
    if not (seed is None):
        qj.set_seed(seed)
        rnd.seed(seed)
        np.random.seed(seed)
    result = [(0, 0) for i in range(38)]  # We have 37 tests with QRegistry
    if QItem == qj.QSystem:
        result += [(0, 0), (0, 0)]  # We have 39 tests with QSystem

    for nq in range(minqubits, maxqubits + 1):
        if verbose:
            print("Testing with " + str(nq) + " qubit registries")
        # Apply U gate tests
        result[0] = map(add, result[0],
                        one_gate_tests(nq, verbose=verbose, QItem=QItem))
        # Simple SWAP gate tests
        result[1] = map(add, result[1],
                        simple_swap_tests(nq, False, False, False,
                                          verbose=verbose, QItem=QItem))
        # Simple iSWAP gate tests
        result[2] = map(add, result[2],
                        simple_swap_tests(nq, True, False, False,
                                          verbose=verbose, QItem=QItem))
        # Simple sqrtSWAP gate tests
        result[3] = map(add, result[3],
                        simple_swap_tests(nq, False, True, False,
                                          verbose=verbose, QItem=QItem))
        # Simple iSWAP-1 gate tests
        result[4] = map(add, result[4],
                        simple_swap_tests(nq, True, False, True,
                                          verbose=verbose, QItem=QItem))
        # Simple sqrtSWAP-1 gate tests
        result[5] = map(add, result[5],
                        simple_swap_tests(nq, False, True, True,
                                          verbose=verbose, QItem=QItem))
        # Sparse SWAP gate tests
        result[6] = map(add, result[6],
                        sparse_swap_tests(nq, False, False, False,
                                          verbose=verbose, QItem=QItem))
        # Sparse iSWAP gate tests
        result[7] = map(add, result[7],
                        sparse_swap_tests(nq, True, False, False,
                                          verbose=verbose, QItem=QItem))
        # Sparse sqrtSWAP gate tests
        result[8] = map(add, result[8],
                        sparse_swap_tests(nq, False, True, False,
                                          verbose=verbose, QItem=QItem))
        # Sparse iSWAP-1 gate tests
        result[9] = map(add, result[9],
                        sparse_swap_tests(nq, True, False, True,
                                          verbose=verbose, QItem=QItem))
        # Sparse sqrtSWAP-1 gate tests
        result[10] = map(add, result[10],
                         sparse_swap_tests(nq, False, True, True,
                                           verbose=verbose, QItem=QItem))
        # Simple XX gate tests
        result[11] = map(add, result[11],
                         simple_ising_tests(nq, 0, False,
                                            verbose=verbose, QItem=QItem))
        # Simple YY gate tests
        result[12] = map(add, result[12],
                         simple_ising_tests(nq, 1, False,
                                            verbose=verbose, QItem=QItem))
        # Simple ZZ gate tests
        result[13] = map(add, result[13],
                         simple_ising_tests(nq, 2, False,
                                            verbose=verbose, QItem=QItem))
        # Simple XX-1 gate tests
        result[14] = map(add, result[14],
                         simple_ising_tests(nq, 0, True,
                                            verbose=verbose, QItem=QItem))
        # Simple YY-1 gate tests
        result[15] = map(add, result[15],
                         simple_ising_tests(nq, 1, True,
                                            verbose=verbose, QItem=QItem))
        # Simple ZZ-1 gate tests
        result[16] = map(add, result[16],
                         simple_ising_tests(nq, 2, True,
                                            verbose=verbose, QItem=QItem))
        # Sparse XX gate tests
        result[17] = map(add, result[17],
                         sparse_ising_tests(nq, 0, False,
                                            verbose=verbose, QItem=QItem))
        # Sparse YY gate tests
        result[18] = map(add, result[18],
                         sparse_ising_tests(nq, 1, False,
                                            verbose=verbose, QItem=QItem))
        # Sparse ZZ gate tests
        result[19] = map(add, result[19],
                         sparse_ising_tests(nq, 2, False,
                                            verbose=verbose, QItem=QItem))
        # Sparse XX-1 gate tests
        result[20] = map(add, result[20],
                         sparse_ising_tests(nq, 0, True,
                                            verbose=verbose, QItem=QItem))
        # Sparse YY-1 gate tests
        result[21] = map(add, result[21],
                         sparse_ising_tests(nq, 1, True,
                                            verbose=verbose, QItem=QItem))
        # Sparse ZZ-1 gate tests
        result[22] = map(add, result[22],
                         sparse_ising_tests(nq, 2, True,
                                            verbose=verbose, QItem=QItem))
        # Controlled U gate tests
        result[23] = map(add, result[23],
                         controlled_gate_tests(nq,
                                               verbose=verbose, QItem=QItem))
        # Controlled SWAP gate tests
        result[24] = map(add, result[24],
                         controlled_swap_tests(nq, False, False, False,
                                               verbose=verbose, QItem=QItem))
        # Controlled iSWAP gate tests
        result[25] = map(add, result[25],
                         controlled_swap_tests(nq, True, False, False,
                                               verbose=verbose, QItem=QItem))
        # Controlled sqrtSWAP gate tests
        result[26] = map(add, result[26],
                         controlled_swap_tests(nq, False, True, False,
                                               verbose=verbose, QItem=QItem))
        # Controlled iSWAP-1 gate tests
        result[27] = map(add, result[27],
                         controlled_swap_tests(nq, True, False, True,
                                               verbose=verbose, QItem=QItem))
        # Controlled sqrtSWAP-1 gate tests
        result[28] = map(add, result[28],
                         controlled_swap_tests(nq, False, True, True,
                                               verbose=verbose, QItem=QItem))
        # Controlled XX gate tests
        result[29] = map(add, result[29],
                         controlled_ising_tests(nq, 0, False,
                                                verbose=verbose, QItem=QItem))
        # Controlled YY gate tests
        result[30] = map(add, result[30],
                         controlled_ising_tests(nq, 1, False,
                                                verbose=verbose, QItem=QItem))
        # Controlled ZZ gate tests
        result[31] = map(add, result[31],
                         controlled_ising_tests(nq, 2, False,
                                                verbose=verbose, QItem=QItem))
        # Controlled XX-1 gate tests
        result[32] = map(add, result[32],
                         controlled_ising_tests(nq, 0, True,
                                                verbose=verbose, QItem=QItem))
        # Controlled YY-1 gate tests
        result[33] = map(add, result[33],
                         controlled_ising_tests(nq, 1, True,
                                                verbose=verbose, QItem=QItem))
        # Controlled ZZ-1 gate tests
        result[34] = map(add, result[34],
                         controlled_ising_tests(nq, 2, True,
                                                verbose=verbose, QItem=QItem))
        if QItem == qj.QRegistry:
            # Registry measurement tests
            result[35] = map(add, result[35],
                             measure_registry_tests(nq, remove=False,
                                                    verbose=verbose))
            # Registry measurement tests
            result[36] = map(add, result[36],
                             measure_registry_tests(nq, remove=True,
                                                    verbose=verbose))
        else:
            # QSystem measurement tests without entanglement
            result[35] = map(add, result[35],
                             measure_system_tests(nq, remove=False,
                                                  entangle=False,
                                                  verbose=verbose))
            # QSystem measurement tests with entanglement
            result[36] = map(add, result[36],
                             measure_system_tests(nq, remove=False,
                                                  entangle=True,
                                                  verbose=verbose))
            # QSystem measurement tests without entanglement
            result[37] = map(add, result[37],
                             measure_system_tests(nq, remove=True,
                                                  entangle=False,
                                                  verbose=verbose))
            # QSystem measurement tests with entanglement
            result[38] = map(add, result[38],
                             measure_system_tests(nq, remove=True,
                                                  entangle=True,
                                                  verbose=verbose))

    if QItem == qj.QRegistry:
        # get_state, density_matrix,
        # reduced_density_matrix and reduced_trace tests
        result[37] = tool_test(verbose=verbose)

    for i in range(len(result)):
        result[i] = tuple(result[i])

    return result


def high_level_tests(minqubits, maxqubits, seed=None, verbose=False):
    """Test high level structures: QGate and QCircuit."""
    if not (seed is None):
        qj.set_seed(seed)
        rnd.seed(seed)
        np.random.seed(seed)
    result = [(0, 0) for i in range(17)]  # We have 11 tests

    if verbose:
        print("Testing QGate inversion and application")
    # Dagger/Inversion QGate tests
    result[0] = inversion_tests(verbose=verbose)
    # Entanglement QGate with QRegistry tests
    result[1] = entangle_tests(verbose=verbose, useSystem=False)
    # Entanglement QGate with QSystem tests
    result[2] = entangle_tests(verbose=verbose, useSystem=True)
    for nq in range(minqubits, maxqubits + 1):
        if verbose:
            print("Testing with " + str(nq) + " qubit circuits")
        # Deutsch-Josza algorithm with QRegistry tests
        result[3] = map(add, result[3],
                        deutschTests(4, verbose=verbose, useSystem=False,
                                     optimize=False))
        # Deutsch-Josza algorithm with QSystem tests
        result[4] = map(add, result[4],
                        deutschTests(4, verbose=verbose, useSystem=True,
                                     optimize=False))
        # Deutsch-Josza algorithm with QRegistry tests optimized
        result[5] = map(add, result[5],
                        deutschTests(4, verbose=verbose, useSystem=False,
                                     optimize=True))
        # Deutsch-Josza algorithm with QSystem tests optimized
        result[6] = map(add, result[6],
                        deutschTests(4, verbose=verbose, useSystem=True,
                                     optimize=True))
    # Teleportation algorithm with QRegistry tests
    result[7] = teleportation_tests(verbose=verbose, useSystem=False,
                                    remove=False, optimize=False)
    # Teleportation algorithm with QRegistry tests and remove option
    result[8] = teleportation_tests(verbose=verbose, useSystem=False,
                                    remove=True, optimize=False)
    # Teleportation algorithm with QSystem tests
    result[9] = teleportation_tests(verbose=verbose, useSystem=True,
                                    remove=False, optimize=False)
    # Teleportation algorithm with QSystem tests and remove option
    result[10] = teleportation_tests(verbose=verbose, useSystem=True,
                                     remove=True, optimize=False)
    # Teleportation algorithm with QRegistry tests optimized
    result[11] = teleportation_tests(verbose=verbose, useSystem=False,
                                     remove=False, optimize=True)
    # Teleportation algorithm with QRegistry tests and remove option optimized
    result[12] = teleportation_tests(verbose=verbose, useSystem=False,
                                     remove=True, optimize=True)
    # Teleportation algorithm with QSystem tests optimized
    result[13] = teleportation_tests(verbose=verbose, useSystem=True,
                                     remove=False, optimize=True)
    # Teleportation algorithm with QSystem tests and remove option optimized
    result[14] = teleportation_tests(verbose=verbose, useSystem=True,
                                     remove=True, optimize=True)
    # Control and anticontrol check for QGate
    result[15] = add_line_tests(qstruct=qj.QGate, verbose=verbose)
    # Control and anticontrol check for QCircuit
    result[16] = add_line_tests(qstruct=qj.QCircuit, verbose=verbose)

    for i in range(3, 7):
        result[i] = tuple(result[i])

    return result


def main():
    """Execute all tests."""
    argv = sys.argv[1:]
    if 2 <= len(argv) <= 4 and int(argv[0]) >= 3:
        results = []
        if len(argv) == 2:
            seed = rnd.randrange(2**32 - 1)
            print("Seed: " + str(seed))
            print("\tTesting Gates...")
            results += all_gate_tests()
            print("\tTesting QRegistry...")
            results += data_structure_tests(int(argv[0]), int(argv[1]),
                                            seed=seed)
            print("\tTesting QSystem...")
            results += data_structure_tests(int(argv[0]), int(argv[1]),
                                            seed=seed, QItem=qj.QSystem)
            print("\tTesting QGate and QCircuit...")
            results += high_level_tests(int(argv[0]), int(argv[1]), seed=seed)
        elif len(argv) == 3:
            print("Seed: " + str(int(argv[2])))
            print("\tTesting Gates...")
            results += all_gate_tests(seed=int(argv[2]))
            print("\tTesting QRegistry...")
            results += data_structure_tests(int(argv[0]), int(argv[1]),
                                            seed=int(argv[2]))
            print("\tTesting QSystem...")
            results += data_structure_tests(int(argv[0]), int(argv[1]),
                                            seed=int(argv[2]),
                                            QItem=qj.QSystem)
            print("\tTesting QGate and QCircuit...")
            results += high_level_tests(int(argv[0]), int(argv[1]),
                                        seed=int(argv[2]))
        else:
            print("Seed: " + str(int(argv[2])))
            print("\tTesting Gates...")
            results += all_gate_tests(seed=int(argv[2]), verbose=bool(argv[3]))
            print("\tTesting QRegistry...")
            results += data_structure_tests(int(argv[0]), int(argv[1]),
                                            seed=int(argv[2]),
                                            verbose=bool(argv[3]))
            print("\tTesting QSystem...")
            results += data_structure_tests(int(argv[0]), int(argv[1]),
                                            seed=int(argv[2]),
                                            verbose=bool(argv[3]),
                                            QItem=qj.QSystem)
            print("\tTesting QGate and QCircuit...")
            results += high_level_tests(int(argv[0]), int(argv[1]),
                                        seed=int(argv[2]),
                                        verbose=bool(argv[3]))
        passed = [int(result[0] == result[1]) for result in results]
        noice = sum(passed)
        total = len(results)
        print("Passed: " + str(noice) + "/" + str(total))
        if noice == total:
            print("PEACE AND TRANQUILITY")
            # wb.open_new_tab("https://youtu.be/SHvhps47Lmc")
        else:
            for testid in range(total):
                if passed[testid] == 0:
                    if testid < 15:
                        print("Gate Test", str(testid), "failed!")
                    elif testid < 15 + 37:
                        print("QRegistry Test", str(testid - 15), "failed!")
                    elif testid < 15 + 37 + 39:
                        print("QSystem Test", str(testid - (15 + 37)),
                              "failed!")
                    else:
                        print("High level Test", str(testid - (15 + 37 + 39)),
                              "failed!")
            print("SORROW")
            # wb.open_new_tab("https://youtu.be/4Js-XbNj6Tk?t=37")
        # We assert so the test fails if we failed something
        assert noice == total
    else:
        print("Syntax: " + sys.argv[0] + " <minimum number of qubits (min 3)>",
              "<maximum number of qubits> <seed (optional)>")


if __name__ == "__main__":
    main()
