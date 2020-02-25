#!/usr/bin/python

import sys
import qsimov as qj
import webbrowser as wb
import random as rnd
import numpy as np
from operator import add
import math as m

def gateToId(numpygate, id, nq):
    nqg = int(np.log2(numpygate.shape[0]))
    if 0 < id < nq - nqg:
        a = np.kron(np.eye(2**(nq - id - nqg)), np.kron(numpygate, np.eye(2**id)))
    elif id == 0:
        a = np.kron(np.eye(2**(nq - nqg)), numpygate)
    else:
        a = np.kron(numpygate, np.eye(2**id))
    #print("    NpGate: " + str(a))
    return a[:,0]

def gateToRegId(numpygate, id, nq, reg):
    nqg = int(np.log2(numpygate.shape[0]))
    if 0 < id < nq - nqg:
        a = np.kron(np.eye(2**(nq - id - nqg)), np.kron(numpygate, np.eye(2**id)))
    elif id == 0:
        a = np.kron(np.eye(2**(nq - nqg)), numpygate)
    else:
        a = np.kron(numpygate, np.eye(2**id))
    #print("    NpGate: " + str(a))
    return np.dot(a, reg.reshape(reg.size, 1)).reshape(reg.size)

def H():
    gate = np.ones(4, dtype=complex).reshape(2, 2)
    aux = 1/np.sqrt(2)
    gate[1,1] = -1
    return gate * aux

def X():
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    gate[0,1] = 1
    gate[1,0] = 1
    return gate

def Y():
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    gate[0,1] = -1j
    gate[1,0] = 1j
    return gate

def Z():
    gate = np.eye(2, dtype=complex)
    gate[1,1] = -1
    return gate

def SqrtX(invert):
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    if not invert:
        gate[0,0] = 0.5 + 0.5j
        gate[0,1] = 0.5 - 0.5j
        gate[1,0] = 0.5 - 0.5j
        gate[1,1] = 0.5 + 0.5j
    else:
        gate[0,0] = 0.5 - 0.5j
        gate[0,1] = 0.5 + 0.5j
        gate[1,0] = 0.5 + 0.5j
        gate[1,1] = 0.5 - 0.5j
    return gate

def R(angle, invert):
    gate = np.eye(2, dtype=complex)
    if not invert:
        gate[1,1] = np.cos(angle) + np.sin(angle) * 1j
    else:
        gate[1,1] = np.cos(angle) - np.sin(angle) * 1j
    return gate

def RX(angle, invert):
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    cosan = np.cos(angle/2)
    sinan = np.sin(angle/2)
    if not invert:
        mult = -1j
    else:
        mult = 1j
    gate[0,0] = cosan
    gate[0,1] = mult * sinan
    gate[1,0] = mult * sinan
    gate[1,1] = cosan
    return gate

def RY(angle, invert):
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    cosan = np.cos(angle/2)
    sinan = np.sin(angle/2)
    gate[0,0] = cosan
    gate[1,1] = cosan
    if not invert:
        gate[0,1] = -sinan
        gate[1,0] = sinan
    else:
        gate[0,1] = sinan
        gate[1,0] = -sinan
    return gate

def RZ(angle, invert):
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    if not invert:
        gate[0,0] = np.cos(-angle/2) + np.sin(-angle/2) * 1j
        gate[1,1] = np.cos(angle/2) + np.sin(angle/2) * 1j
    else:
        gate[0,0] = np.cos(-angle/2) - np.sin(-angle/2) * 1j
        gate[1,1] = np.cos(angle/2) - np.sin(angle/2) * 1j
    return gate

def RUnity(n, invert):
    return R(2*np.pi/(2**n), invert)

def HalfDeutsch(angle, invert):
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    cosan = np.cos(angle) * 1j
    sinan = np.sin(angle)
    if invert:
        cosan = -cosan
    gate[0,0] = cosan
    gate[0,1] = sinan
    gate[1,0] = sinan
    gate[1,1] = cosan
    return gate

def U(angle1, angle2, angle3, invert):
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    cosan = np.cos(angle1/2)
    sinan = np.sin(angle1/2)
    mult = 1
    if invert:
        mult = -1
    gate[0,0] = cosan
    if not invert:
        gate[0,1] = -sinan * np.cos(angle3) - sinan * np.sin(angle3) * 1j
        gate[1,0] = sinan * np.cos(angle2) + sinan * np.sin(angle2) * 1j
    else:
        gate[0,1] = sinan * np.cos(angle2) - sinan * np.sin(angle2) * 1j
        gate[1,0] = -sinan * np.cos(angle3) + sinan * np.sin(angle3) * 1j
    gate[1,1] = cosan * np.cos(angle2+angle3) + mult * cosan * np.sin(angle2+angle3) * 1j
    return gate

U3 = U

def U2(angle1, angle2, invert):
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    mult = 1
    if invert:
        mult = -1
    gate[0,0] = 1
    if not invert:
        gate[0,1] = -np.cos(angle2) - np.sin(angle2) * 1j
        gate[1,0] = np.cos(angle1) + np.sin(angle1) * 1j
    else:
        gate[0,1] = -np.cos(angle2) + np.sin(angle2) * 1j
        gate[1,0] = np.cos(angle1) - np.sin(angle1) * 1j
    gate[1,1] = np.cos(angle1+angle2) + mult * np.sin(angle1+angle2) * 1j
    return gate * (1/np.sqrt(2))

def U1(angle, invert):
    return R(angle, invert)

def SWAP():
    gate = np.zeros(4 * 4, dtype=complex)
    gate = gate.reshape(4, 4)

    gate[0][0] = 1
    gate[1][2] = 1
    gate[2][1] = 1
    gate[3][3] = 1

    return gate

def iSWAP(invert):
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
    gate = np.eye(4, dtype=complex)
    if not invert:
        gate[0,3] = np.sin(angle)-np.cos(angle)*1j
        gate[1,2] = -1j
        gate[2,1] = -1j
        gate[3,0] = np.sin(-angle)-np.cos(-angle)*1j
    else:
        gate[0,3] = np.sin(-angle)+np.cos(-angle)*1j
        gate[1,2] = 1j
        gate[2,1] = 1j
        gate[3,0] = np.sin(angle)+np.cos(angle)*1j
    return gate*(1/np.sqrt(2))

def yy(angle, invert):
    gate = np.eye(4, dtype=complex)
    gate = gate * np.cos(angle)
    ansin = np.sin(angle) * 1j
    if not invert:
        gate[0,3] = ansin
        gate[1,2] = -ansin
        gate[2,1] = -ansin
        gate[3,0] = ansin
    else:
        gate[0,3] = -ansin
        gate[1,2] = ansin
        gate[2,1] = ansin
        gate[3,0] = -ansin
    return gate

def zz(angle, invert):
    gate = np.eye(4, dtype=complex)
    gate = gate * np.cos(angle)
    phi2 = angle/2
    if not invert:
        gate[0,0] = np.cos(phi2) + np.sin(phi2) * 1j
        gate[1,1] = np.cos(-phi2) + np.sin(-phi2) * 1j
        gate[2,2] = gate[1,1]
        gate[3,3] = gate[0,0]
    else:
        gate[0,0] = np.cos(phi2) - np.sin(phi2) * 1j
        gate[1,1] = np.cos(-phi2) - np.sin(-phi2) * 1j
        gate[2,2] = gate[1,1]
        gate[3,3] = gate[0,0]
    return gate

def swapDownstairs(id1, id2, nq, reg):
    swap = SWAP()
    for i in range(id1, id2 + 1):
        reg = gateToRegId(swap, i, nq, reg)
    return reg

def swapUpstairs(id1, id2, nq, reg):
    swap = SWAP()
    for i in range(id2, id1 - 1, -1):
        reg = gateToRegId(swap, i, nq, reg)
    return reg

def swapDownstairsList(id1, id2, li):
    for i in range(id1, id2 + 1):
        li[i], li[i+1] = li[i+1], li[i]
    return li

def swapUpstairsList(id1, id2, li):
    for i in range(id2, id1 - 1, -1):
        li[i], li[i+1] = li[i+1], li[i]
    return li

def sparseTwoGate(gate, id1, id2, nq, reg):
    if id2 < id1:
        id1, id2 = id2, id1
    if id1 < 0 or id2 >= nq:
        reg = None
    else:
        if id2 - id1 > 1:
            reg = swapDownstairs(id1, id2 - 1, nq, reg)
        reg = gateToRegId(gate, id2 - 1, nq, reg)
        if id2 - id1 > 1:
            reg = swapUpstairs(id1, id2 - 1, nq, reg)
    return reg

def CU(gate, ncontrols):
    nqgate = int(np.log2(gate.shape[0]))
    cu = np.eye(2**(nqgate+ncontrols), dtype=complex)
    aux = cu.shape[0] - gate.shape[0]
    for i in range(gate.shape[0]):
        for j in range(gate.shape[1]):
            cu[aux + i, aux + j] = gate[i,j]

    return cu

def negateQubits(qubits, nq, reg):
    for id in qubits:
        reg = gateToRegId(X(), id, nq, reg)
    return reg

def applyCACU(gate, id, controls, anticontrols, nq, reg): # Controlled and Anti-Controlled U
    cset = set(controls)
    acset = set(anticontrols)
    cuac = list(cset.union(acset))
    if type(id) == list:
        extended_cuac = id + cuac
    else:
        extended_cuac = [id] + cuac
    qubitIds = [i for i in range(nq)]
    first = cuac[0]

    reg = negateQubits(acset, nq, reg)
    #print("Ids: " + str(qubitIds))
    for i in range(len(extended_cuac)):
        if qubitIds[i] != extended_cuac[i]:
            indaux = qubitIds.index(extended_cuac[i])
            reg = swapUpstairs(i, indaux - 1, nq, reg)
            qubitIds = swapUpstairsList(i, indaux - 1, qubitIds)
            #print("Ids: " + str(qubitIds))
    #print(reg)
    reg = gateToRegId(CU(gate, len(cuac)), 0, nq, reg)
    #print(reg)
    for i in range(nq):
        if qubitIds[i] != i:
            indaux = qubitIds.index(i)
            reg = swapUpstairs(i, indaux - 1, nq, reg)
            qubitIds = swapUpstairsList(i, indaux - 1, qubitIds)
            #print("Ids: " + str(qubitIds))
    reg = negateQubits(acset, nq, reg)
    #print(reg)
    return reg

def gateTests(gatename, verbose=False, hasInv=False, nArgs=0):
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
        numpygate = qj.getGate(gatename)
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
        numpygate = qj.getGate(gatename)
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
        numpygate = qj.getGate(gatename + "-1")
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
            numpygate = qj.getGate(gatename + argstr)
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
            numpygate = qj.getGate(gatename + argstr + "-1")
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

def oneGateTests(nq, verbose=False, QItem=qj.QRegistry):
    passed = 0
    total = nq
    if verbose:
        print(" One gate tests:")
    for id in range(nq):
        gate = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," + str(rnd.random()) + ")"
        numpygate = qj.getGate(gate)
        if verbose:
            print("  Testing gate " + gate + " to qubit " + str(id) + "...")
        #print("    Gate: " + str(numpygate))
        b = QItem(nq)
        a = gateToId(numpygate, id, nq)
        if verbose:
            print("    Numpy done")
        b.applyGate(gate, id)
        if verbose:
            print("    QSimov done")
        allOk = np.allclose(a, b.getState())
        if not allOk:
            if verbose:
                print(a)
                print(b.getState())
                print(a == b.getState())
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

def simpleSwapTests(nq, imaginary, sqrt, invert, verbose=False, QItem=qj.QRegistry):
    passed = 0
    total = nq - 1
    gate2 = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," + str(rnd.random()) + ")"
    numpygate2 = qj.getGate(gate2)
    a = gateToId(numpygate2, 0, nq)
    b = QItem(nq)
    b.applyGate(gate2)
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
        # print("    Gate: " + str(numpygate2))
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
            print("   " + verbstring + gate + "ping qubits " + str(id) + " and " + str(id + 1) + "...")
        a = gateToRegId(numpygate, id, nq, a)
        if verbose:
            print("    Numpy done")
        b.applyGate(gate + "(" + str(id) + "," + str(id + 1) + ")" + invstring)
        if verbose:
            print("    QSimov done")
        allOk = np.allclose(a, b.getState())
        if not allOk:
            if verbose:
                print(a)
                print(b.getState())
                print(a == b.getState())
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

def sparseSwapTests(nq, imaginary, sqrt, invert, verbose=False, QItem=qj.QRegistry):
    passed = 0
    total = ((nq - 1) * (nq - 2))//2
    # for i in range(2, nq):
    #     total += (nq - i)
    gate2 = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," + str(rnd.random()) + ")"
    numpygate2 = qj.getGate(gate2)
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
            a = gateToId(numpygate2, id1, nq)
            b = QItem(nq)
            b.applyGate(gate2, id1)
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
                print("   " + verbstring + gate + "ping qubits " + str(id1) + " and " + str(id2) + "...")
            a = sparseTwoGate(numpygate, id1, id2, nq, a)
            if verbose:
                print("    Numpy done")
            b.applyGate(gate + "(" + str(id1) + "," + str(id2) + ")" + invstring)
            if verbose:
                print("    QSimov done")
            allOk = np.allclose(a, b.getState())
            if not allOk:
                if verbose:
                    print(a)
                    print(b.getState())
                    print(a == b.getState())
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

def simpleIsingTests(nq, type, invert, verbose=False, QItem=qj.QRegistry):
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

    gate2 = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," + str(rnd.random()) + ")"
    numpygate2 = qj.getGate(gate2)
    a = gateToId(numpygate2, 0, nq)
    b = QItem(nq)
    b.applyGate(gate2)
    if verbose:
        print("   Applied gate " + gate2 + " to qubit 0")
        # print("    Gate: " + str(numpygate2))
    for id in range(nq - 1):
        angle = rnd.random()
        numpygate = gatefun(angle, invert)
        if verbose:
            print("   Coupling qubits " + str(id) + " and " + str(id + 1) + " with angle " + str(angle) + "...")
            #print("    Gate: " + str(numpygate))
        a = gateToRegId(numpygate, id, nq, a)
        if verbose:
            print("    Numpy done")
        b.applyGate(gate + "(" + str(angle) + "," + str(id) + "," + str(id + 1) + ")" + invstring)
        if verbose:
            print("    QSimov done")
        allOk = np.allclose(a, b.getState())
        if not allOk:
            if verbose:
                print(a)
                print(b.getState())
                print(a == b.getState())
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

def sparseIsingTests(nq, type, invert, verbose=False, QItem=qj.QRegistry):
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
            gate2 = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," + str(rnd.random()) + ")"
            numpygate2 = qj.getGate(gate2)
            angle = rnd.random()
            numpygate = gatefun(angle, invert)
            a = gateToId(numpygate2, id1, nq)
            b = QItem(nq)
            b.applyGate(gate2, id1)
            if verbose:
                print("   Applied gate " + gate2 + " to qubit 0")
                # print("    Gate: " + str(numpygate2))
            if verbose:
                print("   Coupling qubits " + str(id1) + " and " + str(id2) + " with angle " + str(angle) + "...")
            a = sparseTwoGate(numpygate, id1, id2, nq, a)
            if verbose:
                print("    Numpy done")
            b.applyGate(gate + "(" + str(angle) + "," + str(id1) + "," + str(id2) + ")" + invstring)
            if verbose:
                print("    QSimov done")
            allOk = np.allclose(a, b.getState())
            if not allOk:
                if verbose:
                    print(a)
                    print(b.getState())
                    print(a == b.getState())
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

# def applyCACU(gate, id, controls, anticontrols, nq, reg)
# cascada aleatoria
def controlledGateTests(nq, verbose=False, QItem=qj.QRegistry):
    total = nq - 1
    passed = 0
    isControl = bool(rnd.randint(0, 1))
    qubitIds = np.random.permutation(nq)
    lastid = qubitIds[0]
    control = []
    anticontrol = []
    gate = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," + str(rnd.random()) + ")"
    if verbose:
        print(" Controlled gate tests:")
        print("  Gate: " + gate + " to qubit " + str(lastid))
    numpygate = qj.getGate(gate)
    a = gateToId(numpygate, lastid, nq)
    b = QItem(nq)
    b.applyGate(gate, lastid)
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
        b.applyGate(gate, qubit=id, control=control, anticontrol=anticontrol)
        isControl = not isControl
        lastid = id
        allOk = np.allclose(a, b.getState())
        if not allOk:
            if verbose:
                print(a)
                print(b.getState())
                print(a == b.getState())
                print("    Michael Bay visited your simulator...")
            break
        passed += 1
        if verbose:
            print("    Noice")
    del a
    del b
    return (passed, total)

def cSWAPTests(nq, imaginary, sqrt, invert, verbose=False, QItem=qj.QRegistry):
    total = nq - 1
    passed = 0
    isControl = bool(rnd.randint(0, 1))
    qubitIds = np.random.permutation(nq)
    lastid = qubitIds[:2]
    control = []
    anticontrol = []
    gate = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," + str(rnd.random()) + ")"
    if verbose:
        print(" Controlled gate tests:")
        print("  Gate: " + gate + " to qubit " + str(lastid))
    numpygate = qj.getGate(gate)
    a = gateToId(numpygate, lastid[0], nq)
    b = QItem(nq)
    b.applyGate(gate, lastid[0])

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
    b.applyGate(gate + "(" + str(lastid[0]) + "," + str(lastid[1]) + ")" + invstring)
    allOk = np.allclose(a, b.getState())
    if not allOk:
        if verbose:
            print(a)
            print(b.getState())
            print(a == b.getState())
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
            b.applyGate(gate + "(" + str(lastid[0]) + "," + str(lastid[1]) + ")" + invstring, control=control, anticontrol=anticontrol)
            allOk = np.allclose(a, b.getState())
            if not allOk:
                if verbose:
                    print(a)
                    print(b.getState())
                    print(a == b.getState())
                    print("    Michael Bay visited your simulator...")
                break
            passed += 1
            if verbose:
                print("    Noice")
    del a
    del b
    return (passed, total)

def cIsingTests(nq, type, invert, verbose=False, QItem=qj.QRegistry):
    total = nq - 1
    passed = 0
    isControl = bool(rnd.randint(0, 1))
    qubitIds = np.random.permutation(nq)
    lastid = qubitIds[:2]
    control = []
    anticontrol = []
    gate = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," + str(rnd.random()) + ")"
    if verbose:
        print(" Controlled gate tests:")
        print("  Gate: " + gate + " to qubit " + str(lastid))
    numpygate = qj.getGate(gate)
    a = gateToId(numpygate, lastid[0], nq)
    b = QItem(nq)
    b.applyGate(gate, lastid[0])

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
    b.applyGate(gate + "(" + str(angle) + "," + str(lastid[0]) + "," + str(lastid[1]) + ")" + invstring)
    allOk = np.allclose(a, b.getState())
    if not allOk:
        if verbose:
            print(a)
            print(b.getState())
            print(a == b.getState())
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
            a = applyCACU(gatefun(angle, invert), lastid, control, anticontrol, nq, a)
            b.applyGate(gate + "(" + str(angle) + "," + str(lastid[0]) + "," + str(lastid[1]) + ")" + invstring, control=control, anticontrol=anticontrol)
            allOk = np.allclose(a, b.getState())
            if not allOk:
                if verbose:
                    print(a)
                    print(b.getState())
                    print(a == b.getState())
                    print("    Michael Bay visited your simulator...")
                break
            passed += 1
            if verbose:
                print("    Noice")
    del a
    del b
    return (passed, total)

def allGateTests(seed=None, verbose=False):
    if not (seed is None):
        qj.setRandomSeed(seed)
        rnd.seed(seed)
        np.random.seed(seed)
    result = [(0, 0) for i in range(15)] # We have 15 tests

    result[0] = gateTests("H", verbose=verbose, hasInv=False, nArgs=0) # H gate tests
    result[1] = gateTests("X", verbose=verbose, hasInv=False, nArgs=0) # X gate tests
    result[2] = gateTests("Y", verbose=verbose, hasInv=False, nArgs=0) # Y gate tests
    result[3] = gateTests("Z", verbose=verbose, hasInv=False, nArgs=0) # Z gate tests
    result[4] = gateTests("SqrtX", verbose=verbose, hasInv=True, nArgs=0) # SqrtX gate tests
    result[5] = gateTests("RX", verbose=verbose, hasInv=True, nArgs=1) # RX gate tests
    result[6] = gateTests("RY", verbose=verbose, hasInv=True, nArgs=1) # RY gate tests
    result[7] = gateTests("RZ", verbose=verbose, hasInv=True, nArgs=1) # RZ gate tests
    result[8] = gateTests("R", verbose=verbose, hasInv=True, nArgs=1) # Phase shift gate tests
    result[9] = gateTests("RUnity", verbose=verbose, hasInv=True, nArgs=1) # Roots of unity gate tests
    result[10] = gateTests("HalfDeutsch", verbose=verbose, hasInv=True, nArgs=1) # Partial Deutsch gate tests
    result[11] = gateTests("U", verbose=verbose, hasInv=True, nArgs=3) # U gate tests
    result[12] = gateTests("U3", verbose=verbose, hasInv=True, nArgs=3) # U3 gate tests
    result[13] = gateTests("U2", verbose=verbose, hasInv=True, nArgs=2) # U2 gate tests
    result[14] = gateTests("U1", verbose=verbose, hasInv=True, nArgs=1) # U1 gate tests

    return result

def dataStructureTests(minqubits, maxqubits, seed=None, verbose=False, QItem=qj.QRegistry):
    if not (seed is None):
        qj.setRandomSeed(seed)
        rnd.seed(seed)
        np.random.seed(seed)
    result = [(0, 0) for i in range(35)] # We have 35 tests

    for nq in range(minqubits, maxqubits + 1):
        if verbose:
            print("Testing with " + str(nq) + " qubit registries")
        result[0] = map(add, result[0], oneGateTests(nq, verbose=verbose, QItem=QItem)) # Apply U gate tests
        result[1] = map(add, result[1], simpleSwapTests(nq, False, False, False, verbose=verbose, QItem=QItem)) # Simple SWAP gate tests
        result[2] = map(add, result[2], simpleSwapTests(nq, True, False, False, verbose=verbose, QItem=QItem)) # Simple iSWAP gate tests
        result[3] = map(add, result[3], simpleSwapTests(nq, False, True, False, verbose=verbose, QItem=QItem)) # Simple sqrtSWAP gate tests
        result[4] = map(add, result[4], simpleSwapTests(nq, True, False, True, verbose=verbose, QItem=QItem)) # Simple iSWAP-1 gate tests
        result[5] = map(add, result[5], simpleSwapTests(nq, False, True, True, verbose=verbose, QItem=QItem)) # Simple sqrtSWAP-1 gate tests
        result[6] = map(add, result[6], sparseSwapTests(nq, False, False, False, verbose=verbose, QItem=QItem)) # Sparse SWAP gate tests
        result[7] = map(add, result[7], sparseSwapTests(nq, True, False, False, verbose=verbose, QItem=QItem)) # Sparse iSWAP gate tests
        result[8] = map(add, result[8], sparseSwapTests(nq, False, True, False, verbose=verbose, QItem=QItem)) # Sparse sqrtSWAP gate tests
        result[9] = map(add, result[9], sparseSwapTests(nq, True, False, True, verbose=verbose, QItem=QItem)) # Sparse iSWAP-1 gate tests
        result[10] = map(add, result[10], sparseSwapTests(nq, False, True, True, verbose=verbose, QItem=QItem)) # Sparse sqrtSWAP-1 gate tests
        result[11] = map(add, result[11], simpleIsingTests(nq, 0, False, verbose=verbose, QItem=QItem)) # Simple XX gate tests
        result[12] = map(add, result[12], simpleIsingTests(nq, 1, False, verbose=verbose, QItem=QItem)) # Simple YY gate tests
        result[13] = map(add, result[13], simpleIsingTests(nq, 2, False, verbose=verbose, QItem=QItem)) # Simple ZZ gate tests
        result[14] = map(add, result[14], simpleIsingTests(nq, 0, True, verbose=verbose, QItem=QItem)) # Simple XX-1 gate tests
        result[15] = map(add, result[15], simpleIsingTests(nq, 1, True, verbose=verbose, QItem=QItem)) # Simple YY-1 gate tests
        result[16] = map(add, result[16], simpleIsingTests(nq, 2, True, verbose=verbose, QItem=QItem)) # Simple ZZ-1 gate tests
        result[17] = map(add, result[17], sparseIsingTests(nq, 0, False, verbose=verbose, QItem=QItem)) # Sparse XX gate tests
        result[18] = map(add, result[18], sparseIsingTests(nq, 1, False, verbose=verbose, QItem=QItem)) # Sparse YY gate tests
        result[19] = map(add, result[19], sparseIsingTests(nq, 2, False, verbose=verbose, QItem=QItem)) # Sparse ZZ gate tests
        result[20] = map(add, result[20], sparseIsingTests(nq, 0, True, verbose=verbose, QItem=QItem)) # Sparse XX-1 gate tests
        result[21] = map(add, result[21], sparseIsingTests(nq, 1, True, verbose=verbose, QItem=QItem)) # Sparse YY-1 gate tests
        result[22] = map(add, result[22], sparseIsingTests(nq, 2, True, verbose=verbose, QItem=QItem)) # Sparse ZZ-1 gate tests
        result[23] = map(add, result[23], controlledGateTests(nq, verbose=verbose, QItem=QItem)) # Controlled U gate tests
        result[24] = map(add, result[24], cSWAPTests(nq, False, False, False, verbose=verbose, QItem=QItem)) # Controlled SWAP gate tests
        result[25] = map(add, result[25], cSWAPTests(nq, True, False, False, verbose=verbose, QItem=QItem)) # Controlled iSWAP gate tests
        result[26] = map(add, result[26], cSWAPTests(nq, False, True, False, verbose=verbose, QItem=QItem)) # Controlled sqrtSWAP gate tests
        result[27] = map(add, result[27], cSWAPTests(nq, True, False, True, verbose=verbose, QItem=QItem)) # Controlled iSWAP-1 gate tests
        result[28] = map(add, result[28], cSWAPTests(nq, False, True, True, verbose=verbose, QItem=QItem)) # Controlled sqrtSWAP-1 gate tests
        result[29] = map(add, result[29], cIsingTests(nq, 0, False, verbose=verbose, QItem=QItem)) # Controlled XX gate tests
        result[30] = map(add, result[30], cIsingTests(nq, 1, False, verbose=verbose, QItem=QItem)) # Controlled YY gate tests
        result[31] = map(add, result[31], cIsingTests(nq, 2, False, verbose=verbose, QItem=QItem)) # Controlled ZZ gate tests
        result[32] = map(add, result[32], cIsingTests(nq, 0, True, verbose=verbose, QItem=QItem)) # Controlled XX-1 gate tests
        result[33] = map(add, result[33], cIsingTests(nq, 1, True, verbose=verbose, QItem=QItem)) # Controlled YY-1 gate tests
        result[34] = map(add, result[34], cIsingTests(nq, 2, True, verbose=verbose, QItem=QItem)) # Controlled ZZ-1 gate tests

    for i in range(len(result)):
        result[i] = tuple(result[i])

    return result

def main():
    argv = sys.argv[1:]
    if 2 <= len(argv) <= 4 and int(argv[0]) >= 3:
        results = []
        if len(argv) == 2:
            seed = rnd.randrange(2**32 - 1)
            print("Seed: " + str(seed))
            print("\tTesting Gates...")
            results += allGateTests()
            print("\tTesting QRegistry...")
            results += dataStructureTests(int(argv[0]), int(argv[1]), seed=seed)
            print("\tTesting QSystem...")
            results += dataStructureTests(int(argv[0]), int(argv[1]), seed=seed, QItem=qj.QSystem)
        elif len(argv) == 3:
            print("Seed: " + str(int(argv[2])))
            print("\tTesting Gates...")
            results += allGateTests(seed=int(argv[2]))
            print("\tTesting QRegistry...")
            results += dataStructureTests(int(argv[0]), int(argv[1]), seed=int(argv[2]))
            print("\tTesting QSystem...")
            results += dataStructureTests(int(argv[0]), int(argv[1]), seed=int(argv[2]), QItem=qj.QSystem)
        else:
            print("Seed: " + str(int(argv[2])))
            print("\tTesting Gates...")
            results += allGateTests(seed=int(argv[2]), verbose=bool(argv[3]))
            print("\tTesting QRegistry...")
            results += dataStructureTests(int(argv[0]), int(argv[1]), seed=int(argv[2]), verbose=bool(argv[3]))
            print("\tTesting QSystem...")
            results += dataStructureTests(int(argv[0]), int(argv[1]), seed=int(argv[2]), verbose=bool(argv[3]), QItem=qj.QSystem)
        passed = [int(result[0] == result[1]) for result in results]
        noice = sum(passed)
        total = len(results)
        print("Passed: " + str(noice) + "/" + str(total))
        if noice == total:
            print("PEACE AND TRANQUILITY")
            #wb.open_new_tab("https://youtu.be/SHvhps47Lmc")
        else:
            for testid in range(total):
                if passed[testid] == 0:
                    if testid < 15:
                        print("Gate Test " + str(testid) + " failed!")
                    elif testid < 50:
                        print("QRegistry Test " + str(testid - 15) + " failed!")
                    else:
                        print("QRegistry Test " + str(testid - 50) + " failed!")
            print("SORROW")
            #wb.open_new_tab("https://youtu.be/4Js-XbNj6Tk?t=37")
        assert noice == total
    else:
        print ("Syntax: " + sys.argv[0] + " <minimum number of qubits (min 3)> <maximum number of qubits> <seed (optional)>")

if __name__ == "__main__":
    main()
