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

def oneGateTests(nq, verbose=False):
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
        b = qj.QRegistry(nq)
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

def simpleSwapTests(nq, imaginary, sqrt, invert, verbose=False):
    passed = 0
    total = nq - 1
    gate2 = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," + str(rnd.random()) + ")"
    numpygate2 = qj.getGate(gate2)
    a = gateToId(numpygate2, 0, nq)
    b = qj.QRegistry(nq)
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

def sparseSwapTests(nq, imaginary, sqrt, invert, verbose=False):
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
            b = qj.QRegistry(nq)
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

'''
def swapTests(nq, imaginary, verbose=False):
    passed = 0
    total = 2 * (nq - 1)
    if verbose:
        if not imaginary:
            print(" SWAP tests:")
            print("  SWAP gate:")
        else:
            print("  ISWAP gate:")
    for id in range(1, nq):
        gate = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," + str(rnd.random()) + ")"
        numpygate = qj.getGate(gate)
        if verbose:
            print("   Testing gate " + gate + " to qubit " + str(id) + "...")
            # print("    Gate: " + str(numpygate))
        b = qj.QRegistry(nq)
        b.applyGate(gate)
        if not imaginary:
            a = gateToId(numpygate, id, nq)
            if verbose:
                print("    Numpy done")
            b.applyGate("SWAP(0, " + str(id) + ")")
        else:
            a = gateToId(numpygate, id, nq)
            a = np.array([a[i] if i&1 == (i>>id)&1 else a[i] * 1j for i in range(a.size)])
            if verbose:
                print("    Numpy done")
            b.applyGate("ISWAP(0, " + str(id) + ")")
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
            print("    Noice, reverting")
        a = gateToId(numpygate, 0, nq)
        if verbose:
            print("    Numpy done")
        if not imaginary:
            b.applyGate("SWAP(0, " + str(id) + ")")
        else:
            b.applyGate("ISWAP(0, " + str(id) + ")-1")
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

def sswapTests(nq, verbose=False):
    passed = 0
    total = 2 * (nq - 1)
    if verbose:
        print("  SqrtSWAP gate:")
    for id in range(1, nq):
        gate = "U(" + str(rnd.random()) + "," + str(rnd.random()) + "," + str(rnd.random()) + ")"
        numpygate = qj.getGate(gate)
        if verbose:
            print("   Testing gate " + gate + " to qubit " + str(id) + "...")
            # print("    Gate: " + str(numpygate))
        b = qj.QRegistry(nq)
        b.applyGate(gate)
        a = gateToId(numpygate, id, nq)
        if verbose:
            print("    Numpy done")
        b.applyGate("SqrtSWAP(0, " + str(id) + ")")
        b.applyGate("SqrtSWAP(0, " + str(id) + ")")
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
        b.applyGate("SWAP(0, " + str(id) + ")")
        if verbose:
            print("    Noice, reverting")
        a = gateToId(numpygate, 0, nq)
        if verbose:
            print("    Numpy done")
        b.applyGate("SqrtSWAP(0, " + str(id) + ")")
        b.applyGate("SqrtSWAP(0, " + str(id) + ")-1")
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

def isingTests(nq, verbose=False):
    passed = 0
    total = nq - 1
    for id in range(nq - 1):
        ranin = rnd.randint(0, 2)
        angle = rnd.random()
        gatename = "XX"
        numpygate = xx(angle, False)
        if ranin == 1:
            gatename = "YY"
            numpygate = yy(angle, False)
        elif ranin == 2:
            gatename = "ZZ"
            numpygate = zz(angle, False)
        gate = gatename + "(" + str(angle) + "," + str(id) + "," + str(id+1) + ")"
        # numpygate = qj.getGate(gate)
        if verbose:
            print("   Testing gate " + gate + " to qubit " + str(id) + "...")
            # print("   Gate: " + str(numpygate))
        b = qj.QRegistry(nq)
        a = gateToId(numpygate, id, nq)
        b.applyGate(gate)
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
'''

def simpleIsingTests(nq, type, invert, verbose=False):
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
    b = qj.QRegistry(nq)
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

def sparseIsingTests(nq, type, invert, verbose=False):
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
            b = qj.QRegistry(nq)
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

def tests(minqubits, maxqubits, seed=None, verbose=False):
    if not (seed is None):
        print("Seed: " + str(seed))
        qj.setRandomSeed(seed)
        rnd.seed(seed)
    result = (0, 0)
    for nq in range(minqubits, maxqubits + 1):
        if verbose:
            print("Testing with " + str(nq) + " qubit registries")
        result = map(add, result, gateTests("H", verbose=verbose, hasInv=False, nArgs=0)) # H gate tests
        result = map(add, result, gateTests("X", verbose=verbose, hasInv=False, nArgs=0)) # X gate tests
        result = map(add, result, gateTests("Y", verbose=verbose, hasInv=False, nArgs=0)) # Y gate tests
        result = map(add, result, gateTests("Z", verbose=verbose, hasInv=False, nArgs=0)) # Z gate tests
        result = map(add, result, gateTests("SqrtX", verbose=verbose, hasInv=True, nArgs=0)) # SqrtX gate tests
        result = map(add, result, gateTests("RX", verbose=verbose, hasInv=True, nArgs=1)) # RX gate tests
        result = map(add, result, gateTests("RY", verbose=verbose, hasInv=True, nArgs=1)) # RY gate tests
        result = map(add, result, gateTests("RZ", verbose=verbose, hasInv=True, nArgs=1)) # RZ gate tests
        result = map(add, result, gateTests("R", verbose=verbose, hasInv=True, nArgs=1)) # Phase shift gate tests
        result = map(add, result, gateTests("RUnity", verbose=verbose, hasInv=True, nArgs=1)) # Roots of unity gate tests
        result = map(add, result, gateTests("HalfDeutsch", verbose=verbose, hasInv=True, nArgs=1)) # Partial Deutsch gate tests
        result = map(add, result, gateTests("U", verbose=verbose, hasInv=True, nArgs=3)) # U gate tests
        result = map(add, result, gateTests("U3", verbose=verbose, hasInv=True, nArgs=3)) # U3 gate tests
        result = map(add, result, gateTests("U2", verbose=verbose, hasInv=True, nArgs=2)) # U2 gate tests
        result = map(add, result, gateTests("U1", verbose=verbose, hasInv=True, nArgs=1)) # U1 gate tests
        result = map(add, result, oneGateTests(nq, verbose=verbose)) # Apply U gate tests
        result = map(add, result, simpleSwapTests(nq, False, False, False, verbose=verbose)) # Simple SWAP gate tests
        result = map(add, result, simpleSwapTests(nq, True, False, False, verbose=verbose)) # Simple iSWAP gate tests
        result = map(add, result, simpleSwapTests(nq, False, True, False, verbose=verbose)) # Simple sqrtSWAP gate tests
        result = map(add, result, simpleSwapTests(nq, True, False, True, verbose=verbose)) # Simple iSWAP-1 gate tests
        result = map(add, result, simpleSwapTests(nq, False, True, True, verbose=verbose)) # Simple sqrtSWAP-1 gate tests
        result = map(add, result, sparseSwapTests(nq, False, False, False, verbose=verbose)) # Sparse SWAP gate tests
        result = map(add, result, sparseSwapTests(nq, True, False, False, verbose=verbose)) # Sparse iSWAP gate tests
        result = map(add, result, sparseSwapTests(nq, False, True, False, verbose=verbose)) # Sparse sqrtSWAP gate tests
        result = map(add, result, sparseSwapTests(nq, True, False, True, verbose=verbose)) # Sparse iSWAP-1 gate tests
        result = map(add, result, sparseSwapTests(nq, False, True, True, verbose=verbose)) # Sparse sqrtSWAP-1 gate tests
        #result = map(add, result, swapTests(nq, False, verbose=verbose))    # Old SWAP gate tests
        #result = map(add, result, swapTests(nq, True, verbose=verbose))     # Old ISWAP gate tests
        #result = map(add, result, sswapTests(nq, verbose=verbose))          # Old SqrtSWAP gate tests
        #result = map(add, result, isingTests(nq, verbose=verbose))          # Old Ising gate tests
        result = map(add, result, simpleIsingTests(nq, 0, False, verbose=verbose)) # Simple XX gate tests
        result = map(add, result, simpleIsingTests(nq, 1, False, verbose=verbose)) # Simple YY gate tests
        result = map(add, result, simpleIsingTests(nq, 2, False, verbose=verbose)) # Simple ZZ gate tests
        result = map(add, result, simpleIsingTests(nq, 0, True, verbose=verbose)) # Simple XX-1 gate tests
        result = map(add, result, simpleIsingTests(nq, 1, True, verbose=verbose)) # Simple YY-1 gate tests
        result = map(add, result, simpleIsingTests(nq, 2, True, verbose=verbose)) # Simple ZZ-1 gate tests
        result = map(add, result, sparseIsingTests(nq, 0, False, verbose=verbose)) # Sparse XX gate tests
        result = map(add, result, sparseIsingTests(nq, 1, False, verbose=verbose)) # Sparse YY gate tests
        result = map(add, result, sparseIsingTests(nq, 2, False, verbose=verbose)) # Sparse ZZ gate tests
        result = map(add, result, sparseIsingTests(nq, 0, True, verbose=verbose)) # Sparse XX-1 gate tests
        result = map(add, result, sparseIsingTests(nq, 1, True, verbose=verbose)) # Sparse YY-1 gate tests
        result = map(add, result, sparseIsingTests(nq, 2, True, verbose=verbose)) # Sparse ZZ-1 gate tests
    return tuple(result)

def main():
    argv = sys.argv[1:]
    if 2 <= len(argv) <= 4:
        if len(argv) == 2:
            passed, total = tests(int(argv[0]), int(argv[1]))
        elif len(argv) == 3:
            passed, total = tests(int(argv[0]), int(argv[1]), seed=int(argv[2]))
        else:
            passed, total = tests(int(argv[0]), int(argv[1]), seed=int(argv[2]), verbose=bool(argv[3]))
        print("Passed: " + str(passed) + "/" + str(total))
        if passed == total:
            print("PEACE AND TRANQUILITY")
            #wb.open_new_tab("https://youtu.be/SHvhps47Lmc")
        else:
            print("SORROW")
            #wb.open_new_tab("https://youtu.be/4Js-XbNj6Tk?t=37")
    else:
        print ("Syntax: " + sys.argv[0] + " <minimum number of qubits> <maximum number of qubits> <seed (optional)>")

if __name__ == "__main__":
    main()
