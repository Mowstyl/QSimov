import numpy as np
import scipy.sparse as sp
from qsimov import SimpleGate, Funmatrix
import doki
import re


_h_regex = re.compile(r"^H[0-9]*$")
_everything = []

def get_gate(gate):
    gate_mat = None
    
    if not _h_regex.match(gate):
        raw_sgate = SimpleGate(gate)
        shape = raw_sgate.matrix.shape
        matrix = [[complex(raw_sgate.matrix[i, j]) for j in range(shape[1])] for i in range(shape[0])]
        gate_mat = doki.funmatrix_create(matrix, False)
    else:
        nq = 1
        if len(gate[1:]) > 0:
            nq = int(gate[1:])
        if nq < 1:
            raise ValueError("Hadamard gate for less than one qubit is impossible")
        gate_mat = doki.funmatrix_hadamard(nq, False)
    _everything.append(gate_mat)
    return gate_mat


def get_identity(nq):
    aux = doki.funmatrix_identity(nq, False)
    _everything.append(aux)
    return aux


def get_system(nq):
    rho = doki.funmatrix_statezero(nq, False)
    _everything.append(rho)
    return rho


def get_projector(nq, qubit, value):
    proj = None
    if not bool(value):
        proj = get_system(nq)
    else:
        proj = doki.funmatrix_create([[0, 0], [0, 1]], False)
        _everything.append(proj)
        if qubit > 0:
            proj = doki.funmatrix_kron(get_identity(qubit, False), proj, False)
            _everything.append(proj)
        if qubit < nq-1:
            proj = doki.funmatrix_kron(proj, get_identity(nq-qubit-1), False)
            _everything.append(proj)
    return proj


def apply_gate(nq, sys, gate, qubit, num_controls=0, verbose=False):
    newgate = get_gate(gate)
    nqg = np.log2(doki.funmatrix_shape(newgate, False)[0])
    if nqg % 1 != 0:
        raise ValueError(f"Wrong gate size: 2^{nqg}")
    nqg = int(nqg)
    while num_controls > 0:
        newgate = doki.funmatrix_addcontrol(newgate, verbose)
        _everything.append(newgate)
        num_controls -= 1
        nqg += 1
    if qubit > 0:
        newgate = doki.funmatrix_kron(get_identity(qubit), newgate, verbose)
        _everything.append(newgate)
    if qubit < nq-nqg:
        newgate = doki.funmatrix_kron(newgate, get_identity(nq-qubit-nqg), verbose)
        _everything.append(newgate)
    gatetrans = doki.funmatrix_dagger(newgate, verbose)
    _everything.append(gatetrans)
    res = doki.funmatrix_matmul(newgate, sys, verbose)
    _everything.append(res)
    res = doki.funmatrix_matmul(res, gatetrans, verbose)
    _everything.append(res)
    return res

def measure(nq, sys, qubit, verbose=False):
    P = get_projector(nq, qubit, 1)
    P_rho = doki.funmatrix_matmul(P, sys, verbose)
    _everything.append(P_rho)
    p = doki.funmatrix_trace(P_rho, verbose).real
    roll = np.random.rand()
    res = roll < p
    if not bool(res):
        P = get_projector(nq, qubit, 0)
        P_rho = doki.funmatrix_matmul(P, sys, verbose)
        _everything.append(P_rho)
        p = doki.funmatrix_trace(P_rho, verbose).real
    numerator = doki.funmatrix_matmul(P_rho, P, verbose)
    _everything.append(numerator)
    aux = doki.funmatrix_scalar_div(numerator, p, verbose)
    return (res, aux)


print("Starting!", flush=True)
sys = get_system(2)
print("rho0 =", Funmatrix(sys)[:], flush=True)
sys2 = apply_gate(2, sys, "H", 0)
print("rho1 =", Funmatrix(sys2)[:], flush=True)
#del sys
sys3 = apply_gate(2, sys2, "X", 0, num_controls=0, verbose=True)
print("Abraracurcix", flush=True)
print("rho2 =", Funmatrix(sys3)[:], flush=True)
#del sys2
res, sys4 = measure(2, sys3, 0)
print("rho3 =", Funmatrix(sys4)[:], flush=True)
#del sys3
print(res)
del _everything
_everything = []

sys = get_system(1)
print("rho0 =", Funmatrix(sys)[:], flush=True)
sys2 = apply_gate(1, sys, "H", 0)
print("rho1 =", Funmatrix(sys2)[:], flush=True)
#del sys
res, sys3 = measure(1, sys2, 0)
#del sys2
print("rho2 =", Funmatrix(sys3)[:], flush=True)
print(res)
#del sys3
