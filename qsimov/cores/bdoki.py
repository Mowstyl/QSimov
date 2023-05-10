import numpy as np
import qsimov.connectors.parser as gt
import sympy as sp

from sympy.functions import conjugate
from sympy.matrices import Matrix
from sympy.physics.quantum import TensorProduct


_pauli = [None for i in range(3)]
_pauli[0] = gt.PauliX() # σx
_pauli[1] = gt.PauliY() # σy
_pauli[2] = gt.PauliZ() # σz

# Unitary vectors
_u = [None for i in range(3)]
_u[0] = Matrix([[1, 0, 0]])
_u[1] = Matrix([[0, 1, 0]])
_u[2] = Matrix([[0, 0, 1]])


def gate_new(num_qubits, matrix, verbose):
    if num_qubits > 1:
        raise ValueError("Gates for bloch spheres can affect at most one single qubit")
    return BRotation(matrix, verbose)


def gate_get(gate, verbose):
    return gate.O


def registry_new(num_qubits, verbose):
    return BRegistry(num_qubits, verbose=verbose)


def registry_clone(registry, num_threads, verbose):
    newr = BRegistry(registry.nq, verbose=verbose)
    for id in range(registry.nq):
        newr.r[id].v = Matrix([registry.r[id].v[0, 0],
                               registry.r[id].v[1, 0],
                               registry.r[id].v[2, 0]])
    return newr


def registry_del(registry, verbose):
    del registry


def registry_apply(registry, gate, target_list, control_set, anticontrol_set, num_threads, verbose):
    if control_set != None and len(control_set) > 0:
        raise ValueError("Gates for bloch spheres can have controls")
    if anticontrol_set != None and len(anticontrol_set) > 0:
        raise ValueError("Gates for bloch spheres can have anticontrols")
    return registry.apply_gate(gate, targets=target_list, verbose=verbose)


def registry_join(most_registry, least_registry, num_threads, verbose):
    newr = BRegistry(most_registry.nq + least_registry.nq, init=False, verbose=verbose)
    newr.r = most_registry.r + least_registry.r


def registry_measure(registry, mask, roll_list, num_threads, verbose):
    targets = [id for id in range(registry.nq) if bool(mask & (1 << id))]
    registry.measure(targets=targets, rands=roll_list, verbose=verbose)


def registry_prob(registry, qubit_id, num_threads, verbose):
    return registry.r[qubit_id].v[2, 0]


class BRotation(object):
    def __init__(self, matrix, verbose=False):
        """https://quantumcomputing.stackexchange.com/questions/16533/can-i-find-the-axis-of-rotation-for-any-single-qubit-gate"""
        # e^iα = sqrt(det(O))
        self.verbose = verbose
        self.O = matrix
        if verbose:
            print("Creating gate for matrix", matrix)
        eia = sp.sqrt(sp.det(self.O))
        alpha = sp.atan2(sp.im(eia), sp.re(eia))
        e_ia = conjugate(eia)
        self.theta = 2 * sp.acos(e_ia * sp.trace(self.O) / 2)
        _2isint2 = -2j * sp.sin(self.theta / 2)
        self.n = Matrix([[(e_ia * sp.trace(self.O @ pw)) / _2isint2 for pw in _pauli]])
    
    def __str__(self):
        return f"Rotation θ = {self.theta}, axis = {self.n}"
    
    def __repr__(self):
        return str(self)

    def rotation_matrix(self):
        cost = sp.cos(self.theta)
        sint = sp.sin(self.theta)
        u = self.n
        auxu = [TensorProduct(u.cross(_u[i]), _u[i]) for i in range(3)]
        ux = (auxu[0] + auxu[1] + auxu[2]).reshape(3, 3)
        r = cost * sp.eye(3) + sint * ux + (1 - cost) * TensorProduct(u, u).reshape(3,3)
        return r


class QuBit(object):
    def __init__(self):
        self.v = Matrix([0, 0, 1])
    
    def __str__(self):
        return str(self.v.transpose())
    
    def __repr__(self):
        return str(self)
    
    def apply_gate(self, rot, verbose=False):
        m = rot.rotation_matrix()
        qb = QuBit()
        qb.v = m @ self.v
        if verbose:
            print("Prev:", self.v)
            print("Post:", qb.v)
        return qb

    def measure(self, rand=None):
        zproj = self.v[2, 0]
        qb = QuBit()
        res = 0
        odds = 1 - ((zproj + 1) / 2)
        if rand is None:
            rand = np.random.rand()
        if rand < odds:
            qb.v = Matrix([0, 0, -1])
            res = 1
        return qb, res


class BRegistry(object):
    def __init__(self, nq, init=True, verbose=False):
        self.nq = nq
        self.r = None
        if init:
            self.r = tuple(QuBit() for i in range(nq))

    def apply_gate(self, gate, targets=0, verbose=False):
        if type(targets) == int:
            targets = [targets]
        t_set = set(targets)
        if len(t_set) != len(targets):
            raise ValueError("Targets has repeated elements")
        rot = BRotation(gate, verbose=verbose)
        if verbose:
            print("Rotation matrix:", rot.rotation_matrix())
        newr = BRegistry(self.nq, False)
        newr.r = tuple(self.r[i].apply_gate(rot, verbose=verbose) if i in t_set else self.r[i] for i in range(self.nq))
        return newr

    def measure(self, targets=0, rands=None, verbose=False):
        if type(targets) == int:
            targets = [targets]
        t_set = set(targets)
        if len(t_set) != len(targets):
            raise ValueError("Targets has repeated elements")
        if rands is None:
            rands = np.random.rand(len(t_set))
        if len(rands) != len(t_set):
            raise ValueError("Rands has to have the same number of elements as targets")
        newr = BRegistry(self.nq, False)
        rands = [rands.pop(0) if i in t_set else None for i in range(self.nq)]
        aux = tuple(self.r[i].measure(rands[i]) if i in t_set else self.r[i] for i in range(self.nq))
        newr.r = tuple(aux[i][0] if i in t_set else aux[i] for i in range(self.nq))
        res = tuple(aux[i][1] if i in t_set else None for i in range(self.nq))
        return newr, res
