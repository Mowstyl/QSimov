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

"""Module that simplify the way of working with Doki gates.

Data Structures:
    SimpleGate: Quantum gate with its matrix

Functions:
    add_gate: adds a gate to the list of available gates.
"""
import doki
import numpy as np
import qsimov.connectors.parser as prs
import sympy as sp
import re

from collections.abc import Iterable
from qsimov.structures.qbase import QBase
from sympy.matrices import Matrix
from sympy.codegen.cfunctions import log2
from sympy.physics.quantum.dagger import Dagger


__rep__ = re.compile(r"^" + prs._gate_name_re + "+$")

class SimpleGate(QBase):
    """Quantum gate with its associated matrix."""

    def __init__(self, gate_string):
        """Load a gate that is in the list of gates."""
        name, args, invert, self_invert = prs.get_gate_data(gate_string)
        self.matrix = _get_gate_matrix(name, args, self_invert)
        if self.matrix.shape[0] != self.matrix.shape[1]:
            raise ValueError("Not a square matrix")
        self.num_qubits = int(log2(self.matrix.shape[0]))
        if 2**self.num_qubits != self.matrix.shape[0]:
            raise ValueError("Shape must be 2^n x 2^n. n = number of qubits")
        inverse = self.matrix
        if not self_invert:
            inverse = self.matrix
            if isinstance(self.matrix, Matrix):
                inverse = Dagger(self.matrix)
            else:
                inverse = self.matrix.conj().T
        _check_inverse(gate_string, self.matrix, inverse)
        if invert:
            self.matrix = inverse
        aux = np.array(self.matrix).astype(complex)
        self.gate = doki.gate_new(self.num_qubits, aux.tolist(), False)
        self._str = name
        if args is not None:
            if prs._gate_data[name][2]:  # has_invert_arg
                args = args[:-1]
            self._str += "(" + ",".join([str(arg) for arg in args]) + ")"
        if invert and not self_invert:
            self._str += "-1"

    def __eq__(self, other):
        if not isinstance(other, SimpleGate):
            return False
        if self is other:
            return True
        return self._str == other._str

    def __str__(self):
        gate_str = self._str
        if gate_str[-2:] == "-1":
            gate_str = gate_str[:-2] + "\u2020"
        return gate_str

    def __repr__(self):
        return str(self)

    def get_doki_gate(self):
        return self.gate

    def get_matrix(self):
        return self.matrix

    def get_num_qubits(self):
        return self.num_qubits

    def invert(self):
        gate_str = self._str
        if gate_str[-2:] == "-1":
            gate_str = gate_str[:-2]
        else:
            gate_str += "-1"
        return SimpleGate(gate_str)


def _check_inverse(name, matrix, inverse):
    got = None
    if isinstance(matrix, Matrix):
        got = np.array(matrix @ inverse).astype(complex)
        '''
        if not ex.equals(got) and not ex.equals(sp.N(got)):
            print("Expected:", ex)
            print("Got:", got)
            print("Got:", sp.simplify(got))
            print(ex.equals(got))
            print([[ex[i, j].equals(got[i, j]) for j in range(ex.shape[1])] for i in range(ex.shape[0])])
            print(sp.N(got))
            print(ex.equals(sp.N(got)))
            raise ValueError("Failed testing inverse of " + name)
        '''
    else:
        got = np.dot(matrix, inverse)
    ex = np.eye(matrix.shape[0], dtype=complex)
    if not np.allclose(got, ex):
        raise ValueError("Failed testing inverse of " + name)


def _get_gate_matrix(name, args, self_invert):
    """Get gate matrix from parser's get_gate_data info."""
    gate_function = prs._gate_func[name]
    matrix = None
    if type(args) != str and isinstance(args, Iterable):
        matrix = gate_function(*args)
    elif args is not None:
        matrix = gate_function(args)
    else:
        matrix = gate_function()

    return Matrix(matrix)


def add_gate(name, funct, min_args, max_args, has_invert_arg=None,
             is_own_inverse=False, aliases=None, overwrite=False):
    """Add specified gate to the list of available gates."""
    if has_invert_arg is not None:
        print("[WARNING] has_invert_arg keyworded argument has been deprecated, and will be removed. Currently it does nothing.")
    if __rep__.match(name) is None:
        raise ValueError("Gate names may only contain letters, numbers and underscores, matching [a-zA-Z0-9_]+ regex.")
    if name in prs._gate_func and not overwrite:
        raise ValueError(name + " gate has already been defined. " +
                         "Set overwrite=True if you want to overwrite it.")
    if aliases is None:
        aliases = {name.lower()}
    else:
        aliases = {alias.lower() for alias in aliases}
        aliases.add(name)
    for alias in aliases:
        if alias in prs._gate_alias and prs._gate_alias[alias] != name:
            raise ValueError("Alias " + alias + " for gate " + name +
                             " is in use by gate " + prs._gate_alias[alias] +
                             ". Either remove this alias or modify the older" +
                             " gate.")
    prs._gate_func[name] = funct
    prs._gate_data[name] = (min_args, max_args, is_own_inverse)
    for alias in aliases:
        prs._gate_alias[alias.lower()] = name
