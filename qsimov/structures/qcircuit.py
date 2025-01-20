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

"""Module that provides a data structure representing a quantum circuit.

Data Structures:
    QCircuit: Quantum Circuit, built from gates and measurements
"""
from qsimov.structures.qregistry import QRegistry
from qsimov.structures.qstructure import _get_op_data
from qsimov.structures.qsystem import QSystem
from qsimov.structures.qdesign import QDesign
import qsimov.connectors.qsimovapi as qapi
import numpy as np


class QCircuit(QDesign):
    """Quantum Circuit, built from gates and measurements."""

    def __init__(self, num_qubits, num_bits, name, ancilla=None):
        """Quantum Circuit constructor.

        Positional arguments:
            num_qubits: maximum number of qubits affected by this circuit
            num_bits: maximum number of classical bits affected by this circuit
            name: name of the gate
        Keyworded arguments:
            ancilla: list of values of the ancilliary qubits (0 or 1)
        """
        super().__init__(num_qubits, num_bits, name)
        if ancilla is None:
            ancilla = []
        aux = [int(val) for val in ancilla]
        if not np.allclose(aux, ancilla):
            raise ValueError("ancilla must be a list of integers")
        ancilla = aux
        if not all(0 <= val <= 1 for val in ancilla):
            raise ValueError("ancilla must be a list of 0 and 1")
        self.ancilla = ancilla

    def draw(self):
        raise NotImplementedError("Draw has not been implemented yet")
