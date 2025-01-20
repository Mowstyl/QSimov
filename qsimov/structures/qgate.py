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

"""Module that provides a data structure representing a quantum gate.

Data Structures:
    QGate: Quantum Gate, built from elemental gates or other QGates
"""
import numpy as np
from qsimov.structures.qdesign import QDesign, verify_op_args
from qsimov.structures.qstructure import _get_op_data


class QGate(QDesign):
    """Quantum Gate, built from elemental gates or other QGates."""
    def add_operation(self, gate, targets=None, c_targets=None, outputs=None,
                      controls=None, anticontrols=None,
                      c_controls=None, c_anticontrols=None,
                      target=None, c_target=None, output=None,
                      control=None, anticontrol=None,
                      c_control=None, c_anticontrol=None):
        if isinstance(gate, str):
            lowgate = gate.lower()
            if lowgate == "measure" or lowgate == "end":
                raise ValueError("MEASURE and END cannot be used in a QGate")
        super().add_operation(gate, targets=targets, c_targets=c_targets,
                              outputs=outputs,
                              controls=controls, anticontrols=anticontrols,
                              c_controls=c_controls, c_anticontrols=c_anticontrols,
                              target=target, c_target=c_target, output=output,
                              control=control, anticontrol=anticontrol,
                              c_control=c_control, c_anticontrol=c_anticontrol)

    def draw(self):
        raise NotImplementedError("Draw has not been implemented yet")

    def dagger(self):
        """Return the Conjugate Transpose of the given matrix."""
        return self.invert()

    def invert(self):
        """Return the Conjugate Transpose of the given matrix."""
        invgate = QGate(self.num_qubits, self.num_bits, self.name + "-1")
        for op_data in self.get_operations()[::-1]:
            aux = op_data["gate"]
            if aux != "BARRIER":
                aux = aux.invert()
            invgate.add_operation(aux, targets=op_data["targets"],
                                  controls=op_data["controls"],
                                  anticontrols=op_data["anticontrols"],
                                  c_controls=op_data["c_controls"],
                                  c_anticontrols=op_data["c_anticontrols"])
        return invgate
