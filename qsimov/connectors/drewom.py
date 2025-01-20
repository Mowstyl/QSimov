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

"""Handles design execution on target machine."""
import qsimov.connectors.qsimovapi as qapi
from qsimov.structures.qcircuit import QCircuit


class Drewom(object):
    """Handle QDesign execution on specified machine."""

    def __init__(self, qmachine="doki",
                 extra=None):
        """Initialize Drewom structure."""
        self.qmachine = qmachine
        if extra is None:
            extra = {}
        self.extra = extra

    def execute(self, qcircuit, shots=1, iterations=None):
        """Execute given design in a machine shots times."""
        if iterations is not None:
            print("[WARNING] iterations argument is deprecated. Please use shots instead")
            shots = iterations
        if not isinstance(qcircuit, QCircuit):
            raise ValueError("qcircuit must be a QCircuit")
        if self.qmachine == "doki":
            return qapi.execute(qcircuit, shots=shots, **self.extra)
        else:
            raise NotImplementedError("Drewom does not support that machine")
