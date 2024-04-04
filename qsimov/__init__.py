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

from qsimov.structures.funmatrix import Funmatrix
from qsimov.structures.qregistry import QRegistry, superposition
from qsimov.structures.qsystem import QSystem, join_systems
from qsimov.structures.qgate import QGate
from qsimov.structures.qcircuit import QCircuit
from qsimov.structures.simple_gate import SimpleGate, add_gate
from qsimov.connectors.drewom import Drewom
from qsimov.connectors.parser import get_available_gates, get_gate_aliases
