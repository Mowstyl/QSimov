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

import numpy as np

from abc import abstractmethod
from numbers import Number
from qsimov.structures.qbase import QBase
from qsimov.structures.qstructure import _get_op_data


class QDesign(QBase):
    def __init__(self, num_qubits: Number, num_bits: Number, name):
        if num_qubits is None or not np.allclose(num_qubits % 1, 0):
            raise ValueError("num_qubits must be an integer")
        num_qubits = int(num_qubits)
        if num_bits is None or not np.allclose(num_bits % 1, 0):
            raise ValueError("num_bits must be an integer")
        num_bits = int(num_bits)
        if num_qubits <= 0:
            raise ValueError("num_qubits must be strictly positive")
        if num_bits < 0:
            raise ValueError("num_bits cannot be negative")
        if name is None or name == "":
            raise ValueError("name can't be None / empty string")
        # Raw operation list
        self.raw_ops = []
        # Flattened operation list
        self.flat_ops = []
        # Name of the gate
        self.name = str(name)
        # Maximum number of qubits affected by this gate
        self.num_qubits = num_qubits
        # Maximum number of bits affected by this gate
        self.num_bits = num_bits

    def add_operation(self, gate, targets=None, c_targets=None, outputs=None,
                      controls=None, anticontrols=None,
                      c_controls=None, c_anticontrols=None,
                      target=None, c_target=None, output=None,
                      control=None, anticontrol=None,
                      c_control=None, c_anticontrol=None):
        targets, c_targets, outputs, \
        controls, anticontrols, \
        c_controls, c_anticontrols = verify_op_args(targets, c_targets, outputs,
                                                    controls, anticontrols,
                                                    c_controls, c_anticontrols,
                                                    target, c_target, output,
                                                    control, anticontrol,
                                                    c_control, c_anticontrol)
        if gate is None:
            raise ValueError("Gate can't be None")
        aux = gate
        if isinstance(gate, str):
            upgate = gate.upper()
            if upgate == "BARRIER" or upgate == "END":
                self.raw_ops.append(upgate)
                self.flat_ops.append(upgate)
                return
            elif upgate == "MEASURE":
                aux = None
        op_data = _get_op_data(self.num_qubits, self.num_bits, aux, targets,
                               c_targets, outputs,
                               controls, anticontrols,
                               c_controls, c_anticontrols)
        if aux is None:
            op_data["gate"] = "MEASURE"
        self.raw_ops.append(op_data)

        if not isinstance(aux, QDesign):
            noct = op_data.copy()
            del noct["c_targets"]
            self.flat_ops.append(noct)
        else:
            subops = aux.get_operations(flatten=True)
            for sub_data in subops:
                sub_data["targets"] = [targets[i] for i in sub_data["targets"]]
                sub_data["outputs"] = [c_targets[i] for i in sub_data["outputs"]]
                sub_data["controls"] = {targets[i] for i in sub_data["controls"]} | op_data["controls"]
                sub_data["c_controls"] = {c_targets[i] for i in sub_data["c_controls"]} | op_data["c_controls"]
                sub_data["anticontrols"] = {targets[i] for i in sub_data["anticontrols"]} | op_data["anticontrols"]
                sub_data["c_anticontrols"] = {c_targets[i] for i in sub_data["c_anticontrols"]} | op_data["c_anticontrols"]
                self.flat_ops.append(sub_data)

    def get_operations(self, flatten=False):
        ops = self.raw_ops
        if flatten:
            ops = self.flat_ops
        return [op.copy() for op in ops]

    @abstractmethod
    def draw(self):
        pass

    def get_num_bits(self):
        return self.num_bits

    def get_num_qubits(self):
        return self.num_qubits

    def __repr__(self):
        """Return string representation of the gate."""
        return self.name

    def __str__(self):
        """Return string representation of the gate."""
        return self.name


def verify_op_args(targets, c_targets, outputs, controls, anticontrols,
                   c_controls, c_anticontrols, target, c_target, output,
                   control, anticontrol, c_control, c_anticontrol):
    if target is not None:
        print("[WARNING] target keyworded argument is deprecated. Please use targets instead")
        if targets is not None:
            raise ValueError("target argument can't be set alongside targets")
        targets = target
    if isinstance(targets, Number):
        targets = [targets]
    elif targets is None:
        targets = []
    else:
        targets = list(targets)
    if c_target is not None:
        print("[WARNING] c_target keyworded argument is deprecated. Please use c_targets instead")
        if targets is not None:
            raise ValueError("c_target argument can't be set alongside c_targets")
        c_targets = c_target
    if isinstance(c_targets, Number):
        c_targets = [c_targets]
    elif c_targets is None:
        c_targets = []
    else:
        c_targets = list(c_targets)
    if output is not None:
        print("[WARNING] output keyworded argument is deprecated. Please use outputs instead")
        if outputs is not None:
            raise ValueError("output argument can't be set alongside outputs")
        outputs = output
    if isinstance(outputs, Number):
        outputs = [outputs]
    elif outputs is None:
        outputs = []
    else:
        outputs = list(outputs)
    if control is not None:
        print("[WARNING] control keyworded argument is deprecated. Please use controls instead")
        if controls is not None:
            raise ValueError("control argument can't be set alongside controls")
        controls = control
    if isinstance(controls, Number):
        controls = {controls}
    elif controls is None:
        controls = set()
    else:
        size = len(controls)
        controls = set(controls)
        if size > len(controls):
            print("[WARNING] repeated controls")
    if anticontrol is not None:
        print("[WARNING] anticontrol keyworded argument is deprecated. Please use anticontrols instead")
        if anticontrols is not None:
            raise ValueError("anticontrol argument can't be set alongside anticontrols")
        anticontrols = anticontrol
    if isinstance(anticontrols, Number):
        anticontrols = {anticontrols}
    elif anticontrols is None:
        anticontrols = set()
    else:
        size = len(anticontrols)
        anticontrols = set(anticontrols)
        if size > len(anticontrols):
            print("[WARNING] repeated anticontrols")
    if c_control is not None:
        print("[WARNING] c_control keyworded argument is deprecated. Please use c_controls instead")
        if c_controls is not None:
            raise ValueError("c_control argument can't be set alongside c_controls")
        c_controls = c_control
    if isinstance(c_controls, Number):
        c_controls = {c_controls}
    elif c_controls is None:
        c_controls = set()
    else:
        size = len(c_controls)
        c_controls = set(c_controls)
        if size > len(c_controls):
            print("[WARNING] repeated classic controls")
    if c_anticontrol is not None:
        print("[WARNING] c_anticontrol keyworded argument is deprecated. Please use c_anticontrols instead")
        if c_anticontrols is not None:
            raise ValueError("c_anticontrol argument can't be set alongside c_anticontrols")
        c_anticontrols = c_anticontrol
    if isinstance(c_anticontrols, Number):
        c_anticontrols = {c_anticontrols}
    elif c_anticontrols is None:
        c_anticontrols = set()
    else:
        size = len(c_anticontrols)
        c_anticontrols = set(c_anticontrols)
        if size > len(c_anticontrols):
            print("[WARNING] repeated classic anticontrols")

    return targets, c_targets, outputs, controls, anticontrols, c_controls, c_anticontrols
