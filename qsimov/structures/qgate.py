"""Module that provides a data structure representing a quantum gate.

Data Structures:
    QGate: Quantum Gate, built from elemental gates or other QGates
"""
import numpy as np
from qsimov.structures.qdesign import QDesign
from qsimov.structures.qstructure import _get_op_data


class QGate(QDesign):
    """Quantum Gate, built from elemental gates or other QGates."""

    def __init__(self, num_qubits, num_bits, name):
        """Quantum Gate constructor.

        Positional arguments:
            num_qubits: maximum number of qubits affected by this gate
            name: name of the gate
        """
        if num_qubits is None or not np.allclose(num_qubits % 1, 0):
            raise ValueError("num_qubits must be an integer")
        num_qubits = int(num_qubits)
        if num_qubits <= 0:
            raise ValueError("num_qubits must be positive")
        if name is None:
            raise ValueError("name can't be None")
        # List of lines. A line is a list of gates that can be executed
        # in parallel. One gate per qubit. None means Identity gate
        self.ops = []
        # Name of the gate
        self.name = str(name)
        # Maximum number of qubits affected by this gate
        self.num_qubits = num_qubits
        # Maximum number of bits affected by this gate
        self.num_bits = num_bits

    def __repr__(self):
        """Return string representation of the gate."""
        return self.name

    def __str__(self):
        """Return string representation of the gate."""
        return self.name

    def draw(self):
        raise NotImplementedError("Draw has not been implemented yet")

    def get_num_qubits(self):
        return self.num_qubits

    def get_num_bits(self):
        return self.num_bits

    def get_operations(self):
        return self.ops

    def add_operation(self, gate, targets=None, c_targets=None, outputs=None,
                      controls=None, anticontrols=None,
                      c_controls=None, c_anticontrols=None,
                      target=None, c_target=None, output=None,
                      control=None, anticontrol=None,
                      c_control=None, c_anticontrol=None):
        """Apply specified gate to specified qubit with specified controls.

        Positional arguments:
            comma separated gates, their sizes must match the number of
                qubits in the system. Sorted by their least significant
                target qubit id.
        Keyworded arguments:
            controls: id or set of ids of the qubits that will
                      act as controls
            anticontrols: id or set of ids of the qubits that will
                          act as anticontrols
            c_controls: id or set of ids of the classic bits that will
                        act as controls
            c_anticontrols: id set list of ids of the classic bits that will
                            act as anticontrols
        """
        if target is not None:
            print("[WARNING] target keyworded argument is deprecated. Please use targets instead")
            if targets is not None:
                raise ValueError("target argument can't be set alongside targets")
            targets = target
        if output is not None:
            print("[WARNING] output keyworded argument is deprecated. Please use outputs instead")
            if outputs is not None:
                raise ValueError("output argument can't be set alongside outputs")
            outputs = output
        if control is not None:
            print("[WARNING] control keyworded argument is deprecated. Please use controls instead")
            if controls is not None:
                raise ValueError("control argument can't be set alongside controls")
            controls = control
        if anticontrol is not None:
            print("[WARNING] anticontrol keyworded argument is deprecated. Please use anticontrols instead")
            if anticontrols is not None:
                raise ValueError("anticontrol argument can't be set alongside anticontrols")
            anticontrols = anticontrol
        if c_control is not None:
            print("[WARNING] c_control keyworded argument is deprecated. Please use c_controls instead")
            if c_controls is not None:
                raise ValueError("c_control argument can't be set alongside c_controls")
            c_controls = c_control
        if c_anticontrol is not None:
            print("[WARNING] c_anticontrol keyworded argument is deprecated. Please use c_anticontrols instead")
            if c_anticontrols is not None:
                raise ValueError("c_anticontrol argument can't be set alongside c_anticontrols")
            c_anticontrols = c_anticontrol
        if isinstance(gate, str) and gate.lower() == "barrier":
            self.ops.append("BARRIER")
            return
        if isinstance(gate, str) and gate.lower() == "measure":
            raise ValueError("A QGate can only contain gates or barriers")
        if outputs is not None:
            print("[WARNING] outputs argument is for measurements. Since they can't be used with QGate, do you mean targets?")
        num_qubits = self.num_qubits
        num_bits = self.num_bits
        op_data = _get_op_data(num_qubits, num_bits, gate, targets, c_targets,
                               outputs, controls, anticontrols,
                               c_controls, c_anticontrols)
        self.ops.append(op_data)

    def dagger(self):
        """Return the Conjugate Transpose of the given matrix."""
        return self.invert()

    def invert(self):
        """Return the Conjugate Transpose of the given matrix."""
        invgate = QGate(self.num_qubits, self.num_bits, self.name + "-1")
        for op_data in self.ops[::-1]:
            aux = op_data["gate"]
            if aux != "BARRIER":
                aux = aux.invert()
            invgate.add_operation(aux, targets=op_data["targets"],
                                  controls=op_data["controls"],
                                  anticontrols=op_data["anticontrols"],
                                  c_controls=op_data["c_controls"],
                                  c_anticontrols=op_data["c_anticontrols"])
        return invgate
