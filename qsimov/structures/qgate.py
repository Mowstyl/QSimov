"""Module that provides a data structure representing a quantum gate.

Data Structures:
    QGate: Quantum Gate, built from elemental gates or other QGates
"""
import numpy as np
from qsimov.structures.qdesign import QDesign
from qsimov.structures.qstructure import _get_op_data


class QGate(QDesign):
    """Quantum Gate, built from elemental gates or other QGates."""

    def __init__(self, num_qubits, name):
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

    def get_operations(self):
        return self.ops

    def add_operation(self, gate, targets=None, controls=None,
                      anticontrols=None):
        """Apply specified gate to specified qubit with specified controls.

        Positional arguments:
            comma separated gates, their sizes must match the number of
                qubits in the system. Sorted by their least significant
                target qubit id.
        Keyworded arguments:
            controls: id or list of ids of the qubit that will act as
                      controls
            anticontrols: id or list of ids of the qubit that will act as
                          anticontrols
        """
        if isinstance(gate, str) and gate.lower() == "barrier":
            self.ops.append("BARRIER")
        if isinstance(gate, str) and gate.lower() == "measure":
            raise ValueError("A QGate can only contain gates or barriers")
        num_qubits = self.num_qubits
        op_data = _get_op_data(num_qubits, gate, targets,
                               controls, anticontrols)
        self.ops.append(op_data)

    def dagger(self):
        """Return the Conjugate Transpose of the given matrix."""
        return self.invert()

    def invert(self):
        """Return the Conjugate Transpose of the given matrix."""
        invgate = QGate(self.num_qubits, self.name + "-1")
        for op_data in self.ops[::-1]:
            aux = op_data["gate"]
            if aux != "BARRIER":
                aux = aux.invert()
            invgate.add_operation(aux, targets=op_data["targets"],
                                  controls=op_data["controls"],
                                  anticontrols=op_data["anticontrols"])
        return invgate
