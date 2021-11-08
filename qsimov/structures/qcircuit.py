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

    def __init__(self, num_qubits, name, ancilla=None):
        """Quantum Circuit constructor.

        Positional arguments:
            num_qubits: maximum number of qubits affected by this gate
            name: name of the gate
        Keyworded arguments:
            ancilla: list of values of the ancilliary qubits (0 or 1)
        """
        if num_qubits is None or not np.allclose(num_qubits % 1, 0):
            raise ValueError("num_qubits must be an integer")
        num_qubits = int(num_qubits)
        if num_qubits <= 0:
            raise ValueError("num_qubits must be positive")
        if name is None:
            raise ValueError("name can't be None")
        if ancilla is None:
            ancilla = []
        aux = [int(val) for val in ancilla]
        if not np.allclose(aux, ancilla):
            raise ValueError("ancilla must be a list of integers")
        ancilla = aux
        if not all(0 <= val <= 1 for val in ancilla):
            raise ValueError("ancilla must be a list of 0 and 1")
        self.name = name
        self.ops = []
        self.ancilla = ancilla
        self.num_qubits = num_qubits

    def draw(self):
        raise NotImplementedError("Draw has not been implemented yet")

    def get_num_qubits(self):
        return self.num_qubits

    def get_operations(self):
        return self.ops

    def add_operation(self, gate, targets=None, controls=None,
                      anticontrols=None):
        """Apply specified operation to this QCircuit.

        Positional arguments:
            gate: string with the name of the gate to apply,
                  a QGate or a Measure
        Keyworded arguments:
            targets: id or list of ids of the qubits the gate will target
            controls: id or set of ids of the qubit that will act as controls
            anticontrols: id or set of ids of the qubit that will act as
                          anticontrols
            num_threads: number of threads to use
        """
        if gate is None:
            raise ValueError("Gate can't be None")
        if isinstance(gate, str) and gate.lower() == "barrier":
            self.ops.append("BARRIER")
        aux = gate
        if isinstance(gate, str) and gate.lower() == "measure":
            aux = None
        op_data = _get_op_data(self.num_qubits, aux, targets,
                               controls, anticontrols)
        if aux is None:
            op_data["gate"] = "MEASURE"
        self.ops.append(op_data)

    def execute(self, qubits, iterations=1, qmachine=None,
                args={"useSystem": True}, optimize=True):
        """Execute the circuit with specified configuration options.

        Positional arguments:
            qubits: list of 0 and 1 with the initial value of each
                        non ancilla qubit of the circuit
                OR  QRegistry or QSystem
        Keyworded arguments:
            iterations: the circuit will be executed this amount of times
            qmachine: string with the name of the machine in which the circuit
                      will be executed. Right now only supported QSimov
                      simulator (None or "local").
            args: dictionary with extra arguments for the specified qmachine.
                  For QSimov (local), "useSystem" is the only extra argument.
                  Defaults to True. Makes QSimov use QSystems instead of
                  QRegistries by default.
            optimize: whether or not to use oplines instead of lines.
                      Default = True
        """
        result = None
        lines = self.lines
        if optimize:
            lines = self.oplines
        if ((isinstance(qubits, QRegistry) or isinstance(qubits, QSystem))
                and iterations > 1):
            raise ValueError("Can not do more than one iteration with a " +
                             "precreated registry or system!")
        elif (qmachine is None or qmachine == "local"):
            result = qapi.execute(qubits, iterations, lines, self.ancilla,
                                  args["useSystem"], optimize)
        else:
            raise ValueError("Unsupported qmachine!")
        return result
