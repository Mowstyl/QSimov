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
        self.num_bits = num_bits

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
        if gate is None:
            raise ValueError("Gate can't be None")
        if isinstance(gate, str) and gate.lower() == "barrier":
            self.ops.append("BARRIER")
            return
        aux = gate
        if isinstance(gate, str) and gate.lower() == "measure":
            aux = None
        op_data = _get_op_data(self.num_qubits, self.num_bits, aux, targets,
                               c_targets, outputs,
                               controls, anticontrols,
                               c_controls, c_anticontrols)
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
