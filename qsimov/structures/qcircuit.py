"""Module that provides a data structure representing a quantum circuit.

Data Structures:
    QCircuit: Quantum Circuit, built from gates and measurements
"""
from qsimov.structures.qregistry import QRegistry
from qsimov.structures.qsystem import QSystem
from qsimov.structures.qgate import _get_line_args, _add_gates
from qsimov.structures.measure import Measure
import qsimov.connectors.qsimovapi as qapi
import gc


class QCircuit(object):
    """Quantum Circuit, built from gates and measurements."""

    def __init__(self, name="UNNAMED", size=None, ancilla=[]):
        """Quantum Circuit constructor.

        Keyworded arguments:
            name: name of the circuit. Default="UNNAMED"
            size: maximum number of qubits affected by this gate
                If nothing specified it will be calculated in the
                first call to add_line
            ancilla: list of values of the ancilliary qubits (0 or 1)
        """
        self.empty = True
        self.name = name
        self.lines = []
        self.oplines = []
        self.freeindexes = None
        self.lastindex = -1
        self.ancilla = ancilla
        self.size = size

    def addLine(self, *args, **kwargs):
        """Use add_line method instead. DEPRECATED."""
        print("Method QCircuit.addLine is deprecated.",
              "Please use add_line if you seek the same functionality")
        return self.add_line(*args, **kwargs)

    def _add_measure(self, args, add_to_lines):
        """Add a Measure object to oplines (and to lines if specified)."""
        if len(args) > 1:
            raise ValueError("You can only use 1 measure per line")
        else:
            size = len(args[0].mask)
            if (self.empty):
                self.size = size
                self.empty = False
                self.freeindexes = [0 for i in range(size)]
            if (self.size != size):
                raise ValueError("This circuit requires a measurement mask " +
                                 "for " + str(self.size) + " QuBits. " +
                                 " Received mask for " +
                                 str(size) + " QuBits.")
            if add_to_lines:
                self.lines += [args]
            freeindex = max([self.freeindexes[i]
                             for i in range(len(args[0].mask))])
            for i in range(len(args[0].mask)):
                self.freeindexes[i] = freeindex + 1
                if freeindex > self.lastindex:
                    self.oplines.append([])
                    self.lastindex += 1
                self.oplines[freeindex].append(args[0])

    def add_line(self, *args, **kwargs):
        """Apply specified gate to specified qubit with specified controls.

        Positional arguments:
            comma separated gates, their sizes must match the number of
                qubits in the system. Sorted by their least significant
                target qubit id.
            OR
            Measure object
        Keyworded arguments:
            add_to_lines: whether to add args to lines or only to oplines
            offset: offset to add to all the qubit ids each arg of this line
            control: id or list of ids of the qubit that will act as
                      controls
            anticontrol: id or list of ids of the qubit that will act as
                      anticontrols
        """
        add_to_lines, offset, controls, anticontrols = _get_line_args(**kwargs)
        if args is None or len(args) == 0:
            print("Here goes nothing.")
        try:
            if any(isinstance(e, Measure) for e in args):
                self._add_measure(args, add_to_lines)
            else:
                _add_gates(self, args, add_to_lines, offset,
                           controls, anticontrols)
        finally:
            gc.collect()

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
