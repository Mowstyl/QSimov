from qsimov.structures.qregistry import QRegistry
from qsimov.structures.qsystem import QSystem
from qsimov.structures.qgate import _rebuildGateName, _getParties
from qsimov.structures.measure import Measure
import qsimov.connectors.qsimovapi as qapi
import gc


class QCircuit(object):
    def __init__(self, name="UNNAMED", size=None, ancilla=[]):
        self.empty = True
        self.name = name
        self.lines = []
        self.ancilla = ancilla
        self.size = size

    def addLine(self, *args):
        try:
            if args is not None and len(args) > 0:
                if any(isinstance(e, Measure) for e in args):
                    if len(args) > 1:
                        raise ValueError("You can only use 1 measure per line")
                    else:
                        size = len(args[0].mask)
                        if (self.empty):
                            self.size = size
                            self.empty = False
                        if (self.size != size):
                            raise ValueError("This circuit requires a measurement mask for " + str(self.size) + " QuBits. Received mask for " + str(size) + " QuBits.")
                        self.lines += [args]
                else:
                    args = [_rebuildGateName(gate) for gate in args]
                    parties = _getParties(args)
                    size = len(parties)
                    if size > 0:
                        if (self.empty):
                            self.size = size
                            self.empty = False
                        if (self.size != size):
                            raise ValueError("This circuit requires gates for " + str(self.size) + " QuBits. Received gates for " + str(size) + " QuBits.")

                        self.lines += [args]
                    else:
                        print("No gates. Just Monika.")
            else:
                print("Here goes nothing.")
        finally:
            gc.collect()

    def execute(self, qubits, iterations=1, qmachine=None, args={"useSystem": True}):
        if (isinstance(qubits, QRegistry) or isinstance(qubits, QSystem)) and iterations > 1:
            raise ValueError("Can not do more than one iteration with a precreated registry or system!")
        elif (qmachine is None or qmachine == "local"):
            return qapi.execute(qubits, iterations, self.lines, self.ancilla, args["useSystem"])
        else:
            raise ValueError("Unsupported qmachine!")
