from structures.qregistry import QRegistry
from structures.qgate import QGate, getGateSize, _rebuildGateName, _getQuBitArg
from structures.measure import Measure
import connectors.qsimovapi as qapi
import numpy as np
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
                    args = [_rebuildGateName(gate) if isinstance(type(gate), str) or (isinstance(gate, Iterable) and isinstance(type(gate[0]), str)) else gate for gate in args]
                    size = sum([getGateSize(gate) for gate in args])
                    if size > 0:
                        if (self.empty):
                            self.size = size
                            self.empty = False
                        if (self.size != size):
                            raise ValueError("This circuit requires gates for " + str(self.size) + " QuBits. Received gates for " + str(size) + " QuBits.")
                        parties = set()
                        for i in range(len(args)):
                            myparty = _getQuBitArg(args[i][0])
                            if myparty is None:
                                myparty = set([i])
                            if args[i][0] is not None:
                                if len(myparty.intersection(args[i][1])) == 0:
                                    myparty = myparty.union(args[i][1])
                                else:
                                    raise ValueError("You can't apply a gate to a qubit and use it as a control: " + str(myparty.intersection(args[i][1])))
                            if args[i][0] is not None:
                                if len(myparty.intersection(args[i][2])) == 0:
                                    myparty = myparty.union(args[i][2])
                                else:
                                    raise ValueError("You can't apply a gate to a qubit and use it as a control, or use it as control and anticontrol at the same time: " + str(myparty.intersection(args[i][2])))
                            if len(parties.intersection(myparty)) == 0:
                                parties = parties.union(myparty)
                            else:
                                raise ValueError("You can't apply two or more gates to the same qubit in the same line: " + str(parties.intersection(myparty)))

                        self.lines += [args]
                    else:
                        print("No gates. Just Monika.")
            else:
                print("Here goes nothing.")
        finally:
            gc.collect()

    def execute(self, qregistry, iterations=1, qmachine=None, args=None, useSystem=True):
        if (isinstance(qregistry, QRegistry) or isinstance(qregistry, QSystem)) and iterations > 1:
            raise ValueError("Can not do more than one iteration with a precreated registry or system!")
        elif (qmachine == None or qmachine == "local"):
            return qapi.execute(qregistry, iterations, self.lines, self.ancilla, useSystem)
        else:
            raise ValueError("Unsupported qmachine!")
