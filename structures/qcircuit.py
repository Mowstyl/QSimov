from structures.qregistry import QRegistry
from structures.qgate import QGate, joinGates, getGateSize
from structures.measure import Measure
import connectors.qsimovapi as qapi
import numpy as np
import gc

class QCircuit(object):
    def __init__(self, name="UNNAMED", size=None, ancilla=[]): # You can choose whether to save the circuit and apply gates separately on each computation (faster circuit creation) or to precompute the matrixes (faster execution)
        self.name = name
        self.lines = []
        self.ancilla = ancilla
        self.size = None

    def addLine(self, *args):
        try:
            if (self.lines == [] and len(args) == 1):
                size = getGateSize(args[0])
                if self.size == None:
                    self.size = size
                if self.size == size:
                    if type(args[0]) == QGate:
                        self.lines = args[0].lines
                    else:
                        self.lines = [[args[0]]]
                else:
                    raise ValueError("This cirquit requires " + str(self.size) + " QuBit gates. Received " + str(size) + " QuBit gate.")
            else:
                size = sum(map(getGateSize, args))
                if self.size == None:
                    self.size = size
                if self.size == size:
                    self.simple = False
                    self.lines += joinGates(args)
                else:
                    raise ValueError("This gate requires a " + str(self.size) + " QuBit matrix. Received " + str(size) + " QuBit matrix.")
        finally:
            gc.collect()

    def execute(self, qregistry, iterations=1, qmachine=None, args=None):
        if type(qregistry) == QRegistry and iterations > 1:
            raise ValueError("Can not do more than one iteration with a precreated registry!")
        elif (qmachine == None or qmachine == "local"):
            return qapi.execute(qregistry, iterations, self.lines, self.ancilla)
        else:
            raise ValueError("Unsupported qmachine!")
