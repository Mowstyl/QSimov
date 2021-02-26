from qsimov.structures.qregistry import QRegistry
from qsimov.structures.qsystem import QSystem
from qsimov.structures.qgate import _rebuildGateName, _getParties, QGate
from qsimov.structures.measure import Measure
import qsimov.connectors.qsimovapi as qapi
import gc


class QCircuit(object):
    def __init__(self, name="UNNAMED", size=None, ancilla=[]):
        self.empty = True
        self.name = name
        self.lines = []
        self.oplines = []
        self.freeindexes = None
        self.lastindex = -1
        self.ancilla = ancilla
        self.size = size

    def addLine(self, *args, **kwargs):
        add_line = True
        offset = 0
        if "add_line" in kwargs:
            add_line = kwargs["add_line"]
        if "offset" in kwargs:
            offset = kwargs["offset"]
        if "controls" in kwargs:
            controls = kwargs["controls"]
        if "anticontrols" in kwargs:
            anticontrols = kwargs["anticontrols"]
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
                            self.freeindexes = [0 for i in range(size)]
                        if (self.size != size):
                            raise ValueError("This circuit requires a measurement mask for " + str(self.size) + " QuBits. Received mask for " + str(size) + " QuBits.")
                        if add_line:
                            self.lines += [args]
                        for i in range(len(args)):
                            arg = args[i]
                            if arg is not None:
                                freeindex = max([self.freeindexes[i] for i in range(self.size)])
                                for i in range(self.size):
                                    self.freeindexes[i] = freeindex + 1
                                if freeindex > self.lastindex:
                                    self.oplines.append([])
                                    self.lastindex += 1
                                self.oplines[freeindex].append(arg)
                else:
                    args = [_rebuildGateName(gate) for gate in args]
                    parties = _getParties(args)
                    size = len(parties)
                    if size > 0:
                        if (self.empty):
                            self.size = size
                            self.empty = False
                            self.freeindexes = [0 for i in range(size)]
                        if (self.size != size):
                            raise ValueError("This circuit requires gates for " + str(self.size) + " QuBits. Received gates for " + str(size) + " QuBits.")
                        if add_line:
                            self.lines += [args]
                        for i in range(len(args)):
                            arg = args[i]
                            if arg is not None:
                                if offset > 0:
                                    arg[1] = controls.union({control + offset for control in arg[1]})
                                    arg[2] = anticontrols.union({acontrol + offset for acontrol in arg[2]})
                                if isinstance(arg[0], QGate):
                                    for line in arg[0].oplines:
                                        self.addLine(*[None for j in range(i)],
                                                     *line,
                                                     *[None for j in range(self.size - arg[0].size - i)],
                                                     add_line=False, offset=i,
                                                     controls=arg[1],
                                                     anticontrols=arg[2])
                                else:
                                    parties = _getParties([arg if j == i else None for j in range(len(args))], ignore_empty=True)
                                    num_targets = len(parties) - len(arg[1]) - len(arg[2])
                                    freeindex = max([self.freeindexes[party] for party in parties])
                                    skip = False
                                    if freeindex > 0:
                                        previousGate = self.oplines[freeindex-1][i]
                                        if previousGate is not None and \
                                           arg[0] == _invertStrGate(previousGate[0]) and \
                                           arg[1] == previousGate[1] and arg[2] == previousGate[2]:
                                            self.oplines[freeindex-1] = self.oplines[freeindex-1][:i] + [None for j in range(num_targets)] + self.oplines[freeindex-1][i+1:]
                                            if len(self.oplines[freeindex-1]) == 0 or \
                                               all([elem is None or elem == "I"
                                                    for elem in self.oplines[freeindex-1]]):
                                                del self.oplines[freeindex-1]
                                            self.recalculate_free()
                                            skip = True
                                    if not skip:
                                        for party in parties:
                                            self.freeindexes[party] = freeindex + 1
                                        if freeindex > self.lastindex:
                                            self.oplines.append([None for i in range(self.size)])
                                            self.lastindex += 1
                                        self.oplines[freeindex][i] = arg
                                        while num_targets > 1:
                                            del self.oplines[freeindex][i+1]
                                            num_targets -= 1
                    else:
                        print("No gates. Just Monika.")
            else:
                print("Here goes nothing.")
        finally:
            gc.collect()

    '''
    def _flatten(self):
        return QGate._flatten(self)

    def getOptimizedLines(self):
        # print("ti")
        flines = self._flatten()
        olines = [[None for j in range(self.size)] for i in range(len(flines))]
        lastFull = [-1 for i in range(self.size)]
        for i in range(len(flines)):
            line = flines[i]
            for j in range(self.size):
                if line[j] is not None:
                    if not isinstance(line[j], Measure):
                        done = False
                        gate, cons, acons = line[j]
                        items = [j] + list(cons) + list(acons)
                        last = max([lastFull[item] for item in items])
                        if last >= 0 and not isinstance(olines[last], Measure) and olines[last][j] is not None:
                            prevGate, prevCons, prevACons = olines[last][j]
                            if cons == prevCons and acons == prevACons:
                                if prevGate == gate:
                                    if gate == "X" or gate == "Y" or gate == "Z":
                                        olines[last][j] = None
                                        for item in items:
                                            for auxi in range(lastFull[item]-1, 0, -1):
                                                if isinstance(olines[auxi], Measure) or olines[auxi][item] is not None:
                                                    lastFull[item] = auxi
                                        done = True
                                    elif gate == "SqrtX":
                                        olines[last][j] = ["X", cons, acons]
                                        done = True
                                    elif gate == "SqrtSWAP":
                                        olines[last][j] = ["SWAP", cons, acons]
                                        done = True
                        if not done:
                            olines[last+1][j] = line[j]
                            for item in items:
                                lastFull[item] = last + 1
                    else:
                        eline = max(lastFull) + 1
                        lastFull = [eline for k in lastFull]
                        olines[eline] = line
                        break
        olines = olines[:max(lastFull)+1]
        # print("mi")
        return olines
    '''

    def recalculate_free(self):
        for party in range(self.size):
            self.freeindexes[party] = 0
        self.lastindex = -1
        for i in range(len(self.oplines)):
            line = self.oplines[i]
            if isinstance(line[0], Measure):
                self.lastindex = i
                for party in range(self.size):
                    self.freeindexes[party] = i+1
            else:
                parties = _getParties(line)
                if len(parties) > 0:
                    self.lastindex = i
                for party in parties:
                    self.freeindexes[party] = i + 1

    def execute(self, qubits, iterations=1, qmachine=None, args={"useSystem": True}, optimize=True):
        result = None
        lines = self.lines
        if optimize:
            # print("op")
            lines = self.oplines
        if (isinstance(qubits, QRegistry) or isinstance(qubits, QSystem)) and iterations > 1:
            raise ValueError("Can not do more than one iteration with a precreated registry or system!")
        elif (qmachine is None or qmachine == "local"):
            result = qapi.execute(qubits, iterations, lines, self.ancilla, args["useSystem"], optimize)
            # if optimize:
            #     print("zed")
        else:
            raise ValueError("Unsupported qmachine!")
        return result
