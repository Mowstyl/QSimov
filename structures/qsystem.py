import numpy as np
from structures.qregistry import QRegistry, superposition
from structures.qgate import QGate, getGateData

class QSystem:
    def __init__(self, nqbits):
        # nqbits -> number of QuBits in the registry.
        self.regs = [[QRegistry(1), [id]] for id in range(nqbits)]
        self.qubitMap = {id: id for id in range(nqbits)}
        self.nqubits = nqbits

    def __del__(self):
        for i in range(len(self.regs)):
            if self.regs[i] != None:
                del self.regs[i][0]
        del self.regs
        del self.qubitMap

    def getRegSize(self):
        return ((reg[0].getSize(), reg[1]) if type(reg[0]) == QRegistry else (1, reg[1]) for reg in self.regs)

    def getSize(self):
        total = 0
        for reg in self.regs:
            if type(reg[0]) == QRegistry:
                total += reg[0].getSize()
            else:
                total += 1
        return total

    def getRegNQubits(self):
        return (reg[0].getNQubits() if type(reg[0]) == QRegistry else 1 for reg in self.regs)

    def getNQubits(self):
        return self.nqubits

    def toString(self):
        return str([[reg[0].toString(), reg[1]] if type(reg[0]) == QRegistry else reg for reg in self.regs])

    def measure(self, msk, remove=False): # List of numbers with the QuBits that should be measured. 0 means not measuring that qubit, 1 otherwise. remove = True if you want to remove a QuBit from the registry after measuring
        nqubits = self.getNQubits()
        if (type(msk) != list or len(msk) != nqubits or \
            not all(type(num) == int and (num == 0 or num == 1) for num in msk)):
            raise ValueError('Not valid mask')
        result = []
        for qubit in range(self.nqubits):
            if msk[qubit] == 1:
                regid = self.qubitMap[qubit]
                if self.regs[regid][0].getNQubits() > 1:
                    result += self.regs[regid][0].measure([0 if self.regs[regid][1][i] != qubit else 1 for i in range(len(self.regs[regid][1]))], remove=True)
                    self.regs[regid][1].remove(qubit)
                    newid = self.regs.index(None)
                    self.regs[newid] = [QRegistry(1), [qubit]]
                    self.qubitMap[qubit] = newid
                    if result[-1] == 1:
                        self.regs[newid][0].applyGate("X")
                else:
                    if not remove:
                        result += self.regs[regid][0].measure([1], remove=False)
                    else:
                        if type(self.regs[regid][0]) == QRegistry:
                            result += self.regs[regid][0].measure([1], remove=False)
                            del self.regs[regid][0]
                            self.regs[regid] = [result[qubit], self.regs[regid][0][0]]
                        else:
                            result.append(self.regs[regid])
        return result

    def applyGate(self, *args, **kwargs): # gate, qubit=0, control=None, anticontrol=None):
        if len(args) == 1 or (len(args) == 2 and "qubit" not in kwargs):
            for key in kwargs:
                if key != "qubit" and key != "control" and key != "anticontrol":
                    print('Apart from the gates, you can only specify "qubit", "control" and/or "anticontrol" (lowercase)')
                    return
            gate = args[0]
            if type(gate) == QGate:
                gate = (gate, gate.getNQubits())
            elif gate == "I" or gate is None:
                return
            else:
                gate = (gate, getGateData(gate)[0])
            qubit = kwargs.get("qubit", 0)
            if len(args) == 2:
                qubit = args[1]
            control = kwargs.get("control", [])
            if type(control) != list:
                control = [control]
            anticontrol = kwargs.get("anticontrol", [])
            if type(anticontrol) != list:
                anticontrol = [anticontrol]
            if type(qubit) == list or type(qubit) == set:
                for qid in qubit:
                    self.applyGate(gate[0], qubit=qid, control=control, anticontrol=anticontrol)
            else:
                qubits = set([qubit + i for i in range(1, gate[1])] + control + anticontrol) # Todos los participantes
                rid = self.qubitMap[qubit]
                for qid in qubits:
                    self.superposition(rid, self.qubitMap[qid])
                reg, idlist = self.regs[rid]
                regmap = {idlist[i]: i for i in range(len(idlist))}
                newqubit = regmap[qubit]
                newcontrol = [regmap[qid] for qid in control]
                newanticontrol = [regmap[qid] for qid in anticontrol]
                if type(gate[0]) == QGate:
                    pass
                else:
                    reg.applyGate(gate[0], qubit=newqubit, control=newcontrol, anticontrol=newanticontrol)
        elif len(kwargs) == 0 and len(args) > 0:
            nq = 0
            gates = []
            for arg in args:
                if type(arg) == QGate:
                    gatenq = (arg, arg.getNQubits())
                elif arg == "I" or arg is None:
                    gatenq = (None, 1)
                else:
                    gatenq = (arg, getGateData(arg)[0])
                nq += gatenq[1]
                gates.append(gatenq)
            if nq == self.getNQubits():
                qid = 0
                for gate in gates:
                    if gate[0] != None:
                        self.applyGate(gate[0], qubit=qid)
                    qid += gate[1]
            else:
                print ("You have to specify a gate for each QuBit (or None if you don't want to operate with it)")
        elif len(args) == 0:
            print("You must specify at least one gate")
        else:
            print("You can't apply more than one gate when using " + '"qubit", "control" or "anticontrol"')

    def hopfCoords(self):
        return [self.regs[self.qubitMap[id]][0].hopfCoords() if type(self.regs[self.qubitMap[id]][0]) == QRegistry else None for id in range(self.nqubits)]

    def blochCoords(self):
        return [self.regs[self.qubitMap[id]][0].blochCoords() if type(self.regs[self.qubitMap[id]][0]) == QRegistry else None for id in range(self.nqubits)]

    def superposition(self, regid1, regid2):
        if regid1 != regid2:
            a = self.regs[regid1]
            b = self.regs[regid2]
            newregdata = joinRegs(a, b)
            self.regs[regid1] = newregdata
            self.regs[regid2] = None
            for id in b[1]:
                self.qubitMap[id] = regid1
            del a[0]
            del b[0]

    def getState(self):
        regs = list(filter(None, self.regs))
        if len(regs) == 1:
            return regs[0][0].getState()
        else:
            reg = joinRegs(regs[0], regs[1])
            for i in range(2, len(regs)):
                if regs[i] != None:
                    aux = reg
                    reg = joinRegs(aux, regs[i])
                    del aux[0]
            state = reg[0].getState()
            del reg[0]
            return state

def joinRegs(a, b):
    newregdata = []
    if b[1][0] > a[1][-1]:
        newregdata = [superposition(b[0], a[0]), a[1] + b[1]]
    elif a[1][0] > b[1][-1]:
        newregdata = [superposition(a[0], b[0]), b[1] + a[1]]
    else:
        newregdata = [superposition(b[0], a[0]), a[1] + b[1]]
        newregdata = sortRegdata(newregdata)
    return newregdata

def sortRegdata(regdata): # regdata[1] = [1,3,2,4] -> el qubit 0 del registro es el 1 del sistema
    relen = len(regdata[1])
    for i in range(relen-1):
        aux = regdata[1][i:]
        minid = min(aux)
        minindex = aux.index(minid)
        if aux[0] != minid:
            regdata[0].applyGate("SWAP(" + str(i) + "," + str(minindex) + ")")
            regdata[1][i], regdata[1][minindex] = minid, regdata[1][i]
    return regdata
