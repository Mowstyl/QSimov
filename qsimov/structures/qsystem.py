"""Module that provides a data structure representing a quantum system.

Data Structures:
    QSystem: Quantum System, preferred over QRegistry (can save a lot of space)

Functions:
    superposition: join two registries into one by calculating tensor product.
"""
import qsimov.connectors.parser as prs
import numpy as np
from qsimov.structures.qregistry import QRegistry, superposition
from qsimov.structures.qgate import QGate, get_gate_data
from collections.abc import Iterable


class QSystem:
    """Quantum System, preferred over QRegistry (can save a lot of space)."""

    def __init__(self, nqbits):
        """Initialize QSystem to state 0.

        nqbits -> number of QuBits in the system.
        """
        self.regs = [[QRegistry(1), [id]] for id in range(nqbits)]
        self.qubitMap = {id: id for id in range(nqbits)}
        self.usable = [id for id in range(nqbits)]
        self.nqubits = nqbits

    def __del__(self):
        """Release memory held by the QSystem."""
        del self.regs
        del self.qubitMap
        del self.usable

    def getRegSize(self):
        """Use get_sizes method instead. DEPRECATED."""
        print("Method QSystem.getRegSize is deprecated.",
              "Please use get_sizes if you seek the same functionality")
        return self.get_sizes()

    def get_sizes(self):
        """Return the number of elements of each registry in the system."""
        return ((reg[0].get_state_size(), reg[1])
                if type(reg[0]) == QRegistry
                else (1, reg[1])
                for reg in self.regs)

    def getSize(self):
        """Use get_state_size method instead. DEPRECATED."""
        print("Method QSystem.getSize is deprecated.",
              "Please use get_state_size if you seek the same functionality")
        return self.get_state_size()

    def get_state_size(self):
        """Return the number of elements in the state vector of the system."""
        total = 0
        for reg in self.regs:
            if type(reg[0]) == QRegistry:
                total += reg[0].get_state_size()
            else:
                total += 1
        return total

    def getRegNQubits(self):
        """Use get_split_num_qubits method instead. DEPRECATED."""
        print("Method QSystem.getRegNQubits is deprecated.",
              "Please use get_split_num_qubits if you",
              "seek the same functionality")
        return self.get_split_num_qubits()

    def get_split_num_qubits(self):
        """Return the number of qubits in each registry of the system."""
        return (reg[0].get_num_qubits()
                if type(reg[0]) == QRegistry
                else 1  # When we measure with remove=True
                for reg in self.regs)

    def getNQubits(self):
        """Use get_num_qubits method instead. DEPRECATED."""
        print("Method QSystem.getNQubits is deprecated.",
              "Please use get_num_qubits if you seek the same functionality")
        return self.get_num_qubits()

    def get_num_qubits(self):
        """Return the number of qubits in this system."""
        return self.nqubits

    def toString(self):
        """Use state_strings method instead. DEPRECATED."""
        print("Method QSystem.toString is deprecated.",
              "Please use state_strings if you seek the same functionality")
        return self.state_strings()

    def state_strings(self):
        """Return string representation of each registry in this system."""
        print("state_strings might be removed in future versions")
        return str([[reg[0].state_string(), reg[1]]
                    if type(reg[0]) == QRegistry
                    else reg  # When we measure with remove=True
                    for reg in self.regs])

    def measure(self, msk, remove=False):
        """Measure specified qubits of this system and collapse.

        Positional arguments:
            msk -> List of numbers with the qubits that should be measured
                0 means not measuring that qubit, 1 otherwise
        Keyworded arguments:
            remove = True if you want to turn the registry containing only
                the measured qubit into an int with the measured value

        Return:
            List with the value obtained after each measure
        """
        nqubits = self.get_num_qubits()
        if (not isinstance(msk, Iterable) or len(msk) != nqubits or
                not all(type(num) == int and (num == 0 or num == 1)
                        for num in msk)):
            raise ValueError('Not valid mask: ' + str(msk))
        result = []
        for qubit in range(self.nqubits):
            if msk[qubit] == 1:
                regid = self.qubitMap[qubit]
                if self.regs[regid][0].get_num_qubits() > 1:
                    aux_mask = [0 if self.regs[regid][1][i] != qubit
                                else 1
                                for i in range(len(self.regs[regid][1]))]
                    result += self.regs[regid][0].measure(aux_mask,
                                                          remove=True)
                    self.regs[regid][1].remove(qubit)
                    newid = self.regs.index(None)
                    self.qubitMap[qubit] = newid
                    if not remove:
                        self.regs[newid] = [QRegistry(1), [qubit]]
                        if result[-1] == 1:
                            self.regs[newid][0].apply_gate("X")
                    else:
                        self.regs[newid] = [result[-1], [qubit]]
                        self.usable.remove(qubit)
                else:
                    if not remove:
                        result += self.regs[regid][0].measure([1],
                                                              remove=False)
                    else:
                        if type(self.regs[regid][0]) == QRegistry:
                            result += self.regs[regid][0].measure([1],
                                                                  remove=False)
                            del self.regs[regid][0]
                            self.regs[regid] = [result[-1],
                                                self.regs[regid][0][0]]
                        else:
                            result.append(self.regs[regid])
                        self.usable.remove(qubit)
        return result

    def applyGate(self, *args, **kwargs):
        """Use apply_gate method instead. DEPRECATED."""
        print("Method QSystem.applyGate is deprecated.",
              "Please use apply_gate if you seek the same functionality")
        return self.apply_gate(*args, **kwargs)

    def apply_gate(self, *args, **kwargs):
        """Apply specified gate to specified qubit with specified controls.

        There are two variants for this method, depending on the number of
        gates you want to apply.
        One gate:
            Positional arguments:
                gate: string with the name of the gate to apply, or a QGate
            Keyworded arguments:
                qubit: id of the least significant qubit the gate will target
                control: id or list of ids of the qubit that will act as
                         controls
                anticontrol: id or list of ids of the qubit that will act as
                             anticontrols
                optimize: only for QGates. Whether to use optimized lines or
                          user defined lines
        Multiple gates:
            Positional arguments:
                comma separated gates, their sizes must match the number of
                    qubits in the system. Sorted by their least significant
                    target qubit id.
        """
        optimize = True
        if "optimize" in kwargs:
            optimize = kwargs["optimize"]
        if (len(args) == 1
                or (len(args) == 2
                    and np.issubdtype(type(args[1]), np.integer)
                    and "qubit" not in kwargs)):
            for key in kwargs:
                if (key != "qubit" and key != "control"
                        and key != "anticontrol" and key != "optimize"):
                    raise ValueError('Apart from the gates, you can only ' +
                                     'specify "qubit", "control", ' +
                                     '"anticontrol" (lowercase) and/or ' +
                                     '"optimize"')
            gate = args[0]
            if type(gate) == QGate:
                gate = (gate, gate.size)
            elif gate == "I" or gate is None:
                return
            else:
                gate = (gate, get_gate_data(gate)[0])
            qubit = kwargs.get("qubit", self.usable[0])
            if len(args) == 2:
                qubit = args[1]
            control = kwargs.get("control", [])
            if control is None:
                control = []
            if not isinstance(control, Iterable):
                control = [control]
            anticontrol = kwargs.get("anticontrol", [])
            if anticontrol is None:
                anticontrol = []
            if not isinstance(anticontrol, Iterable):
                anticontrol = [anticontrol]
            if isinstance(qubit, Iterable):
                for qid in qubit:
                    self.apply_gate(gate[0], qubit=qid, control=control,
                                    anticontrol=anticontrol, optimize=optimize)
            else:
                name = ""
                if type(gate[0]) != QGate:
                    name, arg1, arg2, arg3, invert = prs.get_gate_data(gate[0])
                    invstring = ""
                    if invert:
                        invstring = "-1"
                qubits = set(control).union(anticontrol)  # All affected qubits
                if "SWAP" in name:
                    qubit = arg1
                    qubits.add(arg2)
                elif name == "XX" or name == "YY" or name == "ZZ":
                    qubit = arg2
                    qubits.add(arg3)
                else:
                    qubits.update([qubit + i for i in range(1, gate[1])])
                rid = self.qubitMap[qubit]
                for qid in qubits:
                    self._superposition(rid, self.qubitMap[qid])
                reg, idlist = self.regs[rid]
                regmap = {idlist[i]: i for i in range(len(idlist))}
                newqubit = regmap[qubit]
                newcontrol = [regmap[qid] for qid in control]
                newanticontrol = [regmap[qid] for qid in anticontrol]
                if type(gate[0]) == QGate:
                    reg.apply_gate(gate[0], qubit=newqubit, control=newcontrol,
                                   anticontrol=newanticontrol,
                                   optimize=optimize)
                else:
                    if "SWAP" in name:
                        gate = (name + "(" + str(regmap[arg1]) + "," +
                                str(regmap[arg2]) + ")" + invstring, gate[1])
                    if name == "XX" or name == "YY" or name == "ZZ":
                        gate = (name + "(" + str(arg1) + "," +
                                str(regmap[arg2]) + "," + str(regmap[arg3]) +
                                ")" + invstring, gate[1])
                    reg.apply_gate(gate[0], qubit=newqubit, control=newcontrol,
                                   anticontrol=newanticontrol,
                                   optimize=optimize)
        elif len(kwargs) == 0 and len(args) > 0:
            nq = 0
            gates = []
            for arg in args:
                if type(arg) == QGate:
                    gatenq = (arg, arg.size)
                elif arg == "I" or arg is None:
                    gatenq = (None, 1)
                else:
                    gatenq = (arg, get_gate_data(arg)[0])
                nq += gatenq[1]
                gates.append(gatenq)
            if nq == self.get_num_qubits():
                qid = 0
                for gate in gates:
                    if gate[0] is not None:
                        self.apply_gate(gate[0], qubit=qid, optimize=optimize)
                    qid += gate[1]
            else:
                print("You have to specify a gate for each QuBit",
                      "(or None if you don't want to operate with it)")
        elif len(args) == 0:
            print("You must specify at least one gate")
        else:
            print("You can't apply more than one gate when using",
                  '"qubit", "control" or "anticontrol"')

    def blochCoords(self):
        """Use get_bloch_coords method instead. DEPRECATED."""
        print("Method QRegistry.blochCoords is deprecated.",
              "Please use get_bloch_coords if you seek the same functionality")
        return self.get_bloch_coords()

    def get_bloch_coords(self):
        """Get the polar coordinates of all ONE qubit registries."""
        return [self.regs[self.qubitMap[id]][0].get_bloch_coords()
                if (type(self.regs[self.qubitMap[id]][0]) == QRegistry
                    and len(self.regs[self.qubitMap[id]][1]) == 1)
                else None
                for id in range(self.nqubits)]

    def _superposition(self, regid1, regid2):
        """Join two registries into one by calculating tensor product."""
        if regid1 != regid2:
            a = self.regs[regid1]
            b = self.regs[regid2]
            newregdata = join_regs(a, b)
            self.regs[regid1] = newregdata
            self.regs[regid2] = None
            for id in b[1]:
                self.qubitMap[id] = regid1
            del a[0]
            del b[0]

    def getState(self):
        """Use get_state method instead. DEPRECATED."""
        print("Method QSystem.getState is deprecated.",
              "Please use get_state if you seek the same functionality")
        return self.get_state()

    def get_state(self):
        """Return the state vector of the system."""
        regs = list(filter(lambda reg: reg is not None
                           and type(reg[0]) == QRegistry,
                           self.regs))
        if len(regs) == 1:
            return regs[0][0].get_state()
        else:
            reg = join_regs(regs[0], regs[1])
            for i in range(2, len(regs)):
                if regs[i] is not None:
                    aux = reg
                    reg = join_regs(aux, regs[i])
                    del aux[0]
            state = reg[0].get_state()
            del reg[0]
            return state


def joinSystems(a, b):
    """Use join_systems method instead. DEPRECATED."""
    print("Method QSystem.joinSystems is deprecated.",
          "Please use join_systems if you seek the same functionality")
    return join_systems(a, b)


def join_systems(a, b):
    """Return a system that contains both a and b systems."""
    res = QSystem(a.nqubits + b.nqubits)

    for i in range(len(res.regs)):
        if res.regs[i] is not None:
            del res.regs[i][0]
    res.regs.clear()
    res.qubitMap.clear()
    res.regs = [[reg[0], reg[1][:]] if reg is not None
                else None for reg in a.regs]
    offset = a.nqubits
    res.regs += [[reg[0], [id + offset for id in reg[1]]] if reg is not None
                 else None for reg in b.regs]
    res.qubitMap = a.qubitMap.copy()
    res.qubitMap.update({k + offset: b.qubitMap[k] + offset
                         for k in b.qubitMap})

    return res


def join_regs(a, b):
    """Join two registries and sort them by qubit id (descending order).

    Positional arguments:
        a: pair of QRegistry and list of qubit ids in that registry
        b: pair of QRegistry and list of qubit ids in that registry
    Return:
        pair of QRegistry and list of qubit ids in that registry
            result of joining a and b
    """
    newregdata = []
    # We assume a and b are already sorted
    if b[1][0] > a[1][-1]:
        # If the last qubit of a is less than the first of b
        # we calculate b tensor product a
        newregdata = [superposition(b[0], a[0]), a[1] + b[1]]
    elif a[1][0] > b[1][-1]:
        # If the last qubit of b is less than the first of a
        # we calculate a tensor product b
        newregdata = [superposition(a[0], b[0]), b[1] + a[1]]
    else:  # En caso contrario
        # Otherwise we calculate b tensor product a
        newregdata = [superposition(b[0], a[0]), a[1] + b[1]]
        # and we sort them
        newregdata = sort_reg_data(newregdata)
    return newregdata


def sort_reg_data(reg_data):
    """Sort a registry by qubit id (descending order).

    Positional arguments:
        reg_data: pair of QRegistry and list of qubit ids in that registry
    Return:
        reg_data after sorting
    """
    reg_len = len(reg_data[1])  # Number of qubits in the registry
    for i in range(reg_len-1):
        # With n elements, we don't need the last iteration.
        # The last element will already be in the last position
        unsorted = reg_data[1][i:]  # Elements still not sorted
        min_id = min(unsorted)  # Smaller qubit id still not sorted
        min_index = unsorted.index(min_id) + i  # Index of id in the registry
        if unsorted[0] != min_id:  # If it is not yet in the first position
            # We SWAP it with the qubit in that position
            reg_data[0].apply_gate("SWAP(" + str(i) +
                                   "," + str(min_index) + ")")
            # and we update the list of qubit ids accordingly
            reg_data[1][i], reg_data[1][min_index] = min_id, reg_data[1][i]
    return reg_data
