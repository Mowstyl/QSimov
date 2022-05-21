"""Module that provides a data structure representing a quantum system.

Data Structures:
    QSystem: Quantum System, preferred over QRegistry (can save a lot of space)

Functions:
    superposition: join two registries into one by calculating tensor product.
"""
import numpy as np
from qsimov.structures.qstructure import QStructure, _get_qubit_set, \
                                        _get_op_data, _get_key_with_defaults
from qsimov.structures.qregistry import QRegistry, superposition


class QSystem(QStructure):
    """Quantum System, preferred over QRegistry (can save a lot of space)."""

    def __init__(self, num_qubits, doki=None, verbose=False):
        """Initialize QSystem to state 0.

        num_qubits -> number of QuBits in the system.
        """
        if doki is None:
            import doki
        self.doki = doki
        if num_qubits is None:
            self.regs = None
            self.qubitMap = None
            self.num_qubits = 0
        else:
            self.regs = [[QRegistry(1, doki=self.doki), [id]]
                         for id in range(num_qubits)]
            self.qubitMap = {id: id for id in range(num_qubits)}
            self.num_qubits = num_qubits
        self.verbose = verbose

    def free(self, deep=False):
        """Release memory held by the QSystem."""
        if self.regs is not None:
            if deep:
                for reg, _ in self.regs:
                    if isinstance(reg, QRegistry):
                        reg.free()
            del self.regs
            del self.qubitMap
            self.regs = None
            self.qubitMap = None

    def clone(self, deep=False):
        """Clone this QSystem."""
        new_sys = QSystem(None, doki=self.doki)
        new_sys.num_qubits = self.num_qubits
        new_sys.qubitMap = {}
        for id in self.qubitMap:
            new_sys.qubitMap[id] = self.qubitMap[id]
        if deep:
            new_sys.regs = [[self.regs[id][0].clone(), self.regs[id][1][:]]
                            if isinstance(self.regs[id][0], QRegistry)
                            else [self.regs[id][0], self.regs[id][1][:]]
                            for id in range(new_sys.num_qubits)]
        else:
            new_sys.regs = [[self.regs[id][0], self.regs[id][1][:]]
                            if isinstance(self.regs[id][0], QRegistry)
                            else [self.regs[id][0], self.regs[id][1][:]]
                            for id in range(new_sys.num_qubits)]
        return new_sys

    def __del__(self):
        """Clean after deletion."""
        self.free()

    def prob(self, id, num_threads=-1):
        """Get the odds of getting 1 when measuring specified qubit."""
        id = _get_qubit_set(self.get_num_qubits(), [id], True, "argument")[0]
        reg, ids = self.regs[self.qubitMap[id]]
        new_id = None
        for i in range(len(ids)):
            if ids[i] == id:
                new_id = i
                break
        if new_id is None:
            raise RuntimeError("Couldn't find id in any reg, " +
                               "please report this bug.")
        return reg.prob(new_id, num_threads=num_threads)

    def get_sizes(self):
        """Return the number of elements of each registry in the system."""
        return ((reg[0].get_state_size(), reg[1])
                if type(reg[0]) == QRegistry
                else (1, reg[1])
                for reg in self.regs)

    def get_state_size(self):
        """Return the number of elements in the state vector of the system."""
        total = 0
        for reg in self.regs:
            if type(reg[0]) == QRegistry:
                total += reg[0].get_state_size()
            else:
                total += 1
        return total

    def get_split_num_qubits(self):
        """Return the number of qubits in each registry of the system."""
        return (reg[0].get_num_qubits()
                if type(reg[0]) == QRegistry
                else 1  # When we measure with remove=True
                for reg in self.regs)

    def get_num_qubits(self):
        """Return the number of qubits in this system."""
        return self.num_qubits

    def measure(self, ids, random_generator=np.random.rand,
                num_threads=-1, deep=False):
        """Measure specified qubits of this system and collapse.

        Positional arguments:
            ids -> List of QuBit ids that have to be measured
        Keyworded arguments:
            random_generator -> function without arguments that returns
                                a random real number in [0, 1)
        Return:
            List with the value obtained after each measure
        """
        num_qubits = self.get_num_qubits()
        ids = _get_qubit_set(num_qubits, ids, False, "ids")
        if ids is None:
            raise ValueError("ids cannot be None")
        split_ids = {reg_id: set() for reg_id in range(len(self.regs))}
        for qubit_id in ids:
            reg_id = self.qubitMap[qubit_id]
            split_ids[reg_id].add(qubit_id)
        # In split ids we have reg_id -> set of ids to measure in that reg
        split_ids = {k: v for k, v in split_ids.items() if len(v) > 0}
        result = [None for i in range(num_qubits)]
        # Here we have the registries that have not been used
        untouched_regs = {i for i in range(len(self.regs))
                          if i not in split_ids}
        # We create a new QSystem with the regs that have not been used
        new_sys = QSystem(None, doki=self.doki)
        new_sys.regs = []
        new_sys.qubitMap = {}
        exception = None
        try:
            for reg_id in untouched_regs:
                reggie, reg_ids = self.regs[reg_id]
                if deep:
                    reggie = reggie.clone()
                    reg_ids = reg_ids[:]
                new_sys.regs.append((reggie, reg_ids))
            for reg_id in range(len(untouched_regs)):
                for qubit_id in new_sys.regs[reg_id][1]:
                    new_sys.qubitMap[qubit_id] = reg_id
            new_sys.num_qubits = self.num_qubits

            # We iterate through the registries that have a qubit in ids
            for reg_id in split_ids:
                partial_ids = split_ids[reg_id]  # ids of QSystem to measure
                new_reg = None
                partial_result = None
                reg, reg_ids = self.regs[reg_id]
                # Ids to measure in the QRegistry (not in the whole QSystem)
                # mapped to the id in the QSystem
                new_ids = {i: reg_ids[i] for i in range(len(reg_ids))
                           if reg_ids[i] in partial_ids}
                # Not measured ids in this registry
                new_reg_ids = [id for id in reg_ids if id not in partial_ids]
                # We measure registries
                aux = reg.measure(new_ids.keys(),
                                  random_generator=random_generator,
                                  num_threads=num_threads)
                new_reg, partial_result = aux
                # We add the results to the result list
                for local_id in new_ids:
                    result[new_ids[local_id]] = partial_result[local_id]
                # We add the new registry to the list of regs
                if len(new_reg_ids) > 0:
                    aux_id = len(new_sys.regs)
                    new_sys.regs.append((new_reg, reg_ids))
                    for id in new_reg_ids:
                        new_sys.qubitMap[id] = aux_id
                else:
                    new_reg.free()
                # We add new registries with only the measured qubits
                for id in partial_ids:
                    one_reg = QRegistry(1, doki=self.doki,
                                        verbose=self.verbose)
                    if result[id]:
                        one_aux = one_reg.apply_gate("X")
                        one_reg.free()
                        one_reg = one_aux
                    new_sys.regs.append((one_reg, [id]))
                    new_sys.qubitMap[id] = len(new_sys.regs) - 1
        except Exception as ex:
            exception = ex
        if exception is not None:
            del new_sys
            raise exception
        return (new_sys, result)

    def as_qregistry(self, num_threads=-1, canonical=True):
        """Return this system as a QRegistry."""
        aux_reg = None
        new_reg = None
        new_ids = []
        first = True
        for reg_id in range(len(self.regs)):
            reg, ids = self.regs[reg_id]
            if new_reg is None:
                new_reg = reg
            else:
                aux_reg = superposition(new_reg, reg,
                                        num_threads=num_threads,
                                        verbose=self.verbose)
            new_ids = ids + new_ids
            if aux_reg is not None:
                if not first:
                    del new_reg
                first = False
                new_reg = aux_reg
                aux_reg = None
        # Here we remove the unused ids
        q_ids = [id for id in new_ids if new_reg.get_classic(id) is None]
        swap_ids = np.argsort(np.argsort(q_ids))
        # And we sort the remaining qubits by qubit_id
        for i in range(len(swap_ids)):
            while swap_ids[i] != i:
                swap_targets = [swap_ids[i], swap_ids[swap_ids[i]]]
                swap_ids[swap_targets[0]], swap_ids[i] = swap_targets
                aux_reg = new_reg.apply_gate("SWAP",
                                             targets=[i, swap_targets[0]],
                                             num_threads=num_threads)
                if not first:
                    del new_reg
                new_reg = aux_reg
        return new_reg

    def get_state(self, key=None, canonical=False):
        return self.as_qregistry().get_state(key=key, canonical=canonical)

    def get_classic(self, id):
        """Return classic bit value."""
        return None

    def apply_gate(self, gate, targets=None, controls=None, anticontrols=None,
                   num_threads=-1, deep=False):
        """Apply specified gate to specified qubit with specified controls.

        Positional arguments:
            gate: string with the name of the gate to apply, or a QGate
        Keyworded arguments:
            targets: id of the least significant qubit the gate will target
            controls: id or list of ids of the qubit that will act as
                     controls
            anticontrols: id or list of ids of the qubit that will act as
                         anticontrols
            num_threads: number of threads to use
            optimize: only for QGates. Whether to use optimized lines or
                      user defined lines
        """
        if not np.allclose(num_threads % 1, 0):
            raise ValueError("num_threads must be an integer")
        num_threads = int(num_threads)
        num_qubits = self.get_num_qubits()
        op_data = _get_op_data(num_qubits, 0, gate, targets, None, None,
                               controls, anticontrols, None, None)
        gate = op_data["gate"]
        targets = op_data["targets"]
        controls = op_data["controls"]
        anticontrols = op_data["anticontrols"]
        # We create a new system without the data of the parties
        new_sys = QSystem(None, doki=self.doki)
        new_sys.regs = []
        new_reg = None
        aux_reg = None
        exception = None
        try:
            # If any of the affected qubits is marked as not usable
            if any([self.get_classic(qubit_id) is not None
                    for qubit_id in targets]):
                # we raise an exception
                raise ValueError("Trying to apply gate to classic bit")
            cfail = any([self.get_classic(qubit_id) is False
                         for qubit_id in controls])
            acfail = any([self.get_classic(qubit_id) is True
                          for qubit_id in anticontrols])
            if cfail or acfail:
                if deep:
                    return self.clone(deep=True)
                else:
                    return self
            controls = {qubit_id for qubit_id in controls
                        if self.get_classic(qubit_id) is None}
            anticontrols = {qubit_id for qubit_id in anticontrols
                            if self.get_classic(qubit_id) is None}
            # All affected qubits
            parties = controls.union(anticontrols).union(targets)
            touched_regs = {self.qubitMap[qubit_id]
                            for qubit_id in parties}
            for reg_id in range(len(self.regs)):
                if reg_id not in touched_regs:
                    reggie, reg_ideses = self.regs[reg_id]
                    if deep and isinstance(reggie, QRegistry):
                        reggie = reggie.clone()
                    new_sys.regs.append([reggie, reg_ideses[:]])
            # Create new qubit map
            new_sys.qubitMap = {}
            for reg_id in range(len(new_sys.regs)):
                for qubit_id in new_sys.regs[reg_id][1]:
                    new_sys.qubitMap[qubit_id] = reg_id
            new_sys.num_qubits = self.num_qubits
            new_ids = []
            merged = False
            for reg_id in touched_regs:
                curr_reg, curr_ids = self.regs[reg_id]
                if new_reg is not None:
                    aux_reg = superposition(curr_reg, new_reg,
                                            num_threads=num_threads,
                                            verbose=self.verbose)
                    if merged:
                        del new_reg
                    else:
                        merged = True
                    new_reg = aux_reg
                else:
                    new_reg = curr_reg
                new_ids += curr_ids
            inverse_map = {new_ids[qubit_id]: qubit_id
                           for qubit_id in range(len(new_ids))}
            mapped_targets = [inverse_map[qubit_id]
                              for qubit_id in targets]
            mapped_controls = {inverse_map[qubit_id]
                               for qubit_id in controls}
            mapped_anticontrols = {inverse_map[qubit_id]
                                   for qubit_id in anticontrols}
            aux_reg = new_reg.apply_gate(gate, targets=mapped_targets,
                                         controls=mapped_controls,
                                         anticontrols=mapped_anticontrols,
                                         num_threads=num_threads)
            if merged:
                del new_reg
                new_reg = None
            new_sys.regs.append([aux_reg, new_ids])
            for id in new_ids:
                new_sys.qubitMap[id] = len(new_sys.regs) - 1
        except Exception as ex:
            if new_sys is not None:
                del new_sys
            if new_reg is not None and merged:
                del new_reg
            if aux_reg is not None:
                del aux_reg
            new_sys = None
            exception = ex
        if exception is not None:
            raise exception
        return new_sys

    def get_bloch_coords(self, key=None):
        """Get the polar coordinates of all ONE qubit registries."""
        start, stop, step = _get_key_with_defaults(key, self.num_qubits,
                                                   0, self.num_qubits, 1)
        coords = [None for id in range(start, stop, step)]
        for id in range(start, stop, step):
            try:
                reg, ids = self.regs[self.qubitMap[id]]
                new_id = ids.index(id)
                coords[id] = reg.get_bloch_coords(new_id)
            except Exception:
                pass
        if key is not None and type(key) != slice:
            coords = coords[0]
        return coords

    def bloch(self, key=None):
        """Return matplotlib bloch sphere."""
        all_coords = self.get_bloch_coords(key=key)
        if type(all_coords) != list:
            all_coords = [all_coords]
        from qsimov.utils.bloch import draw_bloch_sphere
        figs = [draw_bloch_sphere(coords[0], coords[1]) if coords is not None
                else None
                for coords in all_coords]
        if key is not None and type(key) != slice:
            figs = figs[0]
        return figs


def join_systems(most, least, deep=False):
    """Return a system that contains both a and b systems."""
    res = QSystem(None, doki=most.doki)
    res.regs = []
    res.qubitMap = {}
    exception = None
    try:
        count = 0
        for reg, ids in least.regs:
            new_reg = reg
            if reg == QRegistry:
                if deep:
                    new_reg = reg.clone()
            count += 1
            res.regs.append([new_reg, ids[:]])
        offset = least.get_num_qubits()
        for reg, ids in most.regs:
            new_reg = reg
            if reg == QRegistry:
                if deep:
                    new_reg = reg.clone()
            count += 1
            res.regs.append([new_reg, [id + offset for id in ids]])
        for i in range(len(res.regs)):
            _, ids = res.regs[i]
            for qubit_id in ids:
                res.qubitMap[qubit_id] = i
    except Exception as ex:
        exception = ex
    if exception is not None:
        del res
        raise exception
    return res
