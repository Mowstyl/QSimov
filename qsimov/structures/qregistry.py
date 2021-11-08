"""Module that provides the lowest level data structure: QRegistry.

Data Structures:
    QRegistry: Quantum Registry, base of all quantum related operations

Functions:
    superposition: join two registries into one by calculating tensor product.
"""
import numpy as np

from qsimov.structures.funmatrix import Funmatrix
from qsimov.structures.simple_gate import SimpleGate
from qsimov.structures.qstructure import QStructure, _get_op_data, \
                                         _get_qubit_set
from qsimov.connectors.qsimovapi import apply_design


class QRegistry(QStructure):
    """Quantum Registry, base of all quantum related operations."""

    def __init__(self, num_qubits, doki=None, verbose=False):
        """Initialize QRegistry to state 0.

        num_qubits -> number of QuBits in the registry.
        """
        if doki is None:
            import doki
        self.doki = doki
        if num_qubits is not None:
            if num_qubits > 0:
                self.reg = self.doki.registry_new(num_qubits, verbose)
                self.qubit_map = {i: i for i in range(num_qubits)}
                self.classic_vals = {}
                self.num_qubits = num_qubits
                self.num_bits = 0
                self.size = 2**num_qubits
            else:
                raise ValueError("num_qubits must be greater than 0")
        else:
            self.reg = None
            self.num_qubits = 0
            self.num_bits = 0
            self.size = 0
            self.qubit_map = {}
            self.classic_vals = {}
        self.verbose = verbose

    def __del__(self):
        """Clean after deletion."""
        self.free()

    def get_classic(self, id):
        """Return classic bit value (if qubit has been measured)."""
        if id in self.classic_vals:
            return self.classic_vals[id]
        return None

    def free(self):
        """Release memory held by the QRegistry."""
        if self.reg is not None:
            self.doki.registry_del(self.reg, self.verbose)
            self.reg = None
            self.num_qubits = 0
            self.size = 0

    def get_state_size(self):
        """Return the number of elements in the state vector."""
        return self.size

    def get_num_qubits(self):
        """Return the number of qubits in this registry."""
        return self.num_qubits

    def measure(self, ids, random_generator=np.random.rand, num_threads=-1):
        """Measure specified qubits of this registry and collapse.

        Positional arguments:
            ids -> List of QuBit ids that have to be measured
        Keyworded arguments:
            random_generator -> function without arguments that returns
                                a random real number in [0, 1)
        Return:
            List with the value obtained after each measure
        """
        ids_set = _get_qubit_set(self.num_qubits + self.num_bits,
                                 ids, False, "ids")
        if len(ids_set) == 0:
            raise ValueError('ids cannot be empty')
        num_measured = len(ids_set)
        mask = 0
        rolls = [random_generator() for i in range(num_measured)]
        for id in ids_set:
            if id not in self.qubit_map:
                raise ValueError(f"Id {id} has already been measured")
            mask += 2**self.qubit_map[id]
        raw_reg, mess = self.doki.registry_measure(self.reg, mask, rolls,
                                                   num_threads, self.verbose)
        mess = mess[::-1]
        final_mess = [mess[self.qubit_map[id]] if id in ids_set
                      else None
                      for id in range(self.num_qubits + self.num_bits)]
        rem_qubits = self.num_qubits - num_measured
        if raw_reg is not None:
            # If rem_qubits is 0 then raw_reg must be None
            if rem_qubits == 0:
                raise RuntimeError("Zero qubit registry found while " +
                                   "measuring! Please report this bug.")
        new_reg = QRegistry(None, doki=self.doki)
        new_reg.num_qubits = rem_qubits
        new_reg.num_bits = self.num_bits + num_measured
        new_reg.size = 2**rem_qubits
        new_reg.reg = raw_reg
        new_reg.verbose = self.verbose
        # Copying maps
        old_keys = [k for k in self.qubit_map if k not in ids_set]
        # print("DEBUG -> old_keys:", old_keys)
        new_reg.qubit_map = {old_keys[i]: i
                             for i in range(len(old_keys))}
        num_total = new_reg.num_qubits + new_reg.num_bits
        new_reg.classic_vals = {id: self.classic_vals[id]
                                if id in self.classic_vals
                                else final_mess[id]
                                for id in range(num_total)
                                if id not in new_reg.qubit_map}

        return (new_reg, final_mess)

    def apply_gate(self, gate, targets=None, controls=None, anticontrols=None,
                   num_threads=-1):
        """Add specified gate to this QGate.

        Positional arguments:
            gate: string with the name of the gate to apply, or a QGate
        Keyworded arguments:
            targets: id or list of ids of the qubits the gate will target
            controls: id or set of ids of the qubit that will act as controls
            anticontrols: id or set of ids of the qubit that will act as
                         anticontrols
            num_threads: number of threads to use
        Return:
            A new QRegistry
        """
        if not np.allclose(num_threads % 1, 0):
            raise ValueError("num_threads must be an integer")
        num_threads = int(num_threads)
        num_qubits = self.num_qubits + self.num_bits
        op_data = _get_op_data(num_qubits, gate, targets,
                               controls, anticontrols)
        gate = op_data["gate"]
        targets = op_data["targets"]
        controls = op_data["controls"]
        anticontrols = op_data["anticontrols"]
        classic_controls = {qubit_id for qubit_id in controls
                            if qubit_id in self.classic_vals}
        classic_anticontrols = {qubit_id for qubit_id in anticontrols
                                if qubit_id in self.classic_vals}
        ccheck = all(self.classic_vals[id] for id in classic_controls)
        accheck = not any(self.classic_vals[id] for id in classic_anticontrols)
        if ((len(classic_controls) > 0 and not ccheck)
                or (len(classic_anticontrols) > 0 and not accheck)):
            return self.clone()
        controls -= classic_controls
        controls -= classic_anticontrols
        new_reg = None
        if type(gate) == SimpleGate:
            aux_targets = [self.qubit_map[id] for id in targets]
            aux_controls = {self.qubit_map[id] for id in controls}
            aux_anticontrols = {self.qubit_map[id] for id in anticontrols}
            doki_reg = self.doki.registry_apply(self.reg, gate.gate,
                                                aux_targets, aux_controls,
                                                aux_anticontrols,
                                                num_threads, self.verbose)
            new_reg = QRegistry(None, doki=self.doki)
            new_reg.reg = doki_reg
            new_reg.num_qubits = self.num_qubits
            new_reg.num_bits = self.num_bits
            new_reg.qubit_map = {k: v for k, v in self.qubit_map.items()}
            new_reg.classic_vals = {k: v for k, v in self.classic_vals.items()}
            new_reg.size = self.size
            new_reg.verbose = self.verbose
        else:
            try:
                new_reg, _ = apply_design(gate, self, targets,
                                          controls, anticontrols,
                                          num_threads=num_threads)
            except Exception:
                raise ValueError("Gate must be a valid string or a QGate")
        return new_reg

    def __getitem__(self, key):
        """Return the selected items of the state."""
        return self.get_state(key=key)

    def clone(self, num_threads=-1):
        if not np.allclose(num_threads % 1, 0):
            raise ValueError("num_threads must be an integer")
        num_threads = int(num_threads)
        new_reg = QRegistry(None, doki=self.doki)
        new_reg.num_qubits = self.num_qubits
        new_reg.num_bits = self.num_bits
        new_reg.size = self.size
        if self.reg is not None:
            new_reg.reg = self.doki.registry_clone(self.reg, num_threads,
                                                   self.verbose)
        else:
            new_reg.reg = None
        new_reg.qubit_map = {k: v for k, v in self.qubit_map.items()}
        new_reg.classic_vals = {k: v for k, v in self.classic_vals.items()}
        return new_reg

    def get_state(self, key=None, canonical=False):
        """Return the selected items of the state."""
        if key is None:
            key = slice(0, self.size, 1)
        if type(key) != slice and np.allclose(key % 1, 0):
            key = int(key)
            if (key < 0):
                key = self.size + key
            if (key >= self.size or key < 0):
                raise IndexError(f"index {key} is out of bounds " +
                                 f"for axis 0 with shape {self.size}")
            key = slice(key, key + 1, 1)
        if type(key) != slice:
            raise ValueError("key must be an index or a slice")

        start = key.start if key.start is not None else 0
        if (start < 0):
            start = self.size + start
            if (start < 0):
                start = 0
        stop = key.stop if key.stop is not None else self.size
        if (stop < 0):
            stop = self.size + stop
            if (stop < 0):
                stop = 0
        step = key.step if key.step is not None else 1
        res = np.array([self.doki.registry_get(self.reg, i, canonical,
                                               self.verbose)
                        for i in range(start, stop, step)])
        return res

    def density_matrix(self, canonical=True):
        """Return functional matrix of the density matrix."""
        return Funmatrix(self.doki.registry_density(self.reg, self.verbose),
                         "Rho")

    def get_entropy(self, **kwargs):
        """Calculate Von Newmann Entropy.

        Keyworded arguments:
            base: the base of the logarithm. Default = 2
                The string "e" is a valid value
        """
        base = kwargs.get('base', 2)
        entropy = 0
        size = self.get_state_size()
        rho = self.density_matrix()
        for i in range(size):
            p = rho[i, i].real
            if p > 0:
                if base == "e":
                    entropy += p * np.log(p)
                elif type(base) == int or type(base) == float:
                    entropy += p * np.log(p)/np.log(base)
        return -entropy

    def get_bloch_coords(self):
        """Get the polar coordinates of ONE qubit in the bloch sphere."""
        if self.get_num_qubits() != 1:
            raise NotImplementedError("Bloch sphere is only supported for " +
                                      "1 qubit registries")
        state = self.get_state(canonical=True)
        if state[0] == 0:
            phase = np.angle(state[1])
            state = np.exp(phase * -1j) * state
        theta = np.arccos(state[0])
        aux = state[1] / np.sin(theta)
        phi = np.arctan2(np.imag(aux), np.real(aux))
        return (theta, phi)

    def bra(self):
        """Get the conjugated row form state vector (bra <v|)."""
        k = np.array(self.get_state())
        k.shape = (1, k.shape[0])
        return np.conjugate(k)

    def ket(self):
        """Get the column form state vector (ket |v>)."""
        k = np.array(self.get_state())
        k.shape = (k.shape[0], 1)
        return k

    def prob(self, id, num_threads=-1):
        """Get the odds of getting 1 when measuring specified qubit."""
        id = _get_qubit_set(self.get_num_qubits(), [id], True, "argument")[0]
        return self.doki.registry_prob(self.reg, id, num_threads, self.verbose)


def superposition(a, b, num_threads=-1, verbose=False):
    """Join two registries into one by calculating tensor product."""
    if not np.allclose(num_threads % 1, 0):
        raise ValueError("num_threads must be an integer")
    num_threads = int(num_threads)
    doki_reg = a.doki.registry_join(a.reg, b.reg, num_threads, verbose)
    new_reg = QRegistry(None, doki=a.doki)
    new_reg.num_qubits = a.num_qubits + b.num_qubits
    new_reg.num_bits = a.num_bits + b.num_bits
    new_reg.size = 2**new_reg.num_qubits
    new_reg.reg = doki_reg
    b_total = b.num_qubits + b.num_bits
    new_keys = [k + b_total for k in a.qubit_map.keys()]
    new_keys += list(b.qubit_map.keys())
    new_keys.sort()
    new_reg.qubit_map = {new_keys[i]: i for i in range(len(new_keys))}
    new_reg.classic_vals = {k: v for k, v in b.classic_vals.items()}
    for k, v in a.classic_vals.items():
        new_reg.classic_vals[k + b_total] = v
    return new_reg
