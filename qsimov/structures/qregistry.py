'''
QSimov: A Quantum Computing Toolkit.
Copyright (C) 2017  Hernán Indíbil de la Cruz Calvo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

"""Module that provides the lowest level data structure: QRegistry.

Data Structures:
    QRegistry: Quantum Registry, base of all quantum related operations

Functions:
    superposition: join two registries into one by calculating tensor product.
"""
import numpy as np
import sympy as sp

from sympy.functions.elementary.complexes import arg
from qsimov.structures.funmatrix import Funmatrix
from qsimov.structures.simple_gate import SimpleGate
from qsimov.structures.qstructure import QStructure, _get_op_data, \
                                         _get_qubit_set, _get_key_with_defaults
from qsimov.connectors.qsimovapi import apply_design


class QRegistry(QStructure):
    """Quantum Registry, base of all quantum related operations."""

    def __init__(self, num_qubits, data=None, doki=None, verbose=False):
        """Initialize QRegistry to state 0.

        num_qubits -> number of QuBits in the registry.
        """
        if doki is None:
            import doki
        self.doki = doki
        if data is None:
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
        else:
            self.reg = self.doki.registry_new_data(data["num_qubits"], data["reg"], verbose)
            self.num_qubits = data["num_qubits"]
            self.num_bits = data["num_bits"]
            self.size = data["size"]
            self.qubit_map = data["qubit_map"]
            self.classic_vals = data["classic_vals"]
        self.verbose = verbose

    def __del__(self):
        """Clean after deletion."""
        self.free()

    def get_data(self):
        return {"reg": self.get_state(),
                "num_qubits": self.num_qubits,
                "num_bits": self.num_bits,
                "size": self.size,
                "qubit_map": self.qubit_map,
                "classic_vals": self.classic_vals}

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
        return self.num_qubits + self.num_bits

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
        ids_set = _get_qubit_set(self.get_num_qubits(),
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
                   num_threads=-1, target=None, control=None, anticontrol=None):
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
        if target is not None:
            print("[WARNING] target keyworded argument is deprecated. Please use targets instead")
            if targets is not None:
                raise ValueError("target argument can't be set alongside targets")
            targets = target
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
        if not np.allclose(num_threads % 1, 0):
            raise ValueError("num_threads must be an integer")
        num_threads = int(num_threads)
        num_qubits = self.num_qubits + self.num_bits
        op_data = _get_op_data(num_qubits, 0, gate, targets, None, None,
                               controls, anticontrols, None, None)
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
            auxg = gate.gate
            if hasattr(self.doki, "BRegistry"):
                auxg = gate.matrix
            doki_reg = self.doki.registry_apply(self.reg, auxg,
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
                res = apply_design(gate, self, [], targets,
                                   controls, anticontrols,
                                   num_threads=num_threads, shots=1,
                                   return_struct=True)
                new_reg, _ = res[0]
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
        start, stop, step = _get_key_with_defaults(key, self.size,
                                                   0, self.size, 1)
        res = np.array([self.doki.registry_get(self.reg, i, canonical,
                                               self.verbose)
                        for i in range(start, stop, step)])
        return res

    def density_matrix(self, canonical=False):
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

    def get_bloch_coords(self, key=None):
        """Get the polar coordinates of ONE qubit in the bloch sphere."""
        if hasattr(self.doki, "BRegistry"):
            start, stop, step = _get_key_with_defaults(key, self.num_qubits,
                                                       0, self.num_qubits, 1)
            ids = [id for id in range(start, stop, step)]
            coords = [(sp.acos(self.reg.r[id].v[2, 0]),
                       sp.atan2(self.reg.r[id].v[1, 0], self.reg.r[id].v[0, 0]))
                      if self.reg.r[id].v[0, 0] != 0
                      else (sp.acos(self.reg.r[id].v[2, 0]), 0)
                      for id in ids]
            return coords
        if self.num_qubits != 1:
            raise NotImplementedError("Bloch sphere is only supported for " +
                                      "1 qubit registries")
        if key is not None and key >= self.get_num_qubits():
            raise ValueError(f"Qubit {key} is not in this registry")
        if key is not None and key in self.classic_vals:
            raise ValueError("That is a classical bit, not a qubit")
        state = self.get_state(canonical=True)
        cos_t2 = state[0].real
        sin_t2 = abs(state[1])
        phi = arg(state[1])
        if phi == sp.nan:
            phi = 0
        theta = sp.atan2(sin_t2, cos_t2) * 2
        if theta < 0:
            theta += 2 * np.pi
        if phi < 0:
            phi += 2 * np.pi
        return (theta, phi)

    def bra(self):
        """Get the conjugated row form state vector (bra <v|)."""
        k = self.get_state()
        k.shape = (1, k.shape[0])
        return np.conjugate(k)

    def ket(self):
        """Get the column form state vector (ket |v>)."""
        k = self.get_state()
        k.shape = (k.shape[0], 1)
        return k

    def prob(self, id, num_threads=-1):
        """Get the odds of getting 1 when measuring specified qubit."""
        id = _get_qubit_set(self.get_num_qubits(), [id], True, "argument")[0]
        if id in self.classic_vals:
            return self.classic_vals[id]
        for i in range(self.get_num_qubits(), id):
            if i in self.classic_vals:
                id -= 1
        return self.doki.registry_prob(self.reg, id, num_threads, self.verbose)


def superposition(a, b, num_threads=-1, verbose=False):
    """Join two registries into one by calculating tensor product."""
    if not np.allclose(num_threads % 1, 0):
        raise ValueError("num_threads must be an integer")
    num_threads = int(num_threads)
    doki_reg = None
    if a.reg is None and b.reg is None:
        doki_reg = None
    elif a.reg is None:
        doki_reg = b.reg
    elif b.reg is None:
        doki_reg = a.reg
    else:
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
