import numpy as np
from abc import abstractmethod
from collections.abc import Iterable
from qsimov.structures.qbase import QBase
from qsimov.structures.simple_gate import SimpleGate


class QStructure(QBase):
    @abstractmethod
    def __init__(self, num_qubits, doki=None, verbose=False):
        pass

    @abstractmethod
    def apply_gate(self, gate, targets=None, controls=None, anticontrols=None,
                   num_threads=-1):
        pass

    @abstractmethod
    def measure(self, ids, random_generator=np.random.rand):
        pass

    @abstractmethod
    def prob(self, id):
        pass

    @abstractmethod
    def get_state(self, key=None, canonical=False):
        pass

    @abstractmethod
    def get_classic(self, id):
        pass

    @abstractmethod
    def clone(self, num_threads=-1):
        pass

    @abstractmethod
    def free(self):
        pass

    @abstractmethod
    def get_bloch_coords(self, key=None):
        pass

    @abstractmethod
    def bloch(self, key=None):
        pass


def _get_op_data(num_qubits, num_bits, gate, targets, c_targets, outputs,
                 controls, anticontrols, c_controls, c_anticontrols,
                 empty=False):
    """Do basic error checking for arguments and return them."""
    targets = _get_qubit_set(num_qubits, targets, True, "targets")
    c_targets = _get_qubit_set(num_bits, c_targets, True, "classic targets")
    if gate is not None:
        num_c_targets = 0
        if type(gate) == str:
            gate = SimpleGate(gate)
        elif type(gate) != SimpleGate:
            num_c_targets = gate.num_bits
        num_targets = gate.num_qubits
        if len(targets) == 0:  # By default we use the least significant qubits
            targets = [i for i in range(num_targets)]
        if len(c_targets) == 0:  # By default we use the least significant bits
            c_targets = [i for i in range(num_c_targets)]
        if len(targets) != num_targets:
            raise ValueError(f"Specified gate is for {num_targets} qubits." +
                             f" {len(targets)} qubit ids given")
        if len(c_targets) != num_c_targets:
            raise ValueError(f"Specified gate is for {num_c_targets} bits." +
                             f" {len(c_targets)} bit ids given")
    controls = _get_qubit_set(num_qubits, controls, False, "controls")
    c_controls = _get_qubit_set(num_bits, c_controls,
                                False, "classic controls")
    anticontrols = _get_qubit_set(num_qubits, anticontrols,
                                  False, "anticontrols")
    c_anticontrols = _get_qubit_set(num_bits, c_anticontrols,
                                    False, "classic anticontrols")
    outputs = _get_qubit_set(num_bits, outputs, True, "outputs")
    _check_no_intersection(targets, controls, anticontrols)
    _check_no_intersection(c_targets, c_controls, c_anticontrols, True)
    if gate is None:
        if not empty:
            if len(outputs) != len(targets):
                raise ValueError(f"Expected {len(targets)} output bits." +
                                 f" {len(outputs)} ids given")
            if len(controls) + len(anticontrols) > 0:
                raise ValueError("Measures can only be controlled by bits")
    elif len(outputs) != 0:
        raise ValueError("Gate applications can't have classical outputs")

    return {"gate": gate,
            "targets": targets, "c_targets": c_targets, "outputs": outputs,
            "controls": controls, "anticontrols": anticontrols,
            "c_controls": c_controls, "c_anticontrols": c_anticontrols}


def _check_no_intersection(targets, controls, anticontrols, classic=False):
    """Raise an exception if any qubit id is used more than once."""
    aux = ""
    if classic:
        aux = "classic "
    if len(controls.intersection(targets)) > 0:
        raise ValueError(f"A {aux}target cannot also be a control")
    if len(anticontrols.intersection(targets)) > 0:
        raise ValueError(f"A {aux}target cannot also be an anticontrol")
    if len(controls.intersection(anticontrols)) > 0:
        raise ValueError(f"A {aux}control cannot also be an anticontrol")


def _get_qubit_set(max_qubits, raw_ids, sorted, name):
    """Get a set or sorted set (list) of qubit ids from raw_ids."""
    if raw_ids is None:
        if sorted:
            return []
        else:
            return set()
    if not isinstance(raw_ids, Iterable):
        raw_ids = [raw_ids]
    num_ids = len(raw_ids)
    ids_check = all([np.allclose(qubit_id % 1, 0)
                     and qubit_id < max_qubits
                     and qubit_id >= 0
                     for qubit_id in raw_ids])
    if not ids_check:
        raise ValueError(f"Invalid id found in {name}")
    if sorted:
        id_list = [int(raw_ids[i]) for i in range(num_ids)]
    else:
        id_list = [raw_id for raw_id in raw_ids]
    id_set = set(id_list)
    if num_ids != len(id_set):  # Check duplicates
        raise ValueError(f"{name} list cannot have duplicated ids")
    if sorted:
        return id_list
    return id_set


def _get_key_with_defaults(key, size, def_start, def_stop, def_step):
    if key is None:
        key = slice(def_start, def_stop, def_step)
    if type(key) != slice and np.allclose(key % 1, 0):
        key = int(key)
        if (key < 0):
            key = size + key
        if (key >= size or key < 0):
            raise IndexError(f"index {key} is out of bounds " +
                             f"for axis 0 with shape {size}")
        key = slice(key, key + 1, 1)
    if type(key) != slice:
        raise ValueError("key must be an index or a slice")

    start = key.start if key.start is not None else def_start
    if (start < 0):
        start = size + start
        if (start < 0):
            start = 0
    stop = key.stop if key.stop is not None else def_stop
    if (stop < 0):
        stop = size + stop
        if (stop < 0):
            stop = 0
    step = key.step if key.step is not None else def_step

    return (start, stop, step)
