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


def _get_op_data(num_qubits, gate, targets, controls, anticontrols):
    """Do basic error checking for arguments and return them."""
    targets = _get_qubit_set(num_qubits, targets, True, "targets")
    if gate is not None:
        if type(gate) == str:
            gate = SimpleGate(gate)
        num_targets = gate.num_qubits
        if len(targets) == 0:  # By default we use the least significant qubits
            targets = [i for i in range(num_targets)]
        if len(targets) != num_targets:
            raise ValueError(f"Specified gate is for {num_targets} qubits." +
                             f" {len(targets)} qubit ids given")
    controls = _get_qubit_set(num_qubits, controls, False, "controls")
    anticontrols = _get_qubit_set(num_qubits, anticontrols,
                                  False, "anticontrols")
    _check_no_intersection(targets, controls, anticontrols)

    return {"gate": gate, "targets": targets,
            "controls": controls, "anticontrols": anticontrols}


def _check_no_intersection(targets, controls, anticontrols):
    """Raise an exception if any qubit id is used more than once."""
    if len(controls.intersection(targets)) > 0:
        raise ValueError("A target cannot also be a control")
    if len(anticontrols.intersection(targets)) > 0:
        raise ValueError("A target cannot also be an anticontrol")
    if len(controls.intersection(anticontrols)) > 0:
        raise ValueError("A control cannot also be an anticontrol")


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
