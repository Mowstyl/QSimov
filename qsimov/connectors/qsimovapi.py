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

"""Module that executes a circuit using Doki (QSimov core).

This module provides functions that execute a given set of quantum gates
on a given quantum system, using Doki (the QSimov core written in C) to
perform the simulation.
"""
import numpy as np
from qsimov.structures.qstructure import QStructure, _get_op_data
from qsimov.structures.qdesign import QDesign


def _check_classical(classical_reg, c_controls, c_anticontrols):
    for k in c_controls:
        if classical_reg[k] is None:
            raise ValueError("Undefined value for classic " +
                             f"control bit {k}.")
    for k in c_anticontrols:
        if classical_reg[k] is None:
            raise ValueError("Undefined value for classic " +
                             f"anticontrol bit {k}.")
    return (all(classical_reg[k] for k in c_controls) and
            not any(classical_reg[k] for k in c_anticontrols))


def to_lines(qdesign):
    pass


def apply_design(qdesign: QDesign, qstruct: QStructure, classical_reg, targets=None,
                 c_targets=None, controls=None, anticontrols=None,
                 c_controls=None, c_anticontrols=None,
                 random_generator=np.random.rand, num_threads=-1):
    """Apply specified gate to specified qubit with specified controls.

    Positional arguments:
        qstruct: QRegistry or QSystem affected by this gate
    Keyworded arguments:
        targets: id or list of ids of the qubits the gate will target
        controls: id or set of ids of the qubit that will act as controls
        anticontrols: id or set of ids of the qubit that will act as
                     anticontrols
        num_threads: number of threads to use
    """
    if not isinstance(qdesign, QDesign):
        raise ValueError("qdesign must be a QGate or a QCircuit")
    if not isinstance(qstruct, QStructure):
        raise ValueError("qstruct must be a QRegistry or a QSystem")
    num_qubits = qstruct.get_num_qubits()
    num_bits = qdesign.get_num_bits()
    op_data = _get_op_data(num_qubits, num_bits, None,
                           targets, c_targets, None,
                           controls, anticontrols, c_controls, c_anticontrols,
                           empty=True)
    targets = op_data["targets"]
    c_targets = op_data["c_targets"]
    controls = op_data["controls"]
    anticontrols = op_data["anticontrols"]
    c_controls = op_data["c_controls"]
    c_anticontrols = op_data["c_anticontrols"]
    num_targets = qdesign.get_num_qubits()
    if len(targets) == 0:  # By default we use the least significant qubits
        targets = [i for i in range(num_targets)]
    if len(c_targets) == 0:  # By default we use the least significant qubits
        c_targets = [i for i in range(num_bits)]
    if len(targets) != num_targets:
        raise ValueError(f"Specified gate is for {num_targets} qubits." +
                         f" {len(targets)} qubit ids given")
    for k in c_controls:
        if (k not in classical_reg):
            raise ValueError("Undefined value for classic " +
                             f"control bit {k}.")
    for k in c_anticontrols:
        if (k not in classical_reg):
            raise ValueError("Undefined value for classic " +
                             f"anticontrol bit {k}.")
    if not _check_classical(classical_reg, c_controls, c_anticontrols):
        return (qstruct, classical_reg)
    new_struct = qstruct.clone()
    aux = None
    exception = None
    try:
        for gate_data in qdesign.get_operations(flatten=True):
            if gate_data == "BARRIER":
                continue
            curr_c_controls = {c_targets[i] for i in gate_data["c_controls"]}
            curr_c_acontrols = {c_targets[i]
                                for i in gate_data["c_anticontrols"]}
            if not _check_classical(classical_reg,
                                    curr_c_controls, curr_c_acontrols):
                continue
            aux = new_struct
            curr_targets = [targets[i] for i in gate_data["targets"]]
            curr_controls = {targets[i] for i in gate_data["controls"]}
            curr_controls = curr_controls.union(controls)
            curr_anticontrols = {targets[i]
                                 for i in gate_data["anticontrols"]}
            curr_anticontrols = curr_anticontrols.union(anticontrols)
            # print(gate_data)
            # print(c_targets)
            curr_outputs = [c_targets[i] for i in gate_data["outputs"]]
            gate = gate_data["gate"]
            if gate == "MEASURE":
                new_struct, m = aux.measure(curr_targets,
                                            random_generator=random_generator)
                for i in range(len(curr_outputs)):
                    classical_reg[curr_outputs[i]] = m[curr_targets[i]]
            else:
                new_struct = aux.apply_gate(gate,
                                            targets=curr_targets,
                                            controls=curr_controls,
                                            anticontrols=curr_anticontrols,
                                            num_threads=num_threads)
            if aux is not qstruct:
                del aux
                aux = None
    except Exception as ex:
        exception = ex
        if aux is not None and aux is not qstruct:
            del aux
        if new_struct is not None:
            del new_struct
    if exception is not None:
        raise exception
    return (new_struct, classical_reg)


def execute(qcircuit, random_generator=np.random.rand, num_threads=-1,
            use_system=True, return_struct=False, core=None, iterations=1):
    """Execute the gates in lines on a qsystem.

    Gets repeated the specified number of iterations.
    Returns the result of each iteration.
    """
    from qsimov.structures.qsystem import QSystem
    from qsimov.structures.qregistry import QRegistry
    num_qubits = qcircuit.get_num_qubits()
    num_bits = qcircuit.get_num_bits()
    num_ancilla = len(qcircuit.ancilla)
    if core is None or core.upper() == "CPU":
        import doki
    elif core.upper() == "GPU":
        import doki_gpu as doki
    elif core.upper() == "MPI":
        import doki_mpi as doki
    elif isinstance(core, str):
        raise ValueError("Unknown Doki core. Valid values: CPU, GPU or MPI")
    else:
        doki = core
    this_struct = QSystem
    if not use_system:
        this_struct = QRegistry
    old_sys = this_struct(num_qubits, doki=doki)
    for i in range(len(qcircuit.ancilla)):
        if qcircuit.ancilla[i] == 1:
            aux = old_sys.apply_gate("X", targets=num_qubits-num_ancilla+i)
            del old_sys
            old_sys = aux
    results = []
    for i in range(iterations):
        classical_reg = [None for i in range(num_bits)]
        new_sys, _ = apply_design(qcircuit, old_sys, classical_reg,
                                  targets=[i for i in range(num_qubits)],
                                  c_targets=[i for i in range(num_bits)],
                                  random_generator=random_generator,
                                  num_threads=num_threads)
        if return_struct:
            results.append([new_sys, classical_reg])
        else:
            results.append(classical_reg)
            del new_sys
    del old_sys
    return results
