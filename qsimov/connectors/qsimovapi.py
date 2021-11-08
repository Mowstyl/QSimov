"""Module that executes a circuit using Doki (QSimov core).

This module provides functions that execute a given set of quantum gates
on a given quantum system, using Doki (the QSimov core written in C) to
perform the simulation.
"""
import numpy as np
from qsimov.structures.qstructure import QStructure, _get_op_data
from qsimov.structures.qdesign import QDesign


def apply_design(qdesign, qstruct, targets=None,
                 controls=None, anticontrols=None,
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
    op_data = _get_op_data(num_qubits, None, targets,
                           controls, anticontrols)
    targets = op_data["targets"]
    controls = op_data["controls"]
    anticontrols = op_data["anticontrols"]
    num_targets = qdesign.get_num_qubits()
    if len(targets) == 0:  # By default we use the least significant qubits
        targets = [i for i in range(num_targets)]
    if len(targets) != num_targets:
        raise ValueError(f"Specified gate is for {num_targets} qubits." +
                         f" {len(targets)} qubit ids given")
    new_struct = qstruct.clone()
    measures = []
    aux = None
    exception = None
    try:
        for gate_data in qdesign.get_operations():
            if gate_data == "BARRIER":
                continue
            aux = new_struct
            curr_targets = [targets[i] for i in gate_data["targets"]]
            curr_controls = {targets[i] for i in gate_data["controls"]}
            curr_controls = curr_controls.union(controls)
            curr_anticontrols = {targets[i]
                                 for i in gate_data["anticontrols"]}
            curr_anticontrols = curr_anticontrols.union(anticontrols)
            gate = gate_data["gate"]
            if isinstance(gate, QDesign):
                new_struct, m = apply_design(gate, aux,
                                             targets=curr_targets,
                                             controls=curr_controls,
                                             anticontrols=curr_anticontrols,
                                             num_threads=num_threads)
                if len(m) > 0:
                    measures.append(m)
            elif gate == "MEASURE":
                cvals = [aux.get_classic(id) for id in curr_controls]
                acvals = [aux.get_classic(id) for id in curr_anticontrols]
                if any(val is None for val in cvals + acvals):
                    raise ValueError("Can't do a measurement controlled " +
                                     "by a qubit, only by classic bits")
                if not all(cvals) or any(acvals):
                    continue
                new_struct, m = aux.measure(curr_targets,
                                            random_generator=random_generator)
                measures.append(m)
            else:
                new_struct = aux.apply_gate(gate,
                                            targets=curr_targets,
                                            controls=curr_controls,
                                            anticontrols=curr_anticontrols,
                                            num_threads=num_threads)
            if aux is not qstruct:
                aux.free()
                aux = None
    except Exception as ex:
        exception = ex
        if aux is not None and aux is not qstruct:
            aux.free()
        if new_struct is not None:
            new_struct.free()
    if exception is not None:
        raise exception
    return (new_struct, measures)


def execute(qcircuit, random_generator=np.random.rand, num_threads=-1,
            use_system=True, return_struct=False, core=None, iterations=1):
    """Execute the gates in lines on a qsystem.

    Gets repeated the specified number of iterations.
    Returns the result of each iteration.
    """
    from qsimov.structures.qsystem import QSystem
    from qsimov.structures.qregistry import QRegistry
    num_qubits = qcircuit.get_num_qubits()
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
            old_sys.free()
            old_sys = aux
    results = []
    for i in range(iterations):
        new_sys, mess = apply_design(qcircuit, old_sys,
                                     targets=[i for i in range(num_qubits)],
                                     random_generator=random_generator,
                                     num_threads=num_threads)
        if return_struct:
            results.append([new_sys, mess])
        else:
            results.append(mess)
            if new_sys is not None:
                new_sys.free()
            del new_sys
    old_sys.free()
    del old_sys
    return results
