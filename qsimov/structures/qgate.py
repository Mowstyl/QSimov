"""Module that provides a data structure representing a quantum gate.

Data Structures:
    QGate: Quantum Gate, built from elemental gates or other QGates

Functions:
    get_gate: join two registries into one by calculating tensor product.
"""
import numpy as np
import ctypes as ct
import platform as plat
import qsimov.connectors.parser as prs
import os
from collections.abc import Iterable
from qsimov.structures.measure import Measure
from os.path import sep

# DLL Load
if plat.system() == "Windows":
    extension = ".dll"
else:
    extension = ".so"
_root_folder = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
_lib_folder = _root_folder + sep + "lib"
_qsimov_path = _lib_folder + sep + "libqsimov" + extension
if hasattr(os, "add_dll_directory"):
    _qsimov_path = "libqsimov" + extension
_qsimov = ct.CDLL(_qsimov_path)

c_double_p = ct.POINTER(ct.c_double)


_cGetQGate = _qsimov.getQGate
_cGetQGate.argtypes = [ct.c_char_p, ct.c_double, ct.c_double, ct.c_double,
                       ct.c_int]
_cGetQGate.restype = ct.c_void_p
"""C QGate getQGate function.
Positional arguments:
    string -> Name of the gate
    double -> First parameter of the gate (optional)
    double -> Second parameter of the gate (optional)
    double -> Third parameter of the gate (optional)
    int -> Boolean. Whether or not to invert the gate. 1 -> Invert
Return:
    pointer to C QGate
"""

_cGetQGateQubits = _qsimov.getQGateQubits
_cGetQGateQubits.argtypes = [ct.c_void_p]
_cGetQGateQubits.restype = ct.c_int
"""C QGate getQGateQubits function.
Positional arguments:
    pointer: C QGate
Return:
    int: Max number of qubits affected by the gate
"""

_cGetQGateSize = _qsimov.getQGateSize
_cGetQGateSize.argtypes = [ct.c_void_p]
_cGetQGateSize.restype = ct.c_uint
"""C QGate getQGateQubits function.
Positional arguments:
    pointer: C QGate
Return:
    unsigned int: number of rows or columns in the matrix of the QGate
"""

_cGet2dGate = _qsimov.get2dGate
_cGet2dGate.argtypes = [ct.c_void_p]
_cGet2dGate.restype = c_double_p
"""C QGate getQGateQubits function.
Positional arguments:
    pointer: C QGate
Return:
    double array: first half -> real part of flattened matrix
                  second half -> imaginary part of flattened matrix
"""


class QGate(object):
    """Quantum Gate, built from elemental gates or other QGates."""

    def __init__(self, name="UNNAMED", size=None):
        """Quantum Gate constructor.

        Keyworded arguments:
            name: name of the gate. Default="UNNAMED"
            size: maximum number of qubits affected by this gate
                If nothing specified it will be calculated in the
                first call to add_line
        """
        # Whether or not add_line has been successfully called
        self.empty = True
        # List of lines. A line is a list of gates that can be executed
        # in parallel. One gate per qubit. None means Identity gate
        self.lines = []
        # List of optimized lines. Like the previous, but only with
        # native gates and optimized by removing unnecesary operations
        self.oplines = []
        # List with the first free index for each qubit in oplines
        self.freeindexes = None
        # The greater index used in freeindexes
        self.lastindex = -1
        # Name of the gate
        self.name = name
        # Maximum number of qubits affected by this gate
        self.size = size

    '''
    def __getitem__(self, key):
        return self.getMatrix()[key]
    '''

    def __repr__(self):
        """Return string representation of the gate."""
        return self.name

    def __str__(self):
        """Return string representation of the gate."""
        return self.name

    def addLine(self, *args, **kwargs):
        """Use add_line method instead. DEPRECATED."""
        print("Method QGate.addLine is deprecated.",
              "Please use add_line if you seek the same functionality")
        return self.add_line(*args, **kwargs)

    def add_line(self, *args, **kwargs):
        """Apply specified gate to specified qubit with specified controls.

        Positional arguments:
            comma separated gates, their sizes must match the number of
                qubits in the system. Sorted by their least significant
                target qubit id.
        Keyworded arguments:
            add_to_lines: whether to add args to lines or only to oplines
            offset: offset to add to all the qubit ids each arg of this line
            control: id or list of ids of the qubit that will act as
                      controls
            anticontrol: id or list of ids of the qubit that will act as
                      anticontrols
        """
        add_to_lines, offset, controls, anticontrols = _get_line_args(**kwargs)
        _add_gates(self, args, add_to_lines, offset, controls, anticontrols)

    def setName(self, name):
        """Use set_name method instead. DEPRECATED."""
        print("Method QGate.setName is deprecated.",
              "Please use set_name if you seek the same functionality")
        return self.set_name(name)

    def set_name(self, name):
        """Set the name of this gate."""
        self.name = name

    def dagger(self):
        """Return the Conjugate Transpose of the given matrix."""
        return self.invert()

    def invert(self):
        """Return the Conjugate Transpose of the given matrix."""
        invgate = QGate(self.name + "-1")
        for line in self.lines[::-1]:
            invgate.add_line(*[_invert_gate(gate) for gate in line])
        return invgate

    def _apply_gate(self, registry, qubit, control, anticontrol,
                    optimize=False):
        """Apply specified gate to specified qubit with specified controls.

        Positional arguments:
            registry: QRegistry or QSystem affected by this gate
            qubit: id of the least significant qubit the gate will target
            control: id or list of ids of the qubit that will act as controls
            anticontrol: id or list of ids of the qubit that will act as
                         anticontrols
        Keyworded arguments:
            optimize: only for QGates. Whether to use optimized lines or
                      user defined lines
        """
        if control is None:
            control = set()
        if anticontrol is None:
            anticontrol = set()
        lines = self.lines
        if optimize:
            lines = self.oplines
        for line in lines:
            currbit = 0
            for i in range(len(line)):
                if (line[i] is not None
                        and line[i][0] is not None
                        and (isinstance(line[i][0], QGate)
                             or line[i][0].lower() != "i")):
                    g_control = set()
                    g_anticontrol = set()
                    if line[i][1] is not None:
                        g_control = set([qubit + c
                                         for c in line[i][1]]).union(control)
                    if line[i][2] is not None:
                        g_anticontrol = set([qubit + ac
                                             for ac in line[i][2]]) \
                                             .union(anticontrol)
                    if type(line[i][0]) == str:
                        registry.apply_gate(_add_qubit_offset(line[i][0],
                                                              qubit),
                                            qubit=currbit+qubit,
                                            control=g_control,
                                            anticontrol=g_anticontrol)
                    else:
                        line[i][0]._apply_gate(registry, qubit+currbit,
                                               g_control, g_anticontrol)
                    currbit += get_gate_qubits(line[i][0])
                else:
                    currbit += 1


def _check_line_size(q_design, size):
    """Check if valid line size and initializes structures if empty."""
    if size <= 0:
        print("No gates. Just Monika.")
        return False
    if (q_design.empty):
        q_design.size = size
        q_design.empty = False
        q_design.freeindexes = [0 for i in range(size)]
    if (q_design.size != size):
        raise ValueError("Required gates for " + str(q_design.size) +
                         " qubits. Received gates for " + str(size) +
                         " qubits.")
    return True


def _add_gates(q_design, args, add_to_lines, offset, controls, anticontrols):
    args = [_rebuild_gate_name(gate) for gate in args]
    parties = _get_parties(args)
    size = len(parties)

    if not _check_line_size(q_design, size):
        return
    if add_to_lines:
        q_design.lines += [args]
    for i in range(len(args)):
        arg = args[i]
        if arg is not None:
            arg[1] = controls.union({control + offset
                                     for control in arg[1]})
            arg[2] = anticontrols.union({acontrol + offset
                                         for acontrol in arg[2]})
            if isinstance(arg[0], QGate):
                for line in arg[0].oplines:
                    none_after_gate = q_design.size - arg[0].size - i
                    q_design.add_line(*[None for j in range(i)], *line,
                                      *[None for j in range(none_after_gate)],
                                      add_to_lines=False, offset=i,
                                      control=arg[1],
                                      anticontrol=arg[2])
            else:
                _update_oplines(q_design, i, args)


def _recalculate_free(q_design):
    """Recalculate the free indexes in oplines for each qubit."""
    for party in range(q_design.size):
        q_design.freeindexes[party] = 0
    q_design.lastindex = -1
    for i in range(len(q_design.oplines)):
        line = q_design.oplines[i]
        if isinstance(line[0], Measure):
            parties = range(q_design.size)
            q_design.lastindex = i
        else:
            parties = _get_parties(line)
            if len(parties) > 0:
                q_design.lastindex = i
        for party in parties:
            q_design.freeindexes[party] = i + 1


def _update_oplines(q_design, i, args):
    """Update the optimized lines.

    Positional arguments:
        i: index of the gate argument
        args: list of gates, their sizes must match the number of
            qubits in the system. Sorted by their least significant
            target qubit id.
    """
    arg = args[i]
    parties = _get_parties([arg if j == i else None
                            for j in range(len(args))], ignore_empty=True)
    num_targets = len(parties) - len(arg[1]) - len(arg[2])
    freeindex = max([q_design.freeindexes[party]
                     for party in parties])
    skip = False
    if freeindex > 0:
        previous_gate = q_design.oplines[freeindex-1][i]
        if previous_gate is not None \
                and not isinstance(previous_gate, Measure) \
                and arg[0] == _invert_str_gate(previous_gate[0]) \
                and arg[1] == previous_gate[1]  \
                and arg[2] == previous_gate[2]:
            aux1 = q_design.oplines[freeindex-1][:i]
            aux2 = [None for j in range(num_targets)]
            aux3 = q_design.oplines[freeindex-1][i+1:]
            q_design.oplines[freeindex-1] = aux1 + aux2 + aux3

            if len(q_design.oplines[freeindex-1]) == 0 or \
               all([elem is None or elem == "I"
                    for elem in q_design.oplines[freeindex-1]]):
                del q_design.oplines[freeindex-1]
            _recalculate_free(q_design)
            skip = True
    if not skip:
        for party in parties:
            q_design.freeindexes[party] = freeindex + 1
        if freeindex > q_design.lastindex:
            q_design.oplines.append([None for i in range(q_design.size)])
            q_design.lastindex += 1
        q_design.oplines[freeindex][i] = arg
        while num_targets > 1:
            del q_design.oplines[freeindex][i+1]
            num_targets -= 1


def _get_line_args(**kwargs):
    """Return given or default values of add_line keyworded arguments."""
    add_to_lines = True
    if "add_line" in kwargs:
        print("Argument add_line is deprecated. Use add_to_lines instead")
        add_to_lines = kwargs["add_line"]
    if "add_to_lines" in kwargs:
        add_to_lines = kwargs["add_to_lines"]

    offset = 0
    if "offset" in kwargs:
        offset = kwargs["offset"]

    controls = set()
    if "controls" in kwargs:
        print("Argument controls is deprecated. Use control instead")
        controls = kwargs["controls"]
    if "control" in kwargs:
        controls = kwargs["control"]

    anticontrols = set()
    if "anticontrols" in kwargs:
        print("Argument anticontrols is deprecated.",
              "Use anticontrol instead")
        anticontrols = kwargs["anticontrols"]
    if "anticontrol" in kwargs:
        anticontrols = kwargs["anticontrol"]
    return add_to_lines, offset, controls, anticontrols


def getGate(gate_name):
    """Use get_gate method instead. DEPRECATED."""
    print("Method QGate.getGate is deprecated.",
          "Please use get_gate if you seek the same functionality")
    return get_gate(gate_name)


def get_gate(gate_name):
    """Return the matrix of the specified gate."""
    name, arg1, arg2, arg3, invert = prs.get_gate_data(gate_name)
    if arg1 is None:
        arg1 = 0
    if arg2 is None:
        arg2 = 0
    if arg3 is None:
        arg3 = 0
    qgate = ct.c_void_p(_cGetQGate(ct.c_char_p(name.encode()),
                                   ct.c_double(arg1), ct.c_double(arg2),
                                   ct.c_double(arg3), ct.c_int(int(invert))))
    # NOTNULLPTR and True = True, NOTNULLPTR or False = NOTNULLPTR
    if (qgate or False) == qgate:
        size = int(_cGetQGateSize(qgate))
        plainmatrix2d = _cGet2dGate(qgate)[:size*size*2]
        rematrix2d = plainmatrix2d[:size*size]
        immatrix2d = plainmatrix2d[size*size:size*size*2]
        matrix = np.array([complex(rematrix2d[i], immatrix2d[i])
                           for i in range(size * size)])
        matrix = matrix.reshape(size, size)
    else:
        matrix = None
        print("Error while getting the specified gate!")

    return matrix


def _rebuild_gate_name(gate):
    """Standarize name of native gate and return tuple with gate data.

    Positional arguments:
        gate: string with name of a native gate
                OR
              list containing
              0: string with name of a native gate
              1 (optional): control qubit or list of control qubits ids
              2 (optional): anticontrol qubit or list of anticontrol qubits ids
    Return:
        list containing
        0: string with standarized name of a native gate
        1: set of control qubits ids (can be empty)
        2: set of anticontrol qubits ids (can be empty)
    """
    cons = set()
    acons = set()
    gatename = gate
    if not isinstance(gate, str) and isinstance(gate, Iterable):
        gatename = gate[0]
        if (len(gate) > 1):
            if gate[1] is not None:
                if isinstance(gate[1], Iterable) and len(gate[1]) > 0:
                    cons = set(gate[1])
                elif not isinstance(gate[1], Iterable):
                    cons = set([gate[1]])
        if (len(gate) > 2):
            if gate[2] is not None:
                if isinstance(gate[2], Iterable) and len(gate[2]) > 0:
                    acons = set(gate[2])
                elif not isinstance(gate[2], Iterable):
                    acons = set([gate[2]])
    if not isinstance(gatename, QGate):
        if gatename is not None and gatename.lower() != "i":
            name, arg1, arg2, arg3, invert = prs.get_gate_data(gatename)
            if arg1 is not None:
                name += "(" + str(arg1)
            if arg2 is not None:
                name += "," + str(arg2)
            if arg3 is not None:
                name += "," + str(arg3) + ")"
            elif arg1 is not None:
                name += ")"
            if invert:
                name += "-1"
        else:
            name = None
    else:
        name = gatename
    return [name, cons, acons] if name is not None else None


def _get_qubit_arg(gatename):
    """Get arguments from string with gate name."""
    args = None
    if isinstance(gatename, str):
        qubitargs = ["XX", "YY", "ZZ"]
        name, arg1, arg2, arg3, invert = prs.get_gate_data(gatename)
        if name in qubitargs or "SWAP" in name:
            if name in qubitargs:
                args = set([arg2, arg3])
            else:
                args = set([arg1, arg2])
    return args


def _add_qubit_offset(gatename, offset):
    """Return the gate name string adding the offset to any id argument."""
    qubitargs = ["XX", "YY", "ZZ"]
    name, arg1, arg2, arg3, invert = prs.get_gate_data(gatename)
    if name in qubitargs or "SWAP" in name:
        arg2 += offset
        if name in qubitargs:
            arg3 += offset
        else:
            arg1 += offset
        if arg1 is not None:
            name += "(" + str(arg1)
        if arg2 is not None:
            name += str(arg2)
        if arg3 is not None:
            name += str(arg3) + ")"
        if invert:
            name += "-1"
        gatename = name
    return gatename


def _invert_gate(gate_data):
    """Invert the given gate."""
    gate = gate_data
    control = None
    anticontrol = None
    is_list = isinstance(gate_data, Iterable)
    if is_list:
        gate = gate_data[0]
        control = gate_data[1]
        anticontrol = gate_data[2]
    if isinstance(gate, QGate):
        gate = gate.dagger()
    elif gate is not None:
        gate = _invert_str_gate(gate)
    if is_list:
        gate = [gate, control, anticontrol]
    return gate


def _invert_str_gate(gatename):
    """Invert the given native gate."""
    selfinvert = ["X", "Y", "Z", "H", "SWAP", "I"]
    name, arg1, arg2, arg3, invert = prs.get_gate_data(gatename)
    if name not in selfinvert:
        if invert:
            gatename = gatename[:-2]
        else:
            gatename += "-1"
    return gatename


def _get_parties(args, ignore_empty=False, offset=0):
    """Return the set of qubit ids affected by the given list of gates.

    Positional arguments:
        args: list of gates
    Keyworded arguments:
        ignore_empty: whether or not to ignore qubits
                      affected by identity gates
        offset: number to add to all the ids
    """
    parties = set()
    empties = set()
    currbit = 0
    for i in range(len(args)):
        if args[i] is not None:
            myparty = _get_qubit_arg(args[i][0])
            if myparty is None:
                gateSize = get_gate_qubits(args[i][0])
                if gateSize is None:
                    continue
                myparty = set([currbit + j for j in range(gateSize)])
            if args[i][0] is not None:
                controls = {control + offset for control in args[i][1]}
                if len(myparty.intersection(controls)) == 0:
                    myparty = myparty.union(controls)
                else:
                    raise ValueError("You can't apply a gate to a qubit and " +
                                     "use it as a control: " +
                                     str(myparty.intersection(controls)))
                anticontrols = {control + offset for control in args[i][2]}
                if len(myparty.intersection(anticontrols)) == 0:
                    myparty = myparty.union(anticontrols)
                else:
                    raise ValueError("You can't apply a gate to a qubit and " +
                                     "use it as a control, or use it as " +
                                     "control and anticontrol at the same " +
                                     "time: " +
                                     str(myparty.intersection(anticontrols)))
            if len(parties.intersection(myparty)) == 0:
                parties = parties.union(myparty)
            else:
                raise ValueError("You can't apply two or more gates to the " +
                                 "same qubit in the same line: " +
                                 str(parties.intersection(myparty)))
            currbit += get_gate_qubits(args[i][0])
        else:
            empties.add(currbit)
            currbit += 1
    if ignore_empty:
        return parties
    return parties.union(empties)


def get_gate_data(gatename):
    """Return number of qubits and number of rows or columns of native gate."""
    if isinstance(gatename, QGate):
        shape = (gatename.size, 2**gatename.size)
    else:
        name, arg1, arg2, arg3, invert = prs.get_gate_data(gatename)
        if arg1 is None:
            arg1 = 0
        if arg2 is None:
            arg2 = 0
        if arg3 is None:
            arg3 = 0
        if "SWAP" in name or name == "XX" or name == "YY" or name == "ZZ":
            shape = (2, 4)
        else:
            qgate = ct.c_void_p(_cGetQGate(ct.c_char_p(name.encode()),
                                           ct.c_double(arg1),
                                           ct.c_double(arg2),
                                           ct.c_double(arg3),
                                           ct.c_int(int(invert))))
            # NOTNULLPTR and True = True, NOTNULLPTR or False = NOTNULLPTR
            if (qgate or False) == qgate:
                nqubits = int(_cGetQGateQubits(qgate))
                size = int(_cGetQGateSize(qgate))
                shape = (nqubits, size)
            else:
                shape = None
                print("Error while getting the specified gate!")

    return shape


'''
def joinGates(gates):
    maxgatelen = max(map(getLines, gates))
    newlines = []
    for i in range(maxgatelen):
        newline = []
        for gate in gates:
            if type(gate) != Measure:
                if i < len(gate.lines):
                    newline += gate.lines[i]
                else:
                    newline += [None for i in range(gate.size)]
            else:
                if i < 1:
                    newline += [gate]
                else:
                    newline += [None for i in range(get_gate_qubits(gate))]
        newlines += [newline]
    return newlines
'''


def get_gate_qubits(gate):
    """Return the number of qubits affected by gate."""
    size = 0
    if isinstance(gate, str):
        size = get_gate_data(gate)[0]
    elif isinstance(gate, Iterable):
        cacs = 0
        if len(gate) > 1:
            if gate[1] is not None:
                cacs += len(gate[1])
        if len(gate) > 2:
            if gate[2] is not None:
                cacs += len(gate[2])
        size = get_gate_data(gate[0])[0] + cacs
    elif isinstance(gate, Measure):
        size = len(gate.mask)
    elif isinstance(gate, QGate):
        size = gate.size
    elif gate is None:
        size = 1
    else:
        raise ValueError(str(gate) + " is not a gate!")
    return size


'''
def get_lines(gate):
    if type(gate) == Measure:
        return 1
    return len(gate.lines)
'''
