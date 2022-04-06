"""Module with gate name parsing stuff.

This module has all name parsing stuff
"""

import numpy as np
import re

__rep__ = re.compile(r"^([a-zA-Z0-9]+)" +
                     r"(\((?:(?:(?:[a-zA-Z]+)|" +
                     r"(?:[\+\-]?[0-9]+(?:\.[0-9]+)?(?:e[\+\-][0-9]+)?))" +
                     r"\,\s*)*(?:(?:(?:[a-zA-Z]+)|" +
                     r"(?:[\+\-]?[0-9]+(?:\.[0-9]+)?" +
                     r"(?:e[\+\-][0-9]+)?)))\))?(\-1)?$")


def parse_groups(groups):
    """Parse the result of get_groups function, passed as parameter."""
    errored = False
    g1 = groups[0]
    g4 = groups[2] is not None
    if groups[1] is not None:
        aux = groups[1][1:-1].split(",")
        g2 = len(aux)
        g3 = []
        for attr in aux:
            attr = attr.strip()
            if len(attr) == 0:
                errored = True
                break
            is_neg = attr[0] == '-'
            is_pos = attr[0] == '+'
            if is_neg or is_pos:
                attr = attr[1:]
            if len(attr) == 0:
                errored = True
                break
            if "." in attr:
                attr = float(attr)
            elif attr[0] in "0123456789":
                attr = int(attr)
            elif attr.lower() == "pi":
                attr = np.pi
            elif attr.lower() == "tau":
                attr = 2 * np.pi
            elif attr.lower() == "e":
                attr = np.e
            else:
                print(attr)
                errored = True
                break
            if is_neg:
                attr = -attr
            g3.append(attr)
    else:
        g2 = 0
        g3 = None

    if not errored:
        return (g1, g2, g3, g4)
    else:
        return None


def get_groups(str_gate):
    """Get matching groups using __rep__ regular expression."""
    res = __rep__.match(str_gate)
    return parse_groups(res.groups()) if res is not None else None


def Hadamard():
    """Return numpy array with Hadamard gate for 1 qubit."""
    gate = np.ones(4, dtype=complex).reshape(2, 2)
    aux = 1/np.sqrt(2)
    gate[1, 1] = -1
    return gate * aux


def PauliX():
    """Return numpy array with Pauli X gate (aka NOT)."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    gate[0, 1] = 1
    gate[1, 0] = 1
    return gate


def PauliY():
    """Return numpy array with PauliY gate."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    gate[0, 1] = -1j
    gate[1, 0] = 1j
    return gate


def PauliZ():
    """Return numpy array with PauliZ gate."""
    gate = np.eye(2, dtype=complex)
    gate[1, 1] = -1
    return gate


def SqrtX(invert):
    """Return numpy array with sqrt(NOT) gate."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    if not invert:
        gate[0, 0] = 0.5 + 0.5j
        gate[0, 1] = 0.5 - 0.5j
        gate[1, 0] = 0.5 - 0.5j
        gate[1, 1] = 0.5 + 0.5j
    else:
        gate[0, 0] = 0.5 - 0.5j
        gate[0, 1] = 0.5 + 0.5j
        gate[1, 0] = 0.5 + 0.5j
        gate[1, 1] = 0.5 - 0.5j
    return gate


def R(angle, invert):
    """Return numpy array with R gate."""
    gate = np.eye(2, dtype=complex)
    if not invert:
        gate[1, 1] = np.cos(angle) + np.sin(angle) * 1j
    else:
        gate[1, 1] = np.cos(angle) - np.sin(angle) * 1j
    return gate


def RX(angle, invert):
    """Return numpy array with rotation gate around X axis."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    cosan = np.cos(angle/2)
    sinan = np.sin(angle/2)
    if not invert:
        mult = -1j
    else:
        mult = 1j
    gate[0, 0] = cosan
    gate[0, 1] = mult * sinan
    gate[1, 0] = mult * sinan
    gate[1, 1] = cosan
    return gate


def RY(angle, invert):
    """Return numpy array with rotation gate around Y axis."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    cosan = np.cos(angle/2)
    sinan = np.sin(angle/2)
    gate[0, 0] = cosan
    gate[1, 1] = cosan
    if not invert:
        gate[0, 1] = -sinan
        gate[1, 0] = sinan
    else:
        gate[0, 1] = sinan
        gate[1, 0] = -sinan
    return gate


def RZ(angle, invert):
    """Return numpy array with rotation gate around Z axis."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    if not invert:
        gate[0, 0] = np.cos(-angle/2) + np.sin(-angle/2) * 1j
        gate[1, 1] = np.cos(angle/2) + np.sin(angle/2) * 1j
    else:
        gate[0, 0] = np.cos(-angle/2) - np.sin(-angle/2) * 1j
        gate[1, 1] = np.cos(angle/2) - np.sin(angle/2) * 1j
    return gate


def RUnity(n, invert):
    """Return numpy array with nth root of unity rotation gate."""
    return R(2*np.pi/(2**n), invert)


def HalfDeutsch(angle, invert):
    """Return numpy array with a portion of the Deutsch gate.

    This gate, when double controlled, is called Deutsch gate.
    """
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    cosan = np.cos(angle) * 1j
    sinan = np.sin(angle)
    if invert:
        cosan = -cosan
    gate[0, 0] = cosan
    gate[0, 1] = sinan
    gate[1, 0] = sinan
    gate[1, 1] = cosan
    return gate


def U(angle1, angle2, angle3, invert):
    """Return numpy array with U gate (IBM)."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    cosan = np.cos(angle1/2)
    sinan = np.sin(angle1/2)
    mult = 1
    if invert:
        mult = -1
    gate[0, 0] = cosan
    if not invert:
        gate[0, 1] = -sinan * np.cos(angle3) - sinan * np.sin(angle3) * 1j
        gate[1, 0] = sinan * np.cos(angle2) + sinan * np.sin(angle2) * 1j
    else:
        gate[0, 1] = sinan * np.cos(angle2) - sinan * np.sin(angle2) * 1j
        gate[1, 0] = -sinan * np.cos(angle3) + sinan * np.sin(angle3) * 1j
    gate[1, 1] = cosan * np.cos(angle2+angle3) \
        + mult * cosan * np.sin(angle2 + angle3) * 1j
    return gate


def U2(angle1, angle2, invert):
    """Return numpy array with U2 gate (IBM)."""
    gate = np.zeros(4, dtype=complex).reshape(2, 2)
    mult = 1
    if invert:
        mult = -1
    gate[0, 0] = 1
    if not invert:
        gate[0, 1] = -np.cos(angle2) - np.sin(angle2) * 1j
        gate[1, 0] = np.cos(angle1) + np.sin(angle1) * 1j
    else:
        gate[0, 1] = -np.cos(angle2) + np.sin(angle2) * 1j
        gate[1, 0] = np.cos(angle1) - np.sin(angle1) * 1j
    gate[1, 1] = np.cos(angle1+angle2) + mult * np.sin(angle1+angle2) * 1j
    return gate * (1/np.sqrt(2))


def U1(angle, invert):
    """Return numpy array with U1 gate (IBM)."""
    return R(angle, invert)


def SWAP():
    """Return numpy array with SWAP gate."""
    gate = np.zeros(4 * 4, dtype=complex)
    gate = gate.reshape(4, 4)

    gate[0][0] = 1
    gate[1][2] = 1
    gate[2][1] = 1
    gate[3][3] = 1

    return gate


def iSWAP(invert):
    """Return numpy array with iSWAP gate."""
    gate = np.zeros(4 * 4, dtype=complex)
    gate = gate.reshape(4, 4)

    gate[0][0] = 1
    if not invert:
        gate[1][2] = 1j
        gate[2][1] = 1j
    else:
        gate[1][2] = -1j
        gate[2][1] = -1j
    gate[3][3] = 1

    return gate


def sqrtSWAP(invert):
    """Return numpy array with sqrt(SWAP) gate."""
    gate = np.zeros(4 * 4, dtype=complex)
    gate = gate.reshape(4, 4)

    gate[0][0] = 1
    if not invert:
        gate[1][1] = 0.5 + 0.5j
        gate[1][2] = 0.5 - 0.5j
        gate[2][1] = 0.5 - 0.5j
        gate[2][2] = 0.5 + 0.5j
    else:
        gate[1][1] = 0.5 - 0.5j
        gate[1][2] = 0.5 + 0.5j
        gate[2][1] = 0.5 + 0.5j
        gate[2][2] = 0.5 - 0.5j
    gate[3][3] = 1

    return gate


def xx(angle, invert):
    """Return numpy array with Ising Coupling XX gate."""
    gate = np.eye(4, dtype=complex)
    if not invert:
        gate[0, 3] = np.sin(angle)-np.cos(angle)*1j
        gate[1, 2] = -1j
        gate[2, 1] = -1j
        gate[3, 0] = np.sin(-angle)-np.cos(-angle)*1j
    else:
        gate[0, 3] = np.sin(-angle)+np.cos(-angle)*1j
        gate[1, 2] = 1j
        gate[2, 1] = 1j
        gate[3, 0] = np.sin(angle)+np.cos(angle)*1j
    return gate*(1/np.sqrt(2))


def yy(angle, invert):
    """Return numpy array with Ising Coupling YY gate."""
    gate = np.eye(4, dtype=complex)
    gate = gate * np.cos(angle)
    ansin = np.sin(angle) * 1j
    if not invert:
        gate[0, 3] = ansin
        gate[1, 2] = -ansin
        gate[2, 1] = -ansin
        gate[3, 0] = ansin
    else:
        gate[0, 3] = -ansin
        gate[1, 2] = ansin
        gate[2, 1] = ansin
        gate[3, 0] = -ansin
    return gate


def zz(angle, invert):
    """Return numpy array with Ising Coupling ZZ gate."""
    gate = np.eye(4, dtype=complex)
    gate = gate * np.cos(angle)
    phi2 = angle/2
    if not invert:
        gate[0, 0] = np.cos(phi2) + np.sin(phi2) * 1j
        gate[1, 1] = np.cos(-phi2) + np.sin(-phi2) * 1j
        gate[2, 2] = gate[1, 1]
        gate[3, 3] = gate[0, 0]
    else:
        gate[0, 0] = np.cos(phi2) - np.sin(phi2) * 1j
        gate[1, 1] = np.cos(-phi2) - np.sin(-phi2) * 1j
        gate[2, 2] = gate[1, 1]
        gate[3, 3] = gate[0, 0]
    return gate


_gate_alias = {}
_gate_alias["x"] = "X"
_gate_alias["not"] = "X"
_gate_alias["sqrtnot"] = "SqrtX"
_gate_alias["sqrtx"] = "SqrtX"
_gate_alias["v"] = "SqrtX"
_gate_alias["y"] = "Y"
_gate_alias["z"] = "Z"
_gate_alias["rx"] = "RX"
_gate_alias["ry"] = "RY"
_gate_alias["rz"] = "RZ"
_gate_alias["r"] = "R"
_gate_alias["phaseshift"] = "R"
_gate_alias["phasechange"] = "R"
_gate_alias["runity"] = "RootPhase"
_gate_alias["rootphase"] = "RootPhase"
_gate_alias["h"] = "H"
_gate_alias["u"] = ("U1", "U2", "U3")
_gate_alias["u3"] = "U3"
_gate_alias["u2"] = "U2"
_gate_alias["u1"] = "U1"
_gate_alias["d"] = "HalfDeutsch"
_gate_alias["deutsch"] = "HalfDeutsch"
_gate_alias["halfdeutsch"] = "HalfDeutsch"
_gate_alias["partialdeutsch"] = "HalfDeutsch"
_gate_alias["xx"] = "XX"
_gate_alias["isingx"] = "XX"
_gate_alias["isingxx"] = "XX"
_gate_alias["yy"] = "YY"
_gate_alias["isingy"] = "YY"
_gate_alias["isingyy"] = "YY"
_gate_alias["zz"] = "ZZ"
_gate_alias["isingz"] = "ZZ"
_gate_alias["isingzz"] = "ZZ"
_gate_alias["swap"] = "SWAP"
_gate_alias["iswap"] = "ISWAP"
_gate_alias["sqrtswap"] = "SqrtSWAP"

# min_args, max_args, has_invert_arg, is_self_invert
_gate_data = {}
_gate_func = {}
_gate_data["X"] = (0, 0, False, True)
_gate_func["X"] = PauliX
_gate_data["Y"] = (0, 0, False, True)
_gate_func["Y"] = PauliY
_gate_data["Z"] = (0, 0, False, True)
_gate_func["Z"] = PauliZ
_gate_data["H"] = (0, 0, False, True)
_gate_func["H"] = Hadamard
_gate_data["SqrtX"] = (0, 0, True, False)
_gate_func["SqrtX"] = SqrtX
_gate_data["RX"] = (1, 1, True, False)
_gate_func["RX"] = RX
_gate_data["RY"] = (1, 1, True, False)
_gate_func["RY"] = RY
_gate_data["RZ"] = (1, 1, True, False)
_gate_func["RZ"] = RZ
_gate_data["R"] = (1, 1, True, False)
_gate_func["R"] = R
_gate_data["RootPhase"] = (1, 1, True, False)
_gate_func["RootPhase"] = RUnity
_gate_data["U3"] = (3, 3, True, False)
_gate_func["U3"] = U
_gate_data["U2"] = (2, 2, True, False)
_gate_func["U2"] = U2
_gate_data["U1"] = (1, 1, True, False)
_gate_func["U1"] = U1
_gate_data["HalfDeutsch"] = (1, 1, True, False)
_gate_func["HalfDeutsch"] = HalfDeutsch
_gate_data["XX"] = (3, 3, True, False)
_gate_func["XX"] = xx
_gate_data["YY"] = (3, 3, True, False)
_gate_func["YY"] = yy
_gate_data["ZZ"] = (3, 3, True, False)
_gate_func["ZZ"] = zz
_gate_data["SWAP"] = (0, 0, False, True)
_gate_func["SWAP"] = SWAP
_gate_data["ISWAP"] = (0, 0, True, False)
_gate_func["ISWAP"] = iSWAP
_gate_data["SqrtSWAP"] = (0, 0, True, False)
_gate_func["SqrtSWAP"] = sqrtSWAP


def get_gate_data(gateraw):
    """Get the data of the gate associated with the given string."""
    gate = None
    if type(gateraw) == str:
        groups = get_groups(gateraw)
        if groups is not None:
            alias, nargs, args, invert = groups
            alias = alias.lower()
            if alias in _gate_alias:
                gatename = _gate_alias[alias]
                if alias == "u":
                    if 1 <= nargs <= 3:
                        gatename = gatename[nargs - 1]
                    else:
                        raise ValueError("U gate needs from 1 to 3 parameters")
                if gatename in _gate_data:
                    minargs, maxargs, has_inv, self_inv = _gate_data[gatename]
                    if self_inv:  # If gate . gate = identity
                        invert = False
                    if has_inv:  # If function has invert as last argument
                        if args is not None:
                            args = args + [invert]
                        else:
                            args = [invert]
                    if minargs <= nargs <= maxargs:  # Adoro Python
                        gate = (gatename, args, invert, self_inv)
                    else:
                        raise ValueError(gatename + " gate number of args" +
                                         " must be between " + str(minargs) +
                                         " and " + str(maxargs))
                else:
                    raise ValueError("Couldn't find data for gate " + gatename)
            else:
                raise ValueError(gateraw + " is not a gate." +
                                 " Have you already added it?")
        else:
            raise ValueError(gateraw + " can't be used with QSimovAPI")
    else:
        raise ValueError("You can only use a string!")
    return gate
