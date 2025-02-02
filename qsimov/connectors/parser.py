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

"""Module with gate name parsing stuff.

This module has all name parsing stuff
"""

import re
import sympy as sp

from sympy.matrices import Matrix
from sympy.parsing.sympy_parser import parse_expr


_gate_name_re = r"[a-zA-Z0-9_]+"
__rep__ = re.compile(r"^(" + _gate_name_re + r")(\(.*\))?(\-1)?$")


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
            attr = parse_expr(attr)
            if not attr.is_number:
                print("Expression", attr, " is not a number")
                errored = True
                break
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
    """Return sympy matrix with Hadamard gate for 1 qubit."""
    return (sp.sqrt(2) / 2) * Matrix([[1, 1], [1, -1]])


def PauliX():
    """Return sympy matrix with Pauli X gate (aka NOT)."""
    return Matrix([[0, 1], [1, 0]])


def PauliY():
    """Return sympy matrix with PauliY gate."""
    return Matrix([[0, -1j], [1j, 0]])


def PauliZ():
    """Return sympy matrix with PauliZ gate."""
    return Matrix([[1, 0], [0, -1]])


def SqrtX():
    """Return sympy matrix with sqrt(NOT) gate."""
    return Matrix([[0.5 + 0.5j, 0.5 - 0.5j], [0.5 - 0.5j, 0.5 + 0.5j]])


def P(angle):
    """Return sympy matrix with PhaseChange gate."""
    _R = sp.eye(2)
    _R[1, 1] = sp.exp(sp.I * angle)
    return _R


def R(theta, phi):
    """Return sympy matrix with R gate."""
    _R = sp.eye(2) * sp.cos(theta / 2)
    _isinth2 = -sp.I * sp.sin(theta / 2)
    _R[0, 1] = sp.exp(-sp.I * phi) * _isinth2
    _R[1, 0] = sp.exp(sp.I * phi) * _isinth2
    return _R


def RX(angle):
    """Return sympy matrix with rotation gate around X axis."""
    gate = Matrix([[0, 0], [0, 0]])
    cosan = sp.cos(angle / 2)
    sinan = -1j * sp.sin(angle / 2)
    gate[0, 0] = cosan
    gate[0, 1] = sinan
    gate[1, 0] = sinan
    gate[1, 1] = cosan
    return gate


def RY(angle):
    """Return sympy matrix with rotation gate around Y axis."""
    gate = Matrix([[0, 0], [0, 0]])
    cosan = sp.cos(angle / 2)
    sinan = sp.sin(angle / 2)
    gate[0, 0] = cosan
    gate[1, 1] = cosan
    gate[0, 1] = -sinan
    gate[1, 0] = sinan
    return gate


def RZ(angle):
    """Return sympy matrix with rotation gate around Z axis."""
    gate = Matrix([[0, 0], [0, 0]])
    gate[0, 0] = sp.cos(-angle / 2) + sp.sin(-angle / 2) * 1j
    gate[1, 1] = sp.cos(angle / 2) + sp.sin(angle / 2) * 1j
    return gate


def RUnity(n):
    """Return sympy matrix with nth root of unity rotation gate."""
    return P(2*sp.pi/(2**n))


def HalfDeutsch(angle):
    """Return sympy matrix with a portion of the Deutsch gate.

    This gate, when double controlled, is called Deutsch gate.
    """
    gate = Matrix([[0, 0], [0, 0]])
    cosan = sp.cos(angle) * sp.I
    sinan = sp.sin(angle)
    gate[0, 0] = cosan
    gate[0, 1] = sinan
    gate[1, 0] = sinan
    gate[1, 1] = cosan
    return gate


def U(th, ph, la):
    """Return sympy matrix with U(θ, φ, λ) gate (IBM)."""
    gate = Matrix([[0, 0], [0, 0]])
    costh2 = sp.cos(th / 2)
    sinth2 = sp.sin(th / 2)
    gate[0, 0] = costh2
    gate[0, 1] = -sp.exp(sp.I * la) * sinth2
    gate[1, 0] = sp.exp(sp.I * ph) * sinth2
    gate[1, 1] = sp.exp(sp.I * (ph + la)) * costh2
    return gate


def U3(th, ph, la):
    """Return sympy matrix with U3(θ, φ, λ) gate (IBM)."""
    print("[WARNING] This gate has been deprecated in OpenQASM standard. Use U(θ, φ, λ) instead.")
    return U(th, ph, la)


def U2(ph, la):
    """Return sympy matrix with U2 gate (IBM)."""
    print("[WARNING] This gate has been deprecated in OpenQASM standard. Use U(π/2, φ, λ) instead.")
    return U(sp.pi / 2, ph, la)


def U1(angle):
    """Return sympy matrix with U1 gate (IBM)."""
    print("[WARNING] This gate has been deprecated in OpenQASM standard. Use U(0, 0, λ) or P(λ) instead.")
    return P(angle)


def SWAP():
    """Return sympy matrix with SWAP gate."""
    gate = sp.zeros(4, 4)
    gate[0, 0] = 1
    gate[1, 2] = 1
    gate[2, 1] = 1
    gate[3, 3] = 1
    return gate


def iSWAP():
    """Return sympy matrix with iSWAP gate."""
    gate = sp.zeros(4, 4)
    gate[0, 0] = 1
    gate[1, 2] = 1j
    gate[2, 1] = 1j
    gate[3, 3] = 1
    return gate


def fSWAP():
    """Return sympy matrix with fermionic swap fSWAP gate."""
    gate = sp.zeros(4, 4)
    gate[0, 0] = 1
    gate[1, 2] = 1
    gate[2, 1] = 1
    gate[3, 3] = -1
    return gate


def sqrtSWAP():
    """Return sympy matrix with sqrt(SWAP) gate."""
    gate = sp.zeros(4, 4)
    gate[0, 0] = 1
    gate[1, 1] = 0.5 + 0.5j
    gate[1, 2] = 0.5 - 0.5j
    gate[2, 1] = 0.5 - 0.5j
    gate[2, 2] = 0.5 + 0.5j
    gate[3, 3] = 1
    return gate


def xx(angle):
    """Return sympy matrix with Ising Coupling XX gate. AKA Mølmer–Sørensen gate"""
    gate = sp.eye(4)
    phi2 = angle / 2
    gate = gate * sp.cos(phi2)
    _isinphi2 = -1j * sp.sin(phi2)
    gate[0, 3] = _isinphi2
    gate[1, 2] = _isinphi2
    gate[2, 1] = _isinphi2
    gate[3, 0] = _isinphi2
    return gate


def yy(angle):
    """Return sympy matrix with Ising Coupling YY gate."""
    gate = sp.eye(4)
    phi2 = angle / 2
    gate = gate * sp.cos(phi2)
    isinphi2 = 1j * sp.sin(phi2)
    gate[0, 3] = isinphi2
    gate[1, 2] = -isinphi2
    gate[2, 1] = -isinphi2
    gate[3, 0] = isinphi2
    return gate


def zz(angle):
    """Return sympy matrix with Ising Coupling ZZ gate."""
    gate = sp.eye(4)
    phi2 = angle / 2
    cosphi2 = sp.cos(phi2)
    isinphi2 = 1j * sp.sin(phi2)
    gate[0, 0] = cosphi2 - isinphi2
    gate[1, 1] = cosphi2 + isinphi2
    gate[2, 2] = cosphi2 + isinphi2
    gate[3, 3] = cosphi2 - isinphi2
    return gate


def xy(angle):
    """Return sympy matrix with Ising Coupling ZZ gate."""
    gate = sp.eye(4)
    phi2 = angle / 2
    cosphi2 = sp.cos(phi2)
    _isinphi2 = -1j * sp.sin(phi2)
    gate[1, 2] = _isinphi2
    gate[1, 1] = cosphi2
    gate[2, 2] = cosphi2
    gate[2, 1] = _isinphi2
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
_gate_alias["p"] = "P"
_gate_alias["phaseshift"] = "P"
_gate_alias["phasechange"] = "P"
_gate_alias["runity"] = "RootPhase"
_gate_alias["rootphase"] = "RootPhase"
_gate_alias["h"] = "H"
_gate_alias["u"] = "U"
_gate_alias["u3"] = "U3"
_gate_alias["u2"] = "U2"
_gate_alias["u1"] = "U1"
_gate_alias["d"] = "HalfDeutsch"
_gate_alias["deutsch"] = "HalfDeutsch"
_gate_alias["halfdeutsch"] = "HalfDeutsch"
_gate_alias["partialdeutsch"] = "HalfDeutsch"
_gate_alias["xx"] = "XX"
_gate_alias["rxx"] = "XX"
_gate_alias["yy"] = "YY"
_gate_alias["ryy"] = "YY"
_gate_alias["zz"] = "ZZ"
_gate_alias["rzz"] = "ZZ"
_gate_alias["xy"] = "XY"
_gate_alias["rxy"] = "XY"
_gate_alias["swap"] = "SWAP"
_gate_alias["iswap"] = "iSWAP"
_gate_alias["fswap"] = "fSWAP"
_gate_alias["sqrtswap"] = "SqrtSWAP"

# min_args, max_args, has_invert_arg, is_self_invert
_gate_data = {}
_gate_func = {}
_gate_data["X"] = (0, 0, True)
_gate_func["X"] = PauliX
_gate_data["Y"] = (0, 0, True)
_gate_func["Y"] = PauliY
_gate_data["Z"] = (0, 0, True)
_gate_func["Z"] = PauliZ
_gate_data["H"] = (0, 0, True)
_gate_func["H"] = Hadamard
_gate_data["SqrtX"] = (0, 0, False)
_gate_func["SqrtX"] = SqrtX
_gate_data["RX"] = (1, 1, False)
_gate_func["RX"] = RX
_gate_data["RY"] = (1, 1, False)
_gate_func["RY"] = RY
_gate_data["RZ"] = (1, 1, False)
_gate_func["RZ"] = RZ
_gate_data["R"] = (2, 2, False)
_gate_func["R"] = R
_gate_data["P"] = (1, 1, False)
_gate_func["P"] = P
_gate_data["RootPhase"] = (1, 1, False)
_gate_func["RootPhase"] = RUnity
_gate_data["U"] = (3, 3, False)
_gate_func["U"] = U
_gate_data["U3"] = (3, 3, False)
_gate_func["U3"] = U3
_gate_data["U2"] = (2, 2, False)
_gate_func["U2"] = U2
_gate_data["U1"] = (1, 1, False)
_gate_func["U1"] = U1
_gate_data["HalfDeutsch"] = (1, 1, False)
_gate_func["HalfDeutsch"] = HalfDeutsch
_gate_data["XX"] = (1, 1, False)
_gate_func["XX"] = xx
_gate_data["YY"] = (1, 1, False)
_gate_func["YY"] = yy
_gate_data["ZZ"] = (1, 1, False)
_gate_func["ZZ"] = zz
_gate_data["XY"] = (1, 1, False)
_gate_func["XY"] = xy
_gate_data["SWAP"] = (0, 0, True)
_gate_func["SWAP"] = SWAP
_gate_data["iSWAP"] = (0, 0, False)
_gate_func["iSWAP"] = iSWAP
_gate_data["fSWAP"] = (0, 0, False)
_gate_func["fSWAP"] = fSWAP
_gate_data["SqrtSWAP"] = (0, 0, False)
_gate_func["SqrtSWAP"] = sqrtSWAP


def get_available_gates():
    return set(_gate_data.keys())


def get_gate_aliases():
    return _gate_alias.copy()


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
                if gatename in _gate_data:
                    minargs, maxargs, self_inv = _gate_data[gatename]
                    if self_inv:  # If gate . gate = identity
                        invert = False
                    if minargs <= nargs <= maxargs:  # Adoro Python
                        gate = (gatename, args, invert, self_inv)
                    else:
                        if gatename == "R" and nargs == 1:
                            print("[ERROR] The old R gate is now called P in OpenQASM standard. Thus, we did the same.")
                        raise ValueError(f"{gatename} gate number of args" +
                                         f" must be between {minargs} and {maxargs}." +
                                         f" Got {nargs} from {gateraw}")
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
