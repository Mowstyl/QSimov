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

"""Deutsch-Jozsa with n qubit."""
from qsimov import QGate

import sympy as sp


def first_column_gate(num_qubits, vector):
    """Create a gate with normalized vector as first matrix column."""
    gates = get_gates(vector, num_qubits)
    # print("Gates:", gates)
    # print("V:", vector)
    vector_oracle = QGate(num_qubits, "Vector oracle")
    if len(gates) > 0:
        for angle, qid, controls, anticontrols in gates:
            vector_oracle.add_operation("U(" + str(angle * 2) + ",0,pi)",
                                        targets=qid, controls=controls,
                                        anticontrols=anticontrols)
    # print("Vector:", vector)
    # print("Gates:", gates)
    return vector_oracle


def get_gates(vector, num_qubits):
    """Get all gates from the vector."""
    gates = []
    for j in range(1, num_qubits + 1):
        for i in range(2**(j - 1)):
            control_string = ""
            format_string = "{:0" + str(j - 1) + "b}"
            if j > 1:
                control_string = format_string.format(i)
            angle = get_angle_from_j(vector, j, num_qubits, control_string)
            controls = {num_qubits - i - 1
                        for i in range(len(control_string))
                        if control_string[i] == "1"}
            anticontrols = {num_qubits - i - 1
                            for i in range(len(control_string))
                            if control_string[i] == "0"}
            if angle != 0:
                gates.append((angle, num_qubits - j, controls, anticontrols))

    return gates


def get_angle_from_j(vector, j, n, control_string):
    """Get the rotation angle."""
    numerator = 0
    denominator = 0
    is_last = j == n
    last_imag = False
    if not is_last:
        for id in get_ids_from_j(j, n, control_string, True):
            numerator += abs(vector[id])**2
        for id in get_ids_from_j(j, n, control_string, False):
            denominator += abs(vector[id])**2
    else:
        numerator = vector[get_ids_from_j(j, n, control_string, True)[0]]
        denominator = vector[get_ids_from_j(j, n, control_string, False)[0]]
        if numerator.imag != 0 or denominator.imag != 0:
            last_imag = True
            numerator = abs(numerator)**2
            denominator = abs(denominator)**2
    angle = 0
    if numerator != 0 and denominator != 0:
        if not is_last or last_imag:
            numerator = sp.sqrt(numerator)
            denominator = sp.sqrt(denominator)
        angle = sp.atan2(numerator, denominator).evalf()
    else:  # 0/0
        subvector = vector
        for i in range(len(control_string)):
            bit_value = int(control_string[i])
            bit_id = n - i - 1
            subvector = get_bit_elements(subvector, bit_id, bit_value)
        subvector = get_bit_elements(subvector, n - j, 0)
        if all([subvector[i] == 0 for i in range(len(subvector))]):
            angle = (sp.pi/2).evalf()
        else:
            angle = 0

    return angle


def get_ids_from_j(j, n, control_string, is_numerator):
    """Get the id of the vector items."""
    static_str = control_string + str(int(is_numerator))
    bits = n - j
    if bits > 0:
        bit_format = "{:0" + str(bits) + "b}"
        ids = [int(static_str + bit_format.format(i), 2)
               for i in range(2**bits)]
    else:
        ids = [int(static_str, 2)]
    return ids


def get_bit_elements(vector, bit, value):
    """Return a list with the elements of vector that has bit with value."""
    subvector = []
    found = 0
    max_found = len(vector)//2
    for i in range(len(vector)):
        bin_i_rev = ("{:0" + str(2**bit) + "b}").format(i)[::-1]
        if int(bin_i_rev[bit]) == value:
            subvector.append(vector[i])
            found += 1
            if found == max_found:
                break
    return subvector
