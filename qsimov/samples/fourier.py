from qsimov import QGate


def get_QFT(num_qubits):
    QFT = QGate(num_qubits, 0, f"QFT{num_qubits}")
    for i in range(num_qubits):
        QFT.add_operation("H", targets=[i])
        for j in range(num_qubits - i - 1):
            QFT.add_operation(f"RUnity({j+2})", targets=[i], controls={j+i+1})
    return QFT


def get_swapped_QFT(num_qubits):
    QFT = QGate(num_qubits, 0, f"QFT{num_qubits}")
    for i in range(num_qubits):
        QFT.add_operation("H", targets=[num_qubits-i-1])
        for j in range(num_qubits - i - 1):
            QFT.add_operation(f"RUnity({j+2})", targets=[num_qubits-i-1],
                              controls={num_qubits-i-j-2})
    return QFT
