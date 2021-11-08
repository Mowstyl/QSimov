from abc import ABC, abstractmethod


class QBase(ABC):
    @abstractmethod
    def get_num_qubits(self):
        pass
