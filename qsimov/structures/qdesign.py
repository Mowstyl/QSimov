from abc import abstractmethod
from qsimov.structures.qbase import QBase


class QDesign(QBase):
    @abstractmethod
    def add_operation(self, gate, targets=None, outputs=None, controls=None,
                      anticontrols=None, c_controls=None, c_anticontrols=None):
        pass

    @abstractmethod
    def get_operations(self):
        pass

    @abstractmethod
    def draw(self):
        pass

    @abstractmethod
    def get_num_bits(self):
        pass
