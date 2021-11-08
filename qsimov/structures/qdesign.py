from abc import abstractmethod
from qsimov.structures.qbase import QBase


class QDesign(QBase):
    @abstractmethod
    def add_operation(self, gate, targets=None, controls=None,
                      anticontrols=None):
        pass

    @abstractmethod
    def get_operations(self):
        pass

    @abstractmethod
    def draw(self):
        pass
