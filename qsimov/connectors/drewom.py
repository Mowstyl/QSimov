"""Handles design execution on target machine."""
import qsimov.connectors.qsimovapi as qapi
from qsimov.structures.qcircuit import QCircuit


class Drewom(object):
    """Handle QDesign execution on specified machine."""

    def __init__(self, qmachine="doki",
                 extra=None):
        """Initialize Drewom structure."""
        self.qmachine = qmachine
        if extra is None:
            extra = {}
        self.extra = extra

    def execute(self, qcircuit, iterations=1):
        """Execute given design in a machine iterations times."""
        if not isinstance(qcircuit, QCircuit):
            raise ValueError("qcircuit must be a QCircuit")
        if self.qmachine == "doki":
            return qapi.execute(qcircuit, iterations=iterations, **self.extra)
        else:
            raise NotImplementedError("Drewom does not support that machine")
