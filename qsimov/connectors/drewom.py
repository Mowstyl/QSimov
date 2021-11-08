"""Handles design execution on target machine."""
import numpy as np
import qsimov.connectors.qsimovapi as qapi
from qsimov.structures.qcircuit import QCircuit


class Drewom(object):
    """Handle QDesign execution on specified machine."""

    def __init__(self, qmachine="doki",
                 extra={"num_threads": -1,
                        "random_generator": np.random.rand,
                        "use_system": True,
                        "return_struct": False,
                        "core": None}):
        """Initialize Drewom structure."""
        self.qmachine = qmachine
        self.extra = extra

    def execute(self, qcircuit, iterations=1):
        """Execute given design in a machine iterations times."""
        if not isinstance(qcircuit, QCircuit):
            raise ValueError("qcircuit must be a QCircuit")
        if self.qmachine == "doki":
            return qapi.execute(qcircuit, iterations=iterations, **self.extra)
        else:
            raise NotImplementedError("Drewom does not support that machine")
