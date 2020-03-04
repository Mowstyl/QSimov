import pyquil as pq
import pyquil.gates as pqg


def execute(qcircuit, processor=None, iterations=1):
    qc = pq.get_qc
