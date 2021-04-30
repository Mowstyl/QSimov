"""Module providing Measure related stuff."""


class Measure(object):
    """Data structure that represents a measurement in a circuit."""

    def __init__(self, mask, conds=[], remove=False):
        """Measure constructor.

        Positional arguments:
            mask: list of 0 and 1, where 0 means not to measure a certain
                  qubit and 1 means that the qubit has to be measured
        Keyworded arguments:
            conds: list of Condition objects to be checked after measuring
            remove: whether or not to remove the measured qubits
        """
        self.mask = mask
        self.conds = conds
        self.remove = remove

    def __repr__(self):
        """Return the string representation of the Measure."""
        return str(["Measure" if i == 1 else None for i in self.mask])

    def __str__(self):
        """Return the string representation of the Measure."""
        return self.__repr__()

    def _mesToList(self, mresults):
        """Build a list with result for measured qubits and None for others."""
        lin = 0
        mres = []
        for m in self.mask:
            tap = None
            if m == 1:
                tap = mresults[lin]
                lin += 1
            mres.append(tap)
        return mres

    def check(self, qregistry, optimize=True):
        """Do de measure and evaluate the conditions. Return the results."""
        res = qregistry.measure(self.mask, remove=self.remove)
        res = self._mesToList(res)
        # print ("Measure result: " + str(res))
        r = (qregistry, [res])
        for cond in self.conds:
            aux = cond.evaluate(r[0], res, optimize=optimize)
            if type(aux) == tuple:
                r = (aux[0], r[1] + aux[1])
            else:
                r = (aux, r[1])
        return r
