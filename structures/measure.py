class Measure(object):
    def __init__(self, mask, conds=[], remove=False):
        # Mask is a list of 0 and 1, where 0 means not to measure a certain QuBit and 1 means that the QuBit has to be measured
        self.mask = mask
        self.conds = conds
        self.remove = remove

    def __repr__(self):
        return str(["Measure" if i == 1 else "I" for i in self.mask])

    def __str__(self):
        return self.__repr__()

    def _mesToList(self, mresults):
        lin = 0
        mres = []
        for m in self.mask:
            tap = None
            if m == 1:
                tap = mresults[lin]
                lin += 1
            mres.append(tap)
        return mres

    def check(self, qregistry):
        res = qregistry.measure(self.mask, remove=self.remove)
        res = self._mesToList(res)
        # print ("Measure result: " + str(res))
        r = (qregistry, [res])
        for cond in self.conds:
            aux = cond.evaluate(r[0], res)
            if type(aux) == tuple:
                r = (aux[0], r[1] + aux[1])
            else:
                r = (aux, r[1])
        return r
