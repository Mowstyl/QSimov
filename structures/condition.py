def _specialCompare(a, b):
    same = len(a) == len(b)
    if (same):
        for i in range(len(a)):
            if a[i] is not None and a[i] != b[i]:
                same = False
                break
    return same

class Condition(object):
    def __init__(self, cond, ifcase, elcase, typeif, typeel):
        # cond is an array of what we expect to have measured in each QuBit. None if we don't care about a certain value. Example: [0, 1, None, None, 1].
        # ifcase and elcase can be Conditions or QCircuits to be applied to the registry. They can also be functions that take the registry and the result as a parameter.
        # typeif and typeel are integers from -1 to 2.
        # -1 means that nothing has to be done in that case.
        # 0 means that ifcase/elcase is another Condition.
        # 1 means that ifcase/elcase is a QGate to be applied.
        # 2 means that ifcase/elcase is a QCircuit.
        self.cond = cond
        self.ifcase = ifcase
        self.elcase = elcase
        self.typeif = typeif
        self.typeel = typeel

    def evaluate(self, qregistry, mresults):
        case = self.elcase
        t = self.typeel
        if _specialCompare(self.cond, mresults):
            case = self.ifcase
            t = self.typeif
        if t == 0: # Condition
            r = case.evaluate(qregistry, mresults)
        elif t == 1: # QGate
            r = qregistry
            r.applyGate(case)
        elif t == 2: # QCircuit
            r = case._executeOnce(qregistry)
        else: # Do nothing
            r = qregistry
        return r
