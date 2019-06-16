import numpy as np
import re

__rep__ = re.compile("^((?:C\-)+)?([a-zA-Z0-9]+)(\-1)?(\((?:(?:(?:[a-zA-Z]+)|(?:[0-9]+(?:\.[0-9]+)?))\,\s*)*(?:(?:(?:[a-zA-Z]+)|(?:[0-9]+(?:\.[0-9]+)?)))\))?$")

def parseGroups(groups):
    errored = False
    g1 = groups[0].count("C") if groups[0] is not None else 0
    g2 = groups[1]
    g3 = groups[2] is not None
    if groups[3] is not None:
        aux = groups[3][1:-1].split(",")
        g4 = len(aux)
        g5 = []
        for attr in aux:
            attr = attr.strip()
            if "." in attr:
                attr = float(attr)
            elif attr[0] in "0123456789":
                attr = int(attr)
            elif attr.lower() == "pi":
                attr = np.pi
            elif attr.lower() == "tau":
                attr = 2 * np.pi
            elif attr.lower() == "e":
                attr = np.e
            else:
                print (attr)
                errored = True
                break
            g5.append(attr)
    else:
        g4 = 0
        g5 = None

    if not errored:
        return (g1, g2, g3, g4, g5)
    else:
        return None

def getGroups(str_gate):
    res = __rep__.match(str_gate)
    return parseGroups(res.groups()) if res is not None else None
