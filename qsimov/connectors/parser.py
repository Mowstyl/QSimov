import numpy as np
import re

__rep__ = re.compile("^([a-zA-Z0-9]+)(\\((?:(?:(?:[a-zA-Z]+)|(?:[0-9]+(?:\\.[0-9]+)?(?:e[\\+\\-][0-9]+)?))\\,\\s*)*(?:(?:(?:[a-zA-Z]+)|(?:[0-9]+(?:\\.[0-9]+)?(?:e[\\+\\-][0-9]+)?)))\\))?(\\-1)?$")


def parseGroups(groups):
    errored = False
    g1 = groups[0]
    g4 = groups[2] is not None
    if groups[1] is not None:
        aux = groups[1][1:-1].split(",")
        g2 = len(aux)
        g3 = []
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
                print(attr)
                errored = True
                break
            g3.append(attr)
    else:
        g2 = 0
        g3 = None

    if not errored:
        return (g1, g2, g3, g4)
    else:
        return None


def getGroups(str_gate):
    res = __rep__.match(str_gate)
    return parseGroups(res.groups()) if res is not None else None


__gateDict__ = {}
__gateDict__["x"] = ("X", 0, 0)
__gateDict__["not"] = ("X", 0, 0)
__gateDict__["sqrtnot"] = ("SqrtX", 0, 0)
__gateDict__["sqrtx"] = ("SqrtX", 0, 0)
__gateDict__["v"] = ("SqrtX", 0, 0)
__gateDict__["y"] = ("Y", 0, 0)
__gateDict__["z"] = ("Z", 0, 0)
__gateDict__["rx"] = ("RX", 1, 1)
__gateDict__["ry"] = ("RY", 1, 1)
__gateDict__["rz"] = ("RZ", 1, 1)
__gateDict__["r"] = ("R", 1, 1)
__gateDict__["phaseshift"] = ("R", 1, 1)
__gateDict__["phasechange"] = ("R", 1, 1)
__gateDict__["runity"] = ("RootPhase", 1, 1)
__gateDict__["rootphase"] = ("RootPhase", 1, 1)
__gateDict__["h"] = ("H", 0, 1)
__gateDict__["u"] = ("U", 1, 3)
__gateDict__["u3"] = ("U", 3, 3)
__gateDict__["u2"] = ("U2", 2, 2)
__gateDict__["u1"] = ("U1", 1, 1)
__gateDict__["d"] = ("HalfDeutsch", 1, 1)
__gateDict__["deutsch"] = ("HalfDeutsch", 1, 1)
__gateDict__["halfdeutsch"] = ("HalfDeutsch", 1, 1)
__gateDict__["partialdeutsch"] = ("HalfDeutsch", 1, 1)
__gateDict__["xx"] = ("XX", 3, 3)
__gateDict__["isingx"] = ("XX", 3, 3)
__gateDict__["isingxx"] = ("XX", 3, 3)
__gateDict__["yy"] = ("YY", 3, 3)
__gateDict__["isingy"] = ("YY", 3, 3)
__gateDict__["isingyy"] = ("YY", 3, 3)
__gateDict__["zz"] = ("ZZ", 3, 3)
__gateDict__["isingz"] = ("ZZ", 3, 3)
__gateDict__["isingzz"] = ("ZZ", 3, 3)
__gateDict__["swap"] = ("SWAP", 2, 2)
__gateDict__["iswap"] = ("ISWAP", 2, 2)
__gateDict__["sqrtswap"] = ("SqrtSWAP", 2, 2)


def getGateData(gateraw):
    gate = None
    if type(gateraw) == str:
        groups = getGroups(gateraw)
        if not (groups is None):
            gatename, nargs, args, invert = groups
            gatename = gatename.lower()
            if gatename in __gateDict__:
                gatemet, minargs, maxargs = __gateDict__[gatename]
                if gatename == "u":
                    if nargs == 3:
                        gatemet = "U"
                        minargs = 3
                    elif nargs == 2:
                        gatemet = "U2"
                        minargs, maxargs = 2, 2
                    elif nargs == 1:
                        gatemet = "U1"
                        minargs, maxargs = 1, 1

                if minargs <= nargs <= maxargs:  # Adoro Python
                    if nargs == 0:
                        gate = (gatemet, None, None, None, invert)
                    elif nargs == 1:
                        gate = (gatemet, args[0], None, None, invert)
                    elif nargs == 2:
                        gate = (gatemet, args[0], args[1], None, invert)
                    else:
                        gate = (gatemet, args[0], args[1], args[2], invert)
                else:
                    #print("Received: " + gateraw)
                    #print("Parsed: " + gate)
                    raise ValueError(gatename + " gate number of args must be between " + str(minargs) + " and " + str(maxargs))
            else:
                raise ValueError(gatename + " can't be used with QSimovAPI")
        else:
            raise ValueError(gateraw + " can't be used with QSimovAPI")
    else:
        raise ValueError("You can only use a string!")
    return gate
