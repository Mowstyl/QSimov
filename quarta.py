import qsimov as qj
import qalg as qa
import sys
import getopt

def nTof(num):
    gate = qj.CNOT()
    if num > 1:
        for i in range(num - 1):
            gate = qj.CU(gate)
    return gate

def main(argv):
    seed = None
    if (len(argv) >= 2):
        nEntradas = int(argv[0])
        if (int(argv[1]) == 1):
            gate = qa.I(nEntradas + 1)
        else:
            gate = nTof(nEntradas)
    else:
        print ("Usage: python quarta.py <Number of Inputs> <0 for sat. function, 1 for insat.> [seed]")
        return -1

    if (len(argv) >= 3):
        seed = int(argv[2])

    qj.setRandomSeed(seed)

    qr = qj.QRegistry(nEntradas + 1)

    for i in range(10):
        print ("Iteracion " + str(i + 1))
        qr.applyGate(qj.H(nEntradas), qj.I(1))
        qr.applyGate(gate)
        res = qr.measure([0 for i in range(nEntradas)] + [1])[0]
        if (res == 1):
            print ("Satisfactible")
            break
        else:
            print (qr.getState())
    print ("Inversión sobre la media")
    qr.applyGate(qj.IAA(nEntradas), qj.I(1))
    print (qr.getState())

    return 0

if __name__ == "__main__":
    main(sys.argv[1:])
