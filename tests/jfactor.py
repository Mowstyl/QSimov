import time
from fractions import gcd
from random import randrange

def factor2(num): # Algoritmo Las Vegas que factoriza un numero N = p*q, p y q primos
    start = time.process_time()
    a = randrange(1, num) # Inicializo a con un valor en Gnum
    print ("choosen a = ", a)
    p = gcd(a, num) # Compruebo si son coprimos. Si el resultado es distinto a 1 habremos obtenido p
    if (p != 1):
        print ("Time elapsed: ", (time.process_time() - start), " seconds")
        return (p, int(num/p))
    print ("coprime")
    for r in range(1, num):
        if pow(a, r, num) == 1: # Si se cumple que a^r = 1 (mod p*q)
            print ("r = ", r)
            if r % 2 == 0: # Comprobamos si r es par
                print ("r is even")
                x = pow(a, int(r/2), num)
                if x != num - 1: # Comprobamos si x + 1 = 0 (mod p*q)
                    p = gcd(x - 1, num)
                    q = int(num/p)
                    print ("Time elapsed: ", (time.process_time() - start), " seconds")
                    return (p, q)
                else:
                    print ("x + 1 = 0 (Mod n)")
                    break
            else:
                print ("r is odd")
                break
    print ("Time elapsed: ", (time.process_time() - start), " seconds")
    return False
