from qsimov import *
import qalg as qa
import structures.funmatrix as fm
import sys, getopt
import math
import time as t

def truncate(number, digits) -> float:
	stepper = pow(10.0, digits)
	return math.trunc(stepper * number) / stepper

def testhMat(n): # Devuelve una matriz que al ser multiplicada por 1/sqrt(2^n) resulta en la puerta Hadamard para n bits
	H = np.ones((2,2), dtype=complex)
	H[1,1] = -1
	if n > 1:
		H = np.kron(H, testhMat(n - 1))
	return H

def testH(n): # Devuelve una puerta Hadamard para n QuBits
	return testhMat(n) / np.sqrt(2**n)

def matrixBuilder(r, c):
	return [[(np.random.rand()*2-1)+1j*(np.random.rand()*2-1) for i in range(c)] for i in range(r)]

def hadamardTest(dec):
	print("\tHadamard test starting...")
	rand = np.append(np.random.permutation(np.arange(2, 8))[:3], 1)
	results = []
	testsize = len(rand)
	iteration = 0;
	start = t.time()
	for n in rand:
		print("\t\tTest completion: " + str(iteration*100//testsize) + "%")
		auxr = np.around(H(n)[:], decimals=dec) == np.around(testH(n), decimals=dec)
		auxr.shape = 2**(2*n)
		results.append(all(auxr))
		iteration += 1
	result = all(results)
	end = t.time()
	print("\t\tTest completion: 100%")
	if not result:
		print ("\tHadamard test failed!\n")
		return 1
	print ("\tHadamard test passed in " + str(end-start) + " seconds!\n")
	return 0

def fmatrixTest(dec):
	result = True
	print("\tFunctionalMatrix test starting...")

	m = [2] + np.random.permutation([3, 5, 8, 13, 15, 20]).tolist() + [1, 1, 12, 2]
	n = [2] + np.random.permutation([3, 5, 8, 13, 15, 20]).tolist() + [17, 1, 1, 2]
	x = [2] + np.random.permutation([3, 5, 8, 13, 15, 20]).tolist() + [12, 2, 1, 1]
	y = [2] + np.random.permutation([3, 5, 8, 13, 15, 20]).tolist() + [1, 2, 17, 1]

	testsize = len(m)
	totalsize = 0

	for i in range(testsize):
		totalsize += 3 * m[i] * n[i] + x[i] * y[i]

		if ((m[i], n[i]) == (x[i], y[i])):
			totalsize += 6 * m[i] * n[i]
		elif (m[i] == 1 and y[i] == 1):
			totalsize += 6 * n[i] * x[i]
		elif (n[i] == 1 and x[i] == 1):
			totalsize += 6 * m[i] * y[i]

		if (n[i] == x[i]):
			totalsize += 6 * m[i] * y[i]

		totalsize += 6 * m[i] * n[i] * x[i] * y[i]

	comsize = 0;
	start = t.time()
	for i in range(testsize):
		print("\t\tTest completion: " + str(truncate(comsize * 100/totalsize, 2)) + "%")
		usize = m[i] * n[i] + x[i] * y[i]

		A = matrixBuilder(m[i], n[i])
		B = matrixBuilder(x[i], y[i])

		nA = np.array(A)
		nB = np.array(B)

		fA = fm.FunctionalMatrix(lambda i, j: nA[:][i,j], (m[i], n[i]))
		fB = fm.FunctionalMatrix(lambda i, j: nB[:][i,j], (x[i], y[i]))

		if ((np.around(A, decimals=dec) != np.around(fA[:], decimals=dec)).any() or (np.around(B, decimals=dec) != np.around(fB[:], decimals=dec)).any()):
			result = False
			print ("\tFunctionalMatrix construction error!")
			break
		comsize += 3 * fA.size + fB.size
		print("\t\tTest completion: " + str(truncate(comsize * 100/totalsize, 2)) + "%")

		try:
			fAefB = np.around((fA * fB)[:], decimals=dec)
			fAeB = np.around((fA * B)[:], decimals=dec)
			AefB = np.around((A * fB)[:], decimals=dec)
			fAenB = np.around((fA * nB)[:], decimals=dec)
			nAefB = np.around((fB.__rmul__(nA))[:], decimals=dec)
			nAenB = np.around(nA * nB, decimals=dec)
			comsize += 6 * fAefB.size
			if (not all([(fAefB == nAenB).all(), (fAeB == nAenB).all(), (AefB == nAenB).all(), (fAenB == nAenB).all(), (nAefB == nAenB).all()])):
				result = False
				print ("\tFunctionalMatrix element-wise multiplication error!")
				break
		except ValueError:
			if ((m[i] == x[i] and n[i] == y[i]) or (m[i] == 1 and y[i] == 1) or (n[i] == 1 and x[i] == 1)):
				result = False
				print ("\tFunctionalMatrix element-wise multiplication raises exception with correct shapes!")
				break
		except:
			result = False
			print ("\tFunctionalMatrix element-wise multiplication raises " + str(sys.exc_info()[0].__name__) + " exception!")
			print (sys.exc_info()[1])
			break
		print("\t\tTest completion: " + str(truncate(comsize * 100/totalsize, 2)) + "%")

		try:
			fAxfB = np.around((fA @ fB)[:], decimals=dec)
			fAxB = np.around((fA @ B)[:], decimals=dec)
			AxfB = np.around((A @ fB)[:], decimals=dec)
			fAxnB = np.around((fA @ nB)[:], decimals=dec)
			nAxfB = np.around((fB.__rmatmul__(nA))[:], decimals=dec)
			nAxnB = np.around(nA @ nB, decimals=dec)
			comsize += 6 * fAxfB.size
			if (not all([(fAxfB == nAxnB).all(), (fAxB == nAxnB).all(), (AxfB == nAxnB).all(), (fAxnB == nAxnB).all(), (nAxfB == nAxnB).all()])):
				result = False
				print ("\tFunctionalMatrix multiplication error!")
				break
		except ValueError:
			if (n[i] == x[i]):
				result = False
				print ("\tFunctionalMatrix multiplication raises exception with correct shapes!")
				break
		except:
			result = False
			print ("\tFunctionalMatrix multiplication raises " + str(sys.exc_info()[0].__name__) + " exception!")
			print (sys.exc_info()[1])
			break
		print("\t\tTest completion: " + str(truncate(comsize * 100/totalsize, 2)) + "%")

		try:
			fAkfB = np.around(fm.kron(fA, fB)[:], decimals=dec)
			fAkB = np.around(fm.kron(fA, B)[:], decimals=dec)
			AkfB = np.around(fm.kron(A, fB)[:], decimals=dec)
			fAknB = np.around(fm.kron(fA, nB)[:], decimals=dec)
			nAkfB = np.around(fm.kron(nA, fB)[:], decimals=dec)
			nAknB = np.around(np.kron(nA, nB), decimals=dec)
			comsize += 6 * fAkfB.size
			if (not all([(fAkfB == nAknB).all(), (fAkB == nAknB).all(), (AkfB == nAknB).all(), (fAknB == nAknB).all(), (nAkfB == nAknB).all()])):
				result = False
				print ("\tFunctionalMatrix Kronecker product error!")
				break
		except:
			result = False
			print ("\tFunctionalMatrix Kronecker product raises " + str(sys.exc_info()[0].__name__) + " exception!")
			print (sys.exc_info()[1])
			print (fA.shape)
			print (fB.shape)
			break
	end = t.time()
	print("\t\tTest completion: 100.00%")
	if not result:
		print ("\tFunctionalMatrix test failed!\n")
		return 1
	print ("\tFunctionalMatrix test passed in " + str(end-start) + " seconds!\n")
	return 0

def teleTest(dec, save=True):
	result = True
	print("\tTeleportation test starting...")

	gateList = [PauliX() @ PauliY() @ PauliZ(), u3(*(np.random.random_sample(3) * 2 * np.pi)), u3(*(np.random.random_sample(3) * 2 * np.pi))]
	gl = len(gateList) * 2
	iteration = 0
	start = t.time()
	for gate in gateList:
		for val in [0, 1]:
			print("\t\tTest completion: " + str(iteration * 100//gl) + "%") # ¿Trombombolicos o tromboembolicos?
			c = qa.TeleportationCircuit(gate)

			(r, mess) = c.execute([val]) # Se ejecuta el circuito
			exr = QRegistry(1)
			if (val == 1):
				exr.applyGate(PauliX())
			exr.applyGate(gate)

			if (not all(np.around(r.getState(), decimals=dec) == np.around(exr.getState(), decimals=dec))):
				print ("\tError!")
				print ("\tExpected result:\n\t", exr.state, "\n\tResult:\n\t", r.state)
				print ("\tMeasured values:\n\t", str(mess))
				result = False
				break
			iteration += 1
	end = t.time()
	print("\t\tTest completion: 100.00%")
	if not result:
		print ("\tTeleportation test failed!\n")
		return 1
	print ("\tTeleportation test passed in " + str(end-start) + " seconds!\n")
	return 0

def djTest(dec, save=True):
	result = True
	print("\tDeutsch-Josza test starting...")

	oracleList = [(qa.Bal(4), 4, False), (I(4), 4, True), (qa.Bal(6), 6, False), (I(6), 6, True)] # (U_f, size, f is constant)

	ol = 0
	for triplet in oracleList:
		ol += 2**triplet[1]

	iteration = 0
	start = t.time()
	for triplet in oracleList:
		print("\t\tTest completion: " + str(truncate(iteration * 100 / ol, 2)) + "%")
		c = qa.DJAlgCircuit(triplet[1], triplet[0])
		(r, mess) = c.execute([0 for i in range(triplet[1] - 1)]) # Los qubits se inicializan a cero (x1..xn) excepto el ultimo (y), inicializado a uno por el circuito tal y como se indicó en su construccion
		if (all(i == 0 for i in mess[0][:-1]) != triplet[2]):
			print ("\tError with U_f " + str(iteration))
			print ("\tMeasured values:\n\t", str(mess))
			result = False
			break
		iteration += 2**triplet[1]
	end = t.time()
	print("\t\tTest completion: 100.00%")
	if not result:
		print ("\tDeutsch-Josza test failed!\n")
		return 1
	print ("\tDeutsch-Josza test passed in " + str(end-start) + " seconds!\n")
	return 0

def main(argv):
	d = 10
	if (len(argv) > 0):
		try:
			opts, args = getopt.getopt(argv,"hd:")
		except getopt.GetoptError:
			print ('Syntax: Qtest.py -d <decimal precision>')
			sys.exit(2)
		for opt, arg in opts:
			if opt == '-h':
				print ('Qtest.py -d <decimal precision>')
				sys.exit()
			elif opt in ("-d"):
				d = int(arg)

	tests = [hadamardTest, teleTest, djTest]
	ntests = len(tests)

	maxint = np.iinfo(np.int_).max + 1
	maxseed = 2**32
	maxseed = maxint if maxint < maxseed else maxseed
	see = np.random.randint(maxseed)
	print ("Seed for RNG: " + str(see))
	np.random.seed(see)

	fails = 0
	tnum = 1
	print ("Running " + str(ntests) + " tests...")
	for test in tests:
		print ("\nTest " + str(tnum) + "/" + str(ntests))
		fails += test(d)
		tnum += 1

	code = 0
	if (fails > 0):
		code = 1

	print ('\nResults:')
	print ('Tests run: ' + str(ntests))
	print ('Tests passed: ' + str(ntests - fails))
	print ('Tests failed: ' + str(fails))

	sys.exit(code)

if __name__ == "__main__":
	main(sys.argv[1:])
