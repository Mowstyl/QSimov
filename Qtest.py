from qlibcj import *
import sys, getopt

def testhMat(n): # Devuelve una matriz que al ser multiplicada por 1/sqrt(2^n) resulta en la puerta Hadamard para n bits
	H = np.ones((2,2), dtype=complex)
	H[1,1] = -1
	if n > 1:
		H = np.kron(H, testhMat(n - 1))
	return H

def testH(n): # Devuelve una puerta Hadamard para n QuBits
	return testhMat(n) / np.sqrt(2**n)

def hadamardTest(dec):
	rand = np.append(np.random.permutation(np.arange(2,13))[:3], 1)
	results = []
	for n in rand:
		auxr = np.around(H(n).m, decimals=dec) == np.around(testH(n), decimals=dec)
		auxr.shape = 2**(2*n)
		results.append(all(auxr))
	result = all(results)
	if not result:
		print ("Hadamard test failed!")
		return 1
	return 0

def main(argv):
	d = 14
	if (len(argv) > 0):
		try:
			opts, args = getopt.getopt(argv,"hd:")
		except getopt.GetoptError:
			print ('Sintax: Qtest.py -d <decimal precision>')
			sys.exit(2)
		for opt, arg in opts:
			if opt == '-h':
				print ('Qtest.py -d <decimal precision>')
				sys.exit()
			elif opt in ("-d"):
				d = int(arg)
	
	tests = [hadamardTest]
	ntests = len(tests)
	fails = sum(test(d) for test in tests)
	
	code = 0
	if (fails > 0):
		code = 1
	
	print ('Tests run: ' + str(ntests))
	print ('Tests passed: ' + str(ntests - fails))
	print ('Tests failed: ' + str(fails))
	
	sys.exit(code)

if __name__ == "__main__":
	main(sys.argv[1:])
