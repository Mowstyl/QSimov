from qlibcj import *

def DJAlg(size, U_f, **kwargs): # U_f es el oraculo, que debe tener x1..xn e y como qubits. Tras aplicarlo el qubit y debe valer f(x1..xn) XOR y. El argumento size es n + 1, donde n es el numero de bits de entrada de f.
	rnd.seed(kwargs.get('seed', None)) # Para asegurar la repetibilidad fijamos la semilla antes del experimento.
	r = QRegistry(([0 for i in range(size - 1)] + [1])) # Los qubits se inicializan a cero (x1..xn) excepto el ultimo (y), inicializado a uno
	r.applyGate(H(size)) # Se aplica una compuerta hadamard a todos los qubits
	r.applyGate(U_f) # Se aplica el oraculo
	r.applyGate(H(size - 1), I(1)) # Se aplica una puerta Hadamard a todos los qubits excepto al ultimo
	return r.measure([1 for i in range(size - 1)] + [0]) # Se miden los qubit x, si es igual a 0 la funcion es constante. En caso contrario no lo es.

def ExampleDJCircuit(size, U_f, **kwargs):
	rnd.seed(kwargs.get('seed', None)) # Para asegurar la repetibilidad fijamos la semilla antes del experimento.
	c = DJAlgCircuit(size, U_f, save=kwargs.get('save', True))
	res = c.execute([0 for i in range(size - 1)]) # Los qubits se inicializan a cero (x1..xn) excepto el ultimo (y), inicializado a uno por el circuito tal y como se indicó en su construccion
	print(all(i == 0 for i in res[1][0][:-1]))

	return res # Los qubits se inicializan a cero (x1..xn) excepto el ultimo (y), inicializado a uno por el circuito tal y como se indicó en su construccion

def DJAlgCircuit(size, U_f, save=True): # U_f es el oraculo, que debe tener x1..xn e y como qubits. Tras aplicarlo el qubit y debe valer f(x1..xn) XOR y. El argumento size es n + 1, donde n es el numero de bits de entrada de f.
	c = QCircuit("Deutsch-Josza Algorithm", save=save, ancilla=[1]) # El ultimo QuBit al ejecutar el algoritmo es de ancilla, con su valor a 1
	c.addLine(H(size)) # Se aplica una compuerta hadamard a todos los qubits
	c.addLine(U_f) # Se aplica el oraculo
	c.addLine(H(size - 1), I(1)) # Se aplica una puerta Hadamard a todos los qubits excepto al ultimo

	# f = lambda _, l: print(all(i == 0 for i in l[:-1])) # Funcion que imprimira cierto tras realizar la medida si la funcion es constante

	# c.addLine(Measure([1 for i in range(size - 1)] + [0], tasks=[f])) # Se miden los qubit x, si es igual a 0 la funcion es constante. En caso contrario no lo es.
	c.addLine(Measure([1 for i in range(size - 1)] + [0])) # Se miden los qubit x, si es igual a 0 la funcion es constante. En caso contrario no lo es.

	return c

'''
Crea un oraculo U_f tal y como viene definido en el algoritmo de Deutsch-Josza para una funcion balanceada f: {0,1}^n ---> {0,1}, f(x) = msb(x) (bit mas significativo de x).
El argumento n no es el numero de bits de la entrada de f, sino dicho numero mas 1 (para el qubit de "salida").
'''
def Bal(n):
	b = I(n)
	'''
	Se invierte el valor del qubit y en los casos en los que el bit mas significativo sea 1.
	Una puerta C-NOT serviria de U_f con la definicion de f dada con n = 2. Bal(2) = CNOT().
	'''
	for i in range(int((2**n)/2), (2**n) - 1, 2):
		t = np.copy(b[i,:])
		b[i], b[i+1] = b[i+1, :], t
	return b
'''
U_f generada con n = 3:
1 0 0 0 0 0 0 0
0 1 0 0 0 0 0 0
0 0 1 0 0 0 0 0
0 0 0 1 0 0 0 0
0 0 0 0 0 1 0 0
0 0 0 0 1 0 0 0
0 0 0 0 0 0 0 1
0 0 0 0 0 0 1 0
Las entradas son, empezando por el qubit mas significativo: x1, x2 e y.
Al aplicar el oraculo lo que hara es intercambiar la probabilidad asociada a |100> con la de |101> y la de |110> con |111>.
De forma mas general, la funcion Bal se observa que devuelve siempre una puerta que al ser aplicada a un conjunto x1, ..., xn, y
de qubits aplicara una C-NOT sobre x1 (control) e y (objetivo), dejando el resto de qubits intactos.
De esta forma el oraculo pondra en el qubit y el valor de x1 XOR y. Como para la mitad de las posibles entradas x1 valdra 0
y para la otra mitad 1, la funcion f es balanceada ya que devuelve 0 para la mitad de las posibles entradas y 1 para la otra mitad.
El oraculo U_f a su vez se comporta como se indica en el algoritmo, teniendo que y <- f(x) XOR y.
'''

def Teleportation(qbit, **kwargs): # El qubit que va a ser teleportado. Aunque en un computador cuantico real no es posible ver el valor de un qubit sin que colapse, al ser un simulador se puede. Puede especificarse semilla con seed = <seed>.
	rnd.seed(kwargs.get('seed', None)) # Se fija la semilla antes de comenzar el experimento. En este caso la tomamos por parametro.
	r = QRegistry([qbit, 0, 0]) # Se crea un registro con el qubit que debe ser enviado a Alice, el qubit de Bob y el de Alice, en adelante Q, B y A. B y A estan inicializados a |0>.
	print ("Original registry:\n", r.state) # Se muestra el estado del registro de qubits.
	r.applyGate(I(1), H(1), I(1)) # Se aplica la puerta Hadamard a B, ahora en una superposicion de los estados |0> y |1>, ambos exactamente con la misma probabilidad.
	r.applyGate(I(1), CNOT()) # Se aplica una puerta C-NOT sobre B (control) y A (objetivo).
	print ("With Bell+ state:\n", r.state) # Tras la aplicacion de las anteriores dos puertas tenemos un estado de Bell +, B y A estan entrelazados. Se muestra el valor del registro.
	# Aqui es donde trabajamos con el qubit Q que queremos enviar posteriormente. En este caso de ejemplo le vamos a aplicar Hadamard y despues un cambio de fase de pi/2
	r.applyGate(H(1), I(2))
	r.applyGate(PhaseShift(np.pi/2), I(2))
	# Una vez terminado todo lo que queremos hacerle al QuBit, procedemos a preparar el envio
	r.applyGate(CNOT(), I(1)) # Se aplica una puerta C-NOT sobre Q (control) y B (objetivo).
	r.applyGate(H(1), I(2)) # Se aplica una puerta Hadamard sobre Q.
	print ("\nBefore measurement:\n", r.state) # Se muestra el valor del registro antes de la medida.
	m = r.measure([1,1,0]) # Se miden los qubits Q y B.
	print ("q0 = ", m[0], "\nq1 = ", m[1]) # Se muestra el resultado de la medida
	q0 = 0 # Se crean para ver que la teleportacion se realiza con exito dos qubits, q0 y q1.
	q1 = 0 # Usandolos crearemos un registro con los valores que debe tener si la teleportacion se ha realizado con exito.
	if (m[1] == 1):
		q1 = 1
		r.applyGate(I(2), PauliX()) # Si al medir B obtuvimos un 1, rotamos A en el eje X (Pauli-X o NOT)
	if (m[0] == 1):
		q0 = 1
		r.applyGate(I(2), PauliZ()) # Si al medir Q obtuvimos un 1, rotamos A en el eje Z (Pauli-Z).
	er = QRegistry([q0, q1, qbit]) # Se crea el registro para testeo mencionado anteriormente.
	# Y aplicamos las mismas operaciones para ver que es lo que se debe recibir, en este caso Hadamard y PhaseShift.
	er.applyGate(I(2), H(1))
	er.applyGate(I(2), PhaseShift(np.pi/2))
	print ("\nExpected result:\n", er.state, "\nResult:\n", r.state) # Se muestra el contenido de los registros, tanto el del resultado esperado como el obtenido.
	print ("Assert: " + str(r.state == er.state))
	return r # Se devuelve el registro obtenido tras aplicar el algoritmo.

def TeleportationCircuit(gate, save=True): # Recibe como argumento lo que se va a ejecutar sobre el primer QuBit despues de hacer el estado de Bell con los dos últimos.
	qc = QCircuit("Teleportation", save=save, ancilla=[0, 0])
	qc.addLine(I(1), H(1), I(1))
	qc.addLine(I(1), CNOT())
	# Aqui es donde trabajamos con el qubit Q que queremos enviar posteriormente. Se le aplica la puerta pasada como parámetro
	qc.addLine(gate, I(2))
	# Una vez terminado todo lo que queremos hacerle al QuBit, procedemos a preparar el envio
	qc.addLine(CNOT(), I(1)) # Se aplica una puerta C-NOT sobre Q (control) y B (objetivo).
	qc.addLine(H(1), I(2)) # Se aplica una puerta Hadamard sobre Q.

	c1 = Condition([None, 1, None], PauliX())
	c2 = Condition([1, None, None], PauliZ())

	m = Measure([1, 1, 0], conds=[c1, c2], remove=True)

	qc.addLine(m)

	return qc # Se devuelve el circuito.

def ExampleTC(value, gate, **kwargs): # El valor debe ser 0 o 1, valor inicial del QuBit a teleportar. Gate es la puerta que se va a aplicar sobre el QuBit a teleportar.
	rnd.seed(kwargs.get('seed', None)) # Para asegurar la repetibilidad fijamos la semilla antes del experimento.

	# Diseñamos la puerta que se va a aplicar sobre el QuBit
	#g = QGate()
	#g.addLine(H(1))
	#g.addLine(PhaseShift(np.pi/2))

	c = TeleportationCircuit(gate, save=kwargs.get('save', True))

	r = c.execute([value]) # Se ejecuta el circuito
	exr = QRegistry([value])
	exr.applyGate(gate)
	print ("Expected result:\n", exr.state, "\nResult:\n", r.state)
	print ("Assert: " + str(all((r.state == exr.state)[0])))

def TwoBitSubstractor(nums, **kwargs): # Se pasa como parametro los dos numeros binarios a restar como [A0, A1, B0, B1]. Devuelve el resultado en los qubit de mayor peso y en el tercer qubit indica si ha habido overflow
	rnd.seed(kwargs.get('seed', None)) # Para asegurar la repetibilidad fijamos la semilla antes del experimento.
	r = QRegistry(nums + [0,0,0,0,0,0,0]) # 7 bits de ancilla a 0 son necesarios en esta implementacion
	r.applyGate(I(1), SWAP(), SWAP(), I(6))
	r.applyGate(I(2), SWAP(), SWAP(), I(5))
	r.applyGate(I(3), SWAP(), SWAP(), I(4))
	r.applyGate(I(4), SWAP(), I(5))
	r.applyGate(I(5), Substractor())
	r.applyGate(I(5), SWAP(), I(4))
	r.applyGate(I(4), SWAP(), I(5))
	r.applyGate(I(3), SWAP(), I(6))
	r.applyGate(I(2), SWAP(), I(7))
	r.applyGate(Substractor(), I(5))
	r.applyGate(I(5), SWAP(), I(4))
	r.applyGate(I(4), SWAP(), I(5))
	r.applyGate(I(3), SWAP(), I(6))
	r.applyGate(I(2), SWAP(), I(7))
	r.applyGate(I(1), SWAP(), I(8))
	return r.measure([1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
