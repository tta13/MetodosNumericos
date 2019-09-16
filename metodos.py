import sympy

def metodo_euler(y0, t0, h, n, f_string, saida):
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = saida)
	print("h = ", h, file = saida)
	f = sympy.sympify(f_string)
	results = []
	for j in range(0, n+1):
		results.append(y0)
		print(j, y0, file = saida)
		fy = f.subs(y, y0)
		k1 = fy.subs(t, t0)
		y0 += h*k1
		t0 += h
	return results

def metodo_euler_inverso(y0, t0, h, n, f_string, saida):
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = saida)
	print("h = ", h, file = saida)
	results = []
	f = sympy.sympify(f_string)
	for j in range(0, n+1):
		results.append(y0)
		print(j, y0, file = saida)
		t0 += h
		func = f.subs(t, t0)
		k1 = sympy.solveset(y0 + func*h - y, y)
		y0 = k1.args[0]
	return results

def metodo_euler_aprimorado(y0, t0, h, n, f_string, saida):
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = saida)
	print("h = ", h, file = saida)
	f = sympy.sympify(f_string)
	results = []
	for j in range(0, n+1):
		results.append(y0)
		print(j, y0, file = saida)
		fy = f.subs(y, y0)
		k1 = fy.subs(t, t0)
		yeuler = y0 + h*k1
		t0 += h
		fy2 = f.subs(y, yeuler)
		k2 = fy2.subs(t, t0)
		y0 += (h/2)*(k1+k2)
	return results	

def metodo_runge_kutta(y0, t0, h, n, f_string, saida):
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = saida)
	print("h = ", h, file = saida)
	f = sympy.sympify(f_string)
	results = []
	for j in range(0, n+1):
		results.append(y0)
		print(j, y0, file = saida)
		fy = f.subs(y, y0)
		k1 = fy.subs(t, t0)
		fy2 = f.subs(y, y0 + (0.5*h*k1))
		k2 = fy2.subs(t, t0 + (0.5*h))
		fy3 = f.subs(y, y0 + (0.5*h*k2))
		k3 = fy3.subs(t, t0 + (0.5*h))
		fy4 = f.subs(y, y0 + (h*k3))
		k4 = fy4.subs(t, t0 + h)
		y0 += (h/6)*(k1 + (2*k2) + (2*k3) + k4)
		t0 += h
	return results

def metodo_adams_bashforth(y0, t0, h, n, f_string, ordem, y_list, saida):
	y, t = sympy.symbols("y t")
	f = sympy.sympify(f_string)	
	if(ordem == 1):
		y0 = y_list[0]
		t0 = t0
		for j in range(0, n+1):
			print(j, y0, file = saida)
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			y0 += h*k1
			t0 += h
	elif(ordem == 2):
		y0 = y_list[1]
		t0 += h
		fy = f.subs(y, y_list[1])
		k1 = fy.subs(t, t0)
		fy2 = f.subs(y, y_list[0])
		k2 = fy2.subs(t, (t0 - h))
		for j in range(2, n+1):
			y0 += (h/2)*(3*k1 - k2)
			t0 += h
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			print(j, y0, file = saida)
	elif(ordem == 3):
		y0 = y_list[2]
		t0 += 2*h
		fy = f.subs(y, y_list[2])
		k1 = fy.subs(t, t0)
		fy2 = f.subs(y, y_list[1])
		k2 = fy2.subs(t, (t0 - h))
		fy3 = f.subs(y, y_list[0])
		k3 = fy3.subs(t, t0 - 2*h)
		for j in range(3, n+1):
			y0 += (h/12)*(23*k1 - 16*k2 + 5*k3)
			t0 += h
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			print(j, y0, file = saida)
	elif(ordem == 4):
		y0 = y_list[3]
		t0 += 3*h
		fy = f.subs(y, y_list[3])
		k1 = fy.subs(t, t0)
		fy2 = f.subs(y, y_list[2])
		k2 = fy2.subs(t, (t0 - h))
		fy3 = f.subs(y, y_list[1])
		k3 = fy3.subs(t, t0 - 2*h)
		fy4 = f.subs(y, y_list[0])
		k4 = fy4.subs(t, t0 - 3*h)
		for j in range(4, n+1):
			y0 += (h/24)*(55*k1 - 59*k2 + 37*k3 - 9*k4)
			t0 += h
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			print(j, y0, file = saida)
	elif(ordem == 5):
		y0 = y_list[4]
		t0 += 4*h
		fy = f.subs(y, y_list[4])
		k1 = fy.subs(t, t0)
		fy2 = f.subs(y, y_list[3])
		k2 = fy2.subs(t, t0 - h)
		fy3 = f.subs(y, y_list[2])
		k3 = fy3.subs(t, t0 - 2*h)
		fy4 = f.subs(y, y_list[1])
		k4 = fy4.subs(t, t0 - 3*h)
		fy5 = f.subs(y, y_list[0])
		k5 = fy5.subs(t, t0 - 4*h)
		for j in range(5, n+1):
			y0 += (h/720)*(1901*k1 - 2774*k2 + 2616*k3 - 1274*k4 + 251*k5)
			t0 += h
			k5 = k4
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			print(j, y0, file = saida) 
	elif(ordem == 6):
		y0 = y_list[5]
		t0 += 5*h
		fy = f.subs(y, y_list[5])
		k1 = fy.subs(t, t0)
		fy2 = f.subs(y, y_list[4])
		k2 = fy2.subs(t, t0 - h)
		fy3 = f.subs(y, y_list[3])
		k3 = fy3.subs(t, t0 - 2*h)
		fy4 = f.subs(y, y_list[2])
		k4 = fy4.subs(t, t0 - 3*h)
		fy5 = f.subs(y, y_list[1])
		k5 = fy5.subs(t, t0 - 4*h)
		fy6 = f.subs(y, y_list[0])
		k6 = fy6.subs(t, t0 - 5*h)
		for j in range(6, n+1):
			y0 += (h/1440)*(4277*k1 - 7923*k2 + 9982*k3 - 7298*k4 + 2877*k5 - 475*k6)
			t0 += h
			k6 = k5
			k5 = k4
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0) 
			print(j, y0, file = saida)
	elif(ordem == 7):
		y0 = y_list[6]
		t0 += 6*h
		fy = f.subs(y, y_list[6])
		k1 = fy.subs(t, t0)
		fy2 = f.subs(y, y_list[5])
		k2 = fy2.subs(t, t0 - h)
		fy3 = f.subs(y, y_list[4])
		k3 = fy3.subs(t, t0 - 2*h)
		fy4 = f.subs(y, y_list[3])
		k4 = fy4.subs(t, t0 - 3*h)
		fy5 = f.subs(y, y_list[2])
		k5 = fy5.subs(t, t0 - 4*h)
		fy6 = f.subs(y, y_list[1])
		k6 = fy6.subs(t, t0 - 5*h)
		fy7 = f.subs(y, y_list[0])
		k7 = fy7.subs(t, t0 - 6*h)
		for j in range(7, n+1):
			y0 += (h/60480)*(198721*k1 - 447288*k2 + 705549*k3 - 688256*k4 + 
				407139*k5 - 134472*k6 + 19087*k7)
			t0 += h
			k7 = k6
			k6 = k5
			k5 = k4
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0) 
			print(j, y0, file = saida)
	elif(ordem == 8):
		y0 = y_list[7]
		t0 += 7*h
		fy = f.subs(y, y_list[7])
		k1 = fy.subs(t, t0)
		fy2 = f.subs(y, y_list[6])
		k2 = fy2.subs(t, t0 - h)
		fy3 = f.subs(y, y_list[5])
		k3 = fy3.subs(t, t0 - 2*h)
		fy4 = f.subs(y, y_list[4])
		k4 = fy4.subs(t, t0 - 3*h)
		fy5 = f.subs(y, y_list[3])
		k5 = fy5.subs(t, t0 - 4*h)
		fy6 = f.subs(y, y_list[2])
		k6 = fy6.subs(t, t0 - 5*h)
		fy7 = f.subs(y, y_list[1])
		k7 = fy7.subs(t, t0 - 6*h)
		fy8 = f.subs(y, y_list[0])
		k8 = fy8.subs(t, t0 - 7*h)
		for j in range(8, n+1):
			y0 += (h/120960)*(434241*k1 - 1152169*k2 + 2183877*k3 - 2664477*k4 + 
				2102243*k5 - 1041723*k6 + 295767*k7 - 36799*k8)
			t0 += h
			k8 = k7
			k7 = k6
			k6 = k5
			k5 = k4
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			print(j, y0, file = saida) 

def metodo_adams_moulton(y0, t0, h, n, f_string, ordem, y_list, saida):
	y, t = sympy.symbols("y t")
	f = sympy.sympify(f_string)	
	if(ordem == 1):
		y0 = y0
		t0 = t0
		for j in range(0, n+1):
			print(j, y0, file = saida)
			t0 += h
			func = f.subs(t, t0)
			k1 = sympy.solveset(y0 + func*h - y, y)
			y0 = k1.args[0]
	elif(ordem == 2):
		y0 = y_list[0]
		t0 += h
		fy = f.subs(y, y_list[0])
		k1 = fy.subs(t, t0)
		#implicito
		k2 = f.subs(t, t0+h)
		for j in range(1, n+1):
			solution = y0 + (h/2)*(k1 + k2) - y
			solution = sympy.solveset(solution, y)
			y0 = solution.args[0]
			t0 += h
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			k2 = f.subs(t, t0+h)
			print(j, y0, file = saida)
	elif(ordem == 3):
		y0 = y_list[1]
		t0 += h
		fy = f.subs(y, y0)
		k1 = fy.subs(t0)
		fy2 = f.subs(y, y_list[0])
		k2 = fy2.subs(t, t0 - h)
		#implicito
		k3 = f.subs(t, t0+h)
		for j in range(2, n+1):
			solution = y0 + (h/12)*(5*k3 + 8*k1 -k2) - y
			solution = sympy.solveset(solution, y)
			y0 = solution.args[0]
			t0 += h
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			k3 = f.subs(t, t0+h)
			print(j, y0, file = saida)
	elif(ordem == 4):
		y0 = y_list[2]
		t0 += 2*h
		fy = f.subs(y, y0)
		k1 = fy.subs(t, (t0))
		fy2 = f.subs(y, y_list[1])
		k2 = fy2.subs(t, t0 - h)
		fy3 = f.subs(y, y_list[0])
		k3 = fy3.subs(t, t0 - 2*h)
		#implicito
		k4 = f.subs(t, t0+h)
		for j in range(3, n+1):
			solution = y0 + (h/24)*(9*k4 + 19*k1 - 5*k2 + k3) - y
			solution = sympy.solveset(solution, y)
			y0 = solution.args[0]
			t0 += h
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			k4 = f.subs(t, t0+h)
			print(j, y0, file = saida)
	elif(ordem == 5):
		y0 = y_list[3]
		t0 += 3*h
		fy = f.subs(y, y0)
		k1 = fy.subs(t, (t0))
		fy2 = f.subs(y, y_list[2])
		k2 = fy2.subs(t, t0 - h)
		fy3 = f.subs(y, y_list[1])
		k3 = fy3.subs(t, t0 - 2*h)
		fy4 = f.subs(y, y_list[0])
		k4 = fy4.subs(t, (t0 - 3*h))
		#implicito
		k5 = f.subs(t, t0+h)
		for j in range(4, n+1):
			solution = y0 + (h/720)*(251*k5 + 646*k1 - 264*k2 + 106*k3 - 19*k4) - y
			solution = sympy.solveset(solution, y)
			y0 = solution.args[0]
			t0 += h
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, (t0))
			k5 = f.subs(t, t0+h)
			print(j, y0, file = saida)
	elif(ordem == 6):
		y0 = y_list[4]
		t0 += 4*h
		fy = f.subs(y, y0)
		k1 = fy.subs(t, (t0))
		fy2 = f.subs(y, y_list[3])
		k2 = fy2.subs(t, t0 - h)
		fy3 = f.subs(y, y_list[2])
		k3 = fy3.subs(t, t0 - 2*h)
		fy4 = f.subs(y, y_list[1])
		k4 = fy4.subs(t, (t0 - 3*h))
		fy5 = f.subs(y, y_list[0])
		k5 = fy5.subs(t, t0 - 4*h)
		#implicito
		k6 = f.subs(t, t0+h)
		for j in range(5, n+1):
			solution = y0 + (h/1440)*(475*k6 + 1427*k1 - 798*k2 + 482*k3 - 173*k4
				+ 27*k5) - y
			solution = sympy.solveset(solution, y)
			y0 = solution.args[0]
			t0 += h
			k5 = k4
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, (t0))
			k6 = f.subs(t, t0+h)
			print(j, y0, file = saida)
	elif(ordem == 7):
		y0 = y_list[5]
		t0 += 5*h
		fy = f.subs(y, y0)
		k1 = fy.subs(t, (t0))
		fy2 = f.subs(y, y_list[4])
		k2 = fy2.subs(t, t0 - h)
		fy3 = f.subs(y, y_list[3])
		k3 = fy3.subs(t, t0 - 2*h)
		fy4 = f.subs(y, y_list[2])
		k4 = fy4.subs(t, (t0 - 3*h))
		fy5 = f.subs(y, y_list[1])
		k5 = fy5.subs(t, t0 - 4*h)
		fy6 = f.subs(y, y_list[0])
		k6 = fy6.subs(t, t0 - 5*h)
		#implicito
		k7 = f.subs(t, t0+h)
		for j in range(6, n+1):
			solution = y0 + (h/60480)*(19087*k7 + 65112*k1 - 46461*k2 + 35704*k3 - 20211*k4
				+ 6312*k5 - 863*k6) - y
			solution = sympy.solveset(solution, y)
			y0 = solution.args[0]
			t0 += h
			k6 = k5
			k5 = k4
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, (t0))
			k7 = f.subs(t, t0+h)
			print(j, y0, file = saida)
	elif(ordem == 8):
		y0 = y_list[6]
		t0 += 6*h
		fy = f.subs(y, y0)
		k1 = fy.subs(t, (t0))
		fy2 = f.subs(y, y_list[5])
		k2 = fy2.subs(t, t0 - h)
		fy3 = f.subs(y, y_list[4])
		k3 = fy3.subs(t, t0 - 2*h)
		fy4 = f.subs(y, y_list[3])
		k4 = fy4.subs(t, (t0 - 3*h))
		fy5 = f.subs(y, y_list[2])
		k5 = fy5.subs(t, t0 - 4*h)
		fy6 = f.subs(y, y_list[1])
		k6 = fy6.subs(t, t0 - 5*h)
		fy7 = f.subs(y, y_list[0])
		k7 = fy7.subs(t, t0-6*h) 
		#implicito
		k8 = f.subs(t, t0+h)
		for j in range(7, n+1):
			solution = y0 + (h/120960)*(36799*k8 + 139849*k1 - 121797*k2 + 123133*k3 - 88547*k4
				+ 41499*k5 - 11351*k6 + 1375*k7) - y
			solution = sympy.solveset(solution, y)
			y0 = solution.args[0]
			t0 += h
			k7 = k6
			k6 = k5
			k5 = k4
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, (t0))
			k8 = f.subs(t, t0+h)
			print(j, y0, file = saida)

entrada = open("entrada.txt")
saida = open("saida.txt", "w")
lines = entrada.readlines()
for i in lines:
	string_parts = i.split()
	if(string_parts[0] == "euler"):
		print("Metodo de Euler", file = saida)
		metodo_euler(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[4]), string_parts[-1], saida)
		saida.write("\n")
	elif(string_parts[0] == "euler_inverso"):
		print("Metodo de Euler Inverso", file = saida)
		metodo_euler_inverso(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[4]), string_parts[-1], saida)
		saida.write("\n")
	elif(string_parts[0] == "euler_aprimorado"):
		print("Metodo de Euler Aprimorado", file = saida)
		metodo_euler_aprimorado(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[4]), string_parts[-1], saida)
		saida.write("\n")
	elif(string_parts[0] == "runge_kutta"):
		print("Metodo de Runge-Kutta", file = saida)
		metodo_runge_kutta(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[4]), string_parts[-1], saida)
		saida.write("\n")
	elif(string_parts[0] == "adam_bashforth"):
		print("Metodo de Adams-Bashforth", file = saida)
		print("y( ", float(string_parts[-5]), " ) = ", string_parts[1], file = saida)
		print("h = ", string_parts[-4], file = saida)
		y_list = []
		for j in range(1, int(string_parts[-1])+1):
			y_list.append(float(string_parts[j]))
			print(j, string_parts[j], file = saida)
		metodo_adams_bashforth(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida)
		saida.write("\n")
	elif(string_parts[0] == "adam_bashforth_by_euler"):
		print("Metodo Adams-Bashforth por Euler", file = saida)
		y_list = []
		y_list = metodo_euler(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 1, string_parts[-2], saida)
		metodo_adams_bashforth(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida)
		saida.write("\n")
	elif(string_parts[0] == "adam_bashforth_by_euler_inverso"):
		print("Metodo Adams-Bashforth por Euler Inverso", file = saida)
		y_list = []
		y_list = metodo_euler_inverso(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 1, string_parts[-2], saida)
		metodo_adams_bashforth(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida)
		saida.write("\n")
	elif(string_parts[0] == "adam_bashforth_by_euler_aprimorado"):
		print("Metodo Adams-Bashforth por Euler Aprimorado", file = saida)
		y_list = []
		y_list = metodo_euler_aprimorado(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 1, string_parts[-2], saida)
		metodo_adams_bashforth(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida)
		saida.write("\n")
	elif(string_parts[0] == "adam_bashforth_by_runge_kutta"):
		print("Metodo Adams-Bashforth por Runge-Kutta (ordem = ", string_parts[-1], ")", file = saida)
		y_list = []
		y_list = metodo_runge_kutta(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 1, string_parts[-2], saida)
		metodo_adams_bashforth(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida)
		saida.write("\n")
	elif(string_parts[0] == "adam_multon"):
		print("Metodo de Adams-Moulton", file = saida)
		print("y( ", float(string_parts[-5]), " ) = ", string_parts[1], file = saida)
		print("h = ", string_parts[-4], file = saida)
		y_list = []
		for j in range(1, int(string_parts[-1])):
			y_list.append(float(string_parts[j]))
			print(j, string_parts[j], file = saida)
		metodo_adams_moulton(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida)
		saida.write("\n")
	elif(string_parts[0] == "adam_multon_by_euler"):
		print("Metodo Adams-Moulton por Euler", file = saida)
		y_list = []
		y_list = metodo_euler(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 2, string_parts[-2], saida)
		metodo_adams_moulton(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida)
		saida.write("\n")
	elif(string_parts[0] == "adam_multon_by_euler_inverso"):
		print("Metodo Adams-Moulton por Euler Inverso", file = saida)
		y_list = []
		y_list = metodo_euler_inverso(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 2, string_parts[-2], saida)
		metodo_adams_moulton(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida)
		saida.write("\n")
	elif(string_parts[0] == "adam_multon_by_euler_aprimorado"):
		print("Metodo Adams-Moulton por Euler Aprimorado", file = saida)
		y_list = []
		y_list = metodo_euler_aprimorado(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 2, string_parts[-2], saida)
		metodo_adams_moulton(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida)
		saida.write("\n")
	elif(string_parts[0] == "adam_multon_by_runge_kutta"):
		print("Metodo Adams-Moulton por Runge-Kutta (ordem = ", string_parts[-1], ")", file = saida)
		y_list = []
		y_list = metodo_runge_kutta(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 2, string_parts[-2], saida)
		metodo_adams_moulton(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida)
		saida.write("\n")
	else: 
		saida.write("Metodo nao reconhecido\n\n")
entrada.close()
saida.close()