import sympy

def metodo_euler(y0, t0, h, n, f_string, arq):
	print("Metodo de Euler", file = arq)
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = arq)
	print("h = ", h, file = arq)
	f = sympy.sympify(f_string)
	results = []
	for j in range(0, n+1):
		results.append(y0)
		print(j, y0, file = arq)
		fy = f.subs(y, y0)
		k1 = fy.subs(t, t0)
		y0 += h*k1
		t0 += h
	return results

def metodo_euler_inverso(y0, t0, h, n, f_string, arq):
	print("Metodo de Euler Inverso", file = arq)
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = arq)
	print("h = ", h, file = arq)
	results = []
	f = sympy.sympify(f_string)
	for j in range(0, n+1):
		results.append(y0)
		print(j, y0, file = arq)
		t0 += h
		func = f.subs(t, t0)
		k1 = sympy.solveset(y0 + func*h - y, y)
		y0 = k1.args[0]
	return results

def metodo_euler_aprimorado(y0, t0, h, n, f_string, arq):
	print("Metodo de Euler Aprimorado", file = arq)
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = arq)
	print("h = ", h, file = arq)
	f = sympy.sympify(f_string)
	results = []
	for j in range(0, n+1):
		results.append(y0)
		print(j, y0, file = arq)
		fy = f.subs(y, y0)
		k1 = fy.subs(t, t0)
		yeuler = y0 + h*k1
		t0 += h
		fy2 = f.subs(y, yeuler)
		k2 = fy2.subs(t, t0)
		y0 += (h/2)*(k1+k2)
	return results	

def metodo_runge_kutta(y0, t0, h, n, f_string, arq):
	print("Metodo de Runge-Kutta", file = arq)
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = arq)
	print("h = ", h, file = arq)
	f = sympy.sympify(f_string)
	results = []
	for j in range(0, n+1):
		results.append(y0)
		print(j, y0, file = arq)
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

def metodo_adams_bashforth(y0, t0, h, n, f_string, ordem, y_list, arq):
	print("Metodo de Adam-Bashforth", file = arq)
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = arq)
	print("h = ", h, file = arq)
	f = sympy.sympify(f_string)	
	if(ordem == 1):
		for j in range(0, n+1):
			print(j, y0, file = arq)
			fy = f.subs(y, y_list[0])
			k1 = fy.subs(t, t0)
			y0 += h*k1
			t0 += h
			y_list[0] = y0
	elif(ordem == 2):
		for j in range(0, n+1):
			print(j, y0, file = arq)
			fy = f.subs(y, y_list[1])
			k1 = fy.subs(t, t0)
			fy2 = f.subs(y, y_list[0])
			k2 = fy2.subs(t, (t0 - h))
			y0 += (h/2)*(3*k1 - k2)
			t0 += h
			y_list[0] = y_list[1]
			y_list[1] = y0
	elif(ordem == 3):
		for j in range(0, n+1):
			print(j, y0, file = arq)
			fy = f.subs(y, y_list[2])
			k1 = fy.subs(t, t0)
			fy2 = f.subs(y, y_list[1])
			k2 = fy2.subs(t, (t0 - h))
			fy3 = f.subs(y, y_list[0])
			k3 = fy3.subs(t, t0 - 2*h)
			y0 += (h/12)*(23*k1 - 16*k2 + 5*k3)
			t0 += h
			y_list[0] = y_list[1]
			y_list[1] = y_list[2]
			y_list[2] = y0
	elif(ordem == 4):
		for j in range(0, n+1):
			print(j, y0, file = arq)
			fy = f.subs(y, y_list[3])
			k1 = fy.subs(t, t0)
			fy2 = f.subs(y, y_list[2])
			k2 = fy2.subs(t, (t0 - h))
			fy3 = f.subs(y, y_list[1])
			k3 = fy3.subs(t, t0 - 2*h)
			fy4 = f.subs(y, y_list[0])
			k4 = fy4.subs(t, t0 - 3*h)
			y0 += (h/24)*(55*k1 - 59*k2 + 37*k3 - 9*k4)
			t0 += h
			y_list[0] = y_list[1]
			y_list[1] = y_list[2]
			y_list[2] = y_list[3]
			y_list[3] = y0
	elif(ordem == 5):
		y0 = y_list[4]
		t0 += 4*h
		fy = f.subs(y, y_list[4])
		k1 = fy.subs(t, t0)
		print(k1)
		fy2 = f.subs(y, y_list[3])
		k2 = fy2.subs(t, t0 - h)
		print(k2)
		fy3 = f.subs(y, y_list[2])
		k3 = fy3.subs(t, t0 - 2*h)
		print(k3)
		fy4 = f.subs(y, y_list[1])
		k4 = fy4.subs(t, t0 - 3*h)
		print(k4)
		fy5 = f.subs(y, y_list[0])
		k5 = fy5.subs(t, t0 - 4*h)
		print(k5)
		for j in range(4, n+1):
			print(j, y0, file = arq)
			y0 += (h/720)*(1901*k1 - 2774*k2 + 2616*k3 - 1274*k4 + 251*k5)
			t0 += h
			k5 = k4
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0) 


entrada = open("entrada.txt")
saida = open("saida.txt", "w")
lines = entrada.readlines()
for i in lines:
	string_parts = i.split()
	if(string_parts[0] == "euler"):
		metodo_euler(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[4]), string_parts[-1], saida)
		saida.write("\n")
	elif(string_parts[0] == "euler_inverso"):
		metodo_euler_inverso(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[4]), string_parts[-1], saida)
		saida.write("\n")
	elif(string_parts[0] == "euler_aprimorado"):
		metodo_euler_aprimorado(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[4]), string_parts[-1], saida)
		saida.write("\n")
	elif(string_parts[0] == "runge_kutta"):
		metodo_runge_kutta(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[4]), string_parts[-1], saida)
		saida.write("\n")
	elif(string_parts[0] == "adam_bashforth"):
		y_list = []
		for j in range(1, int(string_parts[-1])+1):
			y_list.append(float(string_parts[j]))
		metodo_adams_bashforth(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida)
		saida.write("\n")
entrada.close()
saida.close()


