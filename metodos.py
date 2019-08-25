import sympy

def metodo_euler(y0, t0, h, n, f_string, arq):
	print("Metodo de Euler", file = arq)
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = arq)
	print("h = ", h, file = arq)
	f = sympy.sympify(f_string)
	for j in range(0, n+1):
		print(j, y0, file = arq)
		fy = f.subs(y, y0)
		k1 = fy.subs(t, t0)
		y0 += h*k1
		t0 += h

def metodo_euler_inverso(y0, t0, h, n, f_string, arq):
	print("Metodo de Euler Inverso", file = arq)
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = arq)
	print("h = ", h, file = arq)
	f = sympy.sympify(f_string)
	for j in range(0, n+1):
		print(j, y0, file = arq)
		t0 += h
		func = f.subs(t, t0)
		k1 = sympy.solveset(y0 + func*h - y, y)
		y0 = k1.args[0]

def metodo_euler_aprimorado(y0, t0, h, n, f_string, arq):
	print("Metodo de Euler Aprimorado", file = arq)
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = arq)
	print("h = ", h, file = arq)
	f = sympy.sympify(f_string)
	for j in range(0, n+1):
		print(j, y0, file = arq)
		fy = f.subs(y, y0)
		k1 = fy.subs(t, t0)
		yeuler = y0 + h*k1
		t0 += h
		fy2 = f.subs(y, yeuler)
		k2 = fy2.subs(t, t0)
		y0 += (h/2)*(k1+k2)	

def metodo_runge_kutta(y0, t0, h, n, f_string, arq):
	print("Metodo de Runge-Kutta", file = arq)
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = arq)
	print("h = ", h, file = arq)
	f = sympy.sympify(f_string)
	for j in range(0, n+1):
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


entrada = open("entrada.txt")
lines = entrada.readlines()
saida = open("saida.txt", "w")
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
entrada.close()
saida.close()


