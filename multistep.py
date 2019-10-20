#Metodos multi-step para o projeto

import sympy
import singlestep

def runge_kutta(y0, t0, h, n, f_string, saida):
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = saida)
	print("h = ", h, file = saida)
	f = sympy.sympify(f_string)
	lista_y = []
	for j in range(0, n+1):
		lista_y.append(y0)
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
	return lista_y

def adams_bashforth(y0, t0, h, n, f_string, ordem, y_list, saida, printo):
	y, t = sympy.symbols("y t")
	f = sympy.sympify(f_string)	
	lista_y = []
	lista_y = y_list.copy()
	if(ordem == 1):
		y0 = y0
		t0 = t0
		for j in range(1, n+1):
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			y0 += h*k1
			t0 += h
			lista_y.append(y0)
			if(printo != 0):
				print(j, y0, file = saida)
	elif(ordem == 2):
		y0 = y_list[1]
		t0 += h
		fy = f.subs(y, y_list[1])
		k1 = fy.subs(t, t0)
		fy2 = f.subs(y, y_list[0])
		k2 = fy2.subs(t, (t0 - h))
		for j in range(2, n+1):
			y0 += (h/2)*(3*k1 - k2)
			lista_y.append(y0)
			t0 += h
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			if(printo != 0):
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
			lista_y.append(y0)
			t0 += h
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			if(printo != 0):
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
			lista_y.append(y0)
			t0 += h
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			if(printo != 0):
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
			lista_y.append(y0)
			t0 += h
			k5 = k4
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0)
			if(printo != 0):
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
			lista_y.append(y0)
			t0 += h
			k6 = k5
			k5 = k4
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0) 
			if(printo != 0):
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
			lista_y.append(y0)
			t0 += h
			k7 = k6
			k6 = k5
			k5 = k4
			k4 = k3
			k3 = k2
			k2 = k1
			fy = f.subs(y, y0)
			k1 = fy.subs(t, t0) 
			if(printo != 0):
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
			lista_y.append(y0)
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
			if(printo != 0):
				print(j, y0, file = saida) 
	return lista_y		

def adams_moulton(y0, t0, h, n, f_string, ordem, y_list, saida, implicito):
	y, t = sympy.symbols("y t")
	f = sympy.sympify(f_string)	
	lista_y = []
	lista_y = y_list.copy()
	if(ordem == 1):
		lista_y = singlestep.euler_inverso(y0, t0, h, n, f_string, saida, implicito)
	elif(ordem == 2):
		y0 = y_list[0]
		t0 += h
		fy = f.subs(y, y_list[0])
		k1 = fy.subs(t, t0)
		if(implicito != 0):
			k2 = f.subs(t, t0+h)
			for j in range(1, n+1):
				solution = y0 + (h/2)*(k1 + k2) - y
				solution = sympy.solveset(solution, y)
				y0 = solution.args[0]
				lista_y.append(y0)
				t0 += h
				fy = f.subs(y, y0)
				k1 = fy.subs(t, t0)
				k2 = f.subs(t, t0+h)
				print(j, y0, file = saida)
		else:
			for j in range(1, n+1):
				fy2 = f.subs(t, t0+h)
				y_bashforth = adams_bashforth(y0, t0-(ordem-2)*h, h, ordem-1, f_string, ordem - 1, y_list, saida, 0)
				k2 = fy2.subs(y, y_bashforth[-1])
				y0 += (h/2)*(k1 + k2)
				t0 += h
				fy = f.subs(y, y0)
				k1 = fy.subs(t, t0)
				print(j, y0, file = saida)
				lista_y.append(y0)
	elif(ordem == 3):
		y0 = y_list[1]
		t0 += h
		fy = f.subs(y, y0)
		k1 = fy.subs(t, t0)
		fy2 = f.subs(y, y_list[0])
		k2 = fy2.subs(t, t0 - h)
		#implicito
		if(implicito != 0):
			k3 = f.subs(t, t0+h)
			for j in range(2, n+1):
				solution = y0 + (h/12)*(5*k3 + 8*k1 -k2) - y
				solution = sympy.solveset(solution, y)
				y0 = solution.args[0]
				lista_y.append(y0)
				t0 += h
				k2 = k1
				fy = f.subs(y, y0)
				k1 = fy.subs(t, t0)
				k3 = f.subs(t, t0+h)
				print(j, y0, file = saida)
		else:
			for j in range(2, n+1):
				fy3 = f.subs(t, t0+h)
				y_bashforth = adams_bashforth(y0, t0-(ordem-2)*h, h, ordem-1, f_string, ordem - 1, y_list, saida, 0)
				k3 = fy3.subs(y, y_bashforth[-1])
				y0 += (h/12)*(5*k3 + 8*k1 -k2)
				lista_y.append(y0)
				t0 += h
				k2 = k1
				fy = f.subs(y, y0)
				k1 = fy.subs(t, t0)
				y_list[0] = y_list[1]
				y_list[1] = y0
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
		if(implicito != 0):
			k4 = f.subs(t, t0+h)
			for j in range(3, n+1):
				solution = y0 + (h/24)*(9*k4 + 19*k1 - 5*k2 + k3) - y
				solution = sympy.solveset(solution, y)
				y0 = solution.args[0]
				lista_y.append(y0)
				t0 += h
				k3 = k2
				k2 = k1
				fy = f.subs(y, y0)
				k1 = fy.subs(t, t0)
				k4 = f.subs(t, t0+h)
				print(j, y0, file = saida)
		else:
			for j in range(3, n+1):
				fy4 = f.subs(t, t0+h)
				y_bashforth = adams_bashforth(y0, t0-(ordem-2)*h, h, ordem-1, f_string, ordem - 1, y_list, saida, 0)
				k4 = fy4.subs(y, y_bashforth[-1])
				y0 = y0 + (h/24)*(9*k4 + 19*k1 - 5*k2 + k3)
				lista_y.append(y0)
				t0 += h
				k3 = k2
				k2 = k1
				fy = f.subs(y, y0)
				k1 = fy.subs(t, t0)
				y_list[0] = y_list[1]
				y_list[1] = y_list[2]
				y_list[2] = y0
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
		if(implicito != 0):
			k5 = f.subs(t, t0+h)
			for j in range(4, n+1):
				solution = y0 + (h/720)*(251*k5 + 646*k1 - 264*k2 + 106*k3 - 19*k4) - y
				solution = sympy.solveset(solution, y)
				y0 = solution.args[0]
				lista_y.append(y0)
				t0 += h
				k4 = k3
				k3 = k2
				k2 = k1
				fy = f.subs(y, y0)
				k1 = fy.subs(t, (t0))
				k5 = f.subs(t, t0+h)
				print(j, y0, file = saida)
		else:
			for j in range(4, n+1):
				fy5 = f.subs(t, t0+h)
				y_bashforth = adams_bashforth(y0, t0-(ordem-2)*h, h, ordem-1, f_string, ordem - 1, y_list, saida, 0)
				k5 = fy5.subs(y, y_bashforth[-1])
				y0 = y0 + (h/720)*(251*k5 + 646*k1 - 264*k2 + 106*k3 - 19*k4)
				lista_y.append(y0)
				t0 += h
				k4 = k3
				k3 = k2
				k2 = k1
				fy = f.subs(y, y0)
				k1 = fy.subs(t, (t0))
				y_list[0] = y_list[1]
				y_list[1] = y_list[2]
				y_list[2] = y_list[3]
				y_list[3] = y0
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
		if(implicito != 0):
			k6 = f.subs(t, t0+h)
			for j in range(5, n+1):
				solution = y0 + (h/1440)*(475*k6 + 1427*k1 - 798*k2 + 482*k3 - 173*k4
					+ 27*k5) - y
				solution = sympy.solveset(solution, y)
				y0 = solution.args[0]
				lista_y.append(y0)
				t0 += h
				k5 = k4
				k4 = k3
				k3 = k2
				k2 = k1
				fy = f.subs(y, y0)
				k1 = fy.subs(t, (t0))
				k6 = f.subs(t, t0+h)
				print(j, y0, file = saida)
		else:
			for j in range(5, n+1):
				fy6 = f.subs(t, t0+h)
				y_bashforth = adams_bashforth(y0, t0-(ordem-2)*h, h, ordem-1, f_string, ordem - 1, y_list, saida, 0)
				k6 = fy6.subs(y, y_bashforth[-1])
				y0 = y0 + (h/1440)*(475*k6 + 1427*k1 - 798*k2 + 482*k3 - 173*k4
					+ 27*k5)
				lista_y.append(y0)
				t0 += h
				k5 = k4
				k4 = k3
				k3 = k2
				k2 = k1
				fy = f.subs(y, y0)
				k1 = fy.subs(t, (t0))
				y_list[0] = y_list[1]
				y_list[1] = y_list[2]
				y_list[2] = y_list[3]
				y_list[3] = y_list[4]
				y_list[4] = y0
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
		if(implicito != 0):
			k7 = f.subs(t, t0+h)
			for j in range(6, n+1):
				solution = y0 + (h/60480)*(19087*k7 + 65112*k1 - 46461*k2 + 35704*k3 - 20211*k4
					+ 6312*k5 - 863*k6) - y
				solution = sympy.solveset(solution, y)
				y0 = solution.args[0]
				lista_y.append(y0)
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
		else:
			for j in range(6, n+1):
				fy7 = f.subs(t, t0+h)
				y_bashforth = adams_bashforth(y0, t0-(ordem-2)*h, h, ordem-1, f_string, ordem - 1, y_list, saida, 0)
				k7 = fy7.subs(y, y_bashforth[-1])
				y0 = y0 + (h/60480)*(19087*k7 + 65112*k1 - 46461*k2 + 35704*k3 - 20211*k4
					+ 6312*k5 - 863*k6)
				lista_y.append(y0)
				t0 += h
				k6 = k5
				k5 = k4
				k4 = k3
				k3 = k2
				k2 = k1
				fy = f.subs(y, y0)
				k1 = fy.subs(t, (t0))
				y_list[0] = y_list[1]
				y_list[1] = y_list[2]
				y_list[2] = y_list[3]
				y_list[3] = y_list[4]
				y_list[4] = y_list[5]
				y_list[5] = y0
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
		if(implicito != 0):
			k8 = f.subs(t, t0+h)
			for j in range(7, n+1):
				solution = y0 + (h/120960)*(36799*k8 + 139849*k1 - 121797*k2 + 123133*k3 - 88547*k4
					+ 41499*k5 - 11351*k6 + 1375*k7) - y
				solution = sympy.solveset(solution, y)
				y0 = solution.args[0]
				lista_y.append(y0)
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
		else:
			for j in range(7, n+1):
				fy8 = f.subs(t, t0+h)
				y_bashforth = adams_bashforth(y0, t0-(ordem-2)*h, h, ordem-1, f_string, ordem - 1, y_list, saida, 0)
				k8 = fy8.subs(y, y_bashforth[-1])
				y0 = y0 + (h/120960)*(36799*k8 + 139849*k1 - 121797*k2 + 123133*k3 - 88547*k4
					+ 41499*k5 - 11351*k6 + 1375*k7)
				lista_y.append(y0)
				t0 += h
				k7 = k6
				k6 = k5
				k5 = k4
				k4 = k3
				k3 = k2
				k2 = k1
				fy = f.subs(y, y0)
				k1 = fy.subs(t, (t0))
				y_list[0] = y_list[1]
				y_list[1] = y_list[2]
				y_list[2] = y_list[3]
				y_list[3] = y_list[4]
				y_list[4] = y_list[5]
				y_list[5] = y_list[6]
				y_list[6] = y0
				print(j, y0, file = saida)
	return lista_y

def funcao_inversa(y0, t0, h, n, f_string, ordem, y_list, saida, implicito):
	y, t = sympy.symbols("y t")
	f = sympy.sympify(f_string)
	t0 += (ordem-1)*h
	lista_y = []
	lista_y = y_list.copy()
	if(ordem == 1):
		for j in range(1, n+1):
			t0 += h
			func = f.subs(t, t0)
			k1 = sympy.solveset(y0 + func*h - y, y)
			y0 = k1.args[0]
			lista_y.append(y0)
			print(j, y0, file = saida)
	elif(ordem == 2):
		y0 = y_list[1]
		y1 = y_list[0]
		for i in range(2, n+1):
			t0 += h
			k1 = f.subs(t, t0)
			solution = (1/3)*(4*y0 - y1 + 2*h*k1) - y 
			solution = sympy.solveset(solution, y)
			y1 = y0
			y0 = solution.args[0]
			lista_y.append(y0)
			print(j, y0, file = saida)
	elif(ordem == 3):
		y0 = y_list[2]
		y1 = y_list[1]
		y2 = y_list[0]
		for j in range(3, n+1):
			t0 += h
			k1 = f.subs(t, t0)
			solution = (1/11)*(18*y0 - 9*y1 + 2*y2 + 6*h*k1) - y
			solution = sympy.solveset(solution, y)
			y2 = y1
			y1 = y0
			y0 = solution.args[0]
			lista_y.append(y0)
			print(j, y0, file = saida)
	elif(ordem == 4):
		y0 = y_list[3]
		y1 = y_list[2]
		y2 = y_list[1]
		y3 = y_list[0]
		for j in range(4, n+1):
			t0 += h
			k1 = f.subs(t, t0)
			solution = (1/25)*(48*y0 - 36*y1 + 16*y2 - 3*y3 + 12*h*k1) - y
			solution = sympy.solveset(solution, y)
			y3 = y2
			y2 = y1
			y1 = y0
			y0 = solution.args[0]
			lista_y.append(y0)
			print(j, y0, file = saida)
	elif(ordem == 5):
		y0 = y_list[4]
		y1 = y_list[3]
		y2 = y_list[2]
		y3 = y_list[1]
		y4 = y_list[0]
		for j in range(5, n+1):
			t0 += h
			k1 = f.subs(t, t0)
			solution = (1/137)*(300*y0 - 300*y1 + 200*y2 - 75*y3 + 12*y4 + 60*h*k1) - y
			solution = sympy.solveset(solution, y)
			y4 = y3
			y3 = y2
			y2 = y1
			y1 = y0
			y0 = solution.args[0]
			lista_y.append(y0)
			print(j, y0, file = saida)
	elif(ordem == 6):
		y0 = y_list[5]
		y1 = y_list[4]
		y2 = y_list[3]
		y3 = y_list[2]
		y4 = y_list[1]
		y5 = y_list[0]
		for j in range(6, n+1):
			t0 += h
			k1 = f.subs(t, t0)
			solution = (1/147)*(360*y0 - 450*y1 + 400*y2 - 225*y3 + 72*y4 - 10*y5 + 60*h*k1) - y
			solution = sympy.solveset(solution, y)
			y5 = y4
			y4 = y3
			y3 = y2
			y2 = y1
			y1 = y0
			y0 = solution.args[0]
			lista_y.append(y0)
			print(j, y0, file = saida)
	return lista_y