#Metodos single-step para o projeto

import sympy

def euler(y0, t0, h, n, f_string, saida, printo):
	y, t = sympy.symbols("y t")
	if(printo != 0):
		print("y( ", t0, " ) = ", y0, file = saida)
		print("h = ", h, file = saida)
	f = sympy.sympify(f_string)
	lista_y = []
	for j in range(0, n+1):
		lista_y.append(y0)
		if(printo != 0):
			print(j, y0, file = saida)
		fy = f.subs(y, y0)
		k1 = fy.subs(t, t0)
		y0 += h*k1
		t0 += h
	return lista_y

def euler_inverso(y0, t0, h, n, f_string, saida, implicito):
	y, t = sympy.symbols("y t")
	print("y( ", t0, " ) = ", y0, file = saida)
	print("h = ", h, file = saida)
	lista_y = []
	f = sympy.sympify(f_string)
	if(implicito != '0'):
		for j in range(0, n+1):
			lista_y.append(y0)
			print(j, y0, file = saida)
			t0 += h
			fy = f.subs(t, t0)
			k1 = sympy.solveset(y0 + fy*h - y, y)
			y0 = k1.args[0]
	else:
		for j in range(0, n+1):
			lista_y.append(y0)
			print(j, y0, file = saida)
			t0 += h
			fy = f.subs(t, t0)
			yeuler = euler(y0, t0-h, h, 1, f_string, saida, 0)
			k1 = fy.subs(y, yeuler[-1])
			y0 += k1*h
	return lista_y


def euler_aprimorado(y0, t0, h, n, f_string, saida):
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
		yeuler = y0 + h*k1
		t0 += h
		fy2 = f.subs(y, yeuler)
		k2 = fy2.subs(t, t0)
		y0 += (h/2)*(k1+k2)
	return lista_y
