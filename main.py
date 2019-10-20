#função principal do projeto
import singlestep
import multistep
import graficos

entrada = open("entrada.txt")
saida = open("saida.txt", "w")
implicito = int(input("Se voce deseja calcular os metodos a partir de previsao e correcao digite 0: "))
lines = entrada.readlines()
for i in lines:
	string_parts = i.split()
	lista_y = []
	if(string_parts[0] == "euler"):
		print("Metodo de Euler", file = saida)
		lista_y = singlestep.euler(float(string_parts[1]), float(string_parts[2]), 
			float(string_parts[3]), int(string_parts[4]), string_parts[-1], saida, 1) 
		graficos.plotar(string_parts[2], string_parts[3], string_parts[4], lista_y, "Metodo de Euler")
		saida.write("\n")
	elif(string_parts[0] == "euler_inverso"):
		print("Metodo de Euler Inverso", file = saida)
		lista_y = singlestep.euler_inverso(float(string_parts[1]), float(string_parts[2]), 
			float(string_parts[3]), int(string_parts[4]), string_parts[-1], saida, implicito)
		graficos.plotar(string_parts[2], string_parts[3], string_parts[4], lista_y, "Metodo de Euler Inverso")
		saida.write("\n")
	elif(string_parts[0] == "euler_aprimorado"):
		print("Metodo de Euler Aprimorado", file = saida)
		lista_y = singlestep.euler_aprimorado(float(string_parts[1]), float(string_parts[2]), 
			float(string_parts[3]), int(string_parts[4]), string_parts[-1], saida)
		graficos.plotar(string_parts[2], string_parts[3], string_parts[4], lista_y, "Metodo de Euler Aprimorado")
		saida.write("\n")
	elif(string_parts[0] == "runge_kutta"):
		print("Metodo de Runge-Kutta", file = saida)
		lista_y = multistep.runge_kutta(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[4]), string_parts[-1], saida)
		graficos.plotar(string_parts[2], string_parts[3], string_parts[4], lista_y, "Metodo de Runge-Kutta")
		saida.write("\n")
	elif(string_parts[0] == "adam_bashforth"):
		print("Metodo de Adams-Bashforth", file = saida)
		print("y( ", float(string_parts[-5]), " ) = ", string_parts[1], file = saida)
		print("h = ", string_parts[-4], file = saida)
		y_list = []
		for j in range(1, int(string_parts[-1])+1):
			y_list.append(float(string_parts[j]))
			print(j-1, string_parts[j], file = saida)
		lista_y = multistep.adams_bashforth(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, 1)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo de Adams-Bashforth")
		saida.write("\n")
	elif(string_parts[0] == "adam_bashforth_by_euler"):
		print("Metodo Adams-Bashforth por Euler", file = saida)
		y_list = []
		y_list = singlestep.euler(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 1, string_parts[-2], saida, 1)
		lista_y = multistep.adams_bashforth(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, 1)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo de Adams-Bashforth by Euler")
		saida.write("\n")
	elif(string_parts[0] == "adam_bashforth_by_euler_inverso"):
		print("Metodo Adams-Bashforth por Euler Inverso", file = saida)
		y_list = []
		y_list = singlestep.euler_inverso(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 1, string_parts[-2], saida, implicito)
		lista_y = multistep.adams_bashforth(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, 1)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo de Adams-Bashforth by Euler Inverso")
		saida.write("\n")
	elif(string_parts[0] == "adam_bashforth_by_euler_aprimorado"):
		print("Metodo Adams-Bashforth por Euler Aprimorado", file = saida)
		y_list = []
		y_list = singlestep.euler_aprimorado(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 1, string_parts[-2], saida)
		lista_y = multistep.adams_bashforth(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, 1)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo de Adams-Bashforth by Euler Aprimorado")
		saida.write("\n")
	elif(string_parts[0] == "adam_bashforth_by_runge_kutta"):
		print("Metodo Adams-Bashforth por Runge-Kutta (ordem = ", string_parts[-1], ")", file = saida)
		y_list = []
		y_list = multistep.runge_kutta(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 1, string_parts[-2], saida)
		lista_y = multistep.adams_bashforth(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, 1)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo de Adams-Bashforth by Runge-Kutta")
		saida.write("\n")
	elif(string_parts[0] == "adam_multon"):
		print("Metodo de Adams-Moulton", file = saida)
		print("y( ", float(string_parts[-5]), " ) = ", string_parts[1], file = saida)
		print("h = ", string_parts[-4], file = saida)
		y_list = []
		for j in range(1, int(string_parts[-1])):
			y_list.append(float(string_parts[j]))
			print(j-1, string_parts[j], file = saida)
		lista_y = multistep.adams_moulton(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, implicito)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo de Adams-Moulton")
		saida.write("\n")
	elif(string_parts[0] == "adam_multon_by_euler"):
		print("Metodo Adams-Moulton por Euler", file = saida)
		y_list = []
		y_list = singlestep.euler(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 2, string_parts[-2], saida, 1)
		lista_y = multistep.adams_moulton(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, implicito)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo de Adams-Moulton by Euler")
		saida.write("\n")
	elif(string_parts[0] == "adam_multon_by_euler_inverso"):
		print("Metodo Adams-Moulton por Euler Inverso", file = saida)
		y_list = []
		y_list = singlestep.euler_inverso(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 2, string_parts[-2], saida, implicito)
		lista_y = multistep.adams_moulton(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, implicito)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo de Adams-Moulton by Euler Inverso")
		saida.write("\n")
	elif(string_parts[0] == "adam_multon_by_euler_aprimorado"):
		print("Metodo Adams-Moulton por Euler Aprimorado", file = saida)
		y_list = []
		y_list = singlestep.euler_aprimorado(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 2, string_parts[-2], saida)
		lista_y = multistep.adams_moulton(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, implicito)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo de Adams-Moulton by Euler Aprimorado")
		saida.write("\n")
	elif(string_parts[0] == "adam_multon_by_runge_kutta"):
		print("Metodo Adams-Moulton por Runge-Kutta (ordem = ", string_parts[-1], ")", file = saida)
		y_list = []
		y_list = multistep.runge_kutta(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 2, string_parts[-2], saida)
		lista_y = multistep.adams_moulton(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, implicito)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo de Adams-Moulton by Runge-Kutta")
		saida.write("\n")
	elif(string_parts[0] == "formula_inversa"):
		print("Metodo Formula Inversa de Diferenciacao", file = saida)
		print("y( ", float(string_parts[-5]), " ) = ", string_parts[1], file = saida)
		print("h = ", string_parts[-4], file = saida)
		y_list = []
		for j in range(1, int(string_parts[-1])+1):
			y_list.append(float(string_parts[j]))
			print(j-1, string_parts[j], file = saida)
		lista_y = multistep.funcao_inversa(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, implicito)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo Formula Inversa de Diferenciacao")
		saida.write("\n")
	elif(string_parts[0] == "formula_inversa_by_euler"):
		print("Metodo Formula Inversa de Diferenciacao por Euler", file = saida)
		y_list = []
		y_list = singlestep.euler(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 1, string_parts[-2], saida, 1)
		lista_y = multistep.funcao_inversa(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, implicito)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo Formula Inversa de Diferenciacao by Euler")
		saida.write("\n")
	elif(string_parts[0] == "formula_inversa_by_euler_inverso"):
		print("Metodo Formula Inversa de Diferenciacao por Euler Inverso", file = saida)
		y_list = []
		y_list = singlestep.euler_inverso(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 1, string_parts[-2], saida, implicito)
		lista_y = multistep.funcao_inversa(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, implicito)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo Formula Inversa de Diferenciacao by Euler Inversa")
		saida.write("\n")
	elif(string_parts[0] == "formula_inversa_by_euler_aprimorado"):
		print("Metodo Formula Inversa de Diferenciacao por Euler Aprimorado", file = saida)
		y_list = []
		y_list = singlestep.euler_aprimorado(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 1, string_parts[-2], saida)
		lista_y = multistep.funcao_inversa(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, implicito)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo Formula Inversa de Diferenciacao by Euler Aprimorado")
		saida.write("\n")
	elif(string_parts[0] == "formula_inversa_by_runge_kutta"):
		print("Metodo Formula Inversa de Diferenciacao por Runge-Kutta (ordem = ", string_parts[-1], ")", file = saida)
		y_list = []
		y_list = multistep.runge_kutta(float(string_parts[1]), float(string_parts[2]), float(string_parts[3]), 
			int(string_parts[-1]) - 1, string_parts[-2], saida)
		lista_y = multistep.funcao_inversa(float(string_parts[1]), float(string_parts[-5]), float(string_parts[-4]),
			int(string_parts[-3]), string_parts[-2], int(string_parts[-1]), y_list, saida, implicito)
		graficos.plotar(string_parts[-5], string_parts[-4], string_parts[-3], lista_y, "Metodo Formula Inversa de Diferenciacao by Runge-Kutta")
		saida.write("\n")
	else: 
		saida.write("Metodo nao reconhecido\n\n")
		
print("Metodos Calculados no arquivo saida.txt")
entrada.close()
saida.close()