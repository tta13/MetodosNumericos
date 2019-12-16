#Plotar graficos dos metodos
import matplotlib.pyplot as plt

def plotar(t0, h, n, y, titulo):
	t = []
	t0 = float(t0)
	h = float(h)
	n = int(n)
	for i in range(0, n+1):
		t.append(t0)
		t0 += h
	plt.plot(t, y) 
	plt.xlabel('t') 
	plt.ylabel('y') 
	plt.title(titulo) 
	plt.savefig("Graficos\" + titulo) 
	plt.close()