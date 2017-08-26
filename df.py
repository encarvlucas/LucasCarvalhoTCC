# -*- coding: utf-8 -*-
#-------------------------------------- by: LUCAS CARVALHO DE SOUSA --------------------------------------------------
#Solução da equação de Laplace  dT - d²T = Q
#								dt   dx²
import numpy as np
import sys
import os
import math
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.animation import FuncAnimation

def diferencasfinitaspermanente(X, MRE, k, Q, T_0, T_L, dT_0, dT_L):
	#MATRIZ DE SOLUÇÃO
	L = max(X)
	numele = len(X)-1 #NÚMERO DE ESPAÇOS DO DOMÍNIO
	numpnts = numele + 1 #NÚMERO DE PONTOS PARA SEREM CALCULADOS OS VALORES DA FUNÇÃO
	#----------------------MATRIZ GERAL (-K+M)*a = f--------------------------------------------------------------------
	matriz = np.zeros((numpnts,numpnts))
	matrizf = np.copy(Q)
	for elem in MRE:
		matriz[elem[0]][elem[0]] = 2*k*1.0/(X[elem[0]]-X[elem[1]])**2		#\
		matriz[elem[1]][elem[0]] = -k*1.0/(X[elem[0]]-X[elem[1]])**2		# \  MATRIZ	    [ -2  1  0 ...]
		matriz[elem[0]][elem[1]] = -k*1.0/(X[elem[0]]-X[elem[1]])**2		# /        =  k [  1 -2  1 ...]
		# matriz[elem[1]][elem[1]] -= k*1.0/(X[elem[0]]-X[elem[1]])**2		#/				[  0  1 -2 ...]

	if T_0 != None:
		# T(0) = T_0 (CONDIÇÃO DE CONTORNO DE DIRICHLET)
		index = np.where(X==0)[0][0]
		matrizf[index] = T_0 						# VALOR NO PONTO
		matriz[index][index] = 1 							# ÚNICO COEFICIENTE DA LINHA E COLUNA
		a = range(numpnts)
		a.remove(index)
		for i in a:
			matrizf[i] -= matriz[i][index]*matrizf[index]	# SUBTRAÇÂO DA COMPONENTE DA LINHA
			matriz[index][i] = 0 							#{	ELIMINAÇÃO DA LINHA
			matriz[i][index] = 0 							#{	ELIMINAÇÃO DA COLUNA

	if T_L != None:
		# T(L) = T_L (CONDIÇÃO DE CONTORNO DE DIRICHLET)
		index = np.where(X==L)[0][0]
		matrizf[index] = T_L 						# VALOR NO PONTO
		matriz[index][index] = 1 							# ÚNICO COEFICIENTE DA LINHA E COLUNA
		a = range(numpnts)
		a.remove(index)
		for i in a:
			matrizf[i] -= matriz[i][index]*matrizf[index]	# SUBTRAÇÂO DA COMPONENTE DA LINHA
			matriz[index][i] = 0 							#{	ELIMINAÇÃO DA LINHA
			matriz[i][index] = 0 							#{	ELIMINAÇÃO DA COLUNA

	print matriz
	print matrizf
	resp = np.linalg.solve(matriz,matrizf)
	print "Solução:", resp
	return resp

def diferencasfinitasimplicito(X ,MRE ,k ,Q ,Delta_t ,T_inicial ,dT_0 ,dT_L):
	#MATRIZ DE SOLUÇÃO
	L = max(X)
	numele = len(X)-1 #NÚMERO DE ESPAÇOS DO DOMÍNIO
	numpnts = numele + 1 #NÚMERO DE PONTOS PARA SEREM CALCULADOS OS VALORES DA FUNÇÃO
	#----------------------MATRIZ GERAL (-K+M)*a = f--------------------------------------------------------------------
	matriz = np.zeros((numpnts,numpnts))
	matrizf = Q + T_inicial*1.0/Delta_t
	for elem in MRE:
		matriz[elem[0]][elem[0]] = 2*k*1.0/(X[elem[0]]-X[elem[1]])**2 + 1.0/Delta_t		#\  MATRIZ	  [2k/x+1/t  -k/x     0    ...]
		matriz[elem[1]][elem[0]] = -k*1.0/(X[elem[0]]-X[elem[1]])**2					# \        =  [-k/x    2k/x+1/t  -k/x  ...]
		matriz[elem[0]][elem[1]] = -k*1.0/(X[elem[0]]-X[elem[1]])**2					# /			  [  0     -k/x   2k/x+1/t ...]
		# matriz[elem[1]][elem[1]] -= k*1.0/(X[elem[0]]-X[elem[1]])**2					#/

	# T(0) = T_0 (CONDIÇÃO DE CONTORNO DE DIRICHLET)
	index = np.where(X==0)[0][0]
	matrizf[index] = T_inicial[0] 						# VALOR NO PONTO
	matriz[index][index] = 1 							# ÚNICO COEFICIENTE DA LINHA E COLUNA
	a = range(numpnts)
	a.remove(index)
	for i in a:
		matrizf[i] -= matriz[i][index]*matrizf[index]	# SUBTRAÇÂO DA COMPONENTE DA LINHA
		matriz[index][i] = 0 							#{	ELIMINAÇÃO DA LINHA
		matriz[i][index] = 0 							#{	ELIMINAÇÃO DA COLUNA

	# T(L) = T_L (CONDIÇÃO DE CONTORNO DE DIRICHLET)
	index = np.where(X==L)[0][0]
	matrizf[index] = T_inicial[-1] 						# VALOR NO PONTO
	matriz[index][index] = 1 							# ÚNICO COEFICIENTE DA LINHA E COLUNA
	a = range(numpnts)
	a.remove(index)
	for i in a:
		matrizf[i] -= matriz[i][index]*matrizf[index]	# SUBTRAÇÂO DA COMPONENTE DA LINHA
		matriz[index][i] = 0 							#{	ELIMINAÇÃO DA LINHA
		matriz[i][index] = 0 							#{	ELIMINAÇÃO DA COLUNA

	# print matriz
	# print matrizf
	resp = np.linalg.solve(matriz,matrizf)
	# print resp
	return resp

def diferencasfinitasexplicito(X, MRE, k, Q, Delta_t, T_inicial, dT_0, dT_L):
	L = max(X)
	numpnts = len(X) #NÚMERO DE PONTOS PARA SEREM CALCULADOS OS VALORES DA FUNÇÃO
	resp = np.zeros(numpnts)
	resp[0] = T_inicial[0]
	resp[-1] = T_inicial[-1]
	for i in range(1,numpnts-1):
		resp[i] = Delta_t*(Q[i] + (T_inicial[i-1]-2*T_inicial[i]+T_inicial[i+1])/((X[i]-X[i-1])*(X[i+1]-X[i])) + T_inicial[i]/Delta_t) 
	return resp

def grafico(x, y, color, title=None):
	#----------------------------------------GRÁFICO-----------------------------------------------------------
	plt.plot(np.array(x), np.array(y), color=color)
	maior = max(y)
	menor = min(y)
	axes = plt.gca()
	aux = 0.1*abs(maior-menor)
	axes.set_ylim([menor-3*aux, 3*aux+maior])
	plt.grid(True)
	if title:
		plt.title(title)
	# plt.savefig("img/respfinal.jpg")
	# plt.show()
	return

def main():
	args = sys.argv[1:]
	if not args:
		# DEFAULT CHOICE		
		trigger = "ea"
	else:
		trigger = args

	with open("malha.txt", "r") as arq:
		X = []
		for line in arq:
			X.append(float(line))
	with open("MRE.txt", "r") as arq:
		MRE = []
		for line in arq:
			aux = line.split(";")
			aux[0] = int(aux[0])
			aux[1] = int(aux[1])
			MRE.append(tuple(aux))
	X = np.array(X)
	MRE = np.array(MRE, int)
	Q = np.zeros(len(X))
	k = 1.0
	T_inicial = np.zeros(len(X))
	T_inicial[0] = 0.0
	T_inicial[-1] = 1.0
	print "X: ",X
	print "MRE: ",MRE
	print "Q: ",Q

	#Solução Permanente
	perm = diferencasfinitaspermanente(X, MRE, k, Q, T_inicial[0], T_inicial[-1], None, None)
	plt.plot(X, perm, "m", linewidth=4)

	#Solução Inicial
	plt.plot(X, T_inicial, "c--", linewidth=3)
	fig = plt.gcf()
	ax = plt.gca()
	ax.set_xlim([min(X),max(X)])
	aux = 0.1*abs(max(T_inicial)-min(T_inicial))
	ax.set_ylim([min(T_inicial)-1*aux, 1*aux+max(T_inicial)])
	plt.grid(True)

	def plotgif(y):
		startline.set_ydata(np.array(y)) # Use only one line
		# plt.plot(X, np.array(y), color="g") # Keep line history
		return

	#Solução Explicita
	if "e" in trigger:
		t_intervalos = np.zeros(100) + 0.005
		t_acumulado = 0.0
		T_atual = np.copy(T_inicial)
		frames = []

		for i in range(len(t_intervalos)):
			t_acumulado += t_intervalos[i]
			T_atual = diferencasfinitasexplicito(X, MRE, k, Q, t_intervalos[i], T_atual, None, None)
			if "a" in trigger:
				frames.append(T_atual)
			else:
				grafico(X, T_atual, (1-i*1.0/len(t_intervalos),0,0), "Explicit Method t={}".format(t_acumulado))
			# plt.savefig("img/explicito_{}.jpg".format(i)) #FRAME
		if "a" in trigger:
			aux = 0.1*abs(np.amax(frames))-np.amin(frames)
			ax.set_ylim([np.amin(frames)-1*aux, 1*aux+np.amax(frames)]) # Novo Máximo para o maior frame
			startline, = ax.plot(X, frames[0], color="g")
			anim = FuncAnimation(fig, plotgif, frames=frames, interval=50, repeat=False, save_count=0)
		ax.add_artist(AnchoredText("Time step: {0}s\nNumber of steps: {1}".format(np.mean(t_intervalos), len(t_intervalos)), loc=2))
		
		if "s" in trigger:
			if "a" in trigger:
				anim.save("img/explicito_t={}s.gif".format(t_acumulado), dpi=80, writer='imagemagick')
				print "Done!"
			else:
				plt.save("img/explicito_t={}s.jpg".format(t_acumulado))
		plt.show()
	
	#Solução Implicita
	if "i" in trigger:
		t_intervalos = np.zeros(20) + 0.025
		t_acumulado = 0.0
		T_atual = np.copy(T_inicial)
		frames = []

		for i in range(len(t_intervalos)):
			t_acumulado += t_intervalos[i]
			T_atual = diferencasfinitasimplicito(X, MRE, k, Q, t_intervalos[i], T_atual, None, None)
			if "a" in trigger:
				frames.append(T_atual)
			else:
				grafico(X, T_atual, (0,0,1-i*1.0/len(t_intervalos)), "Implicit  Method t={}".format(t_acumulado))
			# plt.savefig("img/explicito_{}.jpg".format(i)) #FRAME
		if "a" in trigger:
			aux = 0.1*abs(np.amax(frames))-np.amin(frames)
			ax.set_ylim([np.amin(frames)-1*aux, 1*aux+np.amax(frames)]) # Novo Máximo para o maior frame
			startline, = ax.plot(X, frames[0], color="g")
			anim = FuncAnimation(fig, plotgif, frames=frames, interval=150, repeat=False, save_count=0)
		ax.add_artist(AnchoredText("Time step: {0}s\nNumber of steps: {1}".format(np.mean(t_intervalos), len(t_intervalos)), loc=2))
		
		if "s" in trigger:
			if "a" in trigger:
				anim.save("img/implicito_t={}s.gif".format(t_acumulado), dpi=80, writer='imagemagick')
				print "Done!"
			else:
				plt.save("img/implicito_t={}s.jpg".format(t_acumulado))
		plt.show()
	return

if __name__ == '__main__':
	main()