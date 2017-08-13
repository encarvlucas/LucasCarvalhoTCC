# -*- coding: utf-8 -*-
#-------------------------------------- by: LUCAS CARVALHO DE SOUSA --------------------------------------------------
#Solução da equação de Laplace  dT - d²T = Q
#								dt   dx²
import numpy as np
import sys
import os
import math
import matplotlib.pyplot as plt

def diferencasfinitasimplicito(X,MRE,k,Q,T_inicial,T_0,T_L,dT_0,dT_L):
	#MATRIZ DE SOLUÇÃO
	L = max(X)
	numele = len(X)-1 #NÚMERO DE ESPAÇOS DO DOMÍNIO
	numpnts = numele + 1 #NÚMERO DE PONTOS PARA SEREM CALCULADOS OS VALORES DA FUNÇÃO
	#----------------------MATRIZ GERAL (-K+M)*a = f--------------------------------------------------------------------
	matriz = np.zeros((numpnts,numpnts))
	matrizf = Q
	for elem in MRE:
		matriz[elem[0]][elem[0]] -= k*1.0/(X[elem[0]]-X[elem[1]])**2		#\
		matriz[elem[1]][elem[0]] -= -k*1.0/(X[elem[0]]-X[elem[1]])**2		# \  MATRIZ	    [ -2  1  0 ...]
		matriz[elem[0]][elem[1]] -= -k*1.0/(X[elem[0]]-X[elem[1]])**2		# /        =  k [  1 -2  1 ...]
		matriz[elem[1]][elem[1]] -= k*1.0/(X[elem[0]]-X[elem[1]])**2		#/				[  0  1 -2 ...]

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
	print resp

	#----------------------------------------GRÁFICO-----------------------------------------------------------
	plt.plot(np.array(X), np.array(resp), "r")
	maior = max(resp)
	menor = min(resp)
	axes = plt.gca()
	axes.set_xlim([0,L])
	aux = 0.1*abs(maior-menor)
	axes.set_ylim([menor-3*aux,3*aux+maior])
	plt.grid(True)
	plt.title("Metodo Implicito")
	plt.savefig("img/respfinal.jpg")
	plt.show()
	return

def diferencasfinitasexplicito(X,MRE,k,Q,T_inicial,Delta_t,T_0,T_L,dT_0,dT_L):
	L = max(X)
	numpnts = len(X) #NÚMERO DE PONTOS PARA SEREM CALCULADOS OS VALORES DA FUNÇÃO
	resp = np.zeros(numpnts)
	resp[0] = T_0
	resp[numpnts-1] = T_L
	for i in range(1,numpnts-1):
		resp[i] = Delta_t*(Q[i] + (T_inicial[i-1]-2*T_inicial[i]+T_inicial[i+1])/((X[i]-X[i-1])*(X[i+1]-X[i])) + T_inicial[i]/Delta_t) 
	
	#----------------------------------------GRÁFICO-----------------------------------------------------------
	plt.plot(np.array(X), np.array(resp), "r")
	maior = max(resp)
	menor = min(resp)
	axes = plt.gca()
	axes.set_xlim([0,L])
	aux = 0.1*abs(maior-menor)
	axes.set_ylim([menor-3*aux,3*aux+maior])
	plt.grid(True)
	return resp

def main():
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
	T_inicial = np.zeros(len(X))
	t_intervalos = np.zeros(150) + 0.001
	print "X: ",X
	print "MRE: ",MRE
	print "Q: ",Q

	t_acumulado = 0.0
	for i in range(len(t_intervalos)):
		t_acumulado += t_intervalos[i]
		T_inicial = diferencasfinitasexplicito(X,MRE,1.0,Q,T_inicial,t_intervalos[i],0,1,None,None)
		plt.title("Metodo Explicito t={}".format(t_acumulado))
		# plt.savefig("img/explicito_{}.jpg".format(i))
		# plt.show()

	plt.savefig("img/explicito_t={}s.jpg".format(t_acumulado))
	return

if __name__ == '__main__':
	main()