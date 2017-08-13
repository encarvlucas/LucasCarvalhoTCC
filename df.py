# -*- coding: utf-8 -*-
#-------------------------------------- by: LUCAS CARVALHO DE SOUSA --------------------------------------------------
#Solução da equação de Laplace  dT - d²T = Q
#								dt   dx²
import numpy as np
import sys
import os
import math
import matplotlib.pyplot as plt

def diferencasfinitas(X,MRE,k,Q,T_0,T_L,dT_0,dT_L):
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
	plt.show()
	return

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
	print "X: ",X
	print "MRE: ",MRE
	print "Q: ",Q
	diferencasfinitas(X,MRE,1.0,Q,0,1,None,None)
	return

if __name__ == '__main__':
	main()