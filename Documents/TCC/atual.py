# -*- coding: utf-8 -*-
#-------------------------------------- by: LUCAS CARVALHO DE SOUSA --------------------------------------------------
import numpy as np
import sys
import os
import math
import matplotlib.pyplot as plt

def pprint(A):
	n = len(A)
	for i in range(0, n):
		line = ""
		for j in range(0, n):
			if A[i][j]>=0:
				line += " "
			line += ("{:.6}".format(str(A[i][j]))) + "\t"
			if (j == n-1):
				line += "| "
		print line
	print "\n"
	return

def openmalha(*args):
	if "default" in args:
		filenames = ("malhadefault.txt","MREdefault.txt")
	else:
		filenames = ("malha.txt","MRE.txt")
	with open(filenames[0], "r") as arq:
		X = []
		for line in arq:
			X.append(float(line))
	with open(filenames[1], "r") as arq:
		MRE = []
		for line in arq:
			aux = line.split(";")
			aux[0] = int(aux[0])
			aux[1] = int(aux[1])
			MRE.append(tuple(aux))
	X = np.array(X)
	MRE = np.array(MRE, int)
	return X, MRE

def geracaodemalha():
	#------------Vetor de posições X-----------------------------------------------------------------------------------
	L = float(raw_input("What is the size of the bar?\n"))
	n = int(raw_input("What is the number of elements?\n"))
	check = raw_input("What is the type of grid distribution:\n (L)inear, (Q)uadratic, (C)ubic, (E)xponential?\n").lower()
	if check not in ("l","1","q","2","c","3","e","4"):
			check = raw_input("Incorrect response, please use choose one of the following:\n (L)inear, (Q)uadratic, (C)ubic, (E)xponential?\n").lower()
	numpnts = n + 1
	X = np.zeros(numpnts)
	#linear
	if check in ("l","1"):
		for x in range(numpnts):
			X[x] = x*(L/(numpnts-1))
	#quadratica
	if check in ("q","2"):
		for x in range(numpnts):
			X[x] = ((x*1.0/(numpnts-1))**2)*L
	#cubica
	if check in ("c","3"):
		for x in range(numpnts):
			X[x] = ((x*1.0/(numpnts-1))**3)*L
	#exponecial
	if check in ("e","4"):
		for x in range(numpnts):
			X[x] = L*math.exp(-x)
		X = X[::-1]
	# print X
	#----------Saving grid to file-------------------------------------
	with open("malha.txt","w") as arq:
		for x in X:
			arq.write("{0}\n".format(x))

	#----------Matriz de relação dos elementos-------------------------------------------------------------
	#---------UNIDIMENSIONAL--------------------------------------------------------------------
	check = raw_input("Use default distribution? (Y/N)\n").lower()
	if check not in ("yes","y","ye","sim","si","s"):
		print "Unfortunately other distributions are not available yet\n Using default distribution"
	MRE = np.zeros([n,2], int)
	for i in range(n):
		MRE[i][0] = i
		MRE[i][1] = i+1
	#----------Saving MRE to file-------------------------------------
	with open("MRE.txt","w") as arq:
		for elem in MRE:
			arq.write("{0};{1}\n".format(elem[0],elem[1]))
	return X, MRE

def questionarCC(X):
	T_0 = None
	T_L = None
	dT_0 = None
	dT_L = None
	#First Boundary condition
	check = raw_input("What is the first type of Boundery Condition? (D/N)\n").lower()
	while check not in ("dirichlet","d","dir","neumann","neu","n"):
		check = raw_input("Inappropiate response, please try again. Answer using D or N -> (D/N)\n").lower()
	if check in ("dirichlet","d","dir"): #Dirichlet OK
		check = raw_input("And does it apply to the first edge (x = {})? (Y/N)\n".format(X[0])).lower()
		while check not in ("yes","y","ye","sim","si","s","no","n","não","nao"):
			check = raw_input("Inappropiate response, please try again. Answer using Y or N -> (Y/N)\n").lower()
		if check in (("yes","y","ye","sim","si","s")): #Dirichlet na primeira borda
			T_0 = float(raw_input("And what's its value?\n"))
		else:											#Dirichlet na segunda borda
			T_L = float(raw_input("And what's its value?\n"))
	else:                                #Neumann OK
		check = raw_input("And does it apply to the first edge (x = {})? (Y/N)\n".format(X[0])).lower()
		while check not in ("yes","y","ye","sim","si","s","no","n","não","nao"):
			check = raw_input("Inappropiate response, please try again. Answer using Y or N -> (Y/N)\n").lower()
		if check in (("yes","y","ye","sim","si","s")): #Neumann na primeira borda
			dT_0 = float(raw_input("And what's its value?\n"))
		else:											#Neumann na segunda borda
			dT_L = float(raw_input("And what's its value?\n"))

	#Second Bondary condition
	check = raw_input("What is the second type of Boundery Condition? (D/N)\n").lower()
	while check not in ("dirichlet","d","dir","neumann","neu","n"):
		check = raw_input("Inappropiate response, please try again. Answer using D or N -> (D/N)\n").lower()
	if check in ("dirichlet","d","dir"): #Dirichlet OK
		if (T_0 != None) or (T_L != None):
			if T_0 != None:
				T_L = float(raw_input("Since the first type of condition was set to Dirichlet for the first edge, please input the value for the second edge\n"))
			else:
				T_0 = float(raw_input("Since the first type of condition was set to Dirichlet for the second edge, please input the value for the first edge\n"))
		else:
			check = raw_input("And does it apply to the first edge (x = {})? (Y/N)\n".format(X[0])).lower()
			while check not in ("yes","y","ye","sim","si","s","no","n","não","nao"):
				check = raw_input("Inappropiate response, please try again. Answer using Y or N -> (Y/N)\n").lower()
			if check in (("yes","y","ye","sim","si","s")): #Dirichlet na primeira borda
				if T_0 != None:
					print
				T_0 = float(raw_input("And what's its value?\n"))
			else:											#Dirichlet na segunda borda
				T_L = float(raw_input("And what's its value?\n"))
	else:                                #Neumann OK
		if (dT_0 != None) or (dT_L != None):
			if dT_0 != None:
				dT_L = float(raw_input("Since the first type of condition was set to Neumann for the first edge, please input the value for the second edge\n"))
			else:
				dT_0 = float(raw_input("Since the first type of condition was set to Neumann for the second edge, please input the value for the first edge\n"))
		else:
			check = raw_input("And does it apply to the first edge (x = {})? (Y/N)\n".format(X[0])).lower()
			while check not in ("yes","y","ye","sim","si","s","no","n","não","nao"):
				check = raw_input("Inappropiate response, please try again. Answer using Y or N -> (Y/N)\n").lower()
			if check in (("yes","y","ye","sim","si","s")): #Neumann na primeira borda
				dT_0 = float(raw_input("And what's its value?\n"))
			else:											#Neumann na segunda borda
				dT_L = float(raw_input("And what's its value?\n"))
	return T_0, T_L, dT_0, dT_L

def elementosfinitoslin(X,MRE,k,T_0,T_L,dT_0,dT_L):
		#MATRIZ DE SOLUÇÃO
	L = max(X)
	numele = len(X)-1 #NÚMERO DE ESPAÇOS DO DOMÍNIO
	numpnts = numele + 1 #NÚMERO DE PONTOS PARA SEREM CALCULADOS OS VALORES DA FUNÇÃO
	#----------------------MATRIZ GERAL (-K+M)*a = f--------------------------------------------------------------------
	matriz = np.zeros((numpnts,numpnts))
	matrizf = np.zeros(numpnts)
	for elem in MRE:
		matriz[elem[0]][elem[0]] -= k*1.0/abs(X[elem[0]]-X[elem[1]])		#\
		matriz[elem[1]][elem[0]] -= -k*1.0/abs(X[elem[0]]-X[elem[1]])		# \  MATRIZ	  k [ 1 -1]
		matriz[elem[0]][elem[1]] -= -k*1.0/abs(X[elem[0]]-X[elem[1]])		# /    K   =  h [-1  1]
		matriz[elem[1]][elem[1]] -= k*1.0/abs(X[elem[0]]-X[elem[1]])		#/				
	print "Before Boundary Conditions:"
	print matriz
	
	#----------------------MATRIZ COM CONDIÇÕES DE CONTORNO (-K+M)*a = f - CC--------------------------------
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

	if dT_0 != None:
		# dT(0)/dx = dT_0 (CONDIÇÃO DE CONTORNO DE NEUMANN)
		index = np.where(X==0.0)[0][0]
		matrizf[index] -= dT_0	# SUBTRAÇÂO DA COMPONENTE DA LINHA

	if dT_L != None:
		# dT(L)/dx = dT_L (CONDIÇÃO DE CONTORNO DE NEUMANN)
		index = np.where(X==L)[0][0]
		matrizf[index] -= dT_L	# SUBTRAÇÂO DA COMPONENTE DA LINHA

	if (numpnts <= 9):
		print "After Boundary Conditions:"
		pprint(matriz)
	resp = np.linalg.solve(matriz,matrizf)

	for i in range(numpnts):
		if (numpnts <= 50):
			print "T({0:.2})= {1}".format(X[i],resp[i])

	#----------------------------------------GRÁFICO-----------------------------------------------------------
	if graficoshow:
		plt.plot(np.array(X), np.array(resp), "r")
		maior = max(resp)
		menor = min(resp)
		axes = plt.gca()
		axes.set_xlim([0,L])
		aux = 0.1*abs(maior-menor)
		axes.set_ylim([menor-3*aux,3*aux+maior])
		plt.grid(True)
	return

def elementosfinitosquad(L,n,k,T_0,T_L,dT_0,dT_L):
	#MATRIZ DE SOLUÇÃO
	numele = n #NÚMERO DE ESPAÇOS DO DOMÍNIO
	numpnts = 2*numele + 1 #NÚMERO DE PONTOS PARA SEREM CALCULADOS OS VALORES DA FUNÇÃO
	h = float(L)/(2*numele)
	
	#----------------------MATRIZ GERAL (-K+M)*a = f---------------------------------------------------
	matriz = [0]*(numpnts)
	for i in range(numpnts):
		matriz[i] = np.zeros(numpnts+1)
		for j in range(numpnts+1): #+1 DEVIDO AO PONTO INICIAL/FINAL DO DOMÍNIO QUE DEVE SER INCLUÍDO
			# matriz[i].append(0)
			if (i == j):
				if (i%2 == 0):
					matriz[i][j] = 2*(-7.0/(6.0*h)+h*4.0/15.0)   #i==j = 0; x2
					if ((i == 0) or (i == numpnts-1)):
						matriz[i][j] /= 2
				else:
					matriz[i][j] = -16.0/(6.0*h)+h*16.0/15.0		#i==j = 1; x1
			if ((j == i+1) or (j == i-1)):
				matriz[i][j] = 8.0/(6.0*h)+h*2.0/15.0 			#|i-j| = 1
			if (((j == i+2) or (j == i-2)) and (j%2 == 0)): # and (j != 0) and (i != 0)):
				matriz[i][j] = -1.0/(6.0*h)-h*1.0/15.0				#|i-j| = 2
			if (j == numpnts):
				if (i%2 == 0):
					matriz[i][j] = 2*(-h*1.0/3)						#i = 0 ou i = 2
					if ((i == 0) or (i == numpnts-1)):
						matriz[i][j] /= 2
				else:
					matriz[i][j] = -h*4.0/3						#i = 1

	#----------------------MATRIZ COM CONDIÇÕES DE CONTORNO (-K+M)*a = f - CC--------------------------------
	if T_0 != None:
		# T(0) = 0 (CONDIÇÃO DE CONTORNO DE DIRICHLET)
		matriz[0][numpnts] = T_0
		matriz[0][0] = 1 										#{
		matriz[0][1] = 0 										#{
		matriz[0][2] = 0 										#{ 
		matriz[1][numpnts] -= matriz[1][0]*matriz[0][numpnts]	#{	ELIMINAÇÃO DA LINHA/COLUNA
		matriz[1][0] = 0 										#{  
		matriz[2][numpnts] -= matriz[2][0]*matriz[0][numpnts]	#{
		matriz[2][0] = 0 										#{

	if T_L != None:
		# T(L) = 0 (CONDIÇÃO DE CONTORNO DE DIRICHLET)
		matriz[numpnts-1][numpnts] = T_L
		matriz[numpnts-1][numpnts-1] = 1						#{
		matriz[numpnts-1][numpnts-2] = 0						#{
		matriz[numpnts-1][numpnts-3] = 0						#{ ELIMINAÇÃO DA LINHA/COLUNA
		matriz[numpnts-2][numpnts] -= matriz[numpnts-2][numpnts-1]*matriz[numpnts-1][numpnts]	#{
		matriz[numpnts-2][numpnts-1] = 0						#{  
		matriz[numpnts-3][numpnts] -= matriz[numpnts-3][numpnts-1]*matriz[numpnts-1][numpnts]	#{
		matriz[numpnts-3][numpnts-1] = 0 

	if dT_0 != None:
		# dT(0)/dx = 1 (CONDIÇÃO DE CONTORNO DE NEUMANN)
		matriz[0][numpnts] -= dT_0

	if dT_L != None:
		# dT(L)/dx = 1 (CONDIÇÃO DE CONTORNO DE NEUMANN)
		matriz[numpnts-1][numpnts] -= dT_L

	if (numpnts <= 10):
		pprint(matriz)
	
	resp = gauss(matriz)

	x = [0]*numpnts
	for i in range(numpnts):
		if (numpnts <= 50):
			print "T({0:.2})= {1}".format(i*L/(numele*2),resp[i])
		x[i] = (i*L/(numele*2))

	#----------------------------------------GRÁFICO-----------------------------------------------------------
	if graficoshow:
		plt.plot(np.array(x), np.array(resp), "r")
		maior = max(resp)
		menor = min(resp)
		axes = plt.gca()
		axes.set_xlim([0,L])
		aux = 0.1*abs(maior-menor)
		axes.set_ylim([menor-3*aux,3*aux+maior])
		plt.grid(True)
	return

def diferencasfinitas(X,MRE,k,T_0,T_L,dT_0,dT_L):
	#MATRIZ DE SOLUÇÃO
	L = max(X)
	numele = len(X)-1 #NÚMERO DE ESPAÇOS DO DOMÍNIO
	numpnts = numele + 1 #NÚMERO DE PONTOS PARA SEREM CALCULADOS OS VALORES DA FUNÇÃO
	#----------------------MATRIZ GERAL (-K+M)*a = f--------------------------------------------------------------------
	matriz = np.zeros((numpnts,numpnts))
	matrizf = np.zeros(numpnts)
	for elem in MRE:
		matriz[elem[0]][elem[0]] -= k*1.0/abs(X[elem[0]]-X[elem[1]])		#\
		matriz[elem[1]][elem[0]] -= -k*1.0/abs(X[elem[0]]-X[elem[1]])		# \  MATRIZ	    [ -2  1  0 ...]
		matriz[elem[0]][elem[1]] -= -k*1.0/abs(X[elem[0]]-X[elem[1]])		# /        =  k [  1 -2  1 ...]
		matriz[elem[1]][elem[1]] -= k*1.0/abs(X[elem[0]]-X[elem[1]])		#/				[  0  1 -2 ...]
	pprint(matriz)
	return

def main():
	args = sys.argv[1:]
	if not args:
		print("usage: show(show graphic)\n By: Lucas Carvalho de Sousa")
	print("Hello! This is a script to solve a diferential equation using the Finite Elements Method, the equation is:\n d²T/dx²=0.")

	#---------------------------------DEBUGGING - VERIFICAÇÂO DE VARIAVEIS---------------------------------------------
	if any("print" in s.lower() for s in args):
		printshow = True
	else:
		printshow = False

	#--------------------------------------------USO DE GRÁFICO?--------------------------------------------------------
	global graficoshow
	if any("show" in s.lower() for s in args):
		graficoshow = True
	else:
		graficoshow = False

	#-------------------------------------------OBTENÇÃO DOS DADOS------------------------------------------------------
	if raw_input("Use Default Values? (Y/N)\n").lower() in ("yes","y","ye","sim","si","s"):
		# L = float (1)		## O TAMANHO DA BARRA
		# n = 3				## O NÚMERO DE DIVISÕES
		X, MRE = openmalha("default") ## OBTÉM OS VALORES DE L E n DIRETO DA MALHA
		L = max(X)			## ASSUMI-SE QUE O REFERENCIAL TEM A ORIGEM EM X=0 e POSIÇÂO MÁXIMA X=L
		n = len(X)-1		## APROXIMAÇÂO LINEAR
		T_0 = 0.0		## O VALOR DA CONDIÇÃO DE CONTORNO DE DIRICHLET NA POSIÇÃO x=0
		T_L = 1.0		## O VALOR DA CONDIÇÃO DE CONTORNO DE DIRICHLET NA POSIÇÃO x=L
		dT_0 = None	## O VALOR DA CONDIÇÃO DE CONTORNO DE NEUMANN NA POSIÇÃO x=0
		dT_L = None	## O VALOR DA CONDIÇÃO DE CONTORNO DE NEUMANN NA POSIÇÃO x=L
		TT_0 = 20.0	## A CONDIÇÃO INICIAL NA BARRA TODA
		k =  1.0		##COEFICIENTE DE CONDUTIVIDADE TÉRMICA
		tempos = (0,20,50,100,500,2000)
	else:
		if raw_input("Use stored grid? (Y/N)\n").lower() in ("yes","y","ye","sim","si","s"):
			X, MRE = openmalha()
		else:
			X, MRE = geracaodemalha()
		k = raw_input("What is the coefficient of thermal conductivity of the bar? (k)\n")
		if k != "skip":
			k = float(k)
			T_0, T_L, dT_0, dT_L = questionarCC(X)
		else:
			k = 1.0
			T_0 = 0.0
			T_L = 1.0
			dT_0 = None
			dT_L = None
	#------------------------------------------------DEBUGGING-------------------------------------------------------
	if printshow:
		print "L = {0}\nn = {1}\nX = ".format(L,n),X
		print "MRE = ",MRE
		if T_0 != None:
			print "T_0 = {}".format(T_0)
		if T_L != None:
			print "T_L = {}".format(T_L)
		if dT_0 != None:
			print "dT_0 = {}".format(dT_0)
		if dT_L != None:
			print "dT_L = {}".format(dT_L)
	#--------------------------------------------FUNÇÕES DE SOLUÇÃO-----------------------------------------------------
	# elementosfinitosquad(L,n,k,T_0=T_0,dT_L=dT_L) #SOLUÇÂO MEF QUADRÀTICA -----OBSOLETA------
	elementosfinitoslin(X,MRE,k,T_0,T_L,dT_0,dT_L)  #SOLUÇÂO MEF LINEAR
	if graficoshow:
		plt.show()
	return

if __name__ == '__main__':
	main()