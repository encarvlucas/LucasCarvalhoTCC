# -*- coding: utf-8 -*-
#-------------------------------------- by: LUCAS CARVALHO DE SOUSA --------------------------------------------------
import numpy as np
import sys
import os
import math
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.animation import FuncAnimation
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

def openmalha(*args):
	if "default" in args:
		filenames = ("malha2ddefault.txt","MRE2ddefault.txt")
	else:
		filenames = ("malha2d.txt","MRE2d.txt")
	with open(filenames[0], "r") as arq:
		XY = []
		for line in arq:
			XY.append(tuple(line.split(";")))
	with open(filenames[1], "r") as arq:
		MRE = []
		for line in arq:
			MRE.append(tuple(line.split(";")))
	return np.array(XY, float), np.array(MRE, int)

def geracaodemalha(Lx,nx,Ly,ny):
	#------------Vetor de posições X-----------------------------------------------------------------------------------
	numpnts_y = ny + 1
	numpnts_x = nx + 1
	X = np.zeros(numpnts_x)
	Y = np.zeros(numpnts_y)

	#----------------Tipos de malha---------------------------------------------------------------------------------------
	def dimlinear(VAL,Lv):
		#linear
		numpnts = len(VAL)
		for val in range(numpnts):
			VAL[val] = val*(Lv/(numpnts-1))
		return 
	def dimquadratica(VAL,Lv):
		#quadratica
		numpnts = len(VAL)
		for val in range(numpnts):
			VAL[val] = ((val*1.0/(numpnts-1))**2)*Lv
		return 
	def dimcubica(VAL,Lv):
		#cubica
		numpnts = len(VAL)
		for val in range(numpnts):
			VAL[val] = ((val*1.0/(numpnts-1))**3)*Lv
		return 
	def dimexponencial(VAL,Lv):
		#exponencial
		numpnts = len(VAL)
		for val in range(numpnts):
			VAL[val] = Lv*math.exp(-val)
		# VAL = VAL[::-1] # BUGGED
		return 

	#---------------------Gerando o IEN---------------------------------------------------------------------------------
	dimlinear(X,Lx)
	dimlinear(Y,Ly)
	xx,yy = np.meshgrid(X,Y)
	# plt.plot(xx, yy, marker=".", color="k", linestyle="none",ms=10)
	xx = np.array(np.reshape(xx,(numpnts_x*numpnts_y,1)))
	yy = np.array(np.reshape(yy,(numpnts_x*numpnts_y,1)))

	MRE = np.zeros((2*((numpnts_x-1)*(numpnts_y-1)),3),dtype=int)
	count = 0
	for i in range(len(MRE)/2):
		MRE[count] = (i,i+numpnts_x+1,i+numpnts_x)
		MRE[count] += i/(numpnts_x-1)
		count += 1
		MRE[count] = (i,i+1,i+numpnts_x+1)
		MRE[count] += i/(numpnts_x-1)
		count += 1

	#----------Saving grid to file------------------------------------------------------------------------------
	with open("malha2d.txt","w") as arq:
		for i in range(len(xx)):
			arq.write("{0};{1}\n".format(xx[i][0],yy[i][0]))

	with open("MRE2d.txt","w") as arq:
		for i in range(len(MRE)):
			arq.write("{0};{1};{2}\n".format(MRE[i][0],MRE[i][1],MRE[i][2]))

	#------------Checagem---------------------------------------------------------------------------------------
	def drawelement(ele_num):
		plt.plot(xx, yy, marker=".", color="k", linestyle="none",ms=10)
		tipodeelem = 3 #triangular
		showelex = np.zeros((tipodeelem+1,2))
		showeley = np.zeros((tipodeelem+1,2))
		for i in range(tipodeelem):
			showelex[i] = xx[MRE[ele_num][i]]
			showeley[i] = yy[MRE[ele_num][i]]
		showelex[tipodeelem] = xx[MRE[ele_num][0]]
		showeley[tipodeelem] = yy[MRE[ele_num][0]]
		plt.plot(showelex,showeley, "r")
		plt.show()
		return

	# drawelement(0)
	# print X, Y, MRE
	# plt.show()
	return X, Y, MRE

def main():
	Lx = 1.0
	Ly = 1.0
	nx = 3
	ny = 3
	geracaodemalha(Lx,nx,Ly,ny)

	k_condx = 1.0 #\ Condutividade térmica
	k_condy = 1.0 #/
	t = 1.0 # Espessura teórica da placa
	xy,IEN = openmalha() # vetor de coordenadas dos pontos, matriz de organização dos elementos
	Q = np.zeros(len(xy)) # Geração de calor

	K = np.zeros((len(xy),len(xy)))
	M = np.copy(K)

	for elem in IEN:
		x = xy[elem,0]
		y = xy[elem,1]
		A = ((x[0]*y[1] - x[1]*y[0]) + 
			 (x[1]*y[2] - x[2]*y[1]) +
			 (x[2]*y[0] - x[0]*y[2]))/2
		b = np.array([y[1]-y[2],
					  y[2]-y[0],
					  y[0]-y[1]])
		c = np.array([x[2]-x[1],
					  x[0]-x[2],
					  x[1]-x[0]])
		k = (t/(4*A))*(k_condx*np.array([[b[0]**2, b[0]*b[1], b[0]*b[2]], 
										 [b[0]*b[1], b[1]**2, b[1]*b[2]],
										 [b[0]*b[2], b[1]*b[2], b[2]**2]])
					 + k_condy*np.array([[c[0]**2, c[0]*c[1], c[0]*c[2]], 
										 [c[0]*c[1], c[1]**2, c[1]*c[2]],
										 [c[0]*c[2], c[1]*c[2], c[2]**2]])) 
		# print k
		# m = (A/12)*np.array([[2,1,1],[1,2,1],[1,1,2]])
		for i in [0,1,2]:
			for j in [0,1,2]:
				K[elem[i]][elem[j]] += k[i][j]
				# M[elem[i]][elem[j]] += m[i][j]

	#---------------------------------Condição de contorno-----------------------------------------------
	# print np.where(xy[:,0]==0)[0]
	for j in np.where(xy[:,0]==0)[0]:
		# print np.where(K[:,j]!=0)[0]
		for i in np.where(K[:,j]!=0)[0]:
			Q[i] -= K[i][j]*(xy[j][1])
			#					/\ - valor da função no ponto
			K[i][j] = 0
		K[j][j] = 1
		Q[j] = xy[j][1]
	for j in np.where(xy[:,0]==max(xy[:,0]))[0]:
		for i in np.where(K[:,j]!=0)[0]:
			Q[i] -= K[i][j]*((xy[j][1])**2+1)
			#					/\ - valor da função no ponto
			K[i][j] = 0
		K[j][j] = 1
		Q[j] = xy[j][1]**2+1
	for j in np.where(xy[:,1]==0)[0]:
		for i in np.where(K[:,j]!=0)[0]:
			Q[i] -= K[i][j]*(xy[j][0])
			#					/\ - valor da função no ponto
			K[i][j] = 0
		K[j][j] = 1
		Q[j] = xy[j][0]
	for j in np.where(xy[:,1]==max(xy[:,1]))[0]:
		for i in np.where(K[:,j]!=0)[0]:
			Q[i] -= K[i][j]*((xy[j][0])**2+1)
			#					/\ - valor da função no ponto
			K[i][j] = 0
		K[j][j] = 1
		Q[j] = xy[j][0]**2+1
	#---------------------------------Solução-------------------------------------------------------------
	T = np.linalg.solve(K,Q)
	# print xy[:,0], "\n", xy[:,1],"\n", T
	print T[5],T[6],T[9],T[10]
	
	# plt.pcolormesh(np.reshape(xy[:,0],(nx+1,ny+1)), np.reshape(xy[:,1],(nx+1,ny+1)),
	# 					 np.reshape(T,(nx+1,ny+1)), cmap='RdBu', vmin=min(T), vmax=max(T))
	# axes = plt.gca()
	# axes.set_xlim([min(xy[:,0])-0.2*abs(np.median(xy[:,0])),max(xy[:,0])+0.2*abs(np.median(xy[:,0]))])
	# axes.set_ylim([min(xy[:,1])-0.2*abs(np.median(xy[:,1])),max(xy[:,1])+0.2*abs(np.median(xy[:,1]))])
	# plt.colorbar()
	# plt.show()

if __name__ == '__main__':
	main()