# -*- coding: utf-8 -*-
#-------------------------------------- by: LUCAS CARVALHO DE SOUSA --------------------------------------------------
import numpy as np
import sys
import os
import math
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.animation import FuncAnimation

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

def geracaodemalha():
	#------------Vetor de posições X-----------------------------------------------------------------------------------
	Lx = 1.0
	Ly = 1.0
	nx = 3
	ny = 3
	
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
		for i in range(len(xx)):
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

	drawelement(0)
	# print X, Y, MRE
	# plt.show()
	return X, Y, MRE

def main():
	geracaodemalha()
	xy,IEN = openmalha()
	for elem in IEN:
		x = xy[elem]
		print x
		A = ((xy[elem[0]][0]*xy[elem[1]][1] - xy[elem[1]][0]*xy[elem[0]][1]) + 
			(xy[elem[1]][0]*xy[elem[2]][1] - xy[elem[2]][0]*xy[elem[1]][1]) +
			(xy[elem[2]][0]*xy[elem[0]][1] - xy[elem[0]][0]*xy[elem[2]][1]))*2
	print A

if __name__ == '__main__':
	main()