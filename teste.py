# -*- coding: utf-8 -*-
#-------------------------------------- by: LUCAS CARVALHO DE SOUSA --------------------------------------------------
import numpy as np
import sys
import os
import math
import matplotlib.pyplot as plt

# def pprint(A):
# 	n = len(A)
# 	for i in range(0, n):
# 		line = ""
# 		for j in range(0, n):
# 			if A[i][j]>=0:
# 				line += " "
# 			line += ("{:.6}".format(str(A[i][j]))) + "\t"
# 			if (j == n-1):
# 				line += "| "
# 		print line
# 	print "\n"
# 	return

tam = 1.0
Lx = 10.0
Ly = 0.5
# n = 5
nx = 100
ny = 120
# blank = np.zeros(n)
x_val = np.zeros(nx)
y_val = np.zeros(ny)
# vecX = np.zeros(n+1)
# M = [np.zeros(n+1) for i in range(n+1)]

# newfile = open("teste.txt","w")
# for i in range(n):
# 	newfile.write("{0};{1}\n".format(i,i+1))
# 	vecX[i+1] = vecX[i]+tam/n  
# newfile.close()
# vecX[1] = 0.25 
# print vecX

# with open("teste.txt","r") as oldfile:
# 	data = []
# 	for line in oldfile:
# 		aux = line.split(";")
# 		aux[0] = int(aux[0])
# 		aux[1] = int(aux[1])
# 		data.append(aux)
# print data

# index = None
# for elem in data:
# 	try:
# 		index = [data.index(elem),elem.index(np.ndarray.tolist(vecX).index(tam))]
# 		index = data[index[0]][index[1]]
# 	except:
# 		pass
# print index

# K = [np.zeros(n+1) for i in range(n+1)]
# for elem in data:
# 	M[elem[0]][elem[0]] = abs(vecX[elem[0]]-vecX[elem[1]])
# 	M[elem[1]][elem[0]] = -abs(vecX[elem[0]]-vecX[elem[1]])
# 	M[elem[0]][elem[1]] = -abs(vecX[elem[0]]-vecX[elem[1]])
# 	M[elem[1]][elem[1]] = abs(vecX[elem[0]]-vecX[elem[1]])
# 	# for i in elem:
# 	# 	for j in elem:
# 	# 		M[i][j] -= K[i][j]#abs(vecX[elem[0]]-vecX[elem[1]])
# pprint(M)

# A = ""
# while A not in ("a","b","c"):
# 	A = raw_input("Bla bla bla? (A/B/C)\n").lower() 

# #linear
for x in range(nx):
	# blank[x] = x*(tam/(n-1))
	x_val[x] = x*(Lx/(nx-1))
# y1 = np.array(blank)
# x_val = np.array(blank)
for y in range(ny):
	# blank[x] = x*(tam/(n-1))
	y_val[y] = y*(Ly/(ny-1))

# #quadratica
# for x in range(n):
# 	blank[x] = ((x*1.0/(n-1))**2)*tam
# y2 = np.array(blank)
# y_val = np.array(blank)

# #cubica
# for x in range(n):
# 	blank[x] = ((x*1.0/(n-1))**3)*tam
# y3 = np.array(blank)

# #exponecial
# for x in range(n):
# 	blank[x] = tam*math.exp(-x)
# y4 = np.array(blank[::-1])
# y_val = np.array(blank)

# x1 = np.linspace(0,tam,n)
# plt.plot(x1,y1,"r")
# plt.plot(x1,y2,"g")
# plt.plot(x1,y4,"b")
# plt.plot(x1,y3,"m")
# plt.show()

# x_val = np.array(blank)
# y_val = np.copy(x_val)
xx,yy = np.meshgrid(x_val,y_val)
# print "x = \n",xx,"\n","y = \n",yy

#Ex 3.2:
A = 0.07
phi = 2.0*np.pi/4.0
lamb = 24.0/6.0
k_ex = 2.0*np.pi/lamb
for i in range(nx):
	yy[ny-1][i] += A*np.sin(2.0*np.pi*xx[ny-1][i]/lamb-phi)
print "x = \n",xx,"\n","y = \n",yy

plt.plot(xx, yy, marker=".", color="k", linestyle="none")
axes = plt.gca()
axes.set_xlim(-Lx*0.1,Lx*1.1)
axes.set_ylim(-Ly*0.1,Ly*1.1)
# plt.show()

# MRE 
x = np.array(np.reshape(xx,(1,nx*ny))[0])
y = np.array(np.reshape(yy,(1,nx*ny))[0])

# Elementos quadrados
# MRE = np.zeros(((nx-1)*(ny-1),4),dtype=int)
# for i in range(len(MRE)):
# 	MRE[i] = (i,i+1,i+nx+1,i+nx)
# 	MRE[i] += i/(nx-1)

# Elementos triangulares
MRE = np.zeros((2*((nx-1)*(ny-1)),3),dtype=int)
count = 0
for i in range(len(MRE)/2):
	MRE[count] = (i,i+nx+1,i+nx)
	MRE[count] += i/(nx-1)
	count += 1
	MRE[count] = (i,i+1,i+nx+1)
	MRE[count] += i/(nx-1)
	count += 1
print len(MRE),"MRE = \n",MRE

#show element ele_num
# ele_num = 3
# tipodeelem = 3 #triangular
# showelex = np.zeros((tipodeelem+1,2))
# showeley = np.zeros((tipodeelem+1,2))
# for i in range(tipodeelem):
# 	showelex[i] = x[MRE[ele_num][i]]
# 	showeley[i] = y[MRE[ele_num][i]]
# showelex[tipodeelem] = x[MRE[ele_num][0]]
# showeley[tipodeelem] = y[MRE[ele_num][0]]
# plt.plot(showelex,showeley, "r")
plt.show()