# -*- coding: utf-8 -*-
#-------------------------------------- by: LUCAS CARVALHO DE SOUSA --------------------------------------------------
import numpy as np
import sys
import os
import math
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

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

def elemfinperm(Lx,Ly,nx,ny,k_condx,k_condy,esp,xy,IEN):
	numpnts_x = nx + 1
	numpnts_y = ny + 1
	#--------------------------Geração das matrizes---------------------------------------
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
		k = -(esp/(4*A))*(k_condx*np.array([[b[0]**2, b[0]*b[1], b[0]*b[2]], 
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
			K[j][i] = 0
		K[j][j] = 1
		Q[j] = xy[j][1]
	for j in np.where(xy[:,0]==max(xy[:,0]))[0]:
		for i in np.where(K[:,j]!=0)[0]:
			Q[i] -= K[i][j]*((xy[j][1])**2+1)
			#					/\ - valor da função no ponto
			K[i][j] = 0
			K[j][i] = 0
		K[j][j] = 1
		Q[j] = xy[j][1]**2+1
	for j in np.where(xy[:,1]==0)[0]:
		for i in np.where(K[:,j]!=0)[0]:
			Q[i] -= K[i][j]*(xy[j][0])
			#					/\ - valor da função no ponto
			K[i][j] = 0
			K[j][i] = 0
		K[j][j] = 1
		Q[j] = xy[j][0]
	for j in np.where(xy[:,1]==max(xy[:,1]))[0]:
		for i in np.where(K[:,j]!=0)[0]:
			Q[i] -= K[i][j]*((xy[j][0])**2+1)
			#					/\ - valor da função no ponto
			K[i][j] = 0
			K[j][i] = 0
		K[j][j] = 1
		Q[j] = xy[j][0]**2+1

	#---------------------------------Solução-------------------------------------------------------------
	T = np.linalg.solve(K,Q)
	# print max(T)
	# print xy[:,0], "\n", xy[:,1],"\n", T
	# print np.reshape(T,(numpnts_x,numpnts_y))

	#---------------------------------Plotting-------------------------------------------------------------
	fig = plt.figure()
	axes = plt.gca()
	def plotmesh2D():
		plt.pcolormesh(np.reshape(xy[:,0],(numpnts_x,numpnts_y)), np.reshape(xy[:,1],(numpnts_x,numpnts_y)),
							 np.reshape(T,(numpnts_x,numpnts_y)), cmap='jet', vmin=min(T), vmax=max(T))
		plt.colorbar()
		return plt.show()
	
	def plot3Dsurf(rotate,*save):
		axes = Axes3D(fig)
		surf = axes.plot_trisurf(xy[:,0],xy[:,1],T, cmap="jet")
		axes.view_init(90,270)
		fig.colorbar(surf, shrink=0.4, aspect=9)
		if rotate:
			def update(angle):
				axes.view_init(50, angle)
				return
			anim = FuncAnimation(fig, update, frames=range(0,360), interval=1, save_count=False)
			if save:
				anim.save("trisurf.gif", dpi=80, writer='imagemagick')
		return plt.show()

	axes.set_xlim([min(xy[:,0])-0.2*abs(np.median(xy[:,0])),max(xy[:,0])+0.2*abs(np.median(xy[:,0]))])
	axes.set_ylim([min(xy[:,1])-0.2*abs(np.median(xy[:,1])),max(xy[:,1])+0.2*abs(np.median(xy[:,1]))])
	# plotmesh2D()
	# plot3Dsurf(False)
	return T

def elemfintrans(Lx,Ly,nx,ny,k_condx,k_condy,esp,xy,IEN,dt,nt,T_0,cc_id):
	numpnts_x = nx + 1
	numpnts_y = ny + 1
	#--------------------------------Geração das matrizes---------------------------------------
	Q = np.zeros(len(xy)) # Geração de calor
	K = np.zeros((len(xy),len(xy)))
	M = np.copy(K)
	T = np.copy(T_0)

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
		k = -(esp/(4.0*A))*(k_condx*np.array([[b[0]**2, b[0]*b[1], b[0]*b[2]], 
										 [b[0]*b[1], b[1]**2, b[1]*b[2]],
										 [b[0]*b[2], b[1]*b[2], b[2]**2]])
					 + k_condy*np.array([[c[0]**2, c[0]*c[1], c[0]*c[2]], 
										 [c[0]*c[1], c[1]**2, c[1]*c[2]],
										 [c[0]*c[2], c[1]*c[2], c[2]**2]])) 
		m = (A/12.0)*np.array([[2,1,1],[1,2,1],[1,1,2]])
		for i in [0,1,2]:
			for j in [0,1,2]:
				K[elem[i]][elem[j]] += k[i][j]
				M[elem[i]][elem[j]] += m[i][j]

	A = M - K*dt

	#---------------------------------Condição de contorno-----------------------------------------------
	def cc1(vec):
		for i in cc_id:
			vec[i] = T_0[i]
		return vec	

	def cc2():
		for i in cc_id:
			A[i,:] = 0
			A[i][i] = 1
			# vec[i] = T_0[i]
		return

	#-----------------------------------Solução no tempo-------------------------------------------
	#PROCURAR método dos gradientes conjugados
	#Gmesh
	#Delauney
	cc2()
	frames = [T_0]
	for t in range(nt):
		B = dt*np.dot(M,Q) + np.dot(M, T)
		B = cc1(B)
		T = np.linalg.solve(A,B)
		frames.append(T)

	#---------------------------------Plotting-------------------------------------------------------------
	fig = plt.figure()
	axes = plt.gca()
	def plotmesh2D():
		plt.pcolormesh(np.reshape(xy[:,0],(numpnts_x,numpnts_y)), np.reshape(xy[:,1],(numpnts_x,numpnts_y)),
							 np.reshape(T,(numpnts_x,numpnts_y)), cmap='jet', vmin=min(T), vmax=max(T))
		plt.colorbar()
		return plt.show()
	
	def plot3Dsurf(gif,*save):
		axes = Axes3D(fig)
		axes.set_zlim3d([np.min(frames),np.max(frames)])
		# axes.view_init(90,270)
		if gif:
			surf = axes.plot_trisurf(xy[:,0],xy[:,1],T_0, cmap="jet", vmin=np.min(frames), vmax=np.max(frames))
			fig.colorbar(surf, shrink=0.4, aspect=9)
			def update(T_atual):
				plt.cla()
				surf = axes.plot_trisurf(xy[:,0],xy[:,1],T_atual, cmap="jet", vmin=np.min(frames), vmax=np.max(frames))
				axes.set_zlim3d([np.min(frames),np.max(frames)])
				return
			anim = FuncAnimation(fig, update, frames=frames[1:], interval=100, save_count=False)
			if save:
				anim.save("trisurf.gif", dpi=80, writer='imagemagick')
		else:
			surf = axes.plot_trisurf(xy[:,0],xy[:,1],T, cmap="jet", vmin=np.min(frames), vmax=np.max(frames))
			fig.colorbar(surf, shrink=0.4, aspect=9)
		return plt.show()

	axes.set_xlim([min(xy[:,0])-0.2*abs(np.median(xy[:,0])),max(xy[:,0])+0.2*abs(np.median(xy[:,0]))])
	axes.set_ylim([min(xy[:,1])-0.2*abs(np.median(xy[:,1])),max(xy[:,1])+0.2*abs(np.median(xy[:,1]))])
	# plotmesh2D()
	# plot3Dsurf(True)
	return frames

def output(VECX,VECY,VECIEN,VECT,ext="VTK",dt=0):
	n = len(VECT)
	num_IEN = len(VECIEN)
	data_name = "Temperature"

	if (ext=="CSV"):
		#-------------------------Saving results to CSV file----------------------------------------------------
		with open("results/resultado2d.csv","w") as arq:
			arq.write("{0}, Points:0, Points:1, Points:2\n".format(data_name))
			for i in range(n):
				arq.write("{0},{1},{2},{3}\n".format(VECT[i], VECX[i], VECY[i], 0))

	if (ext=="VTK"):
		try:
			n = len(VECT[0])
			#---------Saving multiple results to VTK files----------------------------------------------------
			for j in range(len(VECT)):
				with open("results/resultado2d_{}.vtk".format(j),"w") as arq:
					#------------------------------------Header---------------------------------------------
					arq.write("# vtk DataFile Version 3.0\n{0}\n{1}\n\nDATASET {2}\n".format("Cube example","ASCII","POLYDATA"))
					arq.write("FIELD FieldData 1\nTIME 1 1 double\n{}\n".format(dt))
					#------------------------------------Points coordinates----------------------------------------
					arq.write("\nPOINTS {0} {1}\n".format(n, "float"))
					for i in range(n):
						arq.write("{0} {1} 0.0\n".format(VECX[i], VECY[i]))
					#---------------------------------------Cells---------------------------------------------------
					arq.write("\nPOLYGONS {0} {1}\n".format(num_IEN,num_IEN*4))
					for i in range(num_IEN):
						arq.write("{0} {1} {2} {3}\n".format(3, VECIEN[i][0], VECIEN[i][1], VECIEN[i][2]))
					#---------------------------------------Data in each point------------------------------------
					arq.write("\nPOINT_DATA {0}\n\nSCALARS {1} float 1\n".format(n, data_name))
					arq.write("\nLOOKUP_TABLE {0}\n".format(data_name))
					for i in range(n):
						arq.write("{}\n".format(VECT[j][i]))
		except TypeError:
			#-------------------------Saving results to VTK file----------------------------------------------------
			with open("results/resultado2d.vtk","w") as arq:
				#------------------------------------Header---------------------------------------------
				arq.write("# vtk DataFile Version 3.0\n{0}\n{1}\n\nDATASET {2}\n".format("Cube example","ASCII","POLYDATA"))
				#------------------------------------Points coordinates----------------------------------------
				arq.write("\nPOINTS {0} {1}\n".format(n, "float"))
				for i in range(n):
					arq.write("{0} {1} 0.0\n".format(VECX[i], VECY[i]))
				#---------------------------------------Cells---------------------------------------------------
				arq.write("\nPOLYGONS {0} {1}\n".format(num_IEN,num_IEN*4))
				for i in range(num_IEN):
					arq.write("{0} {1} {2} {3}\n".format(3, VECIEN[i][0], VECIEN[i][1], VECIEN[i][2]))
				#---------------------------------------Data in each point------------------------------------
				arq.write("\nPOINT_DATA {0}\n\nSCALARS {1} float 1\n".format(n, data_name))
				arq.write("\nLOOKUP_TABLE {0}\n".format(data_name))
				for i in range(n):
					arq.write("{}\n".format(VECT[i]))
	return

def main():
	# Cadastro dos parâmetros:
	Lx = 1.0
	Ly = 1.0
	nx = 3 # número de divisões no eixo x (número de pontos - 1)
	ny = 3 # número de divisões no eixo y (número de pontos - 1)
	# geracaodemalha(Lx,nx,Ly,ny)

	k_condx = 1.0 #\ Condutividade térmica
	k_condy = 1.0 #/
	esp = 1.0 # Espessura teórica da placa
	xy,IEN = openmalha() # vetor de coordenadas dos pontos, matriz de organização dos elementos

	dt = 0.005
	nt = 100
	T_0 = np.zeros((nx+1)*(ny+1)) # Condição inicial de temperatura
	pontos_cc = []
	for i in range(len(T_0)):
		if (xy[i,0]==0):
			T_0[i] = xy[i,1]
			pontos_cc.append(i)
		if (xy[i,1]==0):
			T_0[i] = xy[i,0]
			pontos_cc.append(i)
		if (xy[i,0]==Lx):
			T_0[i] = xy[i,1]**2+1
			pontos_cc.append(i)
		if (xy[i,1]==Ly):
			T_0[i] = xy[i,0]**2+1
			pontos_cc.append(i)
	# print np.rot90(np.reshape(T_0,(nx+1,ny+1)))

	resp = elemfinperm(Lx, Ly, nx, ny, k_condx, k_condy, esp, xy, IEN)
	print np.rot90(np.reshape(resp,(nx+1,ny+1)))
	output(xy[:,0],xy[:,1],IEN,resp)

	resp = elemfintrans(Lx, Ly, nx, ny, k_condx, k_condy, esp, xy, IEN, dt, nt, T_0, np.unique(pontos_cc))
	print np.rot90(np.reshape(resp[-1],(nx+1,ny+1)))
	output(xy[:,0],xy[:,1],IEN,resp,dt=dt)

	return

if __name__ == '__main__':
	main()