import numpy as np
import time
from math import sqrt

A=np.array([[17,3,10],[3,17,-2],[10,-2,12]],float)
x=np.array([1,3,4])
x0=np.array([0.7,3.2,4.1],float)
b=np.dot(A,x)

def jacobi(A,b,x0):
	start_time=time.time()
	x=np.zeros((2,len(A)),float)
	x[1,:]=x0
	temp1=1
	while sqrt(temp1)>0.000000001:
		x[0,:]=x[1,:]
		for i in range(0,len(A),1):
			temp=0
			for j in range(0,i,1):
				temp=temp+A[i,j]*x[0,j]
			for t in range(i+1,len(A),1):
				temp=temp+A[i,t]*x[0,t]
			x[1,i]=(b[i]-temp)/A[i,i]
		temp1=0
		for l in range(0,len(A),1):
			temp1=temp1+(x[1,l]-x[0,l])*(x[1,l]-x[0,l])	
	print(x)
	print("---%s seconds ---" % (time.time() - start_time))
	print("\n\n")

jacobi(A,b,x0)

def jacobi_vec(A,b,x0):
	start_time=time.time()
	L=np.zeros((len(A),len(A)),float)
	U=np.zeros((len(A),len(A)),float)
	D=np.zeros((len(A),len(A)),float)
	for i in range(0,len(A),1):
		for k in range(0,i,1):
			L[i,k]=A[i,k]
		D[i,i]=A[i,i]
		for t in range(i+1,len(A),1):
			U[i,t]=A[i,t]
	x=np.zeros((2,len(A)),float)
	x[1,:]=x0
	temp1=1
	while sqrt(temp1)>0.000000001:
		x[0,:]=x[1,:]
		x[1]=np.dot(np.linalg.inv(D),(np.dot(-(L+U),x[0])+b))
		temp1=0
		for l in range(0,len(A),1):
			temp1=temp1+(x[1,l]-x[0,l])*(x[1,l]-x[0,l])
	print(x)
	print("---%s seconds ---" % (time.time() - start_time))
	print("\n\n")

jacobi_vec(A,b,x0)