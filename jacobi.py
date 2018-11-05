import numpy as np
import time
from math import sqrt

def jacobi(A,b,x0,acc,it_max):
	t=0
	for i in range(0,len(A),1):
		if 2*abs(A[i][i])>np.sum(abs(A[i])):
			t+=1
		else:
			break
	if t==len(A):
		start_time=time.time()
		x=np.zeros((2,len(A)),float)
		x[1,:]=x0
		temp1=1
		k=0
		while sqrt(temp1)>=acc and k<it_max:
			x[0,:]=x[1,:]
			for i in range(0,len(A),1):
				temp=0
				for j in range(0,i,1):
					temp+=A[i,j]*x[0,j]
				for t in range(i+1,len(A),1):
					temp+=A[i,t]*x[0,t]
				x[1,i]=(b[i]-temp)/A[i,i]
			temp1=0
			for l in range(0,len(A),1):
				temp1+=(x[1,l]-x[0,l])*(x[1,l]-x[0,l])
			k+=1
			exec_time="---%s seconds ---" % (time.time() - start_time)
		if k==it_max:
			print("---Maximum number of iterations is reached!---")
			print("---Solution may not be accurate!---\n\n")
			return x[1],k,exec_time
		else:
			return x[1],k,exec_time
	else:
		print("---Matrix is not correct!---")
		return x0,0,0

def jacobi_vec(A,b,x0,acc,it_max):
	t=0
	for i in range(0,len(A),1):
		if 2*abs(A[i][i])>np.sum(abs(A[i])):
			t+=1
		else:
			break
	if t==len(A):
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
		k=0
		while sqrt(temp1)>acc and k<it_max:
			x[0,:]=x[1,:]
			x[1]=np.dot(np.linalg.inv(D),(np.dot(-(L+U),x[0])+b))
			temp1=0
			for l in range(0,len(A),1):
				temp1+=(x[1,l]-x[0,l])*(x[1,l]-x[0,l])
			k+=1
		exec_time="---%s seconds ---" % (time.time() - start_time)
		if k==it_max:
			print("---Maximum number of iterations is reached!---")
			print("---Solution may not be accurate!---\n\n")
			return x[1],k,exec_time
		else:
			return x[1],k,exec_time
	else:
		print("---Matrix is not correct!---")
		return x0,0,0