#Нахожение максимального по модулю собственного значения
import numpy as np
import pm


def myfunc(N):
	A=np.zeros((N,N),float)
	q0=np.zeros(N,float)
	q0[0]=1
	for i in range(0,N,1):
		for j in range(0,N,1):
			A[i][j]=1/(i+j+1)
	lmbda_max,t=pm.pm(A,q0,0.001,1000)
	print(chr(955)," = ",lmbda_max,"\nExecution time: ",t,"\n\n")

for k in range(2,11,1):
	myfunc(k)