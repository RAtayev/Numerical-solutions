﻿# Решение системы линейных уравнений с трехдиагональной матрицей
import numpy as np

A=np.array([[160,2,0,0],[6,185,5,0],[0,3,193,11],[0,0,8,134]])
b=np.array([10,22,42,72])


def solve_lu3(A,b):
	d=np.array([A[i,i] for i in range(len(A))],float)
	e=np.array([A[i,i+1] for i in range(len(A)-1)],float)
	c=np.array([A[i+1,i] for i in range(len(A)-1)],float)
	alfa=np.zeros(len(A),float)
	alfa[1]=-e[0]/d[0]
	for i in range(2,len(A)):
		alfa[i]=-e[i-1]/(d[i-1]+c[i-1]*alfa[i-1])
	betta=np.zeros(len(A),float)
	betta[1]=b[0]/d[0]
	for i in range(2,len(A)):
		betta[i]=(-c[i-1]*betta[i-1]+b[i-1])/(d[i-1]+c[i-1]*alfa[i-1])
	x=np.zeros(len(A),float)
	x[-1]=(-c[-1]*betta[-1]+b[-1])/(d[-1]+c[-1]*alfa[-1])
	for i in range(len(A)-2,-1,-1):
		x[i]=alfa[i+1]*x[i+1]+betta[i+1]
	return x
    
x = solve_lu3(A,b)
print(np.dot(A,x))