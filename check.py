#coding=utf-8
#ћетоды якоби, «ейдел€, верхней релаксации, сопр€женных градиентов
import matplotlib.pyplot as plt
import numpy as np
import jacobi
import seidel
import sor
import cg

def check(N,a,o):
	alfa=a
	A=np.zeros((N,N),float)
	b=np.zeros(N)
	x0=np.zeros(N)
	A[0,0]=2
	A[0,1]=-1+alfa
	A[N-1,N-1]=2
	A[N-1,N-2]=-1+alfa
	b[0]=1-alfa
	b[N-1]=1+alfa
	x0[0]=0.7
	x0[N-1]=0.8
	for i in range(1,N-1,1):
		A[i,i]=2
		A[i,i+1]=-1+alfa
		A[i,i-1]=-1+alfa
		b[i]=0
		x0[i]=0.1*i
	print(A)
	print(b)
	xj,kj,tj=jacobi.jacobi(A,b,x0,0.00001,1000)
	xjv,kjv,tjv=jacobi.jacobi_vec(A,b,x0,0.00001,1000)
	xs,ks,tjs=seidel.seidel(A,b,x0,0.00001,1000)
	xsv,ksv,tsv=seidel.seidel_vec(A,b,x0,0.00001,1000)
	xsor,ksor,tsor=sor.sor(A,b,x0,o,0.00001,1000)
	xsorv,ksorv,tsorv=sor.sor_vec(A,b,x0,o,0.00001,1000)
	xcg,kcg,tcg=cg.cg(A,b,0.00001,1000)
	return kj,kjv,ks,ksv,ksor,ksorv,kcg

print(check(4,0.3,1.1))