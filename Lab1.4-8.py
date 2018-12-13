#coding=utf-8
#Методы Якоби, Зейделя, верхней релаксации, сопряженных градиентов
import matplotlib.pyplot as plt
import numpy as np
import jacobi
import seidel
import sor
import cg

def check(N,a,o):
	alfa=a
	A=np.zeros((N,N),float)
	AA=np.zeros((N,N),float)
	b=np.zeros(N)
	bb=np.zeros(N,float)
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
	xj,kj,tj=jacobi.jacobi(A,b,x0,0.00001,1000)
	xjv,kjv,tjv=jacobi.jacobi_vec(A,b,x0,0.00001,1000)
	xs,ks,tjs=seidel.seidel(A,b,x0,0.00001,1000)
	xsv,ksv,tsv=seidel.seidel_vec(A,b,x0,0.00001,1000)
	xsor,ksor,tsor=sor.sor(A,b,x0,o,0.00001,1000)
	xsorv,ksorv,tsorv=sor.sor_vec(A,b,x0,o,0.00001,1000)
	xcg,kcg,tcg=cg.cg(A,b,0.00001,1000)
	return kj,kjv,ks,ksv,ksor,ksorv,kcg

T=5
a=0.3
N=[]
kj=[]
kjv=[]
ks=[]
ksv=[]
kcg=[]
for i in range(2,T,1):
	N.append(i+1)
	kj.append(check(i+1,a,1.3)[0])
	kjv.append(check(i+1,a,1.3)[1])
	ks.append(check(i+1,a,1.3)[2])
	ksv.append(check(i+1,a,1.3)[3])
	kcg.append(check(i+1,a,1.3)[6])

al=np.arange(0.1,1,0.01)
kja=[]
kjva=[]
ksa=[]
ksva=[]
for i in range(0,len(al),1):
	kja.append(check(10,al[i],1.3)[0])
	kjva.append(check(10,al[i],1.3)[1])
	ksa.append(check(10,al[i],1.3)[2])
	ksva.append(check(10,al[i],1.3)[3])

o=np.arange(1.1,2,0.1)
NO=5
ksorn=np.zeros((NO,len(o)))
ksorvn=np.zeros((NO,len(o)))
for k in range(3,NO+3,1):
	for i in range(0,len(o),1):
		ksorn[k-3,i]=check(k,a,o[i])[4]
		ksorvn[k-3,i]=check(k,a,o[i])[5]

plt.plot(N,kj,label="jacobi")
plt.plot(N,ks,label="seidel")
plt.plot(N,kcg,label="cg")
plt.legend()
plt.ylabel("Iterations")
plt.xlabel("N")
plt.show()

plt.plot(N,kjv,label="jacobi_vec")
plt.plot(N,ksv,label="seidel_vec")
plt.legend()
plt.ylabel("Iterations")
plt.xlabel("N")
plt.show()

plt.plot(al,kja,label="jacobi")
plt.plot(al,ksa,label="seidel")
plt.legend()
plt.ylabel("Iterations")
plt.xlabel("$\\alpha$")
plt.show()

plt.plot(al,kjva,label="jacobi_vec")
plt.plot(al,ksva,label="seidel_vec")
plt.legend()
plt.ylabel("Iterations")
plt.xlabel("$\\alpha$")
plt.show()

plt.plot(o,ksorn[0],label="sor_N=3")
plt.plot(o,ksorn[1],label="sor_N=4")
plt.plot(o,ksorn[2],label="sor_N=5")
plt.plot(o,ksorn[3],label="sor_N=6")
plt.plot(o,ksorn[4],label="sor_N=7")
plt.legend()
plt.ylabel("Iterations")
plt.xlabel("$\\omega$")
plt.show()

plt.plot(o,ksorvn[0],label="sorvec_N=3")
plt.plot(o,ksorvn[1],label="sorvec_N=4")
plt.plot(o,ksorvn[2],label="sorvec_N=5")
plt.plot(o,ksorvn[3],label="sorvec_N=6")
plt.plot(o,ksorvn[4],label="sorvec_N=7")
plt.legend()
plt.ylabel("Iterations")
plt.xlabel("$\\omega$")
plt.show()