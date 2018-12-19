import numpy as np
from math import pi, sin
import matplotlib.pyplot as plt

def RK(F, tau, T, y0):
	N = int(T/tau)+1
	t_mas = np.array([tau*n for n in range(0, N)])
	y = np.zeros((N, len(y0)), float)
	y[0] = y0
	k = 0
	while tau*(k+1) < T:
		k1 = F(t_mas[k], y[k])
		k11 = [tau*i/2 for i in k1]
		k2 = F(t_mas[k] + tau/2, y[k] + k11)
		k22 = [tau*i/2 for i in k2]
		k3 = F(t_mas[k] + tau/2, y[k] + k22)
		k33 = [tau*i for i in k3]
		k4 = F(t_mas[k] + tau, y[k] + k33)
		k2 = np.array([2*i for i in k2])
		k3 = np.array([2*i for i in k3])
		temp = np.array(k1 + k2 + k3 + k4)
		K = np.array([tau/6*i for i in temp])
		y[k+1] = y[k] + K
		k += 1
	return t_mas, y

def PK(F, tau, T, y0):
	N = int(T/tau)+1
	t_mas = np.array([tau*n for n in range(0, N)])
	y = np.zeros((N, len(y0)), float)
	y[0:4] = RK(F, tau, tau*4, y0)[1][0:4]
	print(y)
	k = 0
	while tau*(k+4) < T:
		k1 = F(t_mas[k], y[k])
		k2 = F(t_mas[k+1], y[k+1])
		k3 = F(t_mas[k+2], y[k+2])
		k4 = F(t_mas[k+3], y[k+3])
		k1 = np.array([-9*i for i in k1])
		k2 = np.array([37*i for i in k2])
		k3 = np.array([-59*i for i in k3])
		k4 = np.array([55*i for i in k4])
		temp = np.array(k1 + k2 + k3 + k4)
		K = np.array([tau/24*i for i in temp])
		y[k+4] = y[k+3] + K
		print(y[k+4])
		k1 = F(t_mas[k+1], y[k+1])
		k2 = F(t_mas[k+2], y[k+2])
		k3 = F(t_mas[k+3], y[k+3])
		k4 = F(t_mas[k+4], y[k+4])
		k2 = np.array([-5*i for i in k2])
		k3 = np.array([19*i for i in k3])
		k4 = np.array([9*i for i in k4])
		temp = np.array(k1 + k2 + k3 + k4)
		K = np.array([tau/24*i for i in temp])
		y[k+4] = y[k+3] + K
		print(y[k+4])
		k += 1
		
	return t_mas, y

def F(tn, yn):
	return np.array([yn[1], -sin(yn[0])])

def F1(tn, yn):
	a = 0.2
	b = 0.2
	c = 2.5
	return np.array([yn[1]-yn[2], yn[0] + a*yn[2], b + yn[2]*(yn[0] - c)])

def F2(tn, yn):
	a = 0.2
	b = 0.2
	c = 5
	return np.array([yn[1]-yn[2], yn[0] + a*yn[2], b + yn[2]*(yn[0] - c)])

t, y = PK(F1, 1, 100, [1, 1, 1])
u = np.array(y[:,0])
print(u)
plt.plot(t, u, label="PK")
plt.legend()
plt.ylabel("y[t]")
plt.xlabel("x")
plt.show()