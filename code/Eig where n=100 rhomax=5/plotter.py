import numpy as np
import matplotlib.pyplot as plt

n = 100
rho = np.linspace(0,5,n)

def plotter(fileval, filevec, n=100):
	l = []
	file = open(fileval)
	l.append((file.readline().split()))
	l = l[0]
	for i in range(n):
		l[i] = float(l[i])

	v = []
	file = open(filevec)
	v.append((file.readline().split()))
	v = v[0]
	for i in range(n):
		v[i] = float(v[i])

	vv = np.zeros((n,n))
	for i in range(n):
		vv[:,i] = v[n*i:n+n*i]

	plt.plot(rho,np.abs(vv[:,0])**2, label = filevec)

plotter("Eigvals_0.txt", "Eigvecs_0.txt")
plotter("Eigvals_1.txt", "Eigvecs_1.txt")
plotter("Eigvals_2_0.01w.txt", "Eigvecs_2_0.01w.txt")
plotter("Eigvals_2_0.5w.txt", "Eigvecs_2_0.5w.txt")
plotter("Eigvals_2_1w.txt", "Eigvecs_2_1w.txt")
plotter("Eigvals_2_5w.txt", "Eigvecs_2_5w.txt")
plt.legend()
plt.ylabel(r"$|u(\rho)|^2$")
plt.xlabel(r"$\rho$")
plt.title(r"$|u(\rho)|^2$ where $\rho_{max}=5$")
plt.show()