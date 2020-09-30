import numpy as np
import matplotlib.pyplot as plt

n = 200
rho = np.linspace(0,5,n)

def plotter(fileval, filevec, n=200):
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

	plt.plot(rho,np.abs(vv[:,0])**2)
	
plotter("Eigvals.txt","Eigvecs.txt")
#plotter("Eigvals_0.txt", "Eigvecs_0.txt")
#plotter("Eigvals_1.txt", "Eigvecs_1.txt")
#plotter("Eigvals_2_0.01w.txt", "Eigvecs_2_0.01w.txt")
#plotter("Eigvals_2_0.5w.txt", "Eigvecs_2_0.5w.txt")
#plotter("Eigvals_2_1w.txt", "Eigvecs_2_1w.txt")
#plotter("Eigvals_2_5w.txt", "Eigvecs_2_5w.txt")
#plt.legend(["Buckling beam","One electron"])
plt.ylabel(r"$|u(\rho)|^2$")
plt.xlabel(r"$\rho$")
plt.title(r"$|u(\rho)|")
plt.show()