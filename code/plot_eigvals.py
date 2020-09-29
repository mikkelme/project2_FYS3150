import numpy as np
import matplotlib.pyplot as plt


txt_file = "Eigvals.txt"


with open(txt_file, "r") as infile:
    eigvals = np.array(infile.readlines()[0].split()).astype(float)


number = np.linspace(1, len(eigvals), len(eigvals))


def analytic_eigvals(N):
    exact_eigvals = np.zeros(N)
    h = 1/N
    for j in range(1, N):
        exact_eigvals[j] = 2/h**2 - 2/h**2*np.cos(j*np.pi/N)
    return exact_eigvals

N = [10, 1e3, 1e4, 1e5]
for i in range(len(N)):
    n = int(N[i])
    exact_eigvals = analytic_eigvals(n)
    plt.subplot(2,2,i+1)
    plt.title(f"Eigenvalues, n = {n}")
    plt.plot(exact_eigvals)
    plt.xlabel("Eigenvalue number")
    plt.ylabel("Eigenvalue")
plt.show()




"""
plt.plot(number[:100] ,eigvals[:100], "o")
plt.title("Eigenvalues ")
plt.xlabel("Eigenvalue number")
plt.ylabel("Eigenvalue")

plt.show()
"""
