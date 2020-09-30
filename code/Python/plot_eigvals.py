import numpy as np
import matplotlib.pyplot as plt
import subprocess

exe_file = "./Eigen_solver.exe"
txt_eigval = "Eigvals.txt"
txt_eigvec = "Eigvecs.txt"



def analytic_eigvals(n):
    exact_eigvals = np.zeros(n)
    h = 1/n
    for j in range(1, n):
        exact_eigvals[j] = 2/h**2 - 2/h**2*np.cos(j*np.pi/n)
    return exact_eigvals

def numerical_eigvals(n, potential):
    exact_eigvals = np.zeros(n)
    print(f"Running: \"{exe_file} {n} {potential}\"")
    subprocess.call([exe_file, str(n), str(potential)])
    with open(txt_eigval, "r") as infile:
        eigvals = np.array(infile.readlines()[0].split()).astype(float)
    with open(txt_eigvec, "r") as infile:
        eigvecs = np.zeros((n,n))
        for i in range(n):
            line = infile.readline()
            eigvecs[i] = np.array(line.split()).astype(float)
    return eigvals, eigvecs.T



N = [10, 50, 100, 200]
fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(7,6))
plt.subplots_adjust(left = 0.12, bottom = 0.09, right = 0.97, top = 0.86, hspace = 0.40, wspace = 0.38)
plt.suptitle("Eigenvalues for Buckling beam problem\n(Found with Jacobi method)")
for i in range(len(N)):
    n = int(N[i])
    eigvals, eigvecs = numerical_eigvals(n,0)
    plt.subplot(2,2,i+1)
    plt.title(f"Eigenvalues, N = {n}")
    plt.plot(eigvals)
    plt.xlabel("Eigenvalue number")
    plt.ylabel("Eigenvalue")
plt.show()


n = 100
fig, axes = plt.subplots(ncols=3, nrows=4, figsize=(10,8))
plt.subplots_adjust(left = 0.08, bottom = 0.07, right = 0.99, top = 0.89, hspace = 0.75, wspace = 0.38)
plt.suptitle(f"Eigenvectors for Buckling beam problem, N = {n}\n(Found with Jacobi method)")
eigvals, eigvecs = numerical_eigvals(n,0)
J = [   0,  1,  2,  \
        9,  39, 51, \
        63, 79, 89, \
        97, 98, 99  ]
for i in range(len(J)):
    j = int(J[i])
    plt.subplot(4,3,i+1)
    plt.title(f"#Eigenvector = {j+1}")
    plt.plot(eigvecs[:,j])
    plt.xlabel("vector element")
    plt.ylabel("vector value")

plt.show()
