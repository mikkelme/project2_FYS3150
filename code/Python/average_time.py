# Usage: Clean the Timeused.txt file and compute N simulations
# to do N simulations in bash: 
# $ for i in {1..N}; do ./Eigen_solver.exe N V rho_max Omega^2; done;
# Then you should have a .txt file with N times. Then run this program.

import numpy as np

n = 10
file = open("Timeused.txt")
t = []
for i in range(n):
	t.append(float(file.readline()))
t = np.array(t)
tavg = np.sum(t)/len(t)
dt = np.std(t)

print(tavg,dt)

