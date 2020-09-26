import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

txt_file = "iter_dim.txt"

#Clear file before run
with open(txt_file, "w") as infile:
    infile.write("")

N_max = 100

#Run program for different dimensions
N = np.linspace(2, N_max, N_max - 1)
print(f"Running program for n: [2:{N_max:d}]\nn: ", end = "")
for n in N: #Run Jacobi_solver for different dimensions n
    subprocess.call(['./Jacobi_solver.exe', str(n)])
    print(f"\r n: {int(n)}/{N_max}", end = "")
print("\nDone")

#Read iterations result
iter = np.zeros(len(N))
with open(txt_file, "r") as infile:
    iter = np.array(infile.readlines()).astype(float)

#Linear regression for logaritmic plot
start = 20
end = N_max

start += -2; end += -2
log_N = np.log(N); log_iter = np.log(iter)
a, b, R, p, std_a = stats.linregress(log_N[start:end], log_iter[start:end])
x = np.linspace(log_N[start], log_N[end], int(1e4))
y = x*a + b

#Plot
fig, axes = plt.subplots(ncols=1, nrows=2, figsize=(7,6))
plt.subplots_adjust(left = 0.09, bottom = 0.09, right = 0.95, top = 0.88, hspace = 0.40)
plt.suptitle("Relation between matrix dimensions N and Iterations I")
plt.subplot(2,1,1)

plt.title("Linear plot")
plt.plot(N, iter, "o")
plt.plot(N[start], iter[start], "o", color = "tab:green", label = "start")
plt.plot(N[end], iter[end], "o", color = "tab:red", label = "end")
plt.xlabel("N")
plt.ylabel("I")
plt.legend()

plt.subplot(2,1,2)
plt.title("Logaritmic plot")
plt.plot(log_N, log_iter, "o")
plt.plot(log_N[start], log_iter[start], "o", color = "tab:green", label = "start")
plt.plot(log_N[end], log_iter[end], "o", color = "tab:red", label = "end")
plt.plot(x,y, label = f"Linear regression (ax + b) \n a={a:.3f}, b={b:0.3f}, $\delta a$={std_a:0.3f}")
plt.xlabel("Log(N)")
plt.ylabel("Log(I)")
plt.legend()
plt.show()
