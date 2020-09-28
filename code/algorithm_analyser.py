import subprocess
import numpy as np
from run_cpp import *

#Variables
exe_file = "./Jacobi_solver.exe"
txt_files = ["iter_dim.txt", "MeanError.txt", "Timeused.txt"]
N_max = 100
potential = 0

#Run exe_file
run_cpp(exe_file, txt_files, N_max, potential)
N = np.linspace(2, N_max, N_max - 1)
#Interval for linear regression
start = 20
end = N_max

#Define title and axis labels
title = [   "Relation between matrix dimension N and Iterations I", \
            "Relation between matrix dimension N\nand mean error $\epsilon$ between numerical and analytical eigenvalues", \
            "Relation between matrix dimension N and time used T"   ]
xlabel = "N"
ylabel = [ "I", "$\epsilon$", "T[s]"]

#Plot results
for i in range(len(txt_files)):
    data = read_data(txt_files[i], N)
    plot(N, data, start, end, title[i], xlabel, ylabel[i])
