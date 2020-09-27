import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def run_cpp(N_max, txt_files):
    """ Run program for different dimensions """
    for txt_file in txt_files: #Clear files before run
        with open(txt_file, "w") as infile:
            infile.write("")

    #Run program for different dimensions
    N = np.linspace(2, N_max, N_max - 1)
    print(f"Running program for n: [2:{N_max:d}]\nn: ", end = "")
    for n in N: #Run Jacobi_solver for different dimensions n
        subprocess.call(['./Jacobi_solver.exe', str(n)])
        print(f"\r n: {int(n)}/{N_max}", end = "")
    print("\nDone")
    return N

def read_data(txt_file, N):
    """ Read data from different txt-files """
    data = np.zeros(len(N))
    with open(txt_file, "r") as infile:
        data = np.array(infile.readlines()).astype(float)
    return data

def plot(N, data, start, end, title, xlabel, ylabel):
    """ Plot linear and logaritmic plot
        with linear regression on the logaritmic plot """
    #Linear regression for logaritmic plot
    start += -2; end += -2
    log_N = np.log(N); log_data = np.log(data)
    a, b, R, p, std_a = stats.linregress(log_N[start:end], log_data[start:end])
    x = np.linspace(log_N[start], log_N[end], int(1e4))
    y = x*a + b

    #Plot
    fig, axes = plt.subplots(ncols=1, nrows=2, figsize=(7,6))
    plt.subplots_adjust(left = 0.12, bottom = 0.09, right = 0.97, top = 0.86, hspace = 0.40)
    plt.suptitle(title)
    plt.subplot(2,1,1)

    plt.title("Linear plot")
    plt.plot(N, data, "o")
    plt.plot(N[start], data[start], "o", color = "tab:green", label = "start")
    plt.plot(N[end], data[end], "o", color = "tab:red", label = "end")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()

    plt.subplot(2,1,2)
    plt.title("Logaritmic plot")
    plt.plot(log_N, log_data, "o")
    plt.plot(log_N[start], log_data[start], "o", color = "tab:green", label = "start")
    plt.plot(log_N[end], log_data[end], "o", color = "tab:red", label = "end")
    plt.plot(x,y, label = f"Linear regression (y = ax + b) \n a={a:.3f}, b={b:0.3f}, $\delta a$={std_a:0.3f}")
    plt.xlabel("Log(" + xlabel + ")")
    plt.ylabel("Log(" + ylabel + ")")
    plt.legend()
    plt.show()
