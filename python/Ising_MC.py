import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def initial_config(size_lat):
    return 2*(np.random.randint(0,2,size=(size_lat, size_lat))-1/2)

def hamiltonian_ising(arr_spin, coupling, ex_field):
    n0, n1 = np.shape(arr_spin)
    hamiltonian = 0.0
    for i in range(-1,n0-1):
        for j in range(-1,n1-1):
            hamiltonian += -coupling*(arr_spin[i][j]*arr_spin[i][j+1] + arr_spin[i][j]*arr_spin[i+1][j])
    hamiltonian += -ex_field*np.sum(arr_spin)
    return hamiltonian

def update_config(arr_spin, temp, coupling, ex_field):
    n0, n1 = np.shape(arr_spin)
    energy = hamiltonian_ising(arr_spin, coupling, ex_field)
    for i in range(0,n0):
        for j in range(0,n1):
            diff_energy = diff_flip(arr_spin, n0, n1, i, j, coupling, ex_field)
            if diff_energy < 0:
                energy += diff_energy
            elif np.random.rand() < np.exp(-diff_energy/temp):
                energy += diff_energy
            else:
                arr_spin[i][j] = -arr_spin[i][j]
    return arr_spin

def diff_flip(arr_spin, n0, n1, i, j, coupling, ex_field):
    arr_spin[i][j] = -arr_spin[i][j]
    fliped = 0.0
    fliped += -coupling*(arr_spin[i][j]*arr_spin[i-1][j] + arr_spin[i][j]*arr_spin[-n0+i+1][j]
                         + arr_spin[i][j]*arr_spin[i][j-1] + arr_spin[i][j]*arr_spin[i][-n1+j+1]) \
              -ex_field*arr_spin[i][j]
    return 2*fliped

def magnetization_plot(size_lat, temp_min, temp_max, coupling, ex_field):
    num_update = 400
    num_update_sub = 400
    sampling_interval = 10
    num_step = 100
    mag_lst =[]
    mag_samples = []
    temp_lst = np.linspace(temp_min, temp_max, num_step)

    arr_spin=initial_config(size_lat)

    for temp_tmp in temp_lst:
        mag_samples = []

        for i in range(0, num_update):
            update_config(arr_spin, temp_tmp, coupling, ex_field)

        for i in range(0, num_update_sub):
            update_config(arr_spin, temp_tmp, coupling, ex_field)
            if (i+1)%sampling_interval==0:
                mag_samples.append(np.sum(arr_spin)/len(arr_spin)**2)

        mag_lst.append(np.sum(mag_samples)/len(mag_samples))

    plt.plot(temp_lst, mag_lst)
    plt.ylabel("magnetization per cite")
    plt.xlabel("T")
    plt.show()
