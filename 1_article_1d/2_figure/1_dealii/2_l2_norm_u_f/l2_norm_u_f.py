
import sys
sys.path.append('../')
import settings

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(5, 4))
ax = fig.add_axes([0,0,1,1])

n_row = 5
n_col = 4

data_l2_norm_u = np.zeros((n_row,n_col))
data_l2_norm_f = np.zeros((n_row,n_col))

vec_case=['2pic_power_minus_2_sin2picx', 'exp', '2pic_power_minus_2_sin2picx_minus_x_square_over_2', '2pic_power_minus_1_sin2picx', 'linear']

for id_case in range(n_col):
    f = open("l2_norm_u_"+vec_case[id_case]+".txt", 'r')
    data_raw = [line.split() for line in f]
    f.close()

    for id_coeff in range(n_row):
        data_l2_norm_u[id_coeff,id_case] = data_raw[id_coeff+1][1]
        data_l2_norm_f[id_coeff,id_case] = data_raw[id_coeff+1][3]
        
    plt.loglog(data_l2_norm_u[:,id_case], data_l2_norm_f[:,id_case], "o-", mfc='none',label="Case "+str(id_case+1))

plt.xticks(fontsize=settings.fontsize_tick)
plt.xlabel(r"$\||u\||_2$",fontsize=settings.fontsize_label)
plt.xlim(1e-6,1e4)

plt.yticks(fontsize=settings.fontsize_tick)
plt.ylabel(r"$\||f\||_2$",fontsize=settings.fontsize_label)
plt.ylim(1e-5,1e5)

plt.legend(fontsize=settings.fontsize_legend, loc=1,ncol=1)

plt.savefig('l2_norm_u_f.pdf', bbox_inches='tight')