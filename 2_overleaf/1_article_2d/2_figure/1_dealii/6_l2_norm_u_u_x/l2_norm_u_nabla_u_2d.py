
import sys
sys.path.append('../../../../3_1d_or_2d')
import settings_fig as set

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(5, 4))
ax = fig.add_axes([0,0,1,1])

n_row = 5
n_col = 3

data_l2_norm_u = np.zeros((n_row,n_col))
data_l2_norm_u_x = np.zeros((n_row,n_col))

vec_case=['1_u_1_over_c_x_m_0p5_square_p_1_over_c_y_m_0p5_square','2_u_x_square_m_c_square','3_u_cx_plus_100_square']

for id_case in range(n_col):
    f = open("l2_norm_"+vec_case[id_case]+".txt", 'r')
    data_raw = [line.split() for line in f]
    f.close()

    for id_coeff in range(n_row):
        data_l2_norm_u[id_coeff,id_case] = data_raw[id_coeff+1][1]
        data_l2_norm_u_x[id_coeff,id_case] = data_raw[id_coeff+1][2]
        
    plt.loglog(data_l2_norm_u[:,id_case], data_l2_norm_u_x[:,id_case], "o-",linestyle='solid',color=set.vec_markeredgecolor[id_case], markeredgecolor=set.vec_markeredgecolor[id_case], mfc='none',label="Case "+str(id_case+1))

plt.xticks(fontsize=set.fontsize_tick)
plt.xlabel(r"$\||u\||_2$",fontsize=set.fontsize_label)
plt.xlim(1e-6,1e6)

plt.yticks(fontsize=set.fontsize_tick)
plt.ylabel(r"$\||\nabla u\||_2$",fontsize=set.fontsize_label)
plt.ylim(1e-6,1e6)

plt.legend(fontsize=set.fontsize_legend, loc=2,ncol=1)

plt.savefig('l2_norm_u_nabla_u_2d.pdf', bbox_inches='tight')