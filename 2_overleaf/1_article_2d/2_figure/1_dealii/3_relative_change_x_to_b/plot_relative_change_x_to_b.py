
import matplotlib.pyplot as plt
import numpy as np

vec_data_src=['Pois_bench_1d_dbc','Helm_d_0_1d_dbc']
id_data_src=1

vec_label=['1e-16','1e-2','1e0']
vec_marker=['o','d','^']

data_raw=[line.strip().split() for line in open("data_relative_change_"+vec_data_src[id_data_src]+".txt")]

n_row=len(data_raw)

data_ndofs=np.zeros(n_row-1)
data_ncond=np.zeros(n_row-1)
data_change=np.zeros((n_row-1,3))

for i in range(0,n_row-1):
    data_ndofs[i]=data_raw[i+1][0]
    data_ncond[i]=data_raw[i+1][1]
    data_change[i][0]=data_raw[i+1][2]
    data_change[i][1]=data_raw[i+1][3]
    data_change[i][2]=data_raw[i+1][4]

f=plt.figure(figsize=(5,4))
axes= f.add_axes([0,0,1,1])
fontsize_label = 18
plt.tick_params(axis = 'both', which = 'major', labelsize = 16)

    
for i in range(0,3):
    plt.loglog(data_ndofs,data_change[:,i],'--',marker=vec_marker[i],mfc="none",label=r'$\delta F=$'+vec_label[i]) 
    
plt.legend()
plt.legend(loc=2, prop={'size': 16})    


plt.xlabel("Number of DoFs", fontsize = fontsize_label)                   # Condition number
plt.ylabel(r"$\Re$", fontsize = fontsize_label) 


axes.set_xlim([1e0, 1e3])
axes.set_ylim([1e0, 1e4])

plt.xticks([1e0, 1e1, 1e2, 1e3])
plt.yticks([1e0,1e1,1e2,1e3,1e4])


f.savefig("py_relative_change_x_to_b_"+vec_data_src[id_data_src]+".pdf", bbox_inches='tight')   
