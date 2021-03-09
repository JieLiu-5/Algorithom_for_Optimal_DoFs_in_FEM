
import sys
sys.path.append('../')
import settings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.ticker as mtick

fig = plt.figure(figsize=(5, 4))
ax = fig.add_axes([0,0,1,1])
ax.xaxis.set_major_locator(MaxNLocator(integer=True))


n_row = 5
n_col = 3

vec_item = ['accuracy', 'cpu_pred_only', 'cpu_pred_plus_cmpt','cpu_brut','cpu_saved','n_opt']

vec_FEM = ['sm','mm']

vec_approach = ['pred','brut']

vec_legend = [r"$E_{\rm min}$",r"CPU time in seconds",r"CPU time in seconds",r"CPU time in seconds",r"Percentage",r"$N_{\rm opt}$"]
vec_marker = ["o","o--"]

coords_x_axis = np.zeros((n_row,1))
data_pure = np.zeros((n_row,n_col))

coords_x_axis = [1e1, 1e2, 1e3, 1e4, 1e5]
labels_x_axis = [1, 2, 3, 4, 5]

id_item = 4                             # '0' for E_min
                                        # '1' for CPU time used for only prediction
                                        # '2' for CPU time used for both prediction and computation
                                        # '3' for CPU time using the brute-force approach
                                        # '4' for CPU time saved
                                        # '5' for N_opt predicted

id_FEM = 1


vec_var=np.zeros((3,1))

range_ylim=np.zeros((2,1))

if id_FEM == 0:
    vec_var = [r"$u$", r"$u_x$", r"$u_{xx}$"]
    if id_item == 0:
        range_ylim=[1e-12,1e-0]
    elif id_item == 1:
        range_ylim=[0,2]
        settings.id_loc_legend = 1
    elif id_item==2 or id_item==3:
        range_ylim=[1e-2,1e6]     
    elif id_item == 4:
        range_ylim=[0,1e2]
        settings.id_loc_legend = 4
        ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    elif id_item==5:
        range_ylim=[1,1e10]

elif id_FEM == 1:
    vec_var = [r"$u$", r"$v$", r"$v_x$"]
    if id_item == 0:
        range_ylim=[1e-14,1e-4]
    elif id_item == 1:
        range_ylim=[0,2]
        settings.id_loc_legend = 1
    elif id_item==2 or id_item==3:
        range_ylim=[1e-2,1e6]
    elif id_item == 4:
        range_ylim=[0,1e2]
        settings.id_loc_legend = 4
        ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    elif id_item==5:
        range_ylim=[1,1e10]        

        

if id_item == 0:
    for id_approach in range(2):         
        f = open("data_accuracy_"+vec_approach[id_approach]+"_"+vec_FEM[id_FEM]+".txt", 'r')
        data_raw = [line.split() for line in f]
        f.close()    
        for id_var in range(3):
            for id_coeff in range(n_row):
                
                data_pure[id_coeff,id_var] = data_raw[id_coeff+1][id_var]
                if(data_pure[id_coeff,id_var]==0):
                    data_pure[id_coeff,id_var] = np.nan  
                    
            if id_approach == 0:
                plt.loglog(coords_x_axis, data_pure[:,id_var], "o-", color = settings.vec_markeredgecolor[id_var], mfc='none', label=vec_var[id_var])
            elif id_approach == 1:
                plt.loglog(coords_x_axis, data_pure[:,id_var], vec_marker[id_item], color = settings.vec_markeredgecolor[id_var])
                
elif id_item==1 or id_item==2 or id_item==3 or id_item==5:
    f = open("data_"+vec_item[id_item]+"_"+vec_FEM[id_FEM]+".txt", 'r')
    data_raw = [line.split() for line in f]
    f.close()    
    for id_var in range(3):
        for id_coeff in range(n_row):
            
            data_pure[id_coeff,id_var] = data_raw[id_coeff+1][id_var]
            if(data_pure[id_coeff,id_var]==0):
                data_pure[id_coeff,id_var] = np.nan  
                
        if id_item==1:
            plt.semilogx(coords_x_axis, data_pure[:,id_var], "o-", color = settings.vec_markeredgecolor[id_var], mfc='none', label=vec_var[id_var])
        elif id_item==2 or id_item==3 or id_item==5:
            plt.loglog(coords_x_axis, data_pure[:,id_var], "o-", color = settings.vec_markeredgecolor[id_var], mfc='none', label=vec_var[id_var])
        
        
#    f = open("data_"+vec_item[1]+"_"+vec_FEM[id_FEM]+".txt", 'r')                   # used for plotting cpu time of only prediction with that of brute force
#    data_raw = [line.split() for line in f]
#    f.close()    
#    for id_var in range(3):
#        for id_coeff in range(n_row):
#            
#            data_pure[id_coeff,id_var] = data_raw[id_coeff+1][id_var]
#            if(data_pure[id_coeff,id_var]==0):
#                data_pure[id_coeff,id_var] = np.nan  
#        plt.semilogx(coords_x_axis, data_pure[:,id_var], "o--", color = settings.vec_markeredgecolor[id_var], mfc='none') 
#             
        
elif id_item == 4:
    f = open("data_cpu_saved_"+vec_FEM[id_FEM]+".txt", 'r')
    data_raw = [line.split() for line in f]
    f.close()      
    for id_var in range(3):
        for id_coeff in range(n_row):
            
            data_pure[id_coeff,id_var] = float(data_raw[id_coeff+1][id_var])*100
            if(data_pure[id_coeff,id_var]==0):
                data_pure[id_coeff,id_var] = np.nan  
        plt.semilogx(coords_x_axis, data_pure[:,id_var], "o-", color = settings.vec_markeredgecolor[id_var], mfc='none', label=vec_var[id_var])
        
plt.minorticks_off()
plt.xticks(coords_x_axis,labels_x_axis,fontsize=settings.fontsize_tick)
plt.xlabel(r"$p$",fontsize=settings.fontsize_label)
#plt.xlim(1,5)
#
plt.yticks(fontsize=settings.fontsize_tick)
plt.ylabel(vec_legend[id_item],fontsize=settings.fontsize_label)
plt.ylim(range_ylim)
#
plt.legend(fontsize=settings.fontsize_legend, loc=settings.id_loc_legend,ncol=1)

plt.savefig('py_comparison_dofinder_'+vec_item[id_item]+'_'+vec_FEM[id_FEM]+'.pdf', bbox_inches='tight')