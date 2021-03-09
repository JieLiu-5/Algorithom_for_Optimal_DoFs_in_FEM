#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 15:21:04 2018

@author: jliu
"""

import matplotlib.pyplot as plt
import numpy as np
#import csv




id_FEM=1        

vec_deg = [5,6,7,8,9]       # [1,2,3,4,5]   [5,6,7,8,9]

vec_legend = [r"1",r"2",r"3",r"4",r"5"]

id_loc_legend = 2
vec_marker=['o','d','^','s','*']  
legend_size=12

f=plt.figure(figsize=(5,4))
axes= f.add_axes([0.1,0.1,0.8,0.8])
fontsize_label = 18

xaxis_up_bound=1e4
yaxis_low_bound=1e0
yaxis_up_bound=1e8



vec_var = ['solu','grad','2ndd']

scalar_slope = 2
scalar_offset = 2e-1

vec_dim=['1d','2d']
id_dim=1

n_err_refer = 28

dof_ref = np.zeros(n_err_refer)
err_round_off_approx = np.zeros(n_err_refer)

for i in range(n_err_refer):
    dof_ref[i] = 2**i 
    

if id_dim==1:
    scalar_slope = 1
    xaxis_up_bound=1e8
    yaxis_up_bound=1e6
    
    id_loc_legend=1

matrix_data_ncond = [line.strip().split() for line in open('data_ncond_'+vec_dim[id_dim]+'.txt','r')]
matrix_data_ndofs = [line.strip().split() for line in open('data_ndofs_'+vec_dim[id_dim]+'.txt','r')]

n_refine = len(matrix_data_ndofs)-1
n_degree = len(matrix_data_ndofs[1])

data_ncond = np.zeros((n_refine,n_degree))
data_ndofs = np.zeros((n_refine,n_degree))


for i in range(n_refine):
    for p in range(n_degree):
        data_ndofs[i][p]=matrix_data_ndofs[i+1][p]
        data_ncond[i][p]=matrix_data_ncond[i+1][p]
            
for i in range(n_refine):
    for p in range(n_degree):
        if(data_ndofs[i][p]==0):
            data_ndofs[i][p] = np.nan                
        if(data_ncond[i][p]==0):
            data_ncond[i][p] = np.nan
            

tria_start = 4
tria_step=2

tria_coeff_x = 1e1
tria_coeff_y = 4e-0

text_offset_coeff_y=2e-1

text_slope_bottom_coeff_x = 0.75
text_slope_bottom_coeff_y = 0.3
text_slope_right_coeff_x = 1.15
text_slope_right_coeff_y = 0.3

text_offset_x = dof_ref[3]

if id_dim==1:
    tria_start = 9
    tria_step=4
    tria_coeff_y = 4e-1
    text_slope_bottom_coeff_x = 0.4
    text_slope_bottom_coeff_y = 0.37
    text_slope_right_coeff_x = 1.2
    text_slope_right_coeff_y = 0.35
    text_offset_x = dof_ref[6]
    text_offset_coeff_y=3e-1
    
    scalar_offset = 1e-1
    
tria_end = tria_start+tria_step

for i in range(n_err_refer):
    err_round_off_approx[i] = scalar_offset* dof_ref[i]**scalar_slope                 

tria_p_1 = [dof_ref[tria_start],err_round_off_approx[tria_start]]
tria_p_2 = [dof_ref[tria_end],err_round_off_approx[tria_start]]
tria_p_3 = [dof_ref[tria_end],err_round_off_approx[tria_end]]

tria_x = [tria_p_1[0],tria_p_2[0],tria_p_3[0],tria_p_1[0]]#*
tria_x = [i*tria_coeff_x for i in tria_x]
tria_y = [tria_p_1[1],tria_p_2[1],tria_p_3[1],tria_p_1[1]]#/tria_coeff_y
tria_y = [i*tria_coeff_y for i in tria_y]

text_slope_bottom_x = (tria_x[0]+tria_x[1])/2*text_slope_bottom_coeff_x
text_slope_bottom_y = tria_y[0]*text_slope_bottom_coeff_y
text_slope_right_x = tria_x[1]*text_slope_right_coeff_x
text_slope_right_y = (tria_y[1]+tria_y[2])/2*text_slope_right_coeff_y    

    
text_offset_y = scalar_offset* text_offset_x**scalar_slope*text_offset_coeff_y       



column_start=0

for p in range(column_start,column_start+5):
    plt.loglog(data_ndofs[:,p], data_ncond[:,p],'k'+vec_marker[p-column_start]+'-', markerfacecolor='none',label='$p=$'+vec_legend[p],linewidth=1.0)            # data_ncond  data_ndofs       


plt.plot(tria_x,tria_y,'k-',linewidth=1.0)
plt.text(text_slope_bottom_x,text_slope_bottom_y,'1',fontsize = 15)
plt.text(text_slope_right_x,text_slope_right_y,str(scalar_slope),fontsize = 15)         
       
plt.legend()
plt.legend(loc=id_loc_legend, prop={'size': legend_size})
  
    
plt.loglog(dof_ref, err_round_off_approx,'k--',linewidth=1.0)

plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(scalar_offset),fontsize = 15)             
           
    

plt.tick_params(axis='both', which='major', labelsize=16)

plt.xlabel('Number of DoFs', fontsize=fontsize_label)             # Condition number   Number of DoFs
plt.ylabel('Condition number', fontsize=fontsize_label)


print('yaxis_low_bound: ', yaxis_low_bound)    
print('yaxis_up_bound: ', yaxis_up_bound)

axes.set_xlim([1, xaxis_up_bound])
axes.set_ylim([yaxis_low_bound, yaxis_up_bound]) 

#plt.yticks([1e0,1e1,1e2,1e3,1e4]) 
    
plt.show()

f.savefig("py_cond_with_dofs_iga_"+vec_dim[id_dim]+".pdf", bbox_inches='tight')    # cond  dofs


