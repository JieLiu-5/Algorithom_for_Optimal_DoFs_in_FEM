#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 15:21:04 2018

@author: jliu
"""      
        
import matplotlib.pyplot as plt
import numpy as np

import math

import sys
sys.path.append('../../../../')
import settings_fig as set


vec_equ=['0_u_1_over_c_x_m_0p5_square_p_y_same','1_u_exp_m_x_m_0p5_square_y_same']
id_equ=1

vec_FEM=['0_sm','1_mm_rt','2_mm_bdm']
id_FEM = 2

degree = 2

print("FEM:", vec_FEM[id_FEM])
print("element degree:",degree)

n_param = 5

vec_offset = [5e-17, 5e-17, 7e-17]
vec_slope = [1, 1, 1]

n_degree = 5
n_err_refer = 28

vec_legend = ['$u_{xx}$', '$u_{xy}$', '$u_{yx}$', '$u_{yy}$', R'$\Delta u$']
if id_FEM == 1 or id_FEM==2:
    vec_legend = ['$v_{1x}$', '$v_{1y}$', '$v_{2x}$', '$v_{2y}$', R'$\nabla \cdot \mathbf{v}$']

vec_markeredgecolor=['k', 'b', 'r', 'k', 'k', 'y', 'g', 'w']

id_loc_legend = 1
vec_marker=['o','d','_','o','1']  
legend_size=16
fontsize_label = 18

xaxis_up_bound=1e8
yaxis_low_bound=1e-20
yaxis_up_bound=1e4
vec_ytick = [1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e0, 1e4]

dof_ref = np.zeros(n_err_refer)
err_round_off_approx = np.zeros(n_err_refer)

for i in range(n_err_refer):
    dof_ref[i] = 2**i
   
        
matrix_data_error = [line.strip().split() for line in open(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_hessian_'+vec_FEM[id_FEM]+'.txt','r')]
n_refine = len(matrix_data_error)-1
    
data_error = np.zeros((n_refine,n_param))
data_xaxis = np.zeros(n_refine)    

id_var_start = 0
id_var_end = 5

f=plt.figure(figsize=(5,4))
axes= f.add_axes([0,0,1,1])

for id_var in range(id_var_start,id_var_end):
        
    print ('var: '+ vec_legend[id_var])      
    
    for i in range(n_refine):
        data_error[i][id_var]=matrix_data_error[i+1][id_var]
        
    for i in range(n_refine):
        if id_FEM==0:
            data_xaxis[i]=(2**(i+1)*degree+1)**2
        elif id_FEM==1:
            dofs_per_quad = 2*degree*(degree+1)
            data_xaxis[i]=2**(i+1)*(degree+1)*((2**(i+1)-1)*2+4) + 4**(i+1)*dofs_per_quad
        elif id_FEM==2:
            dofs_per_quad = degree*(degree-1)
            data_xaxis[i]=2**(i+1)*(degree+1)*((2**(i+1)-1)*2+4) + 4**(i+1)*dofs_per_quad
  
    for i in range(n_refine):
        if(data_xaxis[i]==0):
            data_xaxis[i] = np.nan        
        for p in range(n_param):
            if(data_error[i][p]==0):
                data_error[i][p] = np.nan
                
    plt.loglog(data_xaxis, data_error[:,id_var],'k'+vec_marker[id_var], markerfacecolor='none',markeredgecolor=vec_markeredgecolor[id_var],label=vec_legend[id_var],linewidth=1.0, markersize=8)
    

            
for i in range(n_err_refer):
    err_round_off_approx[i] = vec_offset[id_FEM]* dof_ref[i]**vec_slope[id_FEM]                 

tria_start = 8
tria_end = tria_start+4

tria_coeff_x = 1e1
tria_coeff_y = 1e-2


text_offset_x = dof_ref[1]

text_offset_coeff_y=0.1


text_slope_bottom_coeff_x = 0.4
text_slope_bottom_coeff_y = 0.05
text_slope_right_coeff_x = 1.2
text_slope_right_coeff_y = 0.2
    

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
text_slope_right_y = 0.5*10**((math.log10(tria_y[1])+math.log10(tria_y[2]))/2)
    

text_offset_y = vec_offset[id_FEM]* text_offset_x**vec_slope[id_FEM]*text_offset_coeff_y       


plt.legend()
plt.legend(loc=id_loc_legend, prop={'size': set.fontsize_legend})
   
plt.loglog(dof_ref, err_round_off_approx,'k--',linewidth=1.0)
plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_FEM]),fontsize = 15, color="orange")             
                
plt.plot(tria_x,tria_y,'k-',linewidth=1.0)
plt.text(text_slope_bottom_x,text_slope_bottom_y,'1',fontsize = 15)
plt.text(text_slope_right_x,text_slope_right_y,str(vec_slope[id_FEM]),fontsize = 15)
        

plt.tick_params(axis='both', which='major', labelsize=set.fontsize_tick)

plt.xlabel("Number of DoFs", fontsize=set.fontsize_label)
plt.ylabel('Error', fontsize=set.fontsize_label)

print('yaxis_low_bound: ', yaxis_low_bound)    
print('yaxis_up_bound: ', yaxis_up_bound)   

axes.set_xlim([1, xaxis_up_bound])
axes.set_ylim([yaxis_low_bound, yaxis_up_bound]) 

plt.yticks(vec_ytick) 
    
plt.show()

f.savefig(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+"/py_error_hessian_twod_"+vec_equ[id_equ]+"_"+vec_FEM[id_FEM]+".pdf", bbox_inches='tight')


