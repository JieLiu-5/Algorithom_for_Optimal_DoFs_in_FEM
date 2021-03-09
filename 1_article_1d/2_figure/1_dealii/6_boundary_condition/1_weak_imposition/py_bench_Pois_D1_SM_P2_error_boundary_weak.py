#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 15:21:04 2018

@author: jliu
"""
import sys
sys.path.append('../../')
import settings as set
        
        
import matplotlib.pyplot as plt
import numpy as np
        
n_err_refer = 28

vec_legend = [r"Strong",r"$\rho$=1e6",r"$\rho$=50",r"CG, $10^{-16}$"]

vec_markeredgecolor = ["black","darkcyan","cyan","greenyellow"]

vec_markerfacecolor = ["none","none","none","none"]

vec_var = ['solu','grad','2ndd']


vec_slope = [2, 2, 2]
vec_slope_schur = [1, 1, 2.5]

vec_offset = [2e-17, 5e-17, 5e-16]

id_complex = 0


dof_ref = np.zeros(n_err_refer)
err_round_off_approx = np.zeros(n_err_refer)
err_round_off_approx_schur = np.zeros(n_err_refer)

for i in range(n_err_refer):
    dof_ref[i] = 2**i    
    
id_var_start=0
id_var_end=3

for id_var in range(id_var_start,id_var_end):
        
    print ('var: '+vec_var[id_var])     
    
    matrix_data_error = [line.strip().split() for line in open('data_error'+'_'+vec_var[id_var]+'.txt','r')]
    
    n_row = len(matrix_data_error)-1             # we omit the first row since it is title
    n_column = len(matrix_data_error[0])
    
    data_error = np.zeros((n_row,n_column))
    data_ndofs = np.zeros(n_row) 
    
    for i in range(n_row):
        data_ndofs[i]=2**(i+1)*2+1
        for p in range(n_column):
            data_error[i][p]=float(matrix_data_error[i+1][p])
            

    for i in range(n_row):
        for p in range(n_column):
            
            if(data_error[i][p]==0):
                data_error[i][p] = np.nan
                
            if(data_ndofs[i]==0):
                data_ndofs[i] = np.nan
                
    for i in range(n_err_refer):
        err_round_off_approx[i] = vec_offset[id_var]* dof_ref[i]**vec_slope[id_var]      

    tria_start = 10 
    
    tria_coeff_x = 1e1
    tria_coeff_y = 1e-1

    
    text_slope_bottom_coeff_x = 0.5
    text_slope_bottom_coeff_y = 0.15
    text_slope_right_coeff_x = 1.15
    text_slope_right_coeff_y = 0.12
    
    text_offset_coeff_y=2e-1
    
    text_offset_x = dof_ref[3]
    
    if id_var==2:
        tria_coeff_y = 1e-1
        text_offset_x = dof_ref[1]
        
        
    tria_end = tria_start+3
    
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
            
    
    text_offset_y = vec_offset[id_var]* text_offset_x**vec_slope[id_var]*text_offset_coeff_y       
    
    f=plt.figure(figsize=(5,4))
    axes= f.add_axes([0,0,1,1])
    fontsize_label = 18
    #        fontsize_legend = 20
    
    id_loc_legend = 1
    vec_marker=['o','o','o','s','*']  
    legend_size=14
    
    xaxis_up_bound=1e8
    yaxis_low_bound=1e-16
    yaxis_up_bound=1e0 
    
#    print(data_ndofs)
    
    for i in range(n_column):
        plt.loglog(data_ndofs, data_error[:,i],'kd', markerfacecolor=vec_markerfacecolor[i],markeredgecolor=vec_markeredgecolor[i], color=vec_markeredgecolor[i],linewidth=1.0,label=vec_legend[i])    

    plt.loglog(dof_ref, err_round_off_approx,'--',color="orange",linewidth=1.0)
        
    
    plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15)             
                
    plt.plot(tria_x,tria_y,'k-',linewidth=1.0)
    plt.text(text_slope_bottom_x,text_slope_bottom_y,'1',fontsize = 15)
    plt.text(text_slope_right_x,text_slope_right_y,str(vec_slope[id_var]),fontsize = 15)    
        
    if id_var == 0:    
        plt.legend()
        plt.legend(loc=id_loc_legend, prop={'size': set.fontsize_legend})
    
    plt.tick_params(axis='both', which='major', labelsize=16)
    
    plt.xlabel('Number of DoFs', fontsize=fontsize_label)             # Condition number
    plt.ylabel('Error', fontsize=fontsize_label)
    
    
    print('yaxis_low_bound: ', yaxis_low_bound)    
    print('yaxis_up_bound: ', yaxis_up_bound)
    
    axes.set_xlim([1, xaxis_up_bound])
    axes.set_ylim([yaxis_low_bound, yaxis_up_bound]) 
    
    plt.yticks([1e-16, 1e-12, 1e-8, 1e-4, 1e0]) 
        
    plt.show()
    
    f.savefig("py_bench_Pois_SM_error_boundary_weak_"+vec_var[id_var]+".pdf", bbox_inches='tight')    


