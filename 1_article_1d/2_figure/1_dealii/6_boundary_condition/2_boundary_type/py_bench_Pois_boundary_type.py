#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 15:21:04 2018

@author: jliu
"""

import matplotlib.pyplot as plt
import numpy as np
        
id_loc_legend = 1
vec_marker=['o','d','^','s','*']  

fontsize_label = 18
legend_size=14

xaxis_up_bound=1e8
yaxis_low_bound=1e-16
yaxis_up_bound=1e0 
    
vec_legend = [r"Dirichlet/Dirichlet",r"Dirichlet/Neumann",r"Weak, ($\rho$=1e6)",r"CG, $10^{-16}$"]

vec_markeredgecolor = ["black","green","lightgreen","greenyellow"]

vec_markerfacecolor = ["none","none","none","none"]

vec_FEM=['sm','mm']
vec_boundary=['DBC','MBC']
vec_var = ['solu','grad','2ndd']


vec_slope=np.zeros((3,1))
vec_slope_alter=np.zeros((3,1))
vec_offset=np.zeros((3,1))
vec_offset_alter=np.zeros((3,1))  
        
n_err_refer = 28

dof_ref = np.zeros(n_err_refer)
err_round_off_approx = np.zeros(n_err_refer)
err_round_off_approx_alter = np.zeros(n_err_refer)

for i in range(n_err_refer):
    dof_ref[i] = 2**i    

id_FEM = 0
id_boundary=0

id_var_start=0
id_var_end=3




if id_FEM==0:
    vec_slope = [2, 2, 2]  
    vec_slope_alter = [2, 2, 2]    
    
    vec_offset = [2e-17, 5e-17, 5e-16]
    vec_offset_alter = [7e-17, 1e-16, 5e-16]
    
elif id_FEM==1:
    vec_slope = [1, 1, 1]  
    vec_slope_alter = [1, 1, 1]
    
    vec_offset = [5e-20, 6e-17, 2e-16]
    vec_offset_alter = [3e-17, 1e-17, 2e-16]

    
for id_var in range(id_var_start,id_var_end):

    print ('var: '+vec_var[id_var])   
    
    f=plt.figure(figsize=(5,4))
    
    axes= f.add_axes([0.1,0.1,0.8,0.8])
         
    axes.set_xlim([1, xaxis_up_bound])
    axes.set_ylim([yaxis_low_bound, yaxis_up_bound]) 
    
    matrix_data_error_boundary = [line.strip().split() for line in open('data_error_'+vec_FEM[id_FEM]+'_'+vec_boundary[id_boundary]+'.txt','r')]
    
    n_row = len(matrix_data_error_boundary)-1             
    n_column = len(matrix_data_error_boundary[0])
    
    data_ndofs = np.zeros(n_row) 
    data_error = np.zeros((n_row,n_column))   
    
    for i in range(n_err_refer):
        err_round_off_approx[i] = vec_offset[id_var]* dof_ref[i]**vec_slope[id_var]      
        err_round_off_approx_alter[i] = vec_offset_alter[id_var]* dof_ref[i]**vec_slope_alter[id_var]
                
    
    for i in range(n_row):
        if id_FEM==0:
            data_ndofs[i]=2**(i+1)*2+1
        elif id_FEM==1:
            if id_var == 0:
                data_ndofs[i]=2**(i+1)*4+1
            elif id_var == 1:
                data_ndofs[i]=2**(i+1)*4
        if(data_ndofs[i]==0):
            data_ndofs[i] = np.nan    
            
    for id_boundary in range(2):
        
        matrix_data_error_boundary = [line.strip().split() for line in open('data_error_'+vec_FEM[id_FEM]+'_'+vec_boundary[id_boundary]+'.txt','r')]
                
        for p in range(n_column):
            for i in range(n_row):
                data_error[i][p]=float(matrix_data_error_boundary[i+1][p])
                if(data_error[i][p]==0):
                    data_error[i][p] = np.nan 
        
     
        plt.loglog(data_ndofs, data_error[:,id_var],'k'+vec_marker[id_var], mfc='none',markeredgecolor=vec_markeredgecolor[id_boundary], color=vec_markeredgecolor[id_var],linewidth=1.0,label=vec_legend[id_boundary])    

        plt.loglog(dof_ref, err_round_off_approx,'k--',linewidth=1.0)
        plt.loglog(dof_ref, err_round_off_approx_alter,'g--',linewidth=1.0)
        

    tria_start = 10 
    
    tria_coeff_x = 1e1
    tria_coeff_y = 1e-1
    
    text_slope_bottom_coeff_x = 0.4
    text_slope_bottom_coeff_y = 0.09
    text_slope_right_coeff_x = 1.2
    text_slope_right_coeff_y = 0.05
    

    text_offset_x = dof_ref[3]
    text_offset_x_alter = 1.5*dof_ref[0]
    
    text_offset_coeff_y=2e-1    
    text_offset_coeff_y_alter=5e3

    if id_FEM==0:
        if id_var==1:
            tria_start = 10
    
        elif id_var==2:
            tria_start = 10
            tria_coeff_y = 1e-1
            text_offset_x = dof_ref[1]
    elif id_FEM==1:
        if id_var==0:
            tria_start = 16
            text_offset_x = dof_ref[6]
            text_offset_coeff_y=45
        elif id_var==1:
            text_offset_x = dof_ref[1]
            text_offset_coeff_y=1e2
            text_offset_coeff_y_alter=1
            text_offset_x_alter = dof_ref[8]
        
        
    tria_end = tria_start+4
    
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
    text_offset_y_alter = vec_offset[id_var]* text_offset_x**vec_slope[id_var]*text_offset_coeff_y_alter      


    plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15)             
    plt.text(text_offset_x_alter,text_offset_y_alter,r'$\alpha_{\rm R}=$'+str(vec_offset_alter[id_var]),color=vec_markeredgecolor[1],fontsize = 15)             
       
    plt.plot(tria_x,tria_y,'k-',linewidth=1.0)
    plt.text(text_slope_bottom_x,text_slope_bottom_y,'1',fontsize = 15)
    plt.text(text_slope_right_x,text_slope_right_y,str(vec_slope[id_var]),fontsize = 15)    
            
    if id_var == 0:    
        plt.legend()
        plt.legend(loc=id_loc_legend, prop={'size': legend_size})
        
    plt.tick_params(axis='both', which='major', labelsize=16)
    
    plt.xlabel('Number of DoFs', fontsize=fontsize_label)     
    plt.ylabel('Absolute error', fontsize=fontsize_label)
    

    plt.yticks([1e-16, 1e-12, 1e-8, 1e-4, 1e0]) 
        
    plt.show()
    
    f.savefig("py_bench_Pois_error_boundary_type_"+vec_FEM[id_FEM]+'_'+vec_var[id_var]+".pdf", bbox_inches='tight')    


