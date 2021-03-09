#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 15:21:04 2018

@author: jliu
"""  
        
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append('C:/Users/LIUJI/Dropbox/3_writing/1_1d/3_figure_plot')
import settings
        
n_err_refer = 28

vec_markeredgecolor = ["black","darkcyan","cyan","greenyellow"]

vec_markerfacecolor = ["none","none","none","none"]

vec_var = ['solu','grad','2ndd']

vec_FEM = ["sm", "mm"]

xaxis_up_bound=1e8
yaxis_low_bound=1e-16
yaxis_up_bound=1e0 

id_FEM = 0
id_parameter = 0

id_var_start = 0
id_var_end = 3

id_parameter_start = 0
id_parameter_end = 3

id_loc_legend = 1


dof_ref = np.zeros(n_err_refer)
err_round_off_approx = np.zeros(n_err_refer)
err_round_off_approx_alter = np.zeros(n_err_refer)


for i in range(n_err_refer):
    dof_ref[i] = 2**i    
    

n_var = 3
n_parameter = 4


vec_error_plot_start=[0, 12, 8]
vec_parameter = [r"UMFPACK",r"tol_1em10",r"tol_1em4"]
vec_legend = [r"UMFPACK",r"$tol_{prm}=10^{-10}$",r"$tol_{prm}=10^{-4}$"]
vec_deg = [2,2,2]
vec_slope = [2, 2, 2]  
vec_slope_alter = [1, 1, 2.5]    
vec_offset = [2e-17, 5e-17, 5e-16]
vec_offset_alter = [5e-19, 1e-17, 1e-16]
id_loc_legend = 4
           

for id_var in range(id_var_start, id_var_end):
        
    print ('var: '+vec_var[id_var])     
    
    f=plt.figure(figsize=(5,4))
    axes= f.add_axes([0,0,1,1])
    
    data_raw = [line.split() for line in open('data_error_sm_'+vec_parameter[id_parameter]+'.txt','r')]
    
    n_row = len(data_raw)-1         
    
    data_error = np.zeros((n_row,id_parameter_end-id_parameter_start))
    data_ndofs = np.zeros(n_row) 
    
    for id_parameter in range(id_parameter_start, id_parameter_end):
        
        print(vec_parameter[id_parameter])
        data_raw = [line.split() for line in open('data_error_sm_'+vec_parameter[id_parameter]+'.txt','r')]
        
        for i in range(n_row):
            data_error[i][id_parameter-id_parameter_start]=data_raw[i+1][id_var]    
    
    for i in range(n_row):
        data_ndofs[i]=2**(i+1)*vec_deg[id_var]+1
        
    for i in range(n_row):
        for p in range(id_parameter_start, id_parameter_end):
            
            if(data_error[i][p]==0):
                data_error[i][p] = np.nan
                
            if(data_ndofs[i]==0):
                data_ndofs[i] = np.nan
                
    for i in range(n_err_refer):
        err_round_off_approx[i] = vec_offset[id_var]* dof_ref[i]**vec_slope[id_var]      
        err_round_off_approx_alter[i] = vec_offset_alter[id_var]* dof_ref[i]**vec_slope_alter[id_var]
    
    tria_start = 17
    tria_start_alter = 6
    
    
    tria_coeff_x = 1e1
    tria_coeff_y = 1e-1
    
    tria_coeff_x_alter = 1e1
    tria_coeff_y_alter = 1e-1      
    
    text_slope_bottom_coeff_x = 0.4
    text_slope_bottom_coeff_y = 0.15
    text_slope_right_coeff_x = 1.2
    text_slope_right_coeff_y = 0.06
    
    text_offset_coeff_y=2e-1    
    text_offset_coeff_y_alter=2e-1

    text_slope_right_coeff_y_alter = 0.035
    
    text_offset_x = dof_ref[3]
    text_offset_x_alter = dof_ref[2]
    
    if id_var==1:
        tria_start = 8
    elif id_var==2:
        tria_start = 18
        text_offset_x = dof_ref[14]      
        
    tria_end = tria_start+4
    tria_end_alter = tria_start_alter+4
    
    tria_p_1 = [dof_ref[tria_start],err_round_off_approx[tria_start]]
    tria_p_2 = [dof_ref[tria_end],err_round_off_approx[tria_start]]
    tria_p_3 = [dof_ref[tria_end],err_round_off_approx[tria_end]]
    
    tria_x = [tria_p_1[0],tria_p_2[0],tria_p_3[0],tria_p_1[0]]
    tria_x = [i*tria_coeff_x for i in tria_x]
    tria_y = [tria_p_1[1],tria_p_2[1],tria_p_3[1],tria_p_1[1]]
    tria_y = [i*tria_coeff_y for i in tria_y]
    
    text_slope_bottom_x = (tria_x[0]+tria_x[1])/2*text_slope_bottom_coeff_x
    text_slope_bottom_y = tria_y[0]*text_slope_bottom_coeff_y
    text_slope_right_x = tria_x[1]*text_slope_right_coeff_x
    text_slope_right_y = (tria_y[1]+tria_y[2])/2*text_slope_right_coeff_y    
            

    tria_p_1_alter = [dof_ref[tria_start_alter],err_round_off_approx_alter[tria_start_alter]]
    tria_p_2_alter = [dof_ref[tria_end_alter],err_round_off_approx_alter[tria_start_alter]]
    tria_p_3_alter = [dof_ref[tria_end_alter],err_round_off_approx_alter[tria_end_alter]]
    
    tria_x_alter = [tria_p_1_alter[0],tria_p_2_alter[0],tria_p_3_alter[0],tria_p_1_alter[0]]
    tria_x_alter = [i*tria_coeff_x_alter for i in tria_x_alter]
    tria_y_alter = [tria_p_1_alter[1],tria_p_2_alter[1],tria_p_3_alter[1],tria_p_1_alter[1]]
    tria_y_alter = [i*tria_coeff_y_alter for i in tria_y_alter]
    
    text_slope_bottom_x_alter = (tria_x_alter[0]+tria_x_alter[1])/2*text_slope_bottom_coeff_x
    text_slope_bottom_y_alter = tria_y_alter[0]*text_slope_bottom_coeff_y
    text_slope_right_x_alter = tria_x_alter[1]*text_slope_right_coeff_x
    text_slope_right_y_alter = (tria_y_alter[1]+tria_y_alter[2])/2*text_slope_right_coeff_y_alter  
    
    
    
    text_offset_y = vec_offset[id_var]* text_offset_x**vec_slope[id_var]*text_offset_coeff_y       
    text_offset_y_alter = vec_offset_alter[id_var]* text_offset_x_alter**vec_slope_alter[id_var]*text_offset_coeff_y_alter      
    

    for id_parameter in range(id_parameter_start, id_parameter_end):
        plt.loglog(data_ndofs[vec_error_plot_start[id_parameter]:], data_error[vec_error_plot_start[id_parameter]:,id_parameter],'kd', markerfacecolor=vec_markerfacecolor[id_parameter],markeredgecolor=vec_markeredgecolor[id_parameter], color=vec_markeredgecolor[id_parameter],linewidth=1.0,label=vec_legend[id_parameter])    
        
    plt.loglog(dof_ref, err_round_off_approx,'--',color="orange",linewidth=1.0)
        
    
    plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15)             
                
    plt.plot(tria_x,tria_y,'k-',linewidth=1.0)
    plt.text(text_slope_bottom_x,text_slope_bottom_y,'1',fontsize = 15)
    plt.text(text_slope_right_x,text_slope_right_y,str(vec_slope[id_var]),fontsize = 15)    
        
    if id_var == 0:    
        plt.legend()
        plt.legend(loc=id_loc_legend, prop={'size': settings.fontsize_legend})
    
    elif id_var==2:
        plt.loglog(dof_ref, err_round_off_approx_alter,'k--',linewidth=1.0,color=vec_markeredgecolor[1])
        
        plt.plot(tria_x_alter,tria_y_alter,'k-',linewidth=1.0,color=vec_markeredgecolor[1])
        plt.text(text_slope_bottom_x_alter,text_slope_bottom_y_alter,'1',fontsize = 15,color=vec_markeredgecolor[1])
        plt.text(text_slope_right_x_alter,text_slope_right_y_alter,str(vec_slope_alter[id_var]),fontsize = 15,color=vec_markeredgecolor[1])    
            
        plt.text(text_offset_x_alter,text_offset_y_alter,r'$\alpha_{\rm R}=$'+str(vec_offset_alter[id_var]),fontsize = 15,color=vec_markeredgecolor[1])             
       
    
    plt.tick_params(axis='both', which='major', labelsize=settings.fontsize_tick)
    
    plt.xlabel('Number of DoFs', fontsize=settings.fontsize_label)            
    plt.ylabel('Error', fontsize=settings.fontsize_label)
    
    axes.set_xlim([1, xaxis_up_bound])
    axes.set_ylim([yaxis_low_bound, yaxis_up_bound]) 
    
    plt.yticks([1e-16, 1e-12, 1e-8, 1e-4, 1e0]) 
        
    plt.show()
    
    f.savefig("py_error_comparison_solution_strategy_"+vec_FEM[id_FEM]+"_"+vec_var[id_var]+".pdf", bbox_inches='tight')    


