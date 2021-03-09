#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 15:21:04 2018

@author: jliu
"""    
        
import matplotlib.pyplot as plt
import numpy as np

import sys
sys.path.append('../../../../3_1d_or_2d')
import settings_fig as set

degree=5
degree_for_velocity=4
id_refine=1
dofs_per_quad=1
        
vec_x_axis = ["1_ndofs_as_x", "2_ncond_as_x"]

vec_legend = [r"$P_5$ (deal.II)",r"$RT_5/P_{5}^{\rm disc}$ (deal.II)",r"$BDM_6/Q_{5}^{\rm disc}$ (deal.II)", r"$P_5$ (FEniCS)",r"$BDM_6/Q_{5}^{\rm disc}$ (FEniCS)", r"$S_{5,4}$ (G+Smo)"]
    
vec_markeredgecolor = ["black","darkcyan","cyan","greenyellow", "sandybrown", "brown","r"]

vec_var = ['solu','grad','2ndd']

vec_coord_x = ['ndofs','ncond']
vec_label_coord_x=['Number of DoFs','Condition number']
id_x_axis=0

vec_slope = [2, 2, 2]  
vec_slope_alter = [2, 2, 2]

vec_offset = [2e-17, 5e-17, 5e-16]
vec_offset_alter = [7e-17, 1e-16, 5e-16]

is_direct_results_together=0
    
id_loc_legend = 1

xaxis_up_bound=1e8
yaxis_low_bound=1e-16
yaxis_up_bound=1e0 

n_err_refer = 28
dof_ref = np.zeros(n_err_refer)
err_round_off_approx = np.zeros(n_err_refer)
err_round_off_approx_alter = np.zeros(n_err_refer)

for i in range(n_err_refer):
    dof_ref[i] = 2**i
           
if id_x_axis==0:
    matrix_data_coord_x_iga = [line.strip().split() for line in open(vec_x_axis[id_x_axis]+'/data_ndofs_gismo.txt','r')]
elif id_x_axis==1:
    matrix_data_coord_x_std_dealii = [line.strip().split() for line in open(vec_x_axis[id_x_axis]+'/data_ncond_std_dealii.txt','r')]
    matrix_data_coord_x_PpPpm1 = [line.strip().split() for line in open(vec_x_axis[id_x_axis]+'/data_ncond_mix_PpPpm1.txt','r')]
    matrix_data_coord_x_RTpPp = [line.strip().split() for line in open(vec_x_axis[id_x_axis]+'/data_ncond_mix_RTpPp.txt','r')]
    is_direct_results_together=1
    
if is_direct_results_together==0:
    matrix_data_error = [line.strip().split() for line in open('data_error'+'_'+vec_var[0]+'.txt','r')]
else:
    matrix_data_error = [line.strip().split() for line in open('data_error_direct.txt','r')]    
            
n_row = len(matrix_data_error)-1
n_column = len(matrix_data_error[0])

data_ndofs = np.zeros(n_row)

data_coord_x_analytical = np.zeros((n_row,n_column))
data_error = np.zeros((n_row,n_column))
     
id_var_start=2
id_var_end=3

if(is_direct_results_together==1):
    id_var_end=1

for id_var in range(id_var_start,id_var_end):
        
    print ('var: '+vec_var[id_var])     
    
    if id_x_axis==0:
        if id_var==1 or id_var==2:
            vec_legend = [r"$P_5$ (deal.II)",r"$RT_4/P_{3}^{\rm disc}$ (deal.II)",r"$BDM_4/Q_{3}^{\rm disc}$ (deal.II)", r"$P_5$ (FEniCS)",r"$BDM_4/Q_{3}^{\rm disc}$ (FEniCS)", r"$S_{5,4}$ (G+Smo)"]

    
    if is_direct_results_together==0:
        matrix_data_error = [line.strip().split() for line in open('data_error'+'_'+vec_var[id_var]+'.txt','r')]
    else:
        matrix_data_error = [line.strip().split() for line in open('data_error_direct.txt','r')]        
             
    for i in range(n_row):
        
        id_refine=i+1
        
        for p in range(n_column):
            data_error[i][p]=float(matrix_data_error[i+1][p])
            
        data_coord_x_analytical[i][0]=(2**id_refine*degree+1)**2
        data_coord_x_analytical[i][3]=(2**id_refine*degree+1+2**id_refine*(degree-1))*(2**id_refine+1)+4**id_refine*(4*(degree-1+(degree-2)*(degree-1)/2)+1)
        
        if id_var==0:
            
            data_coord_x_analytical[i][1]=data_coord_x_analytical[i][1]
            data_coord_x_analytical[i][2]=4**id_refine*(degree+1)*(degree+2)/2
            
            data_coord_x_analytical[i][4]=4**id_refine*4*(degree+1)*(degree+1+1)/2
            
        elif id_var==1 or id_var==2:
            
            dofs_per_quad = 2*degree_for_velocity*(degree_for_velocity+1)
            data_coord_x_analytical[i][1]=2**id_refine*(degree_for_velocity+1)*((2**id_refine-1)*2+4) + 4**id_refine*dofs_per_quad      
            
            dofs_per_quad = degree_for_velocity*(degree_for_velocity-1)
            data_coord_x_analytical[i][2]=2**id_refine*(degree_for_velocity+1)*((2**id_refine-1)*2+4) + 4**id_refine*dofs_per_quad      
            
            data_coord_x_analytical[i][4]=2**id_refine*(degree_for_velocity+1)*(2**id_refine+1)*2 + 4**id_refine*(4*(degree_for_velocity+1)+4*(degree_for_velocity**2-1))
    
    for i in range(n_row):
        for p in range(n_column):
            
            if(data_error[i][p]==0):
                data_error[i][p] = np.nan
                
            if(data_coord_x_analytical[i][p]==0):
                data_coord_x_analytical[i][p] = np.nan                
                
    for i in range(n_err_refer):
        err_round_off_approx[i] = vec_offset[id_var]* dof_ref[i]**vec_slope[id_var]      
        err_round_off_approx_alter[i] = vec_offset_alter[id_var]* dof_ref[i]**vec_slope_alter[id_var]
    
    tria_start = 10 
    
    tria_coeff_x = 1e1
    tria_coeff_y = 1e-1

    
    text_slope_bottom_coeff_x = 0.4
    text_slope_bottom_coeff_y = 0.09
    text_slope_right_coeff_x = 1.2
    text_slope_right_coeff_y = 0.05
    
    text_offset_coeff_y=2e-1
    
    
    
    tria_start_alter = 6
    
    tria_coeff_x_alter = 1e1
    tria_coeff_y_alter = 1e-1    
    
    text_slope_right_coeff_y_alter = 0.03
    
    
    text_offset_x = dof_ref[3]
    
    text_offset_x_alter = 1.5*dof_ref[0]
    text_offset_coeff_y_alter=5e4
    
    if id_var==1:
        tria_start = 10
    elif id_var==2:
        tria_start = 10
        tria_coeff_y = 1e-1
        text_offset_x = dof_ref[1]
        id_loc_legend=3
        
        
    tria_end = tria_start+4
    tria_end_alter = tria_start_alter+4
    
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
    
    
    f=plt.figure(figsize=(5,4))
    axes= f.add_axes([0.1,0.1,0.8,0.8])
    fontsize_label = 18
    
    for i in range(n_column):
        plt.loglog(data_coord_x_analytical[:,i], data_error[:,i],'k'+set.vec_marker[i]+'--', markerfacecolor=set.vec_markerfacecolor[i],markeredgecolor=vec_markeredgecolor[i], color=vec_markeredgecolor[i],linewidth=1.0,label=vec_legend[i])    

    plt.legend()
    plt.legend(loc=id_loc_legend, prop={'size': 12})
    
    plt.tick_params(axis='both', which='major', labelsize=16)
    
    plt.xlabel(vec_label_coord_x[id_x_axis], fontsize=fontsize_label)             # Condition number
    plt.ylabel('Error', fontsize=fontsize_label)
    
    
    print('yaxis_low_bound: ', yaxis_low_bound)
    print('yaxis_up_bound: ', yaxis_up_bound)
    
    axes.set_xlim([1, xaxis_up_bound])
    axes.set_ylim([yaxis_low_bound, yaxis_up_bound]) 
    
    plt.yticks([1e-16, 1e-12, 1e-8, 1e-4, 1e0]) 
        
    plt.show()
    
    if is_direct_results_together==1:
        f.savefig(vec_x_axis[id_x_axis]+"/py_bench_pois_twod_pexp_deg_5_UMFPACK_std_vs_mix_"+vec_coord_x[id_x_axis]+"_direct_result.pdf", bbox_inches='tight')    
    else:
        f.savefig(vec_x_axis[id_x_axis]+"/py_bench_pois_twod_pexp_deg_5_UMFPACK_std_vs_mix_"+vec_coord_x[id_x_axis]+"_"+vec_var[id_var]+".pdf", bbox_inches='tight')    


