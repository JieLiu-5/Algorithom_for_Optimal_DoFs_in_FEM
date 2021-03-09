#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 15:21:04 2018

@author: jliu
"""
        
import sys
sys.path.append('../../../../2_overleaf/3_1d_or_2d')
import settings_fig as settings

import matplotlib.pyplot as plt
import numpy as np

vec_marker=['o','d','^','s','*']  

vec_bench_equ = ['bench_Pois','bench_diff','bench_Helm']
vec_legend_bench_sm = ['$P_1$','$P_2$','$P_3$','$P_4$','$P_5$']
vec_legend_bench_mm = [r"$P_1/P_0^{\rm disc}$",r"$P_2/P_1^{\rm disc}$",r"$P_3/P_2^{\rm disc}$",r"$P_4/P_3^{\rm disc}$",r"$P_5/P_4^{\rm disc}$"]

vec_L2_Pois_equ = ['L2_Pois1','L2_Pois2','L2_Pois3','L2_Pois4','L2_Pois5']

vec_infl_d_equ = ['bench_Pois_2_diff_d_1_plus_sincx', 'bench_Pois_2_diff_d_1_plus_cx','bench_Pois_2_diff_d_c']

vec_infl_r_equ = ['bench_Pois_2_Helm_d_1_plus_sinx_r_c', 'bench_Pois_2_Helm_d_1_r_c', 'bench_Pois_2_Helm_d_1_r_1_plus_cx']

xaxis_up_bound = 1e8
yaxis_low_bound = 1e-16
yaxis_up_bound = 1e0

vec_FEM = ['SM','MM']
vec_scaling_scheme_SM = ['scaling_no', 'scaling_S']
vec_scaling_scheme_MM = ['scaling_no', 'scaling_M1', 'scaling_M2']

vec_var = ['solu','grad','2ndd']

is_complex = 0

vec_offset_L2_Pois1_S=[2e-17, 5e-17, 5e-16]
vec_offset_L2_Pois2_S=[2e-17, 5e-17, 5e-16]
vec_offset_L2_Pois3_S=[2e-17, 5e-17, 5e-16]
vec_offset_L2_Pois4_S=[2e-17, 5e-17, 5e-16]
vec_offset_L2_Pois5_S=[2e-17, 5e-17, 5e-16]

vec_offset_L2_Poisx_SM_normal=np.zeros(3)
vec_offset_L2_Pois1_M1 = [5e-18, 5e-17, 5e-16]
vec_offset_L2_Pois1_M2 = [1e-19, 5e-17, 1e-15]
vec_offset_L2_Pois2_M1 = [1e-19, 5e-17, 5e-16]
vec_offset_L2_Pois2_M2 = [1e-19, 1e-16, 1e-15]
vec_offset_L2_Pois3_M1 = [1e-18, 2e-17, 5e-16]
vec_offset_L2_Pois3_M2 = [1e-19, 1e-16, 5e-16]
vec_offset_L2_Pois4_M1 = [1e-18, 2e-17, 5e-16]
vec_offset_L2_Pois4_M2 = [1e-19, 1e-16, 5e-16]
vec_offset_L2_Pois5_M1 = [5e-18, 3e-17, 5e-16]
vec_offset_L2_Pois5_M2 = [1e-18, 1e-16, 5e-16]


vec_offset_L2_Poisx_MM_normal = np.zeros(3)
vec_offset_bench_Pois_SM = [2e-17, 5e-17, 1e-15]
vec_offset_bench_Pois_MM = [1e-19, 5e-17, 2e-16]
vec_offset_bench_diff_SM = [2e-17, 2e-17, 5e-16]
vec_offset_bench_diff_MM = [1e-19, 1e-17, 5e-15]
vec_offset_bench_Helm_SM = [2e-17, 2e-17, 5e-16]
vec_offset_bench_Helm_MM = [1e-17, 1e-17, 2e-16]

vec_L2_u_Pois1 = [9.2e0, 8.8e-1, 1.8e-2, 1.8e-4, 1.8e-6]
vec_L2_u_Pois2 = [1.0, 1.0, 0.92, 0.35, 0.11]
vec_L2_u_Pois3 = [9.2e2, 9.0, 0.23, 0.22, 0.22]
vec_L2_u_Pois4 = [0.58, 0.55, 0.11, 1.1e-2, 1.1e-3]
vec_L2_u_Pois5 = [5.8e3, 5.8e1, 5.8e-1, 5.8e-3, 5.8e-5]

vec_L2_du_Pois1 = [1.6e1, 1.5e0, 1.1e-1, 1.1e-2, 1.1e-3]
vec_L2_du_Pois2 = [5.8e-5, 5.8e-3, 5.0e-1, 3.5, 11]
vec_L2_du_Pois3 = [1.6e3, 15.4, 0.59, 0.58, 0.58]
vec_L2_du_Pois4 = [1.0, 0.94, 0.71, 0.71, 0.71]
vec_L2_du_Pois5 = [1e4, 1e2, 1, 1e-2, 1e-4]

vec_L2_du = np.zeros(3)
vec_L2_u=np.zeros(3)

vec_offset_bench_Pois_2_diff_d_1_plus_sincx_SM = [2e-18, 5e-18, 5e-16]    
vec_offset_bench_Pois_2_diff_d_1_plus_sincx_MM = [1e-19, 5e-18, 2e-16] 

vec_offset_bench_Pois_2_diff_d_1_plus_cx_SM = [2e-18, 1e-16, 5e-16]    
vec_offset_bench_Pois_2_diff_d_1_plus_cx_MM = [2e-18, 7e-15, 1e-12] 

vec_offset_bench_Pois_2_diff_d_c_SM = [5e-20, 2e-18, 5e-16]
vec_offset_bench_Pois_2_diff_d_c_MM = [6e-20, 5e-18, 2e-18]

vec_offset_bench_Pois_2_Helm_d_1_plus_sinx_r_c_SM = [1e-17/1e1, 2e-17/4e0, 2e-16]     
vec_offset_bench_Pois_2_Helm_d_1_plus_sinx_r_c_MM = [1e-20, 5e-18, 2e-16] 

vec_offset_bench_Pois_2_Helm_d_1_r_c_SM = [2e-17, 5e-17, 2e-16]     
vec_offset_bench_Pois_2_Helm_d_1_r_c_MM = [1e-19, 1e-16, 2e-16] 

vec_offset_bench_Pois_2_Helm_d_1_r_1_plus_cx_SM = [1e-17, 2e-17, 2e-16]     
vec_offset_bench_Pois_2_Helm_d_1_r_1_plus_cx_MM = [1e-19, 1e-16, 2e-16] 

n_coeff = 5
vec_coeff = ['1e0', '1e1', '1e2', '1e3', '1e4']

vec_slope = [2, 1]
vec_legend_coeff=['$c_1$','$c_2$','$c_3$','$c_4$','$c_5$']

vec_sensitivity = ["CG_1em16", "CG_1em10", "CG_1em6", "MBC"]

tria_start = 8
tria_end = tria_start+3

tria_coeff_x = 1e1
tria_coeff_y = 1e1

text_slop_bottom_coeff_x = 0.5
text_slop_bottom_coeff_y = 12
text_slop_right_coeff_x = 1.3
text_slop_right_coeff_y = 9

vec_id_equ = vec_bench_equ  

id_loc_legend = 1
fontsize_legend = settings.fontsize_legend
vec_legend=vec_legend_bench_sm

fontsize_label = settings.fontsize_label

fontsize_tick = settings.fontsize_tick

text_offset_loc_x = 4
text_offset_loc_y = 0.1


id_vec_id_equ = 2                            # '0' for benchmark equations
                                             # '1' for Poisson equations with different L2 norms
                                             # '2' for diffusion equations from Poisson
                                             # '3' for Hemlholtz equations from Poisson
                                             # '4' for sensitivity with different element degrees

id_equ = 0
id_FEM = 0

id_sensitivity = 3

id_scaling = 0
id_var_start = 0
id_var_end = 3

n_column = 5
id_case_start = 0
id_case_end = id_case_start+n_column
id_c=2
n_refine = 20
n_err_refer = 28


if(id_FEM==0):
    vec_scaling_scheme = vec_scaling_scheme_SM
elif(id_FEM==1):
    vec_scaling_scheme = vec_scaling_scheme_MM    

if id_vec_id_equ == 0:
    if id_FEM == 0:    
        if id_equ == 0:
            vec_offset = vec_offset_bench_Pois_SM
        elif id_equ == 1:
            vec_offset = vec_offset_bench_diff_SM
        elif id_equ == 2:
            is_complex=1
            vec_offset = vec_offset_bench_Helm_SM    
    elif id_FEM == 1:
        vec_legend = vec_legend_bench_mm
        if id_equ == 0:
            vec_offset = vec_offset_bench_Pois_MM
        elif id_equ == 1:
            vec_offset = vec_offset_bench_diff_MM
        elif id_equ == 2:
            is_complex=1
            vec_offset = vec_offset_bench_Helm_MM            
        
elif id_vec_id_equ ==1:
    vec_id_equ = vec_L2_Pois_equ  
    vec_legend=['1e-4','1e-2','1e0','1e2','1e4']
    if (id_equ == 0 or id_equ == 3):
        vec_legend=['1e-2','1e-1','1e0','1e1','1e2']
    
    yaxis_low_bound = 1e-16
    tria_coeff_x = 1e1
    tria_coeff_y = 1e-1 #*10**float(vec_legend[4])        #1e5 1e1
    text_slop_bottom_coeff_x = 0.5
    text_slop_bottom_coeff_y = 0.05        # 40    20
    text_slop_right_coeff_x = 1.2
    text_slop_right_coeff_y = 0.3         # 15  20

        
    if id_FEM==0:
        if id_equ==0:
            vec_L2_u = vec_L2_u_Pois1
            vec_offset_L2_Poisx_S=vec_offset_L2_Pois1_S
        elif id_equ==1:
            vec_L2_u = vec_L2_u_Pois2
            vec_offset_L2_Poisx_S=vec_offset_L2_Pois2_S        
        elif id_equ==2:
            vec_L2_u = vec_L2_u_Pois3
            vec_offset_L2_Poisx_S=vec_offset_L2_Pois3_S        
        elif id_equ==3:
            vec_L2_u = vec_L2_u_Pois4
            vec_offset_L2_Poisx_S=vec_offset_L2_Pois4_S        
        elif id_equ==4:
            vec_L2_u = vec_L2_u_Pois5
            vec_offset_L2_Poisx_S=vec_offset_L2_Pois5_S      
        
        vec_offset_L2_Poisx_SM_normal=[i*vec_L2_u[id_c] for i in vec_offset_L2_Poisx_S]

        
        if id_scaling==0:
            vec_offset=vec_offset_L2_Poisx_SM_normal
        elif id_scaling==1:
            vec_offset=vec_offset_L2_Pois1_S

    elif id_FEM==1:
        if id_equ==0:
            vec_L2_u = vec_L2_u_Pois1
            vec_L2_du = vec_L2_du_Pois1
            vec_offset_L2_Poisx_M1=vec_offset_L2_Pois1_M1
            vec_offset_L2_Poisx_M2=vec_offset_L2_Pois1_M2
            
        elif id_equ==1:
            vec_L2_u = vec_L2_u_Pois2
            vec_L2_du = vec_L2_du_Pois2    
            vec_offset_L2_Poisx_M1=vec_offset_L2_Pois2_M1
            vec_offset_L2_Poisx_M2=vec_offset_L2_Pois2_M2        
        elif id_equ==2:
            vec_L2_u = vec_L2_u_Pois3
            vec_L2_du = vec_L2_du_Pois3
            vec_offset_L2_Poisx_M1=vec_offset_L2_Pois3_M1
            vec_offset_L2_Poisx_M2=vec_offset_L2_Pois3_M2        
        elif id_equ==3:
            vec_L2_u = vec_L2_u_Pois4
            vec_L2_du = vec_L2_du_Pois4  
            vec_offset_L2_Poisx_M1=vec_offset_L2_Pois4_M1
            vec_offset_L2_Poisx_M2=vec_offset_L2_Pois4_M2        
        elif id_equ==4:
            vec_L2_u = vec_L2_u_Pois5
            vec_L2_du = vec_L2_du_Pois5  
            vec_offset_L2_Poisx_M1=vec_offset_L2_Pois5_M1
            vec_offset_L2_Poisx_M2=vec_offset_L2_Pois5_M2        
            
        vec_offset_L2_Poisx_MM_normal[0] = vec_offset_L2_Poisx_M1[0] * vec_L2_u[id_c] *0.2
        vec_offset_L2_Poisx_MM_normal[1] = vec_offset_L2_Poisx_M2[1] * vec_L2_u[id_c]*1
        vec_offset_L2_Poisx_MM_normal[2] = vec_offset_L2_Poisx_M1[2] * vec_L2_du[id_c]*1
        
        if id_scaling==0:
            vec_offset=vec_offset_L2_Poisx_MM_normal 
        elif (id_scaling==1):
            vec_offset=vec_offset_L2_Poisx_M1
        elif (id_scaling==2):
            vec_offset=vec_offset_L2_Poisx_M2
            
    if (id_scaling==0):
        tria_coeff_x = 1e0
        tria_coeff_y= vec_L2_u[4]*1e0
        
    elif (id_scaling==2):
        tria_start = 5
        tria_end = tria_start+3
        tria_coeff_x = 5e2
        tria_coeff_y = 1e-1             
            
elif id_vec_id_equ == 2:
    vec_id_equ = vec_infl_d_equ
    vec_legend = ['1e0', '1e1', '1e2', '1e3', '1e4']
    vec_legend_coeff=['$c$','$c$','$c$','$c$','$c$']
    
    if id_equ == 2:
        vec_coeff = ['1em2', '1em1', '1e0', '1e1', '1e2']
        vec_legend = ['1e-2', '1e-1', '1e0', '1e1', '1e2']    
    
    if id_FEM == 0:
        if id_equ == 0:
            vec_offset = vec_offset_bench_Pois_2_diff_d_1_plus_sincx_SM
        elif id_equ == 1:
            vec_offset = vec_offset_bench_Pois_2_diff_d_1_plus_cx_SM
        elif id_equ == 2:
            vec_offset = vec_offset_bench_Pois_2_diff_d_c_SM
            
    elif id_FEM == 1:
        if id_equ == 0:
            vec_offset = vec_offset_bench_Pois_2_diff_d_1_plus_sincx_MM
        elif id_equ == 1:
            vec_offset = vec_offset_bench_Pois_2_diff_d_1_plus_cx_MM
        elif id_equ == 2:
            vec_offset = vec_offset_bench_Pois_2_diff_d_c_MM
            
elif id_vec_id_equ == 3:
    
    vec_id_equ = vec_infl_r_equ
    vec_legend = ['1e0', '1e1', '1e2', '1e3', '1e4']
    vec_legend_coeff=['$c$','$c$','$c$','$c$','$c$']
    
    if id_equ == 1:
        vec_coeff = ['1em2', '1em1', '1e0', '1e1', '1e2', '1e3', '1e4']
        vec_legend = ['1e-2', '1e-1', '1e0', '1e1', '1e2', '1e3', '1e4']    
    
    if id_FEM == 0:
        if id_equ == 0:
            vec_offset = vec_offset_bench_Pois_2_Helm_d_1_plus_sinx_r_c_SM
        elif id_equ == 1:
            id_loc_legend = 2
            vec_offset = vec_offset_bench_Pois_2_Helm_d_1_r_c_SM            
        elif id_equ == 2:
            vec_offset = vec_offset_bench_Pois_2_Helm_d_1_r_1_plus_cx_SM
    elif id_FEM == 1:
        if id_equ == 0:
            vec_offset = vec_offset_bench_Pois_2_Helm_d_1_plus_sinx_r_c_MM  
        elif id_equ ==1:
            vec_offset = vec_offset_bench_Pois_2_Helm_d_1_r_c_MM            
        elif id_equ ==2:
            vec_offset = vec_offset_bench_Pois_2_Helm_d_1_r_1_plus_cx_MM
            
elif id_vec_id_equ == 4:
    print("dealing with sensitivity")
    vec_coeff = ['1', '2', '3', '4', '5']
    if id_FEM == 0:
        if id_sensitivity < 3:
            vec_offset = vec_offset_bench_Pois_SM
        elif id_sensitivity == 3:
            vec_offset = [7e-17, 1e-16, 5e-16]
    elif id_FEM == 1:
        vec_legend = vec_legend_bench_mm
        if id_sensitivity < 3:
            vec_offset = vec_offset_bench_Pois_MM
        elif id_sensitivity == 3:
            vec_offset = [3e-17, 1e-17, 2e-16]

print("equation: %s" %(vec_id_equ[id_equ]))
if id_vec_id_equ == 0 & id_equ == 2:    
    print("is_complex: %d" %(is_complex))

print("\n")
print("FEM method: %s" %(vec_FEM[id_FEM]))
print("Scaling scheme: %s" %(vec_scaling_scheme[id_scaling]))  
print("\n")

dof_ref = [0 for i in range(n_err_refer)]
err_round_off_approx = [0 for i in range(n_err_refer)]

for i in range(n_err_refer):
    dof_ref[i] = 2**i      
text_offset_x = dof_ref[2]

for id_var in range(id_var_start, id_var_end):
        
    print ('var: '+vec_var[id_var])   
 
    print('slope: ', vec_slope[id_FEM])
    
    if id_vec_id_equ == 1 or id_vec_id_equ == 2 or id_vec_id_equ == 3 or id_vec_id_equ == 4:
        print('legend: ', vec_legend)

    
    if id_vec_id_equ == 1:
        if id_FEM==0:
            print('L2_u: ',vec_L2_u[id_c])  
            print('scale-independent offset: ', vec_offset_L2_Pois1_S[id_var])     
        elif id_FEM==1:
            if id_var == 0 or id_var==2:
                print('L2_u: ',vec_L2_u[id_c])  
                print('scale-independent offset: ', vec_offset_L2_Poisx_M1[id_var])     
            elif id_var==1:
                print('L2_du: ',vec_L2_du[id_c])  
                print('scale-independent offset: ', vec_offset_L2_Poisx_M2[id_var])     
    
    print('offset in use: ', vec_offset[id_var])  
    
    if id_vec_id_equ == 0: 
        
        dataraw = [line.strip().split() for line in open('data_error_'+vec_id_equ[id_equ]+'_'+vec_FEM[id_FEM]+'_'+vec_var[id_var]+'.txt','r')]
        n_refine = len(dataraw)-1
        
        data_n_dofs_multi_p = np.zeros((n_refine,n_column))     #[[0 for i in range(n_refine)] for y in range(n_column)] 
        data_n_cond_for_degree = np.zeros((n_refine,n_column))
        
        if(id_FEM == 0):
            
            tria_start = 8
            tria_end = tria_start+3
            
            tria_coeff_x = 1e1
            tria_coeff_y = 1e-1
            
            text_slop_bottom_coeff_x = 0.5
            text_slop_bottom_coeff_y = 0.1
            text_slop_right_coeff_x = 1.2
            text_slop_right_coeff_y = 0.12
            
            text_offset_loc_x = 3
            text_offset_loc_y = 0.2 

            if id_equ==0:
                if id_var==0:
                    text_offset_loc_y = 0.15
                    id_loc_legend=4
                    text_offset_loc_x = 1
                elif id_var==1:
                    text_offset_loc_x=2
                elif id_var==2:
                    text_offset_loc_x=1
            elif id_equ==1:
                if id_var==0:
                    text_offset_loc_x=3
                elif id_var==1:
                    text_offset_loc_x=3
                elif id_var==2:
                    text_offset_loc_x=1
            elif id_equ==2:
                if id_var==2:
                    text_offset_loc_x=1         
            
            
            text_offset_x = dof_ref[text_offset_loc_x]          
                            
            for i in range(n_refine):
                for j in range(n_column):
                    data_n_dofs_multi_p[i][j] = (2**(i+1)*(j+1)+1)*2**is_complex          # we use '+1' here because both i and j start from 0   
                    
                    if id_equ==0:
                        data_n_cond_for_degree[i][j]=data_n_dofs_multi_p[i][j]**2
                    elif id_equ==1:
                        data_n_cond_for_degree[i][j]=5*data_n_dofs_multi_p[i][j]**2
                    elif id_equ==2:
                        data_n_cond_for_degree[i][j]=2*data_n_dofs_multi_p[i][j]**2
                        

        elif(id_FEM == 1):
            
            tria_start = 6
            tria_end = tria_start+4
            
            tria_coeff_x = 5e1
            tria_coeff_y = 0.1
            
            text_slop_bottom_coeff_x = 0.4
            text_slop_bottom_coeff_y = 0.09
            text_slop_right_coeff_x = 1.1          # 1.3   0.8
            text_slop_right_coeff_y = 0.2   
            
            
            text_offset_loc_x = 1
            text_offset_loc_y = 0.1

            if id_equ==0:
                if id_var==0:
                    text_offset_loc_x=2
                    text_offset_loc_y = 0.15
                elif id_var==1:
                    tria_coeff_y = 0.01
            elif id_equ==1:
                if id_var==0:
                    text_offset_loc_x=2
            
            
            text_offset_x = dof_ref[text_offset_loc_x]            
            
            yaxis_low_bound = 1e-20
            xaxis_up_bound = 1e8            
            
            for i in range(n_refine):
                for j in range(n_column):            
                    if (id_var == 0):
                        data_n_dofs_multi_p[i][j] = (2**(i+1)*(j+1))*2**is_complex
                    elif (id_var == 1 or id_var == 2):
                        data_n_dofs_multi_p[i][j] = (2**(i+1)*(j+1)+1)*2**is_complex
                        
                    data_n_cond_for_degree[i][j]=((2**(i+1)*(j+1))*2+1)*2**is_complex                            
               
            
    elif id_vec_id_equ == 1:
        
        tria_coeff_x = 1e1
        tria_coeff_y = 1
        
        text_slop_bottom_coeff_y = 0.1
        text_slop_right_coeff_y = 0.1
        text_offset_loc_y = 0.1       # 0.1 0.05         
        
        dataraw = [line.strip().split() for line in open('data_error_'+vec_id_equ[id_equ]+'_'+vec_FEM[id_FEM]+'_'+vec_scaling_scheme[id_scaling]+'_'+vec_var[id_var]+'.txt','r')]
        n_refine = len(dataraw)-1              # omit the title
        data_n_dofs_one_p =  np.zeros((n_refine,1))     # [0]* n_refine # [[0] for y in range(n_refine)]   
            
        if (id_FEM == 0):
            text_offset_loc_x = 4
            text_offset_x =dof_ref[text_offset_loc_x-id_var]
            for i in range(n_refine):
                data_n_dofs_one_p[i] = 2**(i+1)*2+1                 # P2 elements are used for the standard FEM
                
            if (id_scaling==0): 
                
#                for i in range(n_err_refer):
#                    err_round_off_approx[i] = err_round_off_approx[i]*vec_L2[4]                
                
                tria_coeff_x = 1e3             
                tria_coeff_y = 1e-1
                
                if (id_equ==0):
                    id_loc_legend = 3
                    tria_coeff_x = 1e3
                    tria_coeff_y = 1e0
                    text_slop_bottom_coeff_y = 0.05
#                    text_offset_loc_y = 0.005
                elif id_equ==1:
                    tria_coeff_x = 1e1           
                    tria_coeff_y = 5e-2
#                    text_slop_bottom_coeff_y = 0.05
                elif (id_equ==2):
                    tria_coeff_x = 1e0
                    tria_coeff_y=1e-1
                    text_slop_right_coeff_y=0.1
                    text_slop_bottom_coeff_y = 0.05
                elif (id_equ==3):
                    id_loc_legend = 3

                    tria_coeff_x=1e2
                    tria_coeff_y=1e0
                elif (id_equ==4):
                    id_loc_legend = 2
                    text_slop_bottom_coeff_y = 0.05

                    
            elif (id_scaling == 1):
                tria_coeff_y = 1e-1
                if (id_equ==0):
                    id_loc_legend = 4
                    tria_coeff_x = 1e1
                    yaxis_up_bound = 1e4
                    text_slop_bottom_coeff_y = 0.025        
                    text_slop_right_coeff_y = 0.08  
                    text_offset_loc_x=3
                    text_offset_loc_y = 0.05
                elif id_equ==2:
                    text_slop_bottom_coeff_y = 0.05
                elif id_equ==4:
                    text_slop_bottom_coeff_y = 0.05
                    
        elif (id_FEM == 1):
            
            text_slop_bottom_coeff_y=0.1                       
            text_slop_right_coeff_y=0.2              
            
            if (id_var==0):
                for i in range(n_refine):
                    data_n_dofs_one_p[i] = 2**(i+1)*4
            elif (id_var == 1 or id_var == 2):
                for i in range(n_refine):
                    data_n_dofs_one_p[i] = 2**(i+1)*4+1  
    
                    
            if (id_scaling==0):
                
#                for i in range(n_err_refer):
#                    err_round_off_approx[i] = err_round_off_approx[i]*vec_L2[4]                
                
                tria_coeff_y= 0.01
                text_slop_bottom_coeff_y = 0.05
                id_loc_legend = 1
                
                if (id_equ==0):
                    tria_coeff_x = 1e0  
                    if id_var==0 or id_var==1:
                        tria_coeff_y = 1e2            #  uxx 1e-2
                elif(id_equ==1):
                    if(id_var==2):
                        tria_coeff_y=1e-6
                elif (id_equ==2):                       
                    text_slop_right_coeff_y=0.2
                elif (id_equ==3):
                    tria_coeff_y=1e-1
                elif(id_equ==4):
                    tria_coeff_y= 0.5
                    text_slop_bottom_coeff_y=0.1
                    text_slop_right_coeff_y=0.25
                    yaxis_up_bound=1e-4
                    if(id_var==0):
                        tria_coeff_x=1e3
                    
                    
            elif (id_scaling==1):
                
                text_slop_bottom_coeff_y=0.05  
                id_loc_legend = 1
                
                if (id_equ==0):
                    tria_coeff_x = 1e1               
                    tria_coeff_y = 1e-2 
                    text_slop_right_coeff_x = 1.1
                    if (id_var==2):
                        id_loc_legend=4
                elif(id_equ==1):
                    tria_coeff_x=1e2
                    tria_coeff_y=1e-1
                elif(id_equ==2):  
                    tria_coeff_x = 3e1
                    tria_coeff_y = 1e-2
                elif(id_equ==3):
                    tria_coeff_y = 1e-2
                elif(id_equ==4):
                    yaxis_up_bound=1e-4
                    tria_coeff_y = 1e-2
                    text_slop_right_coeff_y = 0.3
                    text_slop_bottom_coeff_y=0.1
                    
            elif (id_scaling==2):
                
                text_slop_bottom_coeff_y=0.05 
                id_loc_legend = 1
                
                if (id_equ==0):
                    tria_coeff_x = 1e2
                elif(id_equ==1):
                    text_slop_bottom_coeff_y=0.05
                    tria_coeff_x = 1e3
                    text_offset_loc_y=0.05
                    tria_coeff_y=1
                    if (id_var==2):
                        tria_coeff_y=1e-4
                elif(id_equ==2):
                    tria_coeff_x=2e2
                    tria_coeff_y=1e-1
                    if(id_var==1):
                        text_offset_loc_y=0.05
                elif(id_equ==3):
                    tria_coeff_x = 1e2
                    tria_coeff_y = 1e-1            
                elif(id_equ==4):                        
                    tria_coeff_x = 1e2
                    tria_coeff_y = 1e-1                    
                    text_slop_bottom_coeff_y=0.1
                    yaxis_up_bound=1e-4
                    
    elif id_vec_id_equ == 2 or id_vec_id_equ == 3:
        
        dataraw = [line.strip().split() for line in open('data_error_'+vec_id_equ[id_equ]+'_'+vec_FEM[id_FEM]+'_coeff_'+vec_coeff[0]+'.txt','r')]
        
        n_refine = len(dataraw)-1
        data_n_dofs_one_p = np.zeros((n_refine,1))
        
        if (id_FEM == 0):
            text_offset_loc_x = 6
            text_offset_x =dof_ref[text_offset_loc_x-id_var]   
            
            for i in range(n_refine):
                data_n_dofs_one_p[i] = 2**(i+1)*2+1                 # P2 elements are used for the standard FEM
        elif id_FEM == 1:
            
            text_offset_loc_x = 10
            text_offset_x =dof_ref[text_offset_loc_x-id_var] 
            
            if (id_var==0):
                for i in range(n_refine):
                    data_n_dofs_one_p[i] = 2**(i+1)*4
            elif (id_var == 1 or id_var == 2):
                for i in range(n_refine):
                    data_n_dofs_one_p[i] = 2**(i+1)*4+1  
                    
    elif id_vec_id_equ == 4:
        dataraw = [line.strip().split() for line in open('data_error_'+vec_id_equ[id_equ]+'_'+vec_sensitivity[id_sensitivity]+'_'+vec_FEM[id_FEM]+'_deg_'+vec_coeff[0]+'.txt','r')]
        n_refine = len(dataraw)-1        
        
        data_n_dofs_multi_p = np.zeros((n_refine,n_column))
        
        if id_FEM == 0:
            for i in range(n_refine):
                for j in range(n_column):
                    data_n_dofs_multi_p[i][j] = (2**(i+1)*(j+1)+1)*2**is_complex     
        elif id_FEM == 1:
            for i in range(n_refine):
                for j in range(n_column):            
                    if (id_var == 0):
                        data_n_dofs_multi_p[i][j] = (2**(i+1)*(j+1))*2**is_complex
                    elif (id_var == 1 or id_var == 2):
                        data_n_dofs_multi_p[i][j] = (2**(i+1)*(j+1)+1)*2**is_complex
                                    
        
    for i in range(n_err_refer):
        err_round_off_approx[i] = vec_offset[id_var]* dof_ref[i]**vec_slope[id_FEM] 
        
        
    data_error = np.zeros((n_refine,n_column))
    
    
    if id_vec_id_equ == 0 or id_vec_id_equ ==1:
        
        for i in range(n_refine):
    #        print (i, end =': ')
            for p in range(n_column):
    #            print(dataraw[i+1][p+1], end = ',')
                if(id_vec_id_equ==0):
                    data_error[i][p]=float(dataraw[i+1][p+1])
                elif(id_vec_id_equ == 1):
                    data_error[i][n_column-1-p]=float(dataraw[i+1][p+1])
                    
    elif id_vec_id_equ == 2 or id_vec_id_equ == 3:
        for p in range(id_case_start, id_case_end):
            
            dataraw = [line.strip().split() for line in open('data_error_'+vec_id_equ[id_equ]+'_'+vec_FEM[id_FEM]+'_coeff_'+vec_coeff[p]+'.txt','r')]
            
            for i in range(n_refine):
                data_error[i][p-id_case_start]=dataraw[i+1][id_var]
                
    elif id_vec_id_equ == 4:
        for p in range(id_case_start, id_case_end):
            dataraw = [line.strip().split() for line in open('data_error_'+vec_id_equ[id_equ]+'_'+vec_sensitivity[id_sensitivity]+'_'+vec_FEM[id_FEM]+'_deg_'+vec_coeff[p]+'.txt','r')]
            
            for i in range(n_refine):
                data_error[i][p-id_case_start]=dataraw[i+1][id_var]
                    
    for i in range(n_refine):
        for p in range(n_column):
            if(data_error[i][p]==0):
                data_error[i][p] = np.nan
                
                
#        dataraw = [line.strip().split() for line in open('data_error_'+vec_id_equ[id_equ]+'_'+vec_FEM[id_FEM]+'_'+vec_scaling_scheme[id_scaling]+'_'+vec_var[id_var]+'.txt','r')]
#                    
#        print ('data_n_dofs_one_p: ',end='')
#        print (data_n_dofs_one_p)            
                    
                    
    text_offset_y = vec_offset[id_var]* text_offset_x**vec_slope[id_FEM]*text_offset_loc_y       
            
          
#    print ('tria_coeff_x: ',tria_coeff_x)  
#    print ('tria_coeff_y: ',tria_coeff_y) 
#    print ('text_slop_bottom_coeff_y: ',text_slop_bottom_coeff_y)                      
#    print ('text_slop_right_coeff_y: ',text_slop_right_coeff_y) 
#    
#    print ('text_offset_loc_x: ', text_offset_loc_x) 
#    print('text_offset_x: ', text_offset_x)  
#    print ('text_offset_loc_y: ', text_offset_loc_y) 
#    print('text_offset_y: ', text_offset_y)    
        
    tria_p_1 = [dof_ref[tria_start],err_round_off_approx[tria_start]]
    tria_p_2 = [dof_ref[tria_end],err_round_off_approx[tria_start]]
    tria_p_3 = [dof_ref[tria_end],err_round_off_approx[tria_end]]
    
    tria_x = [tria_p_1[0],tria_p_2[0],tria_p_3[0],tria_p_1[0]]#*
    tria_x = [i*tria_coeff_x for i in tria_x]
    tria_y = [tria_p_1[1],tria_p_2[1],tria_p_3[1],tria_p_1[1]]#/tria_coeff_y
    tria_y = [i*tria_coeff_y for i in tria_y]
    
    text_bottom_x = (tria_x[0]+tria_x[1])/2*text_slop_bottom_coeff_x
    text_bottom_y = tria_y[0]*text_slop_bottom_coeff_y
    text_right_x = tria_x[1]*text_slop_right_coeff_x
    text_right_y = (tria_y[1]+tria_y[2])/2*text_slop_right_coeff_y
    
          
#    print('\n')
#    print('xaxis_up_bound: ', xaxis_up_bound)
#    print('yaxis_low_bound: ', yaxis_low_bound)    
#    print('yaxis_up_bound: ', yaxis_up_bound)
    
    
    f=plt.figure(figsize=(5,4))
    axes= f.add_axes([0,0,1,1])
    
            
    if id_vec_id_equ == 0 or id_vec_id_equ == 4:
        
        plt.loglog(dof_ref, err_round_off_approx,'--',color='orange',linewidth=1.0)
        plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15)             
        
        for p in range(n_column):           # n_column      
            plt.loglog(data_n_dofs_multi_p[:,p], data_error[:,p],'k'+vec_marker[p], markerfacecolor='none',label=vec_legend[p],linewidth=1.0)    
            
    elif id_vec_id_equ == 1 or id_vec_id_equ ==2 or id_vec_id_equ == 3:
        
        if (id_FEM==0):
            if id_scaling==1:
                plt.loglog(dof_ref, err_round_off_approx,'k--',linewidth=1.0)
                plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15)       
                
        elif (id_FEM==1):
            if id_scaling == 1:
                if (id_var == 0 or id_var == 2):
                    plt.loglog(dof_ref, err_round_off_approx,'k--',linewidth=1.0)
                    plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15) 
            elif id_scaling == 2:
                if (id_var == 0 or id_var == 1):
                    plt.loglog(dof_ref, err_round_off_approx,'k--',linewidth=1.0)
                    plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15)            
        
        for p in range(id_case_start, id_case_end):          # n_column
            plt.loglog(data_n_dofs_one_p, data_error[:,p-id_case_start],'k'+vec_marker[p-id_case_start]+'-', markerfacecolor='none',label=vec_legend_coeff[id_equ]+'='+str(vec_legend[p]),linewidth=1.0)    

    plt.plot(tria_x,tria_y,'k-',linewidth=1.0)
    plt.text(text_bottom_x,text_bottom_y,'1',fontsize = 15)
    plt.text(text_right_x,text_right_y,str(vec_slope[id_FEM]),fontsize = 15)
    
    if id_var == 0:    
        plt.legend()
        plt.legend(loc=id_loc_legend, prop={'size': fontsize_legend})
        
    plt.xlabel('Number of DoFs', fontsize=fontsize_label)             # Condition number
    plt.ylabel('Error', fontsize=fontsize_label)
    
    axes.set_xlim([1, xaxis_up_bound])
    axes.set_ylim([yaxis_low_bound, yaxis_up_bound])           
    
    if id_vec_id_equ == 0 or id_vec_id_equ == 2 or id_vec_id_equ == 3:
        if id_FEM == 0:
            plt.xticks([1e0, 1e2, 1e4, 1e6, 1e8])      
            plt.yticks([1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e0])
                        
        elif id_FEM == 1:
            plt.xticks([1e0, 1e2, 1e4, 1e6, 1e8])      
            plt.yticks([1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e0]) 
            
    elif id_vec_id_equ == 1:
        plt.xticks([1e0, 1e2, 1e4, 1e6, 1e8])    
            
        if (id_FEM==0):
            if(id_equ==0):
                if id_scaling==0:
                    plt.yticks([1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e0]) 
                elif id_scaling==1:
                    plt.yticks([1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e0, 1e4]) 
            elif (id_equ==1 or id_equ==3):
                plt.yticks([1e-16, 1e-12, 1e-8, 1e-4, 1e0]) 
            else:
                plt.yticks([1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e0]) 
        elif (id_FEM==1):
            if(id_equ==4):
                plt.yticks([1e-20, 1e-16, 1e-12, 1e-8, 1e-4])
            else:
                plt.yticks([1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e0]) 
        
            
    plt.tick_params(axis='both', which='major', labelsize=fontsize_tick)
    
    plt.show()
    
    if id_vec_id_equ == 0 or id_vec_id_equ == 2 or id_vec_id_equ == 3:
        f.savefig("py_"+vec_id_equ[id_equ]+"_"+vec_FEM[id_FEM]+"_"+vec_var[id_var]+".pdf", bbox_inches='tight')    
    elif id_vec_id_equ == 1:
        f.savefig("py_"+vec_id_equ[id_equ]+"_"+vec_FEM[id_FEM]+"_"+vec_scaling_scheme[id_scaling]+"_"+vec_var[id_var]+".pdf", bbox_inches='tight')




#plot_error()

#for id_FEM in range(0,2):
#    for id_equ in range(0,3):
#        plot_error()   