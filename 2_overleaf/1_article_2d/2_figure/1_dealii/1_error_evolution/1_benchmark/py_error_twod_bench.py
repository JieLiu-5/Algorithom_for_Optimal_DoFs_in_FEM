#!/usr/bin/env python2
# -*- coding: utf-8 -*-
    
import matplotlib.pyplot as plt
import numpy as np

import math

import sys
sys.path.append('../../../../../3_1d_or_2d/')
import settings_fig as set

vec_dimension=['1d','2d']
id_dimension=0

vec_problem_type=['0_pois','1_diff','2_helm']
id_problem_type=2

id_equ = 10
is_distorted = 0                        # only for Poisson problems


id_equ_segmentation_pois = 3
id_equ_segmentation_diff = 6
id_equ_segmentation_helm = 5


if id_dimension==0:
    if id_problem_type==0:    
        vec_equ=['0p0_u_x_m_0p5_square_uniform']
        vec_equ.append('0p1_u_x_m_0p5_square_randomly_distorted')
        vec_equ.append('0p2_u_x_m_0p5_square_regularly_1_distorted')
        vec_equ.append('0p3_u_x_m_0p5_square_regularly_2_distorted')            # id_equ=3
        vec_equ.append('1p0_u_exp_m_x_m_0p5_square_uniform')
        vec_equ.append('1p1_u_exp_m_x_m_0p5_square_randomly_distorted')
        vec_equ.append('1p2_u_exp_m_x_m_0p5_square_regularly_1_distorted')
        vec_equ.append('1p3_u_exp_m_x_m_0p5_square_regularly_2_distorted')
        vec_equ.append('2_u_2o3_t_x_cubic_p_1o2_t_x_square')     
        vec_equ.append('3_u_x_pow_4_over12_m_x_over4')
    elif id_problem_type==1:
        vec_equ=['0p00_u_x_m_0p5_square_d_1px']
        vec_equ.append('0p01_u_x_m_0p5_square_d_1px_distorted')
        vec_equ.append('0p1_u_x_m_0p5_square_d_1p1em3x')
        vec_equ.append('0p2_u_x_m_0p5_square_d_1p1e3x')
        vec_equ.append('0p3_u_x_m_0p5_square_d_1pxpxsquare')
        vec_equ.append('0p4_u_x_m_0p5_square_d_exp_m_x_m_0p5_square')
        vec_equ.append('0p5_u_x_m_0p5_square_d_0p5_p_cosx_square')         # id_equ=5
        vec_equ.append('1p00_u_exp_m_x_m_0p5_square_d_1px')
        vec_equ.append('1p01_u_exp_m_x_m_0p5_square_d_1px_distorted')
        vec_equ.append('1p1_u_exp_m_x_m_0p5_square_d_1pxpxsquare')
        vec_equ.append('1p2_u_exp_m_x_m_0p5_square_d_exp_m_x_m_0p5_square')
        vec_equ.append('1p3_u_exp_m_x_m_0p5_square_d_0p5_p_cosx_square')
    elif id_problem_type==2:
        vec_equ=['0p0_u_x_m_0p5_square_d_1px_r_1']
        vec_equ.append('0p1_u_x_m_0p5_square_d_1px_r_1px')
        vec_equ.append('0p2_u_x_m_0p5_square_d_1px_r_exp_m_x_m_0p5_square')
        vec_equ.append('0p3_u_x_m_0p5_square_d_exp_m_x_m_0p5_square_r_1')         # id_equ=3
        vec_equ.append('0p4_u_x_m_0p5_square_d_exp_m_x_m_0p5_square_r_1px')
        vec_equ.append('0p5_u_x_m_0p5_square_d_exp_m_x_m_0p5_square_r_exp_m_x_m_0p5_square')
        vec_equ.append('2p0_u_exp_m_x_m_0p5_square_d_1px_r_1')
        vec_equ.append('2p1_u_exp_m_x_m_0p5_square_d_1px_r_1px')    
        vec_equ.append('2p2_u_exp_m_x_m_0p5_square_d_1px_r_exp_m_x_m_0p5_square')
        vec_equ.append('2p3_u_exp_m_x_m_0p5_square_d_exp_m_x_m_0p5_square_r_1')
        vec_equ.append('2p4_u_exp_m_x_m_0p5_square_d_exp_m_x_m_0p5_square_r_1px')
        vec_equ.append('2p5_u_exp_m_x_m_0p5_square_d_exp_m_x_m_0p5_square_r_exp_m_x_m_0p5_square')       

elif id_dimension==1:
    if id_problem_type==0:
        vec_equ=['0_u_1_over_c_x_m_0p5_square_p_y_same']
        vec_equ.append('1_u_exp_m_x_m_0p5_square')
        vec_equ.append('2_u_x_m_0p5_square_plus_1','3_u_x_m_0p5_square')
    elif id_problem_type==1:
        vec_equ=['0_u_x_m_0p5_square_p_y_same_d_1_0p5_0p5_2']
        vec_equ.append('1_u_x_m_0p5_square_p_y_same_d_1pxpy_0_0_1pxpy')
        vec_equ.append('2_u_x_m_0p5_square_p_y_same_d_1px_square_py_square_0_0_same')
        vec_equ.append('3_u_x_m_0p5_square_p_y_same_d_exp_0_0_exp')
        vec_equ.append('4_u_x_m_0p5_square_p_y_same_d_1pxpy_0p5_0p5_1pxpy')
        vec_equ.append('5_u_x_m_0p5_square_p_y_same_d_1pxpy_xy_xy_1pxpy')
        vec_equ.append('6_u_x_m_0p5_square_p_y_same_d_cos_0p5_0p5_sin')


if id_dimension==0:
    vec_FEM=['0_sm']    
elif id_dimension==1:
    vec_FEM=['0_sm','1_mm_RT', '2_mm_BDM']
id_FEM=0

vec_deg = [1,2,3,4,5,6,7]
n_degree_data = 5
if (id_FEM==1 or id_FEM==2) and id_problem_type==0 and id_equ==0:
    n_degree_data = 7

vec_var = ['solu','grad','2ndd']
vec_offset = [2e-17, 2e-16, 2e-15]
vec_slope = [1, 1, 1]

vec_offset_auxiliary = [1e-18, 5e-18, 1e-16]
vec_slope_auxiliary = [2.0,2.0,2.0]
    
if id_problem_type == 0:
    if id_equ > id_equ_segmentation_pois:
        vec_offset_auxiliary = [2e-17, 1e-16, 5e-16]
elif id_problem_type == 1:
    if id_equ > id_equ_segmentation_diff:
        vec_offset_auxiliary = [2e-17, 1e-16, 5e-16]
elif id_problem_type == 2:

    if id_equ > id_equ_segmentation_helm:
        vec_offset_auxiliary = [2e-17, 1e-16, 5e-16]
            

vec_xaxis = ['ndofs', 'ncond']
vec_xlabel = ['Number of DoFs', 'Condition number']
id_xaxis=0


n_err_refer = 28
dof_ref = np.zeros(n_err_refer)
err_round_off_approx = np.zeros(n_err_refer)
err_round_off_approx_auxiliary = np.zeros(n_err_refer)
for i in range(n_err_refer):
    dof_ref[i] = 2**i
    
color_round_off_line_main = 'orange'
color_round_off_line_auxiliary = 'tomato'

if id_problem_type==0 and is_distorted == 0:
    color_round_off_line_main = 'tomato' 
    

vec_legend = ['$Q_1$','$Q_2$','$Q_3$','$Q_4$','$Q_5$','$Q_6$','$Q_7$']
id_loc_legend = 4
fontsize_legend = set.fontsize_legend

if id_FEM==1 or id_FEM==2:
    fontsize_legend = 16

vec_marker=['o','d','^','s','*']

is_solu_grad_together=1

id_error_format = 1

if id_dimension==0:
    if id_problem_type==0:
        if id_equ==0:
            if id_FEM==0:
                vec_offset = [1e-18, 5e-18, 1e-16]
                vec_slope=[2.0,2.0,2.0]   
        elif id_equ==1:
            if id_FEM==0:
                vec_offset = [1e-18, 5e-18, 2e-16]
                vec_slope=[2.0,2.0,2.0]      
        elif id_equ==2:
            if id_FEM==0:
                vec_offset = [1e-18, 5e-18, 2e-16]
                vec_slope=[2.0,2.0,2.0]   
        elif id_equ==3:
            if id_FEM==0:
                vec_offset = [2e-18, 1e-17, 1e-16]
                vec_slope=[1.5,1.5,2.0]                  
        elif id_equ==4:                                     # starting exp
            if id_FEM==0:
                vec_offset = [2e-17, 1e-16, 5e-16]
                vec_slope=[2.0,2.0,2.0]
        elif id_equ==5:
            if id_FEM==0:
                vec_offset = [2e-17, 1e-16, 5e-15]
                vec_slope=[2.0,2.0,2.0]    
        elif id_equ==6:
            if id_FEM==0:
                vec_offset = [2e-17, 1e-16, 2e-15]
                vec_slope=[2.0,2.0,2.0]    
        elif id_equ==7:
            if id_FEM==0:
                vec_offset = [2e-18, 5e-18, 1e-15]
                vec_slope=[2.0,2.0,2.0]
        elif id_equ==8:
            if id_FEM==0:
                vec_offset = [5e-18, 2e-17, 2e-16]
                vec_slope=[2.0,2.0,2.0]
        elif id_equ==9:
            if id_FEM==0:
                vec_offset = [5e-18, 1e-17, 1e-16]
                vec_slope=[2.0,2.0,2.0]                  
    elif id_problem_type==1:
        if id_equ==0:
            if id_FEM==0:
                vec_offset = [2e-18, 1e-17, 1e-16]
                vec_slope=[1.5,1.5,2.0]
        elif id_equ==1:
            if id_FEM==0:
                vec_offset = [2e-18, 1e-17, 1e-16]
                vec_slope=[1.5,1.5,2.0]                
        elif id_equ==2:
            if id_FEM==0:
                vec_offset = [2e-18, 1e-17, 1e-16]
                vec_slope=[1.5,1.5,2.0]    
        elif id_equ==3:
            if id_FEM==0:
                vec_offset = [2e-18, 1e-17, 1e-16]
                vec_slope=[1.5,1.5,2.0]                  
        if id_equ==4:
            if id_FEM==0:
                vec_offset = [2e-18, 1e-17, 1e-16]
                vec_slope=[1.5,1.5,2.0]                        
        elif id_equ==5:
            if id_FEM==0:
                vec_offset = [2e-18, 1e-17, 1e-16]
                vec_slope=[1.5,1.5,2.0]
        elif id_equ==6:
            if id_FEM==0:
                vec_offset = [2e-18, 1e-17, 1e-16]
                vec_slope=[1.5,1.5,2.0]                
        elif id_equ==id_equ_segmentation_diff+1:
            if id_FEM==0:
                vec_offset = [2e-18, 2e-18, 5e-16]
                vec_slope=[2.0,2.0,2.0]
        elif id_equ==id_equ_segmentation_diff+2:
            if id_FEM==0:
                vec_offset = [2e-17, 2e-16, 5e-16]
                vec_slope=[1.5,1.5,2.0]                
        elif id_equ==8:
            if id_FEM==id_equ_segmentation_diff+3:
                vec_offset = [5e-19, 2e-18, 5e-16]
                vec_slope=[2.0,2.0,2.0]                  
        elif id_equ==id_equ_segmentation_diff+4:
            if id_FEM==0:
                vec_offset = [2e-17, 1e-16, 5e-16]
                vec_slope=[1.5,1.5,2.0]              
        elif id_equ==id_equ_segmentation_diff+5:
            if id_FEM==0:
                vec_offset = [2e-17, 1e-16, 5e-16]
                vec_slope=[1.5,1.5,2.0]                
    elif id_problem_type==2:
        if id_equ==0:
            if id_FEM==0:
                vec_offset = [5e-19, 1e-18, 1e-16]
                vec_slope=[2.0,2.0,2.0]
        elif id_equ==1:
            if id_FEM==0:
                vec_offset = [2e-18, 1e-17, 1e-16]
                vec_slope=[1.5,1.5,2.0]                
        elif id_equ==2:
            if id_FEM==0:
                vec_offset = [2e-18, 1e-17, 1e-16]
                vec_slope=[1.5,1.5,2.0]                 
        elif id_equ==3:
            if id_FEM==0:
                vec_offset = [5e-19, 1e-18, 1e-16]
                vec_slope=[2.0,2.0,2.0]
        elif id_equ==4:
            if id_FEM==0:
                vec_offset = [2e-18, 1e-17, 1e-16]
                vec_slope=[1.5,1.5,2.0]
        elif id_equ==id_equ_segmentation_helm:
            if id_FEM==0:
                vec_offset = [2e-18, 1e-17, 1e-16]
                vec_slope=[1.5,1.5,2.0]
        elif id_equ==id_equ_segmentation_helm+1:
            if id_FEM==0:
                vec_offset = [2e-18, 2e-18, 5e-16]
                vec_slope=[2.0,2.0,2.0]                  
        elif id_equ==id_equ_segmentation_helm+2:
            if id_FEM==0:
                vec_offset = [2e-18, 2e-18, 5e-16]
                vec_slope=[2.0,2.0,2.0]                    
        elif id_equ==id_equ_segmentation_helm+3:
            if id_FEM==0:
                vec_offset = [2e-17, 1e-16, 5e-16]
                vec_slope=[1.5,1.5,2.0]                  
        elif id_equ==id_equ_segmentation_helm+4:
            if id_FEM==0:
                vec_offset = [1e-17, 5e-17, 5e-16]
                vec_slope=[2.0,2.0,2.0]
        elif id_equ==id_equ_segmentation_helm+5:
            if id_FEM==0:
                vec_offset = [2e-17, 1e-16, 5e-16]
                vec_slope=[1.5,1.5,2.0]           
        elif id_equ==id_equ_segmentation_helm+6:
            if id_FEM==0:
                vec_offset = [2e-17, 1e-16, 5e-16]
                vec_slope=[1.5,1.5,2.0]
   
               
elif id_dimension==1:
    if id_problem_type==0:
        if id_equ==0:
            if id_FEM==0:
                vec_offset = [1e-17, 2e-17, 2e-16]
            elif id_FEM==1:
                vec_offset=[5e-17, 1e-16, 5e-16]
                vec_slope=[0.25,0.5,0.5]
                vec_legend = [r"$RT_1/Q_1^{\rm disc}$", r"$RT_2/Q_2^{\rm disc}$",r"$RT_3/Q_3^{\rm disc}$",r"$RT_4/Q_4^{\rm disc}$",r"$RT_5/Q_5^{\rm disc}$",r"$RT_6/Q_6^{\rm disc}$",r"$RT_7/Q_7^{\rm disc}$"] 
            elif id_FEM==2:
                vec_offset=[1e-16, 1e-16, 2e-16]
                vec_slope=[0.25,0.5,0.5]    
                vec_legend = [r"$BDM_1/P_0^{\rm disc}$",r"$BDM_2/P_1^{\rm disc}$",r"$BDM_3/P_2^{\rm disc}$",r"$BDM_4/P_3^{\rm disc}$",r"$BDM_5/P_4^{\rm disc}$",r"$BDM_6/P_5^{\rm disc}$",r"$BDM_7/P_6^{\rm disc}$"]
        if id_equ==1:
            if id_FEM==1:
                vec_deg = [5,6,7,8,9]
                vec_legend = [r"$RT_5/Q_5^{\rm disc}$", r"$RT_6/Q_6^{\rm disc}$",r"$RT_7/Q_7^{\rm disc}$",r"$RT_8/Q_8^{\rm disc}$",r"$RT_9/Q_9^{\rm disc}$"] 
                vec_slope=[0.5,0.5,1]
                if id_xaxis==0:
                    vec_offset=[1e-16, 1e-15, 5e-15]
                elif id_xaxis==1:
                    vec_offset=[2e-15, 2e-15, 1e-14]
            elif id_FEM==2:
                vec_deg = [3,4,5,6,7]
                vec_legend = [r"$BDM_3/P_2^{\rm disc}$",r"$BDM_4/P_3^{\rm disc}$",r"$BDM_5/P_4^{\rm disc}$",r"$BDM_6/P_5^{\rm disc}$",r"$BDM_7/P_6^{\rm disc}$"]
                vec_offset=[3.6e-16, 4e-15, 1e-14]
                vec_slope=[0.5,0.5,1]
    elif id_problem_type==1:
        if id_equ==0:
            if id_FEM==0:
                vec_offset = [1e-17, 2e-17, 2e-16]
                vec_slope=[1.0,1.0,1.0]
            elif id_FEM==1:
                vec_offset=[1e-16, 2e-16, 1e-15]
                vec_slope=[0.25,0.5,0.5]
                vec_legend = [r"$RT_1/Q_1^{\rm disc}$", r"$RT_2/Q_2^{\rm disc}$",r"$RT_3/Q_3^{\rm disc}$",r"$RT_4/Q_4^{\rm disc}$",r"$RT_5/Q_5^{\rm disc}$"]             
            elif id_FEM==2:
                vec_offset=[1e-16, 1e-16, 2e-16]
                vec_slope=[0.25,0.5,0.5]    
                vec_legend = [r"$BDM_1/P_0^{\rm disc}$",r"$BDM_2/P_1^{\rm disc}$",r"$BDM_3/P_2^{\rm disc}$",r"$BDM_4/P_3^{\rm disc}$",r"$BDM_5/P_4^{\rm disc}$",r"$BDM_6/P_5^{\rm disc}$",r"$BDM_7/P_6^{\rm disc}$"]            
        elif id_equ==1 or id_equ==2 or id_equ==3 or id_equ==4:
            if id_FEM==0:
                vec_offset = [1e-17, 2e-17, 2e-16]
                vec_slope=[0.75,0.75,1.0]       
        elif id_equ==5:
            if id_FEM==0:
                vec_offset = [1e-17, 2e-17, 2e-16]
                vec_slope=[0.75,0.75,1.0]
            elif id_FEM==1:
                vec_offset=[5e-17, 2e-16, 1e-15]
                vec_slope=[0.25,0.5,0.5]
                vec_legend = [r"$RT_1/Q_1^{\rm disc}$", r"$RT_2/Q_2^{\rm disc}$",r"$RT_3/Q_3^{\rm disc}$",r"$RT_4/Q_4^{\rm disc}$",r"$RT_5/Q_5^{\rm disc}$"] 
        elif id_equ==6:
            if id_FEM==0:
                n_degree_data=7
                vec_offset = [1e-17, 2e-16, 2e-15]
                vec_slope=[1.0,1.0,1.5]    
            
id_deg_start=0
n_deg_plot=n_degree_data
            
        
if id_xaxis==1:
    matrix_data_xaxis = [line.strip().split() for line in open(vec_problem_type[id_problem_type]+'/'+vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_ncond.txt','r')]
  
if id_error_format==1:
    matrix_data_error = [line.strip().split() for line in open(vec_problem_type[id_problem_type]+'/'+vec_dimension[id_dimension]+'/'+vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_deg_1.txt','r')]
elif id_error_format==2:
    matrix_data_error = [line.strip().split() for line in open(vec_problem_type[id_problem_type]+'/'+vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_'+vec_var[0]+'.txt','r')]
    
n_refine_initial = len(matrix_data_error)-1
n_refine = 0

data_error = np.zeros((n_refine_initial,n_degree_data)) 
data_xaxis = np.zeros((n_refine_initial,n_degree_data))
    
id_var_start=0
id_var_end=3

if((id_FEM==1 or id_FEM==2) and is_solu_grad_together==1 and id_xaxis==1):
    id_var_start=1

for id_var in range(id_var_start,id_var_end):
    
    print('id_var:', id_var)
    print ('var: '+vec_var[id_var])

    if id_error_format==1:
        for id_degree in range(0,n_degree_data):
            matrix_data_error = [line.strip().split() for line in open(vec_problem_type[id_problem_type]+'/'+vec_dimension[id_dimension]+'/'+vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_deg_'+str(vec_deg[id_degree])+'.txt','r')]
            n_refine = len(matrix_data_error)-1
            for i in range(n_refine):
                data_error[i][id_degree]=matrix_data_error[i+1][id_var]
                if(data_error[i][id_degree]==0):
                    data_error[i][id_degree] = np.nan               
                if(data_xaxis[i][id_degree]==0):
                    data_xaxis[i][id_degree] = np.nan                
    elif id_error_format==2:
        if ((id_FEM==1 or id_FEM==2) and id_xaxis==1 and is_solu_grad_together==1 and id_var==1):         
            print('solu and grad together')
            matrix_data_error = [line.strip().split() for line in open(vec_problem_type[id_problem_type]+'/'+vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_solu_plus_grad.txt','r')]    
        else:
            matrix_data_error = [line.strip().split() for line in open(vec_problem_type[id_problem_type]+'/'+vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_'+vec_var[id_var]+'.txt','r')]
        for i in range(n_refine):
            for p in range(n_degree_data):
                data_error[i][p]=matrix_data_error[i+1][p]
        
      
    for i in range(n_refine_initial):
        for p in range(n_degree_data):
            if id_xaxis == 0:
                if id_dimension==0:
                    if id_FEM==0:
                        data_xaxis[i][p]=2**(i+1)*vec_deg[p]+1                    
                elif id_dimension==1:
                    if id_FEM==0:
                        data_xaxis[i][p]=(2**(i+1)*vec_deg[p]+1)**2
                    elif id_FEM==1:
                        if id_var==0:
                            data_xaxis[i][p]=(2**(i+1)*(vec_deg[p]+1))**2
                        elif id_var==1 or id_var==2:
                            dofs_per_quad = 2*vec_deg[p]*(vec_deg[p]+1)
                            data_xaxis[i][p]=2**(i+1)*(vec_deg[p]+1)*((2**(i+1)-1)*2+4) + 4**(i+1)*dofs_per_quad
                    elif id_FEM==2:
                        if id_var==0:
                            dofs_per_cell_dgq = (vec_deg[p]+1)*(vec_deg[p]+2)/2
                            data_xaxis[i][p]=4**(i+1)*dofs_per_cell_dgq
                        elif id_var==1 or id_var==2:
                            dofs_per_quad = vec_deg[p]*(vec_deg[p]-1)
                            data_xaxis[i][p]=2**(i+1)*(vec_deg[p]+1)*((2**(i+1)-1)*2+4) + 4**(i+1)*dofs_per_quad
                        
            elif id_xaxis == 1:
                data_xaxis[i][p]=matrix_data_xaxis[i+1][p]
    
            if(data_error[i][p]==0):
                data_error[i][p] = np.nan               
            if(data_xaxis[i][p]==0):
                data_xaxis[i][p] = np.nan
                
    for i in range(n_err_refer):
        err_round_off_approx[i] = vec_offset[id_var]* dof_ref[i]**vec_slope[id_var]             
        err_round_off_approx_auxiliary[i] = vec_offset_auxiliary[id_var]* dof_ref[i]**vec_slope_auxiliary[id_var]
        
    tria_start = 8
    tria_end = tria_start+4
    
    tria_coeff_x = 1e1
    tria_coeff_y = 1e-2
    
    text_offset_coeff_y=0.1

    text_slope_bottom_coeff_x = 0.9
    text_slope_bottom_coeff_y = 0.08
    text_slope_right_coeff_x = 1.2
    text_slope_right_coeff_y = 0.18
    
    
    if id_var == 0: 
        text_offset_coeff_y=0.1
    elif id_var == 1:
        text_slope_right_coeff_y = 0.16 
    elif id_var == 2:
        text_slope_right_coeff_y = 0.18
        
    
    tria_p_1 = [dof_ref[tria_start],err_round_off_approx[tria_start]]
    tria_p_2 = [dof_ref[tria_end],err_round_off_approx[tria_start]]
    tria_p_3 = [dof_ref[tria_end],err_round_off_approx[tria_end]]
    
    tria_x = [tria_p_1[0],tria_p_2[0],tria_p_3[0],tria_p_1[0]]#*
    tria_x = [i*tria_coeff_x for i in tria_x]
    tria_y = [tria_p_1[1],tria_p_2[1],tria_p_3[1],tria_p_1[1]]#/tria_coeff_y
    tria_y = [i*tria_coeff_y for i in tria_y]
    
    text_slope_bottom_x = 10**((math.log10(tria_x[0])+math.log10(tria_x[1]))/2)*text_slope_bottom_coeff_x
    text_slope_bottom_y = tria_y[0]*text_slope_bottom_coeff_y
    text_slope_right_x = tria_x[1]*text_slope_right_coeff_x
    text_slope_right_y = 0.5*10**((math.log10(tria_y[1])+math.log10(tria_y[2]))/2)


    f=plt.figure(figsize=(5,4))
    axes= f.add_axes([0,0,1,1])
    
    xaxis_up_bound=1e8
    yaxis_low_bound=1e-20
    yaxis_up_bound=1e0
    
    for p in range(id_deg_start,id_deg_start+n_deg_plot):
        plt.loglog(data_xaxis[:,p], data_error[:,p],'k'+set.vec_marker[p], markerfacecolor='none',label=vec_legend[p],linewidth=1.0)
    
    
    plt.plot(tria_x,tria_y,'k-',linewidth=1.0)
    plt.text(text_slope_bottom_x,text_slope_bottom_y,'1',fontsize = 15)
    plt.text(text_slope_right_x,text_slope_right_y,str(vec_slope[id_var]),fontsize = 15)         
        
    if (id_var == 0) or ((id_FEM==1 or id_FEM==2) and id_xaxis==1 and is_solu_grad_together==1 and id_var==1):    
        plt.legend()
        plt.legend(loc=id_loc_legend, prop={'size': fontsize_legend})
  
        
    plt.loglog(dof_ref, err_round_off_approx,'k--',linewidth=1.0, color=color_round_off_line_main)
    if (id_problem_type==0 and is_distorted==1) or id_problem_type==1 or id_problem_type==2:              #  and (id_equ<4 or id_equ>5)
        plt.loglog(dof_ref, err_round_off_approx_auxiliary,'k--',linewidth=1.0, color=color_round_off_line_auxiliary)
    
    
    text_offset_x = dof_ref[1]
    
    if (id_FEM==1 and id_xaxis==1 and is_solu_grad_together==1 and id_var==2):
        text_offset_x = dof_ref[2]
    
    text_offset_y = vec_offset[id_var]* text_offset_x**vec_slope[id_var]*text_offset_coeff_y       
    
    plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15, color="black")             
    
    plt.tick_params(axis='both', which='major', labelsize=set.fontsize_tick)
    
    plt.xlabel(vec_xlabel[id_xaxis], fontsize=set.fontsize_label)
    plt.ylabel('Error', fontsize=set.fontsize_label)
    
    print('yaxis_low_bound: ', yaxis_low_bound)    
    print('yaxis_up_bound: ', yaxis_up_bound)
    
    axes.set_xlim([1, xaxis_up_bound])
    axes.set_ylim([yaxis_low_bound, yaxis_up_bound])
    
    plt.yticks([1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e0])
        
    plt.show()
    
    if ((id_FEM==1 or id_FEM==2) and id_xaxis==1 and is_solu_grad_together==1 and id_var==1):
        f.savefig(vec_problem_type[id_problem_type]+'/'+vec_dimension[id_dimension]+'/'+vec_equ[id_equ]+"/"+vec_FEM[id_FEM]+"/py_TwoD_Pois_"+vec_equ[id_equ]+"_"+vec_FEM[id_FEM]+"_UMF_"+vec_xaxis[id_xaxis]+"_solu_plus_grad.pdf", bbox_inches='tight')
    else:
        f.savefig(vec_problem_type[id_problem_type]+'/'+vec_dimension[id_dimension]+'/'+vec_equ[id_equ]+"/"+vec_FEM[id_FEM]+"/py_error_"+vec_var[id_var]+".pdf", bbox_inches='tight')


