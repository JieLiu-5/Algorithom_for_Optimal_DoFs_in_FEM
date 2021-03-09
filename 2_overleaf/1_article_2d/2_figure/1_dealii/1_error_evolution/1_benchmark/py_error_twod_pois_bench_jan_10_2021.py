#!/usr/bin/env python2
# -*- coding: utf-8 -*-
    
import matplotlib.pyplot as plt
import numpy as np

import math

import sys
sys.path.append('../../../../../3_1d_or_2d/')
import settings_fig as set

vec_equ=['0_u_1_over_c_x_m_0p5_square_p_y_same','1_u_exp_m_x_m_0p5_square','2_u_x_m_0p5_square_plus_1','3_u_x_m_0p5_square']
id_equ=0

vec_FEM=['1_sm','2_mm_PPdisc','3_mm_RT', '4_mm_BDM']
id_FEM=3

vec_xaxis = ['ndofs', 'ncond']
vec_xlabel = ['Number of DoFs', 'Condition number']
id_xaxis=0

vec_var = ['solu','grad','2ndd']
vec_offset = [2e-17, 2e-16, 2e-15]
vec_slope = [1, 1, 1]
        

vec_deg = [1,2,3,4,5,6,7]
n_degree_data = 7

id_deg_start=0
n_deg_plot=n_degree_data

n_err_refer = 28
dof_ref = np.zeros(n_err_refer)
err_round_off_approx = np.zeros(n_err_refer)
for i in range(n_err_refer):
    dof_ref[i] = 2**i 
    

dofs_per_cell_dgq = [1,3,6,10,15,21,28]    

vec_legend = ['$P_1$','$P_2$','$P_3$','$P_4$','$P_5$']
id_loc_legend = 1
fontsize_legend = set.fontsize_legend

if id_FEM==2 or id_FEM==3:
    fontsize_legend = 16

vec_marker=['o','d','^','s','*']

is_solu_grad_together=1

id_error_format = 1


if id_equ==0:
    if id_FEM==0:
        vec_offset = [1e-17, 2e-17, 2e-16]
    elif id_FEM==2:
        vec_offset=[1e-17, 1e-16, 2e-16]
        vec_slope=[0.5,0.5,0.5]
        vec_legend = [r"$RT_1/P_1^{\rm disc}$", r"$RT_2/P_2^{\rm disc}$",r"$RT_3/P_3^{\rm disc}$",r"$RT_4/P_4^{\rm disc}$",r"$RT_5/P_5^{\rm disc}$",r"$RT_6/P_6^{\rm disc}$",r"$RT_7/P_7^{\rm disc}$"] 
    elif id_FEM==3:
        vec_offset=[1e-17, 1e-16, 2e-16]
        vec_slope=[0.5,0.5,0.5]    
        vec_legend = [r"$BDM_1/Q_0^{\rm disc}$",r"$BDM_2/Q_1^{\rm disc}$",r"$BDM_3/Q_2^{\rm disc}$",r"$BDM_4/Q_3^{\rm disc}$",r"$BDM_5/Q_4^{\rm disc}$",r"$BDM_6/Q_5^{\rm disc}$",r"$BDM_7/Q_6^{\rm disc}$"]
if id_equ==1:
    if id_FEM==1:
        vec_deg = [5,6,7,8,9]
        vec_legend = [r"$P_5/P_4^{\rm disc}$", r"$P_6/P_5^{\rm disc}$",r"$P_7/P_6^{\rm disc}$",r"$P_8/P_7^{\rm disc}$",r"$P_9/P_8^{\rm disc}$"]
        if id_xaxis==0:
            vec_offset=[1e-17, 1e-16, 5e-16]
            vec_slope=[0.5,0.75,1.25]        
        elif id_xaxis==1:
            vec_offset = [5e-20, 1e-19, 1e-20]
            vec_slope = [1, 1.5, 2.5]
    elif id_FEM==2:
        vec_deg = [5,6,7,8,9]
        vec_legend = [r"$RT_5/P_5^{\rm disc}$", r"$RT_6/P_6^{\rm disc}$",r"$RT_7/P_7^{\rm disc}$",r"$RT_8/P_8^{\rm disc}$",r"$RT_9/P_9^{\rm disc}$"] 
        vec_slope=[0.5,0.5,1]
        if id_xaxis==0:
            vec_offset=[1e-16, 1e-15, 5e-15]
        elif id_xaxis==1:
            vec_offset=[2e-15, 2e-15, 1e-14]
    elif id_FEM==3:
        vec_deg = [3,4,5,6,7]
        vec_legend = [r"$BDM_3/Q_2^{\rm disc}$",r"$BDM_4/Q_3^{\rm disc}$",r"$BDM_5/Q_4^{\rm disc}$",r"$BDM_6/Q_5^{\rm disc}$",r"$BDM_7/Q_6^{\rm disc}$"]
        vec_offset=[3.6e-16, 4e-15, 1e-14]
        vec_slope=[0.5,0.5,1]
elif id_equ==2:
    if id_FEM==1:
        vec_legend = [r"$P_1/P_0^{\rm disc}$",r"$P_2/P_1^{\rm disc}$",r"$P_3/P_2^{\rm disc}$",r"$P_4/P_3^{\rm disc}$",r"$P_5/P_4^{\rm disc}$"]
        if id_xaxis==0:
            vec_offset=[1e-17, 1e-16, 5e-16]
            vec_slope=[0.5,0.75,1.25]        
        elif id_xaxis==1:
            vec_offset = [1e-19, 1e-18, 1e-19]
            vec_slope = [1, 1.5, 2.5]
elif id_equ==3:
    vec_legend = [r"$P_1/P_0^{\rm disc}$",r"$P_2/P_1^{\rm disc}$",r"$P_3/P_2^{\rm disc}$",r"$P_4/P_3^{\rm disc}$",r"$P_5/P_4^{\rm disc}$"]
    if id_FEM==0:
        vec_offset = [2e-18, 1e-17, 1e-16]
    elif id_FEM==1:
        if id_xaxis==0:
            vec_offset = [1e-18, 1e-17, 5e-17]
            vec_slope=[0.5,0.75,1.25]
        elif id_xaxis==1:
            vec_offset = [1e-19, 1e-19, 1e-20]
            vec_slope = [1, 1.5, 2.5]

        
if id_xaxis==1:
    matrix_data_xaxis = [line.strip().split() for line in open(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_ncond.txt','r')]

n_refine = 0
matrix_data_error=np.zeros((n_refine,n_degree_data))
  
if id_error_format==1:
    matrix_data_error = [line.strip().split() for line in open(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_deg_1.txt','r')]
elif id_error_format==2:
    matrix_data_error = [line.strip().split() for line in open(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_'+vec_var[0]+'.txt','r')]
    
n_refine = len(matrix_data_error)-1

data_error = np.zeros((n_refine,n_degree_data)) 
data_xaxis = np.zeros((n_refine,n_degree_data))
    
id_var_start=0
id_var_end=3

if((id_FEM==1 or id_FEM==2) and is_solu_grad_together==1 and id_xaxis==1):
    id_var_start=1

for id_var in range(id_var_start,id_var_end):
    
    print('id_var:', id_var)
    print ('var: '+vec_var[id_var])

    if id_error_format==1:
        for id_degree in range(0,n_degree_data):
            matrix_data_error = [line.strip().split() for line in open(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_deg_'+str(vec_deg[id_degree])+'.txt','r')]
            for i in range(n_refine):
                data_error[i][id_degree]=matrix_data_error[i+1][id_var]
                if(data_error[i][id_degree]==0):
                    data_error[i][id_degree] = np.nan               
                if(data_xaxis[i][id_degree]==0):
                    data_xaxis[i][id_degree] = np.nan                
    elif id_error_format==2:
        if ((id_FEM==1 or id_FEM==2) and id_xaxis==1 and is_solu_grad_together==1 and id_var==1):         
            print('solu and grad together')
            matrix_data_error = [line.strip().split() for line in open(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_solu_plus_grad.txt','r')]    
        else:
            matrix_data_error = [line.strip().split() for line in open(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_'+vec_var[id_var]+'.txt','r')]
        for i in range(n_refine):
            for p in range(n_degree_data):
                data_error[i][p]=matrix_data_error[i+1][p]
        
      
    for i in range(n_refine):
        for p in range(n_degree_data):
            if id_xaxis == 0:
                if id_FEM==0:
                    data_xaxis[i][p]=(2**(i+1)*vec_deg[p]+1)**2
                elif id_FEM==1:
                    if id_var==0:
                        data_xaxis[i][p]=(2**(i+1)*vec_deg[p])**2
    #            data_ndofs[i][p]=data_ndofs[i][p]-((2**(i+1)*vec_deg[p]+1)*4-4)    # ignoring ndofs on the boundary
                    elif id_var==1 or id_var==2:
                        data_xaxis[i][p]=(2**(i+1)*vec_deg[p]+1)**2*2
                elif id_FEM==2:
                    if id_var==0:
                        data_xaxis[i][p]=(2**(i+1)*(vec_deg[p]+1))**2
                    elif id_var==1 or id_var==2:
                        dofs_per_quad = 2*vec_deg[p]*(vec_deg[p]+1)
                        data_xaxis[i][p]=2**(i+1)*(vec_deg[p]+1)*((2**(i+1)-1)*2+4) + 4**(i+1)*dofs_per_quad
                elif id_FEM==3:
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

    text_offset_x = dof_ref[1]
    
    if (id_FEM==1 and id_xaxis==1 and is_solu_grad_together==1 and id_var==2):
        print('dksk')
        text_offset_x = dof_ref[2]
    
    text_offset_y = vec_offset[id_var]* text_offset_x**vec_slope[id_var]*text_offset_coeff_y       
               
          
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
  
        
    plt.loglog(dof_ref, err_round_off_approx,'k--',linewidth=1.0)
    
    plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15)             
    
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
        f.savefig(vec_equ[id_equ]+"/"+vec_FEM[id_FEM]+"/py_TwoD_Pois_"+vec_equ[id_equ]+"_"+vec_FEM[id_FEM]+"_UMF_"+vec_xaxis[id_xaxis]+"_solu_plus_grad.pdf", bbox_inches='tight')
    else:
        f.savefig(vec_equ[id_equ]+"/"+vec_FEM[id_FEM]+"/py_TwoD_Pois_"+vec_equ[id_equ]+"_"+vec_FEM[id_FEM]+"_UMF_"+vec_xaxis[id_xaxis]+"_"+vec_var[id_var]+".pdf", bbox_inches='tight')    # cond  dofs


