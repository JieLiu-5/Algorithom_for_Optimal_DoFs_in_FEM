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
sys.path.append('../../../../1_1d/2_figure')
import settings as set


vec_equ=['1_u_e_m_c_x_m_0p5_square','2_u_1_over_c_x_m_0p5_square','3_u_cx_plus_100_square','4_u_x_square_m_c_square']
id_equ=2

vec_FEM=['1_sm','2_mm_PPdisc','3_mm_RT']
id_FEM=1

degree = 4
if id_FEM==2:
    degree=2


vec_xaxis=['ndofs','ncond']
id_xaxis=0



n_param = 5
vec_param = ['1em2', '1em1', '1e0', '1e1', '1e2']
id_param=0

vec_var = ['solu','grad','2ndd']

vec_offset = [2e-20, 1e-19, 2e-20]
vec_slope = [1, 1.5, 2.5]

vec_label_xaxis=['Number of DoFs','Condition number']

n_degree = 5
n_err_refer = 28

vec_legend = ['1e-2', '1e-1', '1e0', '1e1', '1e2']

id_loc_legend = 1
vec_marker=['o','d','^','s','*']  
legend_size=16
fontsize_label = 18

xaxis_up_bound=1e6
yaxis_low_bound=1e-20
yaxis_up_bound=1e0
vec_ytick = [1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e0]

dof_ref = np.zeros(n_err_refer)
err_round_off_approx = np.zeros(n_err_refer)

for i in range(n_err_refer):
    dof_ref[i] = 2**i
    
    
if id_FEM==0:
    vec_slope=[1,1,1]    
    
if id_xaxis==0:    
    if id_FEM==1:
        vec_slope=[1.0,1.0,1.0]
    elif id_FEM==2:
        vec_slope=[0.5,0.5,1]
elif id_xaxis==1:
    if id_FEM==1:
        vec_slope=[1,1.5,2.5]
    elif id_FEM==2:
        vec_slope=[1,1,2]            
  
if id_equ==0:
    if id_FEM==0:
        vec_offset = [2e-17, 2e-16, 2e-15]
elif id_equ==1:
    if id_FEM==0:
        yaxis_low_bound=1e-20
        yaxis_up_bound=1e-4
        vec_ytick = [1e-20, 1e-16, 1e-12, 1e-8, 1e-4]         
    elif id_FEM==1 or id_FEM==2:
        yaxis_low_bound=1e-24
        yaxis_up_bound=1e-4
        vec_ytick = [1e-24, 1e-20, 1e-16, 1e-12, 1e-8, 1e-4]        
elif id_equ==2:
    if id_FEM==1:
        vec_offset=[1e-14, 1e-12, 1e-15]
    elif id_FEM==2:
        vec_offset=[4e-13, 2e-12, 5e-12]
elif id_equ==3:
    vec_offset=[6e-13, 3e-12, 4e-12]
        
if id_xaxis==1:    
    matrix_data_xaxis = [line.strip().split() for line in open(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_ncond.txt','r')]
    
matrix_data_error = [line.strip().split() for line in open(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_deg_'+str(degree)+'_c_'+vec_param[0]+'.txt','r')]
n_refine = len(matrix_data_error)-1
    
data_error = np.zeros((n_refine,n_param))
data_xaxis = np.zeros(n_refine)    
        

id_var_start = 0
id_var_end = 3

for id_var in range(id_var_start,id_var_end):
        
    print ('var: '+vec_var[id_var])      
          
    f=plt.figure(figsize=(5,4))
    axes= f.add_axes([0,0,1,1])
    
    for id_param in range(0,n_param):      
        matrix_data_error = [line.strip().split() for line in open(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/data_error_deg_'+str(degree)+'_c_'+vec_param[id_param]+'.txt','r')]
        for i in range(n_refine):
            data_error[i][id_param]=matrix_data_error[i+1][id_var]
        
    for i in range(n_refine):
        if id_xaxis==0:
            if id_FEM==0:
               data_xaxis[i]=2**(i+1)*degree+1
            elif id_FEM==1:
                if id_var==0:
                    data_xaxis[i]=2**(i+1)*degree
                elif id_var==1 or id_var==2:
                    data_xaxis[i]=2**(i+1)*degree+1
        elif id_xaxis==1:
            data_xaxis[i]=matrix_data_xaxis[i+1][0]
                        
    for i in range(n_refine):
        if(data_xaxis[i]==0):
            data_xaxis[i] = np.nan        
        for p in range(n_param):
            if(data_error[i][p]==0):
                data_error[i][p] = np.nan


                
    for i in range(n_err_refer):
        err_round_off_approx[i] = vec_offset[id_var]* dof_ref[i]**vec_slope[id_var]                 
    
    tria_start = 8
    tria_end = tria_start+3
    
    tria_coeff_x = 1e1
    tria_coeff_y = 1e-2
    
    
    text_offset_x = dof_ref[1]
    
    text_offset_coeff_y=0.1
    
    
    text_slope_bottom_coeff_x = 0.4
    text_slope_bottom_coeff_y = 0.07
    text_slope_right_coeff_x = 1.2
    text_slope_right_coeff_y = 0.2
    
    

    if id_xaxis==0:
        if id_var == 1:
            text_slope_right_coeff_y = 0.2
        elif id_var == 2:
            text_offset_x = dof_ref[1]
            text_slope_right_coeff_y = 0.13        
    elif id_xaxis==1:
        if id_var == 1:
            text_slope_right_coeff_y = 0.08
        elif id_var == 2:
            text_offset_x = dof_ref[2]
            text_slope_right_coeff_y = 0.02            
        
    
    tria_p_1 = [dof_ref[tria_start],err_round_off_approx[tria_start]]
    tria_p_2 = [dof_ref[tria_end],err_round_off_approx[tria_start]]
    tria_p_3 = [dof_ref[tria_end],err_round_off_approx[tria_end]]
    
    tria_x = [tria_p_1[0],tria_p_2[0],tria_p_3[0],tria_p_1[0]]#*
    tria_x = [i*tria_coeff_x for i in tria_x]
    tria_y = [tria_p_1[1],tria_p_2[1],tria_p_3[1],tria_p_1[1]]#/tria_coeff_y
    tria_y = [i*tria_coeff_y for i in tria_y]
    
    text_slope_bottom_x = 10**((math.log10(tria_x[0])+math.log10(tria_x[1]))/2)
    text_slope_bottom_y = tria_y[0]*text_slope_bottom_coeff_y
    text_slope_right_x = tria_x[1]*text_slope_right_coeff_x
    text_slope_right_y = 0.5*10**((math.log10(tria_y[1])+math.log10(tria_y[2]))/2)
        
    
    text_offset_y = vec_offset[id_var]* text_offset_x**vec_slope[id_var]*text_offset_coeff_y       
    
    for id_param in range(n_param):
        plt.loglog(data_xaxis, data_error[:,id_param],'k'+vec_marker[id_param]+'-', markerfacecolor='none',label='c='+vec_legend[id_param],linewidth=1.0)
    
    
   
    if id_var == 0:    
        plt.legend()
        plt.legend(loc=id_loc_legend, prop={'size': set.fontsize_legend})
    
        if degree==4 or (degree==2 and id_FEM==2):
            
            plt.loglog(dof_ref, err_round_off_approx,'k--',linewidth=1.0)
            plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15)             
                            
            plt.plot(tria_x,tria_y,'k-',linewidth=1.0)
            plt.text(text_slope_bottom_x,text_slope_bottom_y,'1',fontsize = 15)
            plt.text(text_slope_right_x,text_slope_right_y,str(vec_slope[id_var]),fontsize = 15)         
                     
    elif id_var == 1 or id_var == 2:
        plt.loglog(dof_ref, err_round_off_approx,'k--',linewidth=1.0)
        plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15)             
               
        plt.plot(tria_x,tria_y,'k-',linewidth=1.0)
        plt.text(text_slope_bottom_x,text_slope_bottom_y,'1',fontsize = 15)
        plt.text(text_slope_right_x,text_slope_right_y,str(vec_slope[id_var]),fontsize = 15)         
         

        
    
    plt.tick_params(axis='both', which='major', labelsize=set.fontsize_tick)
    
    plt.xlabel(vec_label_xaxis[id_xaxis], fontsize=set.fontsize_label)
    plt.ylabel('Error', fontsize=set.fontsize_label)
    
    
    print('yaxis_low_bound: ', yaxis_low_bound)    
    print('yaxis_up_bound: ', yaxis_up_bound)   
    
    axes.set_xlim([1, xaxis_up_bound])
    axes.set_ylim([yaxis_low_bound, yaxis_up_bound]) 
    
    plt.yticks(vec_ytick) 
        
    plt.show()
    
    f.savefig(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+"/py_error_oned_"+vec_equ[id_equ]+"_"+vec_FEM[id_FEM]+'_deg_'+str(degree)+'_'+vec_xaxis[id_xaxis]+"_"+vec_var[id_var]+".pdf", bbox_inches='tight')
    

