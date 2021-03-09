#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 15:21:04 2018

@author: jliu
"""
import sys
sys.path.append('../')
import settings as set
        
import matplotlib.pyplot as plt
import numpy as np                       
        

vec_marker=['o', 'd', '^', 's', 'p']
vec_markeredgecolor = ['k', 'b', 'r', 'c', 'm', 'y', 'g', 'w']
vec_mfc = ["none","none","none","none","none"]

vec_type_equ = ["bench", "Pois", "diff", "Helm"]
vec_FEM = ['sm','mm']
vec_var = ['solu','grad','2ndd']
vec_legend = [r'$u$', r'$u_x$', r'$u_{xx}$']

vec_xlabel = [r"Equation type", r"$\||u\||_2$", r"$\||d\||_2$", r"$\||r\||_2$"]


vec_offset = np.zeros((3,1))
vec_slope = np.zeros((3,1))
vec_formular_base = ["" for i in range(3)]


trend_line_x_coord = np.logspace(-7, 5, 50)

data_L2_norm_all_cases = np.zeros((5,5))
data_L2_norm_all_cases_u = np.zeros((5,5))
data_L2_norm_all_cases_du = np.zeros((5,5))
data_offset_all_cases = np.zeros((5,5))          

        
f=plt.figure(figsize=(5,4))
axes= f.add_axes([0,0,1,1])
fontsize_label = 18
plt.tick_params(axis = 'both', which = 'major', labelsize = 16)
    
id_loc_legend = 1
fontsize_legend = 15

xaxis_low_bound=1e-7
xaxis_up_bound=1e5
yaxis_low_bound=1e-25
yaxis_up_bound=1e-10

coord_x_legend = 1.2e-3
vec_coord_y_text = np.zeros((3,1))
coeff_legend_adjust_y_axis = 1 

line_auxiliary_x = [1,1,1e-7]
line_auxiliary_y = [1e-25,2e-17,2e-17]


id_vec_type_equ = 1         # '0' for benchmark equations
                            # '1' for Poisson equations
                            # '2' for diffusion equations
                            # '3' for Helmholtz equations

id_FEM = 1

n_case = 3

if id_FEM == 1:
    vec_legend = [r'$u$', r'$v$', r'$v_{x}$']

id_var_start = 0
id_var_end = 3

print ('plotting')

if id_vec_type_equ == 0:
    print("benchmark equations")
    
    data_pure = np.zeros((3,n_case))
    
    coords_x_axis = [1, 2, 3]
    labels_x_axis = ["Poisson", "diffusion", "Helmholtz"]
        
    yaxis_low_bound=1e-20
    yaxis_up_bound=1e-14
    
    id_loc_legend=4
    
    vec_offset_independent = [2e-17, 5e-17, 1e-15]   
    if id_FEM==1:
        vec_offset_independent[1] = 2e-16
                
    
    vec_coord_y_text = [i*0.5 for i in vec_offset_independent]
    vec_coord_y_text[2]=vec_coord_y_text[2]*2.4
    
#    if id_FEM==0:
#        vec_coord_y_text[1]=vec_coord_y_text[1]*3
    vec_l2_d=[1, 1.53, 2.63]

    
    with open('offset_benchmark_'+vec_FEM[id_FEM]+'.txt','r') as file:
        data_raw = [line.split() for line in file]
        
    for id_var in range(3):
        for id_case in range(n_case):
            data_pure[id_var][id_case]=data_raw[id_var+1][id_case]

        coord_x_text=0.92
        
#        if id_FEM==0:
#            if id_var==1:
#                coord_x_text=1.1
        if id_FEM==1:
            
#            vec_offset_independent[1] = 2e-16*vec_l2_d[id_var]
            
            if id_var==0:
                coord_x_text=1.1
                vec_coord_y_text[0]=vec_coord_y_text[0]*2.4
                
        plt.semilogy(coords_x_axis[:n_case],data_pure[id_var],"o--",mfc=vec_mfc[id_var], markeredgecolor=vec_markeredgecolor[id_var],color=vec_markeredgecolor[id_var],label=vec_legend[id_var])
    
        axes.axhline(y=vec_offset_independent[id_var],color=vec_markeredgecolor[id_var])
        
        plt.text(coord_x_text, vec_coord_y_text[id_var], str(vec_offset_independent[id_var]), fontsize=fontsize_legend, color = vec_markeredgecolor[id_var])
            
    

elif id_vec_type_equ == 1:
    
    print('Poisson equations') 
    
    id_case_start = 0
    id_case_end = 5
    
    vec_slope = [1, 1, 1]
    
    id_loc_legend = 4
    
    coord_x_legend=1.5e-7
    
    if id_FEM==0:
        vec_offset = [2e-17, 5e-17, 5e-16]        
        vec_coord_y_text = [3e-18, 7e-17, 8e-16]
        
    elif id_FEM==1:        
        vec_xlabel[id_vec_type_equ] = '$\||u\||_2$ or $\||v\||_2$, with the variable on the\n horizontal axis $\||u\||_2$ for $u$ and $v$, and $\||v\||_2$ for $v_x$'
        vec_offset = [1e-18, 1e-16, 5e-16] 
        vec_coord_y_text = [1.4e-19,1.5e-17,8e-16]
        
    
    for id_var in range(id_var_start, id_var_end):
            
        print (vec_var[id_var] + "...", end = ' ')   
                
        line_auxiliary_y[1]=vec_offset[id_var]
        line_auxiliary_y[2]=line_auxiliary_y[1]
        
        for id_case in range(id_case_start, id_case_end):
            
            text_file_l2_norm = open('l2_norm_case_'+str(id_case+1)+'.txt','r')
            text_file_offset = open('offset_case_'+str(id_case+1)+"_"+vec_FEM[id_FEM]+'.txt','r')
            
            data_L2_norm_per_case = [line.strip().split() for line in text_file_l2_norm]
            data_offset_per_case = [line.strip().split() for line in text_file_offset]
            
            text_file_l2_norm.close()
            text_file_offset.close()
    
            for id_coeff in range(0,5): 
                data_L2_norm_all_cases_u[id_coeff, id_case]=data_L2_norm_per_case[id_coeff+1][0]
                data_offset_all_cases[id_coeff,id_case]=data_offset_per_case[id_coeff+1][id_var]
                if id_FEM == 1:    
                    data_L2_norm_all_cases_du[id_coeff,id_case]=data_L2_norm_per_case[id_coeff+1][1]
                
                if id_FEM==0:
                    data_L2_norm_all_cases = data_L2_norm_all_cases_u
                elif id_FEM==1:                    
                    if id_var==0 or id_var==1:
                        data_L2_norm_all_cases = data_L2_norm_all_cases_u
                    elif id_var==2:
                        data_L2_norm_all_cases = data_L2_norm_all_cases_du
                
            plt.loglog(data_L2_norm_all_cases[:,id_case], data_offset_all_cases[:,id_case], 'o', mfc=vec_mfc[id_case], markeredgecolor=vec_markeredgecolor[id_var], color=vec_markeredgecolor[id_case], linewidth=1.0)    
        
        
        trend_line_rounding = vec_offset[id_var] * trend_line_x_coord**vec_slope[id_var]     
        plt.loglog(trend_line_x_coord, trend_line_rounding, '-' + vec_markeredgecolor[id_var],label=vec_legend[id_var])
            
        plt.plot(line_auxiliary_x,line_auxiliary_y,'--',linewidth=1.0,color='grey')
        plt.text(coord_x_legend, vec_coord_y_text[id_var]*coeff_legend_adjust_y_axis, str(vec_offset[id_var]), fontsize=fontsize_legend, color = vec_markeredgecolor[id_var])

elif id_vec_type_equ == 2:

    print('diffusion equations')   
    
    vec_coord_y_text = [2e-18, 3e-19, 5e-20]       
    
    xaxis_low_bound=1e-3
    xaxis_up_bound=1e3
    yaxis_low_bound=1e-20
    yaxis_up_bound=1e-12   

    line_auxiliary_x[2]=1e-3
    
    with open("l2_norm_diff.txt", "r") as text_file_L2_norm:
        data_L2_norm = [per_row.split() for per_row in text_file_L2_norm.readlines()]
        
        data_id_d = [data_L2_norm[i][0] for i in range(1, len(data_L2_norm))]
        data_L2_norm_d = [float(data_L2_norm[i][1]) for i in range(1, len(data_L2_norm))]
        data_L2_norm_u = [float(data_L2_norm[i][2]) for i in range(1, len(data_L2_norm))]
        data_L2_norm_ux = [float(data_L2_norm[i][3]) for i in range(1, len(data_L2_norm))]
        data_L2_norm_dux = [float(data_L2_norm[i][4]) for i in range(1, len(data_L2_norm))]
        
        data_L2_norm_d_times_L2_norm_u = [data_L2_norm_d[i]*data_L2_norm_u[i] for i in range(0, len(data_L2_norm_d))]
            
        if id_FEM == 0:
            vec_offset = [2e-17, 5e-17, 1e-15]
            vec_slope = [0, 0, 0]
            
            text_file_offset = open("offset_diff_sm.txt", "r")
            data_L2_norm = data_L2_norm_d
            
            vec_coord_y_text = [i*0.4 for i in vec_offset]
            vec_coord_y_text[1]=vec_coord_y_text[1]*3.5
            vec_coord_y_text[2]=vec_coord_y_text[2]*3.5
            
        elif id_FEM == 1:
            vec_offset = [5e-17, 2e-16, 5e-16]
            vec_slope = [0, 1, 1]
            
            text_file_offset = open("offset_diff_mm.txt", "r")
            data_L2_norm = data_L2_norm_d
            
            vec_coord_y_text = [i*0.35 for i in vec_offset]
            id_loc_legend=2
                
            
        data_offset = [per_row.split() for per_row in text_file_offset.readlines()]
        text_file_offset.close()
        
        for id_var in range(id_var_start, id_var_end):
            data_offset_per_var = [float(data_offset[i][id_var]) for i in range(1, len(data_offset))]
            plt.loglog(data_L2_norm, data_offset_per_var, 'o', color = vec_markeredgecolor[id_var], mfc = "none")    
            plt.loglog(trend_line_x_coord, vec_offset[id_var] * trend_line_x_coord**vec_slope[id_var], '-' + vec_markeredgecolor[id_var], label=vec_legend[id_var])
            
            if ((id_FEM==0) or (id_FEM==1 and id_var==0)):
                plt.text(coord_x_legend, vec_coord_y_text[id_var], str(vec_offset[id_var]), fontsize=fontsize_legend, color = vec_markeredgecolor[id_var])
            else:
                vec_coord_y_text = [1.5e-19,7e-17,7e-16]
                line_auxiliary_y[1]=vec_offset[id_var]
                line_auxiliary_y[2]=line_auxiliary_y[1]
                plt.plot(line_auxiliary_x,line_auxiliary_y,'--',linewidth=1.0,color='grey')
                plt.text(coord_x_legend, vec_coord_y_text[id_var]*coeff_legend_adjust_y_axis, str(vec_offset[id_var]), fontsize=fontsize_legend, color = vec_markeredgecolor[id_var])



elif id_vec_type_equ == 3:
    
    print('Helmholtz equations') 
    
    xaxis_low_bound=1e-3
    xaxis_up_bound=1e3
    yaxis_low_bound=1e-20
    yaxis_up_bound=1e-15
    
    
    with open("l2_norm_Helm.txt", "r") as text_file_L2_norm:
        data_L2_norm = [per_row.split() for per_row in text_file_L2_norm.readlines()]
        
        data_id_d = [data_L2_norm[i][0] for i in range(1, len(data_L2_norm))]
        data_L2_norm_d = [float(data_L2_norm[i][1]) for i in range(1, len(data_L2_norm))]
        
        if id_FEM == 0:
            vec_offset = [2e-17, 5e-17, 2e-16]
            vec_slope = [0, 0, 0]
            
            text_file_offset = open("offset_Helm_sm.txt", "r")
            data_L2_norm = data_L2_norm_d
            
            vec_coord_y_text = [i*0.55 for i in vec_offset]
            vec_coord_y_text[1] = vec_coord_y_text[1]*2.3
            vec_coord_y_text[2] = vec_coord_y_text[2]*2.3
            
            id_loc_legend=4
            
        elif id_FEM == 1:
            vec_offset = [1e-19, 1e-16, 2e-16]
            vec_slope = [0, 0, 0]
            vec_formular_base = ["", "", ""]
            
            text_file_offset = open("offset_Helm_mm.txt", "r")
            data_L2_norm = data_L2_norm_d
            
            vec_coord_y_text = [i*0.52 for i in vec_offset]
            vec_coord_y_text[2] = vec_coord_y_text[2]*2.5
                
            id_loc_legend=5
            
        data_offset = [per_row.split() for per_row in text_file_offset.readlines()]
        text_file_offset.close()
        
        for id_var in range(id_var_start, id_var_end):
            data_offset_per_var = [float(data_offset[i][id_var]) for i in range(1, len(data_offset))]
            plt.loglog(data_L2_norm, data_offset_per_var, 'o', color = vec_markeredgecolor[id_var], mfc = "none")    
            plt.loglog(trend_line_x_coord, vec_offset[id_var] * trend_line_x_coord**vec_slope[id_var], '-' + vec_markeredgecolor[id_var], label=vec_legend[id_var])
            plt.text(coord_x_legend, vec_coord_y_text[id_var], str(vec_offset[id_var])+vec_formular_base[id_var], fontsize=fontsize_legend, color = vec_markeredgecolor[id_var])
            

if id_vec_type_equ == 1 or (id_vec_type_equ==2 and id_FEM==1):
    
    tria_start = 35
    
    tria_coeff_y = 1e-4

    coeff_text_slope_bottom_x = 0.25
    coeff_text_slope_bottom_y = 0.14
    coeff_text_slope_right_x = 1.3
    coeff_text_slope_right_y = 0.17
     
    if (id_vec_type_equ == 1):
        if (id_FEM==1):
            tria_start=40
            tria_coeff_y = 1e-6
            
        tria_end = tria_start+6
    
            
    elif (id_vec_type_equ==2 and id_FEM==1):
        tria_start = 35
        tria_end = tria_start+3
        tria_coeff_y = 3e-2
        
        coeff_text_slope_bottom_x = 0.6
        coeff_text_slope_bottom_y = 0.4
        coeff_text_slope_right_x = 1.1
        coeff_text_slope_right_y = 0.54
        
    tria_p_1 = [trend_line_x_coord[tria_start],trend_line_rounding[tria_start]]
    tria_p_2 = [trend_line_x_coord[tria_end],trend_line_rounding[tria_start]]
    tria_p_3 = [trend_line_x_coord[tria_end],trend_line_rounding[tria_end]]
    
    tria_x_all_control = [tria_p_1[0],tria_p_2[0],tria_p_3[0],tria_p_1[0]]
    tria_y_all_control = [tria_p_1[1],tria_p_2[1],tria_p_3[1],tria_p_1[1]]
    tria_y_all_control = [i*tria_coeff_y for i in tria_y_all_control]
    
    text_slope_bottom_x = (tria_x_all_control[0]+tria_x_all_control[1])/2*coeff_text_slope_bottom_x
    text_slope_bottom_y = tria_y_all_control[0]*coeff_text_slope_bottom_y
    text_slope_right_x = tria_x_all_control[1]*coeff_text_slope_right_x
    text_slope_right_y = (tria_y_all_control[1]+tria_y_all_control[2])/2*coeff_text_slope_right_y   
    
    plt.plot(tria_x_all_control,tria_y_all_control,'k-',linewidth=1.0)
    plt.text(text_slope_bottom_x,text_slope_bottom_y,'1',fontsize = 15)
    plt.text(text_slope_right_x,text_slope_right_y,str(vec_slope[id_var]),fontsize = 15)    
    
    
    
plt.xlabel(vec_xlabel[id_vec_type_equ], fontsize = fontsize_label)  

plt.ylabel(r"$\alpha_{\rm R}$", fontsize = fontsize_label)        
               
if id_vec_type_equ != 0:
    axes.set_xlim([xaxis_low_bound, xaxis_up_bound])
    
axes.set_ylim([yaxis_low_bound, yaxis_up_bound])

if id_vec_type_equ==0:
    
    plt.xticks(coords_x_axis,labels_x_axis,fontsize=set.fontsize_tick)
elif id_vec_type_equ==1:
    plt.yticks([1e-25, 1e-22, 1e-19, 1e-16, 1e-13, 1e-10])
elif id_vec_type_equ==2:
    plt.yticks([1e-20, 1e-18, 1e-16, 1e-14, 1e-12])

plt.legend(fontsize=set.fontsize_legend, loc=id_loc_legend)
#            
f.savefig("py_offset_summary_" + vec_type_equ[id_vec_type_equ] + "_" +vec_FEM[id_FEM]+".pdf", bbox_inches='tight')   
        
            
    

