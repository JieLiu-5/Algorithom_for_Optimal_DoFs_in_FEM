

import sys
sys.path.append('../../../../3_1d_or_2d')
import settings_fig as set


import matplotlib.pyplot as plt
import numpy as np


vec_dimension=['1d','2d']
id_dimension=0

vec_problem_type=['0_pois','1_diff','2_helm']
id_problem_type=2

vec_equ=['0_u_x_m_0p5_square', '1_u_exp_m_x_m_0p5_square']
id_equ = 1

vec_markeredgecolor = ['k', 'b', 'r', 'c', 'y', 'm', 'g', 'w']

vec_FEM=['0_sm']
id_FEM=0

vec_legend=['Uniformly',r'Randomly','Regularly-linearly','Regularly-sine']
vec_line_type = ['-','--','--','--','--','--','--']


fontsize_legend=12
id_loc_legend=2

bound_low_xaxis=1
bound_up_xaxis=1e8
bound_low_yaxis=1e-20
bound_up_yaxis=1e0
vec_ytick = [1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e0]

data_x_axis=['$u$',r'$u_{x}$','$u_{xx}$']

file_alpha_R_var_mesh=open(vec_dimension[id_dimension]+'/'+vec_problem_type[id_problem_type]+'/'+vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/alpha_R_var_mesh.txt','r')
data_alpha_R_var_mesh=[line.strip().split() for line in file_alpha_R_var_mesh]
file_alpha_R_var_mesh.close()

file_beta_R_var_mesh=open(vec_dimension[id_dimension]+'/'+vec_problem_type[id_problem_type]+'/'+vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/beta_R_var_mesh.txt','r')
data_beta_R_var_mesh=[line.strip().split() for line in file_beta_R_var_mesh]
file_beta_R_var_mesh.close()


n_cases = len(data_alpha_R_var_mesh)-1

data_alpha_R=np.zeros((n_cases,3))
data_beta_R=np.zeros((n_cases,3))

for i in range(n_cases):
    data_alpha_R[i][:] = data_alpha_R_var_mesh[i+1][:]
    data_beta_R[i][:] = data_beta_R_var_mesh[i+1][:]
    
vec_offset=[1e-18,5e-18,1e-16]
vec_slope=[2.0, 2.0, 2.0]

if id_problem_type==0:
    if id_equ==1:
        vec_offset=[2e-17,1e-16,5e-16]
elif id_problem_type==1:
    vec_legend=['Poisson, Uniformly',r'$D=1+x$, Uniformly',r'$D=1+x$, Distorted',r'$D=e^{-(x-0.5)^2}$, Uniformly',r'$D=0.5+\cos^2(x)$, Uniformly']
    if id_equ==1:
        vec_offset=[2e-17,1e-16,5e-16]
elif id_problem_type==2:
    if id_equ==0:
        vec_legend=['Poisson', r'$D=1+x, r=1$', r'$D=1+x, r=1+x$', r'$D=1+x, r=e^{-(x-0.5)^2}$', r'$D=e^{-(x-0.5)^2}, r=1$', r'$D=e^{-(x-0.5)^2}, r=1+x$', r'$D=e^{-(x-0.5)^2}, r=e^{-(x-0.5)^2}$']
    elif id_equ==1:
        vec_offset=[2e-17,1e-16,5e-16]
        vec_legend=['Poisson', r'$D=1+x, r=1$', r'$D=1+x, r=1+x$', r'$D=1+x, r=e^{-(x-0.5)^2}$', r'$D=e^{-(x-0.5)^2}, r=1$', r'$D=e^{-(x-0.5)^2}, r=1+x$', r'$D=e^{-(x-0.5)^2}, r=e^{-(x-0.5)^2}$']
    
trend_line_x_coord = np.logspace(0, 8, 50)

color_danger_zone = "yellow"
coeff_strip_line_second = 0.02
    
vec_var = ['solu','grad','2ndd']

id_var_start = 0
id_var_end = 3

for id_var in range(id_var_start,id_var_end):
    
    print ('variable: '+vec_var[id_var])
    
    
    f=plt.figure(figsize=(5,4))
    axes= f.add_axes([0,0,1,1])
    fontsize_label = 18
    plt.tick_params(axis = 'both', which = 'major', labelsize = 16)
    
 
    for i in range(n_cases):
        trend_line_rounding = data_alpha_R[i][id_var] * trend_line_x_coord**data_beta_R[i][id_var]
        plt.loglog(trend_line_x_coord, trend_line_rounding, vec_line_type[i] + vec_markeredgecolor[i], mfc='none',label=vec_legend[i])


    if id_var==2:
        coeff_strip_line_second = 20
        
    strip_line_first = data_alpha_R[0][id_var] * trend_line_x_coord**data_beta_R[0][id_var]
    strip_line_second = strip_line_first*coeff_strip_line_second
       
    plt.fill_between(trend_line_x_coord,strip_line_first,strip_line_second,facecolor=color_danger_zone,alpha=1.0)        # facecolor='yellow', 
    
    tria_start = 20
    tria_coeff_y = 1e-4
    
    coeff_text_slope_bottom_x = 0.5
    coeff_text_slope_bottom_y = 0.10
    coeff_text_slope_right_x = 1.3
    coeff_text_slope_right_y = 0.12
    

    if id_problem_type==1:
        coeff_text_slope_right_x = 1.15
        coeff_text_slope_right_y = 0.12
        
        tria_coeff_y = 1e-4

            
    tria_end = tria_start+6
    
    tria_p_1 = [trend_line_x_coord[tria_start],strip_line_first[tria_start]]
    tria_p_2 = [trend_line_x_coord[tria_end],strip_line_first[tria_start]]
    tria_p_3 = [trend_line_x_coord[tria_end],strip_line_first[tria_end]]
    
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


    text_offset_x = trend_line_x_coord[1]
    
    text_offset_coeff_y=0.2
    
    text_offset_y = vec_offset[id_var]* text_offset_x**vec_slope[id_var]*text_offset_coeff_y       
    
    plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15, color="black")             
    
    
    #for i in range(n_cases):
    #    plt.semilogy(data_x_axis, data_alpha_R[i,:],'--o',color = vec_markeredgecolor[i],mfc = "none",label=vec_legend[i])
    #    
    plt.xlabel(r"Number of DoFs", fontsize = fontsize_label)  
    plt.ylabel(r"Error", fontsize = fontsize_label) 
    
    axes.set_xlim([bound_low_xaxis, bound_up_xaxis])
    axes.set_ylim([bound_low_yaxis, bound_up_yaxis])
    
    plt.yticks(vec_ytick)
    
    if id_var==0:
        plt.legend(fontsize=fontsize_legend, loc=id_loc_legend)
    
    f.savefig(vec_dimension[id_dimension]+'/'+vec_problem_type[id_problem_type]+'/'+vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+"/py_alpha_R_beta_R_summary_mesh_"+vec_dimension[id_dimension]+'_'+vec_problem_type[id_problem_type]+'_'+vec_equ[id_equ]+"_"+vec_FEM[id_FEM]+ "_" + vec_var[id_var] +".pdf", bbox_inches='tight')   
