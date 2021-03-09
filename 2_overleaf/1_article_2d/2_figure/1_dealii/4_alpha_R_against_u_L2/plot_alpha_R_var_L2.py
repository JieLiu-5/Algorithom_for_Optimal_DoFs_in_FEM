

import sys
sys.path.append('../../../../3_1d_or_2d')
import settings_fig as set


import matplotlib.pyplot as plt
import numpy as np



vec_equ=['0_all_equations','1_u_e_m_c_x_m_0p5_square','2_u_1_over_c_x_m_0p5_square','3_u_cx_plus_100_square','4_u_x_square_m_c_square','5_u_1_over_c_x_m_0p5_square_p_same_for_y']
id_equ=5

vec_markeredgecolor = ['k', 'b', 'r', 'c', 'm', 'y', 'g', 'w']

vec_FEM=['1_sm','2_mm_PPdisc','3_mm_RT','4_mm_BDM']
id_FEM=3
id_error_space=1


vec_xaxis=['ndofs','ncond']
id_xaxis=0


vec_legend=['$u$',r'${\nabla} u$','$\Delta u$']
vec_offset=[2e-17,1e-16,2e-15]
vec_slope=[1,1,1,1]

id_loc_legend=4

line_auxiliary_x = [1,1,1e-6]
line_auxiliary_y = [1e-25,0,0]

coord_x_legend=1.5e-6
vec_coord_y_text = [3e-18, 2e-16, 4e-15]

xaxis_low_bound=1e-6
xaxis_up_bound=1e4
yaxis_low_bound=1e-25
yaxis_up_bound=1e-10
vec_xtick = [1e-6, 1e-4, 1e-2, 1e0, 1e2, 1e4]
vec_ytick = [1e-25, 1e-22, 1e-19, 1e-16, 1e-13, 1e-10]

trend_line_x_coord = np.logspace(-7, 6, 50)


if id_FEM==1:
    id_loc_legend=2
    vec_legend=['$u$',r'$\mathbf{v}$',r'$\nabla \cdot \mathbf{v}$']
    id_error_space=2
    
    if id_equ==0:
        id_loc_legend=4
        vec_offset=[2e-17,1e-16,3e-16]
        vec_coord_y_text = [ 0.2e-17,0.2e-16,7e-16]  
        xaxis_up_bound=1e6
        vec_xtick = [1e-6, 1e-4, 1e-2, 1e0, 1e2, 1e4, 1e6]         
    elif id_equ==1:
        if id_xaxis==0:
            vec_offset=[1e-17,1e-16,5e-16]
            vec_coord_y_text = [ 0.15e-17,0.15e-16,1e-15]        
        elif id_xaxis==1:
            vec_offset=[1e-19,1e-18,5e-20]
            vec_coord_y_text = [ 1.5e-19,1.5e-18,0.8e-20]
    elif id_equ==3 or id_equ==0:
        id_loc_legend=4
        vec_offset=[3e-17,1.5e-16,3e-16]
        vec_coord_y_text = [ 0.2e-17,0.2e-16,7e-16]  
        xaxis_up_bound=1e6
        vec_xtick = [1e-6, 1e-4, 1e-2, 1e0, 1e2, 1e4, 1e6]        
elif id_FEM==2:
    id_loc_legend=4
    vec_legend=['$u$',r'$\mathbf{v}$',r'$\nabla \cdot \mathbf{v}$']
    id_error_space=2
    if id_equ==2:
        vec_offset=[2e-17,5e-16,5e-16]
        vec_coord_y_text = [ 0.2e-17,7e-16,7e-16]
    elif id_equ==4 or id_equ==0:
        vec_offset=[5e-17,3e-16,3e-16]
        vec_coord_y_text = [ 0.7e-17,7e-16,7e-16] 
        xaxis_up_bound=1e6
        vec_xtick = [1e-6, 1e-4, 1e-2, 1e0, 1e2, 1e4, 1e6]
    elif id_equ==5:
        vec_offset=[7e-17,3e-16,3e-16]
        vec_coord_y_text = [ 1.2e-17,7e-16,7e-16]
elif id_FEM==3:
    id_loc_legend=4
    vec_legend=['$u$',r'$\mathbf{v}$',r'$\nabla \cdot \mathbf{v}$']
    id_error_space=2
    if id_equ==5:
        vec_offset=[4e-16,4e-15,1e-14]
        vec_coord_y_text = [ 5e-17,7e-16,2e-14]
           

txt_alpha_R_var_L2=open(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+'/alpha_R_var_L2_'+vec_xaxis[id_xaxis]+'.txt','r')

data_alpha_R_var_L2=[line.strip().split() for line in txt_alpha_R_var_L2]

txt_alpha_R_var_L2.close()

n_param = len(data_alpha_R_var_L2)-1

data_u_L2=np.zeros(n_param)
data_u_x_L2=np.zeros(n_param)
data_alpha_R=np.zeros((n_param,3))

for id_param in range(0,n_param):
    data_u_L2[id_param]=data_alpha_R_var_L2[id_param+1][0]
    if id_FEM==1 or id_FEM==2 or id_FEM==2:
        data_u_x_L2[id_param]=data_alpha_R_var_L2[id_param+1][1]

    for id_var in range(0,3):
        data_alpha_R[id_param][id_var]=data_alpha_R_var_L2[id_param+1][id_var+id_error_space]


f=plt.figure(figsize=(5,4))
axes= f.add_axes([0,0,1,1])
fontsize_label = 18
plt.tick_params(axis = 'both', which = 'major', labelsize = 16)


id_var_start = 0
id_var_end = 3

for id_var in range(id_var_start,id_var_end):
    line_auxiliary_y[1]=vec_offset[id_var]
    line_auxiliary_y[2]=line_auxiliary_y[1]
        
    if id_FEM==0:
        plt.loglog(data_u_L2,data_alpha_R[:,id_var],'o',color = vec_markeredgecolor[id_var],mfc = "none",label=vec_legend[id_var])
    elif id_FEM==1 or id_FEM==2 or id_FEM==3:
        if id_var==0 or id_var==1:
            plt.loglog(data_u_L2,data_alpha_R[:,id_var],'o',color = vec_markeredgecolor[id_var],mfc = "none",label=vec_legend[id_var])
        elif id_var==2:
            plt.loglog(data_u_L2,data_alpha_R[:,id_var],'o',color = vec_markeredgecolor[id_var],mfc = "none",label=vec_legend[id_var])
    
    plt.plot(line_auxiliary_x,line_auxiliary_y,'--',linewidth=1.0,color='grey')
    
    plt.text(coord_x_legend, vec_coord_y_text[id_var], str(vec_offset[id_var]), fontsize=set.fontsize_legend, color = vec_markeredgecolor[id_var])
    
    trend_line_rounding = vec_offset[id_var] * trend_line_x_coord**vec_slope[id_var]     
    plt.loglog(trend_line_x_coord, trend_line_rounding, '-' + vec_markeredgecolor[id_var])
             
    
tria_start = 32
tria_coeff_y = 1e-4

coeff_text_slope_bottom_x = 0.25
coeff_text_slope_bottom_y = 0.14
coeff_text_slope_right_x = 1.3
coeff_text_slope_right_y = 0.17

if id_FEM==1 or id_FEM==2 or id_FEM==3:
    tria_start = 32
    tria_coeff_y = 1e-3
    if id_equ==2:
        tria_coeff_y = 1e-3

tria_end = tria_start+6

trend_line_rounding = 1e-16 * trend_line_x_coord**vec_slope[id_FEM]
 
    
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
plt.text(text_slope_right_x,text_slope_right_y,str(vec_slope[id_FEM]),fontsize = 15)    
    


    
plt.xlabel(r"$\||u\||_2$", fontsize = fontsize_label)  
plt.ylabel(r"$\alpha_{\rm R}$", fontsize = fontsize_label) 

axes.set_xlim([xaxis_low_bound, xaxis_up_bound])
axes.set_ylim([yaxis_low_bound, yaxis_up_bound])

plt.xticks(vec_xtick)
plt.yticks(vec_ytick)

plt.legend(fontsize=set.fontsize_legend, loc=id_loc_legend)

f.savefig(vec_equ[id_equ]+'/'+vec_FEM[id_FEM]+"/py_offset_summary_"+vec_equ[id_equ]+"_"+vec_FEM[id_FEM]+'_' +vec_xaxis[id_xaxis]+".pdf", bbox_inches='tight')   
