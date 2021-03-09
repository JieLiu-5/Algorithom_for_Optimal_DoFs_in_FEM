

import sys
sys.path.append('../../3_1d_or_2d/')
import settings_fig as set

import matplotlib.pyplot as plt
import numpy as np


vec_package = ['1_dealii', '2_fenics', '3_gismo']
id_package = 1

print("results of", vec_package[id_package])

vec_dimension = ['1d', '2d']
id_dimension = 0

vec_problem = ["0_u_x_m_0p5_square"]
vec_problem.append("1_u_x_m_0p5_square_d_1px")
vec_problem.append("2_u_x_m_0p5_square_d_1pxpxsquare")
vec_problem.append("3_u_x_m_0p5_square_d_exp_m_x_m_0p5_square")
vec_problem.append("4_u_x_m_0p5_square_d_0p5_p_cosx_square")

if id_dimension==1:
    vec_problem = ["0_u_x_m_0p5_square_p_y_same"]
    vec_problem.append("1_u_exp_m_x_m_0p5_square")
    
id_problem = 0

vec_FEM = ["0_sm", "1_bdm"]
id_FEM = 0

color_round_off_line_main = 'orange'
color_round_off_line_auxiliary = 'tomato'

if id_problem==0:
    color_round_off_line_main = 'tomato'


xaxis_up_bound = 1e8
yaxis_low_bound = 1e-20
yaxis_up_bound = 1e0


vec_xticks=[1e0, 1e2, 1e4, 1e6, 1e8]
vec_yticks=[1e-20,1e-16, 1e-12, 1e-8, 1e-4, 1e0]

id_loc_legend = 4


vec_var = ['solu','grad','2ndd']
vec_offset = [2e-17, 1e-16, 2e-15]
vec_slope = [1.0, 1.0, 1.0]

vec_offset_auxiliary = [5e-18, 1e-17, 1e-16]
vec_slope_auxiliary = [2.0, 2.0, 2.0]


if id_dimension==0:
    if id_problem==0:
        if id_FEM==0:
            vec_offset = [5e-18, 1e-17, 1e-16]
            vec_slope = [2.0, 2.0, 2.0]
    elif id_problem==1:
        if id_FEM==0:
            vec_offset = [5e-18, 2e-17, 1e-16]
            vec_slope = [1.5, 1.5, 2.0] 
    elif id_problem==2:
        if id_FEM==0:
            vec_offset = [5e-18, 1e-17, 1e-16]
            vec_slope = [1.5, 1.5, 2.0]  
    elif id_problem==3:
        if id_FEM==0:
            vec_offset = [5e-18, 1e-17, 1e-16]
            vec_slope = [1.5, 1.5, 2.0]  
    elif id_problem==4:
        if id_FEM==0:
            vec_offset = [5e-18, 1e-17, 1e-16]
            vec_slope = [1.5, 1.5, 2.0]              
elif id_dimension==1:
    
    yaxis_low_bound = 1e-16
    vec_yticks=[1e-16, 1e-12, 1e-8, 1e-4, 1e0]
    
    if id_problem==0:
        if id_FEM==0:
            vec_offset = [2e-17, 1e-16, 2e-15]
            vec_slope = [1.0, 1.0, 1.0]            
        elif id_FEM==1:
            vec_offset = [5e-15, 2e-15, 5e-15]
            vec_slope = [0.0, 0.5, 1.0]
    elif id_problem==1:
        if id_FEM==0:
            vec_offset = [2e-16, 5e-16, 1e-14]
        elif id_FEM==1:
            vec_offset = [5e-15, 2e-15, 5e-15]
            vec_slope = [0.0, 0.5, 1.0]

vec_coeff = [1, 2, 3, 4, 5, 6, 7, 8]

id_deg_start = 0
n_column = 5

n_err_refer = 28

dof_ref = np.zeros(n_err_refer)
trend_line_rounding = np.zeros(n_err_refer)
trend_line_rounding_auxiliary = np.zeros(n_err_refer)

for i in range(n_err_refer):
    dof_ref[i] = 2**i  


data_raw = [line.strip().split() for line in open(vec_package[id_package]+'/2_error/'+vec_dimension[id_dimension]+'/'+vec_problem[id_problem]+'/'+vec_FEM[id_FEM]+'/data_error_deg_'+str(vec_coeff[id_deg_start])+'.txt','r')]

n_refine_initial = len(data_raw)-1

data_ndofs = np.zeros((n_refine_initial, n_column)) 
data_error = np.zeros((n_refine_initial, n_column))


trend_line_x_coord = np.logspace(0, 8, 50)
    

id_case_start = 0
id_case_end = id_case_start+n_column

id_refine=1

id_var_start = 0
id_var_end = 3


for id_var in range(id_var_start, id_var_end):
        
    print ('var: '+vec_var[id_var])  
    
    f=plt.figure(figsize=(5,4))
    axes= f.add_axes([0,0,1,1])
    
    for id_case in range(id_case_start, id_case_end):
        
        print('deg:',vec_coeff[id_deg_start+id_case])

        data_raw = [line.strip().split() for line in open(vec_package[id_package]+'/2_error/'+vec_dimension[id_dimension]+'/'+vec_problem[id_problem]+'/'+vec_FEM[id_FEM]+'/data_error_deg_'+str(vec_coeff[id_deg_start+id_case])+'.txt','r')]
        
        n_refine = len(data_raw)-1
        for i in range(n_refine):
            
            id_refine=i+1
            
            if id_FEM == 0:
                if id_dimension==0:
                    data_ndofs[i][id_case-id_case_start]=2**id_refine*vec_coeff[id_case]+1
                elif id_dimension==1:
                    data_ndofs[i][id_case-id_case_start]=((2**id_refine*vec_coeff[id_case]+1)+2**id_refine*(vec_coeff[id_case]-1))*(2**id_refine+1)+4**id_refine*(4*(vec_coeff[id_case]-1)+4*((vec_coeff[id_case]-2)*(vec_coeff[id_case]-1)/2)+1)
                
                
                data_error[i][id_case-id_case_start]=data_raw[i+1][id_var]
            elif id_FEM == 1:
                if id_var==0:
                    data_ndofs[i][id_case-id_case_start]=4**(id_refine+1)*(vec_coeff[id_case]*(vec_coeff[id_case]+1)/2)
                elif id_var==1 or id_var==2:
                    data_ndofs[i][id_case-id_case_start]=(2**id_refine*(vec_coeff[id_case]+1))*(2**id_refine+1)*2+4**id_refine*(4*(vec_coeff[id_case]+1)+4*(vec_coeff[id_case]**2-1))
                data_error[i][id_case-id_case_start]=data_raw[i+1][id_var]
    
    
            if(data_ndofs[i][id_case]==0):
                data_ndofs[i][id_case] = np.nan
                
            if(data_error[i][id_case]==0):
                data_error[i][id_case] = np.nan
                    
    
        plt.loglog(data_ndofs[:,id_case-id_case_start], data_error[:,id_case-id_case_start],'k'+set.vec_marker[id_case-id_case_start], markerfacecolor='none',label='$p$='+str(vec_coeff[id_deg_start+id_case]),linewidth=1.0)    
                
    trend_line_rounding = vec_offset[id_var] * trend_line_x_coord**vec_slope[id_var]
    trend_line_rounding_auxiliary = vec_offset_auxiliary[id_var] * trend_line_x_coord**vec_slope_auxiliary[id_var]
    
    plt.loglog(trend_line_x_coord, trend_line_rounding, '--',color=color_round_off_line_main)
    if id_problem > 0:
        plt.loglog(trend_line_x_coord, trend_line_rounding_auxiliary, '--',color=color_round_off_line_auxiliary)
    
    id_loc_offset=3
    
    tria_start = 8
    
    tria_coeff_x = 1e1
    tria_coeff_y = 1e0
    
    text_slope_bottom_coeff_x = 0.5
    text_slope_bottom_coeff_y = 0.1
    text_slope_right_coeff_x = 1.15
    text_slope_right_coeff_y = 0.2
    
    coeff_text_offset_y=2e-1
    
    
    if id_dimension==0:
        id_loc_offset=1
    elif id_dimension==1:
        
        text_slope_bottom_coeff_y = 0.15
        text_slope_right_coeff_y = 0.5
        
        if id_problem == 0:
            id_loc_offset=8
            if id_FEM==0:
                if id_var==0:
                    tria_start = 14
                if id_var==1:
                    id_loc_offset=5
                elif id_var==2:
                    id_loc_offset=1
            elif id_FEM==1:
                if id_var==0:
                    id_loc_offset=0
                    tria_start = 16
                    coeff_text_offset_y = 2
                elif id_var==1:
                    id_loc_offset=1
                    coeff_text_offset_y = 2e1
                elif id_var==2:
                    id_loc_offset=1            
        elif id_problem == 1:
            id_loc_offset=1
            if id_package==1:
                if id_var==0:
                    id_loc_offset=2
            
    tria_end = tria_start+4
    
    
                
    text_offset_x = dof_ref[id_loc_offset]
    text_offset_y = vec_offset[id_var]* text_offset_x**vec_slope[id_var]*coeff_text_offset_y  
    plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15) 
                



    tria_p_1 = [dof_ref[tria_start],trend_line_rounding[tria_start]]
    tria_p_2 = [dof_ref[tria_end],trend_line_rounding[tria_start]]
    tria_p_3 = [dof_ref[tria_end],trend_line_rounding[tria_end]]
    
    tria_x = [tria_p_1[0],tria_p_2[0],tria_p_3[0],tria_p_1[0]]#*
    tria_x = [i*tria_coeff_x for i in tria_x]
    tria_y = [tria_p_1[1],tria_p_2[1],tria_p_3[1],tria_p_1[1]]#/tria_coeff_y
    tria_y = [i*tria_coeff_y for i in tria_y]
    
    text_slope_bottom_x = (tria_x[0]+tria_x[1])/2*text_slope_bottom_coeff_x
    text_slope_bottom_y = tria_y[0]*text_slope_bottom_coeff_y
    text_slope_right_x = tria_x[1]*text_slope_right_coeff_x
    text_slope_right_y = (tria_y[1]+tria_y[2])/2*text_slope_right_coeff_y         
    
    plt.plot(tria_x,tria_y,'k-',linewidth=1.0)
    plt.text(text_slope_bottom_x,text_slope_bottom_y,'1',fontsize = 15)
    plt.text(text_slope_right_x,text_slope_right_y,str(vec_slope[id_var]),fontsize = 15)  
    
    
    if id_var == 0:    
        plt.legend()
        plt.legend(loc=id_loc_legend, prop={'size': set.fontsize_legend})
        
    plt.xlabel('Number of DoFs', fontsize=set.fontsize_label)
    plt.ylabel('Error', fontsize=set.fontsize_label)
    
    axes.set_xlim([1, xaxis_up_bound])
    axes.set_ylim([yaxis_low_bound, yaxis_up_bound])           
    
    plt.xticks(vec_xticks)      
    plt.yticks(vec_yticks)
                

    plt.tick_params(axis='both', which='major', labelsize=set.fontsize_tick)
    
    plt.show()
    
    f.savefig(vec_package[id_package]+"/2_error/"+vec_dimension[id_dimension]+'/'+vec_problem[id_problem]+'/'+vec_FEM[id_FEM]+"/py_error_"+vec_var[id_var]+".pdf", bbox_inches='tight')    
    
