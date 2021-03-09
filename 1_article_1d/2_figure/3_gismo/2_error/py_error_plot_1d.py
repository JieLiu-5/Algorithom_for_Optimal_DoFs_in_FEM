

import sys
sys.path.append('../../../../')
import settings as set

import matplotlib.pyplot as plt
import numpy as np


xaxis_up_bound = 1e6
yaxis_low_bound = 1e-16
yaxis_up_bound = 1e0


vec_marker=['o','d','^','s','*']  
vec_markeredgecolor = ['k', 'b', 'r', 'c', 'm', 'y', 'g', 'w']

vec_FEM = ['SM','MM']

vec_var = ['solu','grad','2ndd']
vec_offset = [5e-18, 1e-17, 5e-16]


n_coeff = 5
vec_coeff = ['1', '2', '3', '4', '5']


vec_slope = [2, 2]
vec_legend=['1','2','3','4','5']


n_err_refer = 28
dof_ref = np.zeros(n_err_refer)
err_round_off_approx = np.zeros(n_err_refer)

for i in range(n_err_refer):
    dof_ref[i] = 2**i  


id_FEM = 0

id_sensitivity = 1

id_scaling = 0

id_var_start = 0
id_var_end = 2



n_column = 5

id_case_start = 0
id_case_end = id_case_start+n_column

id_c=2

n_refine = 20

n_err_refer = 28


dataraw = [line.strip().split() for line in open('data_error_deg_'+vec_coeff[0]+'.txt','r')]

n_refine = len(dataraw)-1



data_ndofs = np.zeros((n_refine, id_case_end-id_case_start))            
data_error = np.zeros((n_refine, id_case_end-id_case_start))


trend_line_x_coord = np.logspace(0, 6, 50)

id_loc_offset=10

text_offset_x = trend_line_x_coord[id_loc_offset]
text_offset_y = 1e-16

    

for id_var in range(id_var_start, id_var_end):
        
    print ('var: '+vec_var[id_var])  
    
    f=plt.figure(figsize=(5,4))
    axes= f.add_axes([0,0,1,1])
          
    for i in range(n_err_refer):
        err_round_off_approx[i] = vec_offset[id_var]* dof_ref[i]**vec_slope[id_var]                

    for id_case in range(id_case_start, id_case_end):
        
        print('deg: '+vec_coeff[id_case])

        dataraw = [line.strip().split() for line in open('data_error_deg_'+vec_coeff[id_case]+'.txt','r')]
               
        
        for i in range(n_refine):
            data_ndofs[i][id_case-id_case_start]=dataraw[i+1][0]    
            data_error[i][id_case-id_case_start]=dataraw[i+1][id_var+1]    
    
            if(data_error[i][id_case]==0):
                data_error[i][id_case] = np.nan
                    
    
        plt.loglog(data_ndofs[:,id_case-id_case_start], data_error[:,id_case-id_case_start],'k'+vec_marker[id_case-id_case_start], markerfacecolor='none',label='$p$='+(vec_legend[id_case]),linewidth=1.0)    
                
    trend_line_rounding = vec_offset[id_var] * trend_line_x_coord**vec_slope[id_FEM]     
    plt.loglog(trend_line_x_coord, trend_line_rounding, '--',color='orange')
    
    text_offset_y = trend_line_rounding[id_loc_offset]*0.2
    plt.text(text_offset_x,text_offset_y,r'$\alpha_{\rm R}=$'+str(vec_offset[id_var]),fontsize = 15) 
                
    tria_start = 8
    
    tria_coeff_x = 1e1
    tria_coeff_y = 1e-2

    
    text_slope_bottom_coeff_x = 0.5
    text_slope_bottom_coeff_y = 0.15
    text_slope_right_coeff_x = 1.15
    text_slope_right_coeff_y = 0.12
    
    text_offset_coeff_y=2e-1
    

    tria_end = tria_start+3
    
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
    
    
    plt.plot(tria_x,tria_y,'k-',linewidth=1.0)
    plt.text(text_slope_bottom_x,text_slope_bottom_y,'1',fontsize = 15)
    plt.text(text_slope_right_x,text_slope_right_y,str(vec_slope[id_var]),fontsize = 15)       
    
    if id_var == 0:    
        plt.legend()
        plt.legend(loc=set.id_loc_legend, prop={'size': set.fontsize_legend})
        
    plt.xlabel('Number of DoFs', fontsize=set.fontsize_label)
    plt.ylabel('Error', fontsize=set.fontsize_label)
    
    axes.set_xlim([1, xaxis_up_bound])
    axes.set_ylim([yaxis_low_bound, yaxis_up_bound])           
    

    plt.xticks([1e0, 1e2, 1e4, 1e6])      
    plt.yticks([1e-16, 1e-12, 1e-8, 1e-4, 1e0])
                

    plt.tick_params(axis='both', which='major', labelsize=set.fontsize_tick)
    
    plt.show()
    
    f.savefig("py_"+vec_var[id_var]+".pdf", bbox_inches='tight')    
    