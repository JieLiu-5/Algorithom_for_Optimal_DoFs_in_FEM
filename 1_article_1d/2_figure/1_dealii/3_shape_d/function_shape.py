#clear all;
#close all;

import matplotlib.pyplot as plt
import numpy as np

#import math


vec_color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

coord = np.linspace(0, 1, 1000)

vec_func = ["sin_2pi_c_x", "exp_minus_c_x_minus_0p5_square"]               # they represent types of functions

def func_u(x,c):                                                           # initial definition for functions
    return c*x**2 

def func_du(x,c):
    return 2*c*x

def func_ddu(x,c):
    return 1+0.5*np.sin(c*x)

vec_coeff = [1e-2, 1e-1, 1, 1e1, 1e2, 1e4]                                      # initial vector of coeff

vec_legend = ["0.01", "0.1", "1", "10", "100", "10000"]                   # initial vector of legend

vec_var = ["Solution", "First derivative", "Second derivative"]

id_func = 1            # '1' for u=(1/(2*pi*c))^2*sin(2*pi*c*x)

id_var = 1

id_var_start = 2
id_var_end = 3

id_coeff_start = 2
id_coeff_end = 5

fontsize_label = 18

id_loc_legend = 2
legend_size = 14


print("plotting started\n")
        
for id_var in range(id_var_start, id_var_end):
    
    print(vec_var[id_var])

    f=plt.figure(figsize=(5,4))
    axes= f.add_axes([0.1, 0.1, 0.8, 0.8])
    
#    plt.xticks([1e0, 1e2, 1e4, 1e6, 1e8])      
    plt.yticks([0.2, 0.6, 1.0, 1.4, 1.8])    
    
    xaxis_up_bound = 1
    yaxis_low_bound = 0.2
    yaxis_up_bound = 1.8
    
    
    axes.set_xlim([0, xaxis_up_bound])
    axes.set_ylim([yaxis_low_bound, yaxis_up_bound])       
    
    for id_coeff in range(id_coeff_start, id_coeff_end):
        
        print('id_coeff:', id_coeff, end =", ")
        print('value:', vec_coeff[id_coeff])
        
        if id_var == 0:
            plt.plot(coord, func_u(coord, vec_coeff[id_coeff]), label="c"+vec_legend[id_coeff])
        elif id_var == 1:
            plt.plot(coord, func_du(coord, vec_coeff[id_coeff]), 'k', label=vec_legend[id_coeff])
        elif id_var == 2:
            plt.plot(coord, func_ddu(coord, vec_coeff[id_coeff]), vec_color[id_coeff], label="c="+str(vec_legend[id_coeff]))
        else:
            print("Invalid variable\n")


    plt.xlabel('x', fontsize=fontsize_label)
    plt.ylabel('y', fontsize=fontsize_label)                       # vec_var[id_var]
                        
    plt.legend()
    plt.legend(loc=id_loc_legend, prop={'size': legend_size})

    plt.tick_params(axis='both', which='major', labelsize=16)
    
    plt.show()
#    
    f.savefig("py_"+"var_1_plus_sincx"+".pdf", bbox_inches='tight')                  # vec_func[id_func]+"_"+vec_var[id_var]
#    
    print("plotting finished")
