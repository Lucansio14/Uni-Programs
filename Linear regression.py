# -*- coding: utf-8 -*-
"""
Created on Tue Oct 21 19:10:17 2025
Algorithm for the utilization of a simple linear regression statistical method (without using SciPy), which attempts to find a linear relationship between
a dependent variable (y) and an independent variable (x) in this manner: y = bx + a, where b is the slope of the line and a is the y-intercept.
More general information on these websites: https://en.wikipedia.org/wiki/Linear_regression
                                            https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
@author: Lucas Romero Fernández
"""
import time
import numpy as np
import matplotlib.pyplot as plt
#main_program
start_time_program=time.process_time()#To calculate the program execution time
#Definition of variables, constants, arrays and functions
x_array=np.array([0.79,1.05,1.57,1.83,2.09,2.36])#Example data for x
err_x=0.09#Error bar value for x
y_array=np.array([0.022,0.031,0.039,0.048,0.056,0.062])#Example data for y
err_y=0.003#Error bar value for y
#Linear regresion method
sample_size=x_array.size
mean_x=x_array.mean()
mean_y=y_array.mean()
sum_dif_xy=((x_array-mean_x)*(y_array-mean_y)).sum()
sum_dif_x2=((x_array-mean_x)**2).sum()
sum_dif_y2=((y_array-mean_y)**2).sum()
r=sum_dif_xy/np.sqrt(sum_dif_x2*sum_dif_y2)#Pearson correlation coefficient
b=sum_dif_xy/sum_dif_x2
a=mean_y-b*mean_x
delta=np.sqrt(((a+b*x_array-y_array)**2).sum()/(sample_size-2))
err_a=delta*np.sqrt(1/sample_size+mean_x*mean_x/sum_dif_x2)
err_b=delta*np.sqrt(1/sum_dif_x2)
print("Results of the sample data fitting to the form of y = a + bx:",sep="")
print("r = ",r,", r**2 = ",r**2,",",sep="")
print("b = ",b,", a = ",a,",",sep="")
print("err_b = ±",err_b,", err_a = ±",err_a,",",sep="")
#Plot of the results
plt.figure(figsize=(9,5))
plt.plot(x_array,a+b*x_array,color="blue",label="$y=({b:7.3f}±{err_b:7.3f})x + ({a:7.3f}±{err_a:7.3f})$".format(b=b,err_b=err_b,a=a,err_a=err_a))
plt.errorbar(x_array,y_array,xerr=err_x,yerr=err_y,fmt="g.",label="Sample data",markersize=12,capsize=5)
plt.plot(x_array,(a+err_a)+(b+err_b)*x_array,"r--")
plt.plot(x_array,(a-err_a)+(b-err_b)*x_array,"r--")
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("Linear regression of the sample data")
plt.legend(fontsize=14)
plt.grid()
plt.tight_layout()
plt.show()
print("Program execution time:",time.process_time()-start_time_program,"seconds.")