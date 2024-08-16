# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 18:29:42 2023
Is *blank* a Prime Number? (Method 1)
@author: Lucas Romero Fernández
"""
import time
import sys
#main_program
start_time_program=time.process_time()#To calculate the program execution time
reclimit=5000#This value increases the recursion limit (default: 1000). For bigger numbers, this value has to increase (with caution, may lead to crash).
sys.setrecursionlimit(reclimit)#Recursion limit
def Prime(num,n=2):
    if n>=num:
        print(num,"is a prime number.")
        return True
    elif num%n!=0:
        return Prime(num,n+1)
    else:
        print(num," is not a prime number, ",n," is a divisor of ",num,".",sep="")
        return False
number=abs(int(input("What number/integer do you want to check if it is a prime number or not? ")))
print("")
Prime(number)
print("")
print("Program execution time:",time.process_time()-start_time_program,"seconds.")
print("")
print("Goodbye.")
