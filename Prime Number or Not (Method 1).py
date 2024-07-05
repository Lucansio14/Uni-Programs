# -*- coding: utf-8 -*-
"""
Created on Tue Aug 22 18:29:42 2023
Is *blank* a Prime Number? (Method 1)
@author: Lucas Romero Fernández
"""
from timeit import timeit
import sys
#main_program
reclimit=5000#This value increases the recursion limit (default: 1000). For bigger numbers, this value has to increase (with caution, may lead to crash).
sys.setrecursionlimit(reclimit)
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
print("Program execution time:",timeit(number=1),"seconds.")
print("")
print("Goodbye.")