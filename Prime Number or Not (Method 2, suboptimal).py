# -*- coding: utf-8 -*-
"""
Created on Tue Dec 3 17:50:03 2019
Is *blank* a Prime Number? (Method 2, suboptimal)
@author: Lucas Romero Fernández
"""
from timeit import timeit
import sys
#main_program
n=0
question="What number/integer do you want to check if it is a prime number or not? "
number=int(input(question))
divisors=[2,3,5,7,11,13,17,19]#Knowing, beforehand, that this numbers are prime numbers...
for i in divisors:
    if number<divisors[n]:
        n+=1
        continue
    if number>=divisors[n]:
        res_div=number/i
        if (res_div.is_integer()) and (res_div!=float(1)):
            print("The number ",number," is not a prime number, it is divisible by ", divisors[n],".",sep="")
            print("")
            print("Program execution time:",timeit(number=1),"seconds.")
            print("")
            sys.exit("Goodbye.")
        else:
            n+=1
print("")
print("The number",number,"is a prime number!")
print("")
print("Program execution time:",timeit(number=1),"seconds.")
print("")
print("Goodbye.")