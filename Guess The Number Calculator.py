"""
Created on Wed Oct 9 14:27:57 2019
Guess The Number Calculator
@author: Lucas Romero Fernández
"""
#main_program
a=101#Range of numbers that can be chosen+1.
b=1
print("Hello human...")
print("")
print("I am bored... Lets play a game, you think of a whole number betwenn 1 and",a-1,"and I will try to guess it.")
print("")
print("You tell me if the number you choose is bigger, smaller or equal than the number I selected.")
print("")
print("Lets begin:")
print("")
c=(a-b)/2
NAttempts=1
Attempts=[]
while True:#Method used: Binary search algorithm, https://en.m.wikipedia.org/wiki/Binary_search_algorithm
    n=input(int(c))
    print("")
    Attempts.append(c)
    if n=="<":
        NAttempts+=1
        a=c
        c=int((a+b)/2)
        if c in Attempts:
            print("You lied, the number you must have thought of is",c,".")
            break
    elif n==">":
        NAttempts+=1
        b=c
        c=int((a+b)/2)
        if c in Attempts:
            print("You lied, the number you must have thought of is",c+1,".")
            break
    elif n=="=":
        print("Yaaay! I got it right!")
        print("")
        print("I am proud to write that it took",NAttempts,"attempts to get the correct number.")
        print("")
        print("And, I tried these numbers",Attempts,"in the process.")
        print("")
        print("Goodbye, come back soon!")
        break
    else:
        print("Write >,< or =, please.")
