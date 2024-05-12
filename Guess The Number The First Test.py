"""
Created on Wed Oct 9 13:57:44 2019

@author: Lucas Romero Fernández
"""
import random
import numpy as np
#main_program
a=100#Range of numbers that can be chosen.
subject_number=random.randint(0,1000)
secret_number=random.randint(1,a)
NAttempts=0
Attempts=[]
print("Hello Subject",subject_number,"...")
print("")
print("You must have questions about your current situation, do not worry, everything will reveal itself in time...")
print("")
print("What matters now is that you have been brought here to take part in a series of tests, the objective of these cannot be disclosed but if you succeed in all of them, there will be a reward for you...")
print("")
print("If you fail any of them... Better to not describe what happens...")
print("")
if subject_number<=200:
    print("You are one of the first subjects to try these tests, because of that, we have high hopes in you, so, try to not disapoint, okay?")
if subject_number>200 and subject_number<=600:
    print("There has been quite a few test subjects before you by now and they all have failed, what a tragedy...")
    print("")
    print("We do not expect much from you, but, your success could be a welcome surprise.")
if subject_number>600:
    print("You are one of the many test subjects that we have had, we expect failure by this point and you are not an exception...")
    print("")
    print("To be honest, I do not know why these tests keep being performed, we should get rid of all of you...")
print("")
print("Lets not waste any more time, the first test is some sort of a game, just to warm up and see if you have not been too much 'damaged' by the Incident, a random whole number between 1 and",a,"has been selected.")
print("")
print("Try to guess it in as few attempts as possible.")
print("")
print("I will tell you if the random number chosen is bigger or smaller than the number you presented, to help you out.")
print("")
print("Lets begin, shall we?")
print("")
while True:
    l=int(input("What whole number it is? "))
    print("")
    if (l<1) or (l>100):
        print("Write a whole number within the given range, please.")
        print("")
        Attempts.append(l)
        NAttempts+=1
        continue
    if secret_number==l:
        Attempts.append(l)
        NAttempts+=1
        print("Congratulations,",secret_number,"was the selected number.")
        print("")
        print("Lets see, the number of attempts has been",NAttempts,"...")
        print("")
        if NAttempts<=np.log2(a)+1:#Maximum number of attempts using the binary search algorithm, https://en.m.wikipedia.org/wiki/Binary_search_algorithm
            print("In the range of optimal number of attempts, good.")
        if NAttempts>np.log2(a)+1:
            print("It did not go as well as it could have gone... Maybe the Incident did affect you more than expected or you are not that bright after all...")
        print("")
        print("Also, you have tried the numbers",Attempts,", interesting...")
        break
    elif secret_number>l:
        print("Try again, it is bigger.")
        print("")
        NAttempts+=1
        Attempts.append(l)
    elif secret_number<l:
        print("Try again, it is smaller.")
        print("")
        NAttempts+=1
        Attempts.append(l)
print("")
print("As you can see, failure was not possible in this test, do not get overconfident, the rest will be more difficult and you can fail in them...")
print("")
print("Rest for now, tomorrow we will continue...")