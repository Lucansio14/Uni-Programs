# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 09:25:06 2019
Space Commando (Text-based Game)
@author: Lucas Romero Fernández
"""
from random import random
from random import randrange
import sys
import time
import matplotlib.pyplot as plt
NumOptions=0
def Persuasion(p):
    x=random()
    y="You are good at getting what you want with diplomacy, you are a very eloquent person (+30 Persuasion)."
    z="You are terrible at persuading people, you never get what you want by talking (-30 Persuasion)."
    v="You are not known for your persuasion skills (your Persuasion skill is not affected)."
    if x>0.66:
        p=0.3
        return [y,p]
    elif x<0.33:
        p=-0.3
        return [z,p]
    elif x<0.66 and x>0.33:
        p=0
        return [v,p]
def Combat(c):
    x=random()
    y="You are a master of hand-to-hand combat and all types of guns, including spacecraft weaponry, a complete killing machine (+30 Combat)."
    z="You are weak, clumsy and leaden, you are not made for any type of physical confrontation or learning about combat technology (-30 Combat)."
    v="Your combat prowess in all fields is like that of any beginner special agent, it does not stand out (your Combat skill is not affected)."
    if x>0.66:
        c=0.3
        return [y,c]
    elif x<0.33:
        c=-0.3
        return [z,c]
    elif x<0.66 and x>0.33:
        c=0
        return [v,c]
def Stealth(s):
    x=random()
    y="You are very good at hiding and going unnoticed, you are like a shadow in the night (+30 Stealth)."
    z="They hear you from miles away and even the blind can see you, it is clear that stealth is not your forte (-30 Stealth)."
    v="You're not an ace at stealth but you're not terrible either (your Stealth skill is not affected)."
    if x>0.66:
        s=0.3
        return [y,s]
    elif x<0.33:
        s=-0.3
        return [z,s]
    elif x<0.66 and x>0.33:
        s=0
        return [v,s]
def InfiltrationEncounter(c,p,s):
    while True:
        question="How do you want to do inflitrate (*f*,*p* or *s*)? "
        answer=input(question)
        print("")
        if answer=="f":
            x=random()
            if NumOptions==1:
                print("Enemies roll:",x)
                print("")
                print("Player roll:",float(0.4+c))
                print("")
            if x<float(0.4+c):
                print("You decide to confront things head-on, disabling the comms of the ship with the exterior before dispatching quickly one by one every small attack ship of the USS Horizon. With the coast clear, you enter the USS Horizon ship itself.")
                break
            else:
                print("You fail to disable the ship's comms in a hasty attack to the ship. You are quickly surrounded and destroyed.")
                print("")
                sys.exit("GAME OVER. Reason: Subject has died.")
        elif answer=="p":
            x=random()
            if NumOptions==1:
                print("Enemies roll:",x)
                print("")
                print("Player roll:",float(0.4+p))
                print("")
            if x<float(0.4+p):
                print("You approach the ship in a normal manner, getting contacted via intercom by the, especially lazy and uninterested, ship's personnel. You masterfully persuade them that you are a new technician for the ship (fake ID and all) and that you were late due to family business and problems with this personal spacecraft, they are convinced by that and let you through. You enter the ship unimpeded.")
                break
            else:
                print("They are not convinced by your lies, you try to escape but small fighter ships destroy your propulsion systems. With an evil intent, they let you to die in your now disabled and helpless ship.")
                print("")
                sys.exit("GAME OVER. Reason: Subject lost in space.")
        elif answer=="s":
            x=random()
            if NumOptions==1:
                print("Enemies roll:",x)
                print("")
                print("Player roll:",float(0.4+s))
                print("")
            if x<float(0.4+s):
                print("You activate the invisibility shield and radar inhibitors in your ship, skillfully dodging the small fighter ships patrols detection fields. After you surpass all the patrols, you enter the USS Horizon ship itself.")
                break
            else:
                print("You activate the invisibility shield and radar inhibitors in your ship, but you are not able to dodge the detection fields of the patrols, they discover you and get shoot down.")
                print("")
                sys.exit("GAME OVER. Reason: Subject has died.")
        else:
            print("You have to type *f*,*p* or *s* to choose between fighting, persuading or sneaking past, respectively: ")
            print("")
            continue
def LandingZone(x):
    x=random()
    v="one of the ship's corridors."
    z="a storeroom."
    y="the soldiers' barracks."
    l="the ship's dining room."
    if x>0.75:
        return z
    elif x<0.25:
        return v
    elif x<0.75 and x>0.50:
        return y
    elif x<0.50 and x>0.25:
        return l
def Encounters(c,p,s,Meds):
    print("Tarkanian soldiers have appeared!")
    print("")
    while True:
        question="What do you want to do (*f*,*p* or *s*)? "
        answer=input(question)
        print("")
        if answer=="f":
            x=random()
            if NumOptions==1:
                print("Enemies roll:",x)
                print("")
                print("Player roll:",float(0.3+c))
                print("")
            if x<float(0.3+c):
                print("You successfully neutralize the threats, they never stood a chance.")
                break
            else:
                print("You fail to incapacitate them in time before they discover you, they shoot you and wound you.")
                print("")
                if Meds==0:
                    print("You have no medicine to patch you up, you bleed to death aboard the USS Horizon, without having completed your mission.")
                    print("")
                    sys.exit("GAME OVER. Reason: Subject has died.")
                else:
                    print("You manage to hide and heal yourself with the medicine you have and, luckily, escape the situation altogether.")
                    Meds+=-1
                    break
        elif answer=="p":
            x=random()
            if NumOptions==1:
                print("Enemies roll:",x)
                print("")
                print("Player roll:",float(0.3+p))
                print("")
            if x<float(0.3+p):
                print("You talk your way through and convince them that you are the ship's new superluminal engine technician, you do not know if that position really exists, but they believe it and let you pass.")
                break
            else:
                print("They see through your bluff, they start shooting, wounding you in the process.")
                print("")
                if Meds==0:
                    print("You have no medicine to patch you up, you bleed to death aboard the USS Horizon, without having completed your mission.")
                    print("")
                    sys.exit("GAME OVER. Reason: Subject has died.")
                else:
                    print("You manage to hide and heal yourself with the medicine you have and, luckily, escape the situation altogether.")
                    Meds+=-1
                    break
        elif answer=="s":
            x=random()
            if NumOptions==1:
                print("Enemies roll:",x)
                print("")
                print("Player roll:",float(0.3+s))
                print("")
            if x<float(0.3+s):
                print("You go unnoticed and escape from the room without them noticing your presence, like a good space spy.")
                break
            else:
                print("They discover you while you are trying to get away from the situation, lasers rain down on you. You pray that one does not hit you, but you get shot anyway.")
                print("")
                if Meds==0:
                    print("You have no medicine to patch you up, you bleed to death aboard the USS Horizon, without having completed your mission.")
                    print("")
                    sys.exit("GAME OVER. Reason: Subject has died.")
                else:
                    print("You manage to hide and heal yourself with the medicine you have and, luckily, escape the situation altogether.")
                    Meds+=-1
                    break
        else:
            print("You have to type *f*,*p* or *s* to choose between fighting, persuading or sneaking past, respectively: ")
            print("")
            continue
    return Meds
def DoorCommandRoomEncounter(c,p,s,Meds):
    while True:
        question="What do you want to do (*f*,*p* or *s*)? "
        answer=input(question)
        print("")
        if answer=="f":
            x=random()
            if NumOptions==1:
                print("Enemies roll:",x)
                print("")
                print("Player roll:",float(0.3+c))
                print("")
            if x<float(0.3+c):
                print("You cause a small noise near them and while they are distracted, you rush them and quickly incapacitate the three of them, taking them by surprise. You are a little bit surprised by how well that went.")
                break
            else:
                print("You try to rush them while they seem specially distracted bickering between each other, but the numbers advantage quickly makes its presence known and you are shoot and wounded.")
                print("")
                if Meds==0:
                    print("You have no medicine to patch you up, you bleed to death aboard the USS Horizon, without having completed your mission.")
                    print("")
                    sys.exit("GAME OVER. Reason: Subject has died.")
                else:
                    print("You manage to hide and heal yourself with the medicine you have and, luckily, escape the situation altogether.")
                    Meds+=-1
                    break
        elif answer=="p":
            x=random()
            if NumOptions==1:
                print("Enemies roll:",x)
                print("")
                print("Player roll:",float(0.3+p))
                print("")
            if x<float(0.3+p):
                print("You discover a ship's crew suit left behind nearby, you put it on and you approach the guards in a hurry, convincing them that a catastrophic failure in the black hole engine has occurred and that assistance is needed. In the turmoil of the moment, you stay behind to 'warn' the captain...")
                break
            else:
                print("They see through your bluff by recognizing that your face and way of talking does not ressemble any of the crew members of the ship, you try to escape the situation but they start shooting, wounding you in the process.")
                print("")
                if Meds==0:
                    print("You have no medicine to patch you up, you bleed to death aboard the USS Horizon, without having completed your mission.")
                    print("")
                    sys.exit("GAME OVER. Reason: Subject has died.")
                else:
                    print("You manage to hide and heal yourself with the medicine you have and, luckily, escape the situation altogether.")
                    Meds+=-1
                    break
        elif answer=="s":
            x=random()
            if NumOptions==1:
                print("Enemies roll:",x)
                print("")
                print("Player roll:",float(0.3+s))
                print("")
            if x<float(0.3+s):
                print("You sabotage the wiring and machinery of a nearby room, causing a huge racket and gaining their attention in the process. You sneak past them while they are investigating the noise and arrive undetected to the door.")
                break
            else:
                print("You produce a considerable noise distraction nearby but you are discovered while you are trying to sneak past to the door, you make a run to the door but you get shot before you get to it.")
                print("")
                if Meds==0:
                    print("You have no medicine to patch you up, you bleed to death aboard the USS Horizon, without having completed your mission.")
                    print("")
                    sys.exit("GAME OVER. Reason: Subject has died.")
                else:
                    print("You manage to hide and heal yourself with the medicine you have and, luckily, escape the situation altogether.")
                    Meds+=-1
                    break
        else:
            print("You have to type *f*,*p* or *s* to choose between fighting, persuading or sneaking past, respectively: ")
            print("")
            continue
    return Meds
def RedFoxEncounter(c,p,s,Meds):
    print("You encounter Captain Red Fox!")
    print("")
    print("Red Fox seems to not have noticed your presence yet...")
    print("")
    while True:
        question="What do you want to do (*f*,*p* or *s*)? "
        answer=input(question)
        print("")
        if answer=="f":
            x=random()+0.15
            if NumOptions==1:
                print("Red Fox roll:",x)
                print("")
                print("Player roll:",float(0.4+c))
                print("")
            if x<float(0.4+c):
                print("While Red Fox seems distracted, you vault over the command table and attack Red Fox, iniciating a fight to the death that lasts minutes, ending in a struggle where you emerge victorious by desarming and stabbing Red Fox in the neck. After that, you pause for a second to recompose yourself.")
                break
            else:
                print("While Red Fox seems distracted, you vault over the command table and attack Red Fox, iniciating a fight that lasts minutes, ending in a struggle that you lose, getting severely stabbed. ")
                print("")
                if Meds==0:
                    print("You have no medicine to patch you up, you bleed to death aboard the USS Horizon, without having completed your mission. The last thing you see and hear is Red Fox laughing evilly.")
                    print("")
                    sys.exit("GAME OVER. Reason: Subject has died.")
                else:
                    print("When you are about to heal yourself with medicine, Red Fox kicks it away. You bleed to death aboard the USS Horizon, without having completed your mission. The last thing you see and hear is Red Fox laughing evilly.")
                    print("")
                    sys.exit("GAME OVER. Reason: Subject has died.")
        elif answer=="p":
            x=random()+0.2
            if NumOptions==1:
                print("Red Fox roll:",x)
                print("")
                print("Player roll:",float(0.4+p))
                print("")
            if x<float(0.4+p):
                print("You miraculously convince Red Fox that you are an ally, betraying the ESA by joining the crew in the USS Horizon. After a while, when you are given your own combat knife and Red Fox is distracted, you murder him in cold blood.")
                break
            else:
                print("Red Fox sees through your lies and persuasion without you noticing it. Nevertheless, Red Fox goes along with it and when you think you are safe, Red Fox stabs you in the stomach.")
                print("")
                if Meds==0:
                    print("You have no medicine to patch you up, you bleed to death aboard the USS Horizon, without having completed your mission. The last thing you see and hear is Red Fox laughing evilly.")
                    print("")
                    sys.exit("GAME OVER. Reason: Subject has died.")
                else:
                    print("When you are about to heal yourself with medicine, Red Fox kicks it away. You bleed to death aboard the USS Horizon, without having completed your mission. The last thing you see and hear is Red Fox laughing evilly.")
                    print("")
                    sys.exit("GAME OVER. Reason: Subject has died.")
        elif answer=="s":
            x=random()-0.1
            if NumOptions==1:
                print("Red Fox roll:",x)
                print("")
                print("Player roll:",float(0.4+s))
                print("")
            if x<float(0.4+s):
                print("You sneak behind and take the knife of Red Fox out of the sheath. Then, you cut the calf muscles and stab beneath the chin, killing Red Fox in the process.")
                break
            else:
                print("You sneak behind and try to take the knife that Red Fox has. Unfortunately, Red Fox hears you and unsheathes the knife at lightning speed, trying to end your life with a single cut to the neck. you try to block it but it is too fast...")
                print("")
                if Meds==0:
                    print("You have no medicine to patch you up, you bleed to death aboard the USS Horizon, without having completed your mission. The last thing you see and hear is Red Fox laughing evilly.")
                    print("")
                    sys.exit("GAME OVER. Reason: Subject has died.")
                else:
                    print("When you are about to heal yourself with medicine, Red Fox kicks it away. You bleed to death aboard the USS Horizon, without having completed your mission. The last thing you see and hear is Red Fox laughing evilly.")
                    print("")
                    sys.exit("GAME OVER. Reason: Subject has died.")
        else:
            print("You have to type *f*,*p* or *s* to choose between fighting, persuading or sneaking past, respectively: ")
            print("")
            continue
def Places1(x):
    x=random()
    v="a ship's bathroom."
    z="the ship's armory."
    y="the commanders' bedrooms."
    l="the janitor's closet."
    if x>0.75:
        return z
    elif x<0.25:
        return v
    elif x<0.75 and x>0.50:
        return y
    elif x<0.50 and x>0.25:
        return l
def Places2(x):
    x=random()
    y="the ship's gym."
    z="Captain Red Fox's bedroom."
    v="the engine room."
    if x>0.66:
        return z
    elif x<0.33:
        return y
    elif x<0.66 and x>0.33:
        return v
def Objects1(x):
    x=random()
    y="an untagged GunTek1000 plasma pistol (+15 Combat)."
    z="a cutting-edge Tarkanian combat knife (+10 Combat)."
    v="a book on the art of diplomacy (+10 Persuasion)."
    l="nothing particularly useful."
    if x>0.90:
        r=0.15
        return y,r,0,0
    elif x<0.40:
        return l,0,0,0
    elif x<0.55 and x>0.40:
        r=0.1
        return z,r,0,0
    elif x<0.90 and x>0.55:
        d=0.1
        return v,0,d,0
def Objects2(x):
    x=random()
    y="a drug called -100% Space Politician-, you've heard of it, it's supposed to make you more elocuent (+10 Persuasion)."
    z="a ship's guard suit (+10 Stealth)."
    v="a small device called -Stealth Boy-, apparently it makes you invisible for a short time (+20 Stealth)."
    l="nothing particularly useful."
    if x>0.90:
        a=0.2
        return v,0,0,a
    elif x<0.30:
        return l,0,0,0
    elif x<0.55 and x>0.30:
        a=0.10
        return z,0,0,a
    elif x<0.90 and x>0.55:
        b=0.10
        return y,0,b,0
def Objects3(x):
    x=random()
    y="a bunch of untagged plasma grenades, very effective in combat scenarios (+10 Combat)."
    z="a book on how to go unnoticed by the naked eye (+10 Stealth)."
    v="a wad of space euros, useful for bribing guards (+15 Persuasion)."
    l="nothing particularly useful."
    if x>0.90:
        n=0.15
        return v,0,n,0
    elif x<0.30:
        return l,0,0,0
    elif x<0.55 and x>0.30:
        k=0.1
        return z,0,0,k
    elif x<0.90 and x>0.55:
        i=0.1
        return y,i,0,0
#main program
print("Hello! Welcome to my game SPACE COMMANDO, I am writing this to give you a few instructions that can help you during the game:")
print("")
print("1-SPACE COMMANDO is a game where luck and randomness are important factors, so there will be playthroughs where everything seems lost and unfair.")
print("")
print("2-Each playthrough you play will be different because randomness is a fundamental part of the game, every time you play the map will change, your objects will change, etc... All this to always give an -unique- experience each time you play.")
print("")
print("3-The game mechanics are simple, when you encounter a tricky situation, you will have three options, fight, persuade or sneak past. To choose them, you just have to write *f*, *p* or *s*, respectively.")
print("")
print("4-This game is not permissive and there are no checkpoints, it ends when you are injured and you don't have medicine or when you complete your mission and escape alive.")
print("")
print("5-Different actions in the game can unlock achievements. Try to get them all!")
print("")
question="First of all, type *1* in the terminal if you want to see the 'dice rolls'/numbers generated of the player and enemies in encounters: "
answer=input(question)
if answer=="1":
    NumOptions=1
print("")
question="Now, when you are ready to start the game, just type *i* in the terminal: "
answer=input(question)
print("")
lw=0
DumbMeter=0
while True:
    if answer=="i":
        w=0
        if DumbMeter>4 and DumbMeter<10:
            print("Achievement unlocked: Wake-up call -> You have changed your attitude.")
            print("")
            lw+=1
        print("Context: You are a rookie special agent recently hired by the Espionage Space Agency (ESA) to infiltrate the USS Horizon to steal documents of vital importance to the ESA and the human colony on Mars as your first real mission.")
        print("")
        print("Your mission will not be easy, the USS Horizon is well guarded, there are armed Tarkanian guards spread through out the ship, but the most dangerous is its captain, Red Fox, a ruthless being capable of anything to reach his plans for universal domination.")
        print("")
        print("How you complete your mission is up to you, do whatever is necessary to steal those documents and if possible, take down Red Fox in a subtle and inadvertent way, Red Fox is a nuisance for the ESA but you cannot risk compromising the agency.")
        print("")
        print("You will be sent alone and unarmed, they cannot provide you with equipment due to lack of time and personnel, you will have to use your ingenuity and what you find to achieve your objective. Remember that the ranged weapons of the Tarkanian guards while on patrol are DNA-tagged and can only be used by themselves, you will have to find untagged ones.")
        print("")
        question="Ready to begin the mission? "
        answer=input(question)
        print("")
        if answer=="yes" or answer=="Yes" or answer=="YES":
            print("That's the spirit!")
            print("")
            q=1
            print("Achievement unlocked: Positivity -> Always yes to everything.")
            print("")
            lw+=1
        else:
            q=0
            print("Ready or not, it does not matter, there is no time to spare.")
            print("")
        print("Good luck, agent.")
        print("")
        print("**********************COMMENCING MISSION***************************")
        print("")
        break
    if answer=="p" or answer=="f" or answer=="s":
        if DumbMeter>4 and DumbMeter<10:
            print("Achievement unlocked: Wake-up call -> You have changed your attitude.")
            print("")
            lw+=1
        print("You are really eager to start, good...")
        w=1
        print("")
        print("Achievement unlocked: Eager to play -> Be complimented for your impatience at the beginning of the game.")
        print("")
        lw+=1
        print("Context: You are a rookie special agent recently hired by the Espionage Space Agency (ESA) to infiltrate the USS Horizon to steal documents of vital importance to the ESA and the human colony on Mars as your first real mission.")
        print("")
        print("Your mission will not be easy, the USS Horizon is well guarded, there are armed Tarkanian guards spread through out the ship, but the most dangerous is its captain, Red Fox, a ruthless being capable of anything to reach his plans for universal domination.")
        print("")
        print("How you complete your mission is up to you, do whatever is necessary to steal those documents and if possible, take down Red Fox in a subtle and inadvertent way, Red Fox is a nuisance for the ESA but you cannot risk compromising the agency.")
        print("")
        print("You will be sent alone and unarmed, they cannot provide you with equipment due to lack of time and personnel, you will have to use your ingenuity and what you find to achieve your objective.")
        print("")
        question="Ready to begin the game? "
        answer=input(question)
        print("")
        if answer=="yes" or answer=="Yes" or answer=="YES":
            print("That's the spirit!")
            print("")
            q=1
            print("Achievement unlocked: Positivity -> Always yes to everything.")
            print("")
            lw+=1
        else:
            print("Ready or not, it does not matter, there is no time to spare.")
            print("")
        print("Good luck, agent.")
        print("")
        print("**********************COMMENCING MISSION****************************")
        print("")
        break
    else:
        answer=input("You have to type *i* to start the game: ")
        print("")
        DumbMeter+=1
        if DumbMeter==5:
            print("Are you dumb?")
            print("")
        if DumbMeter==10:
            print("Understood, if you are going to be a nuisance, it is better for you to not proceed.")
            print("")
            sys.exit("GAME OVER. Reason: Subject stupidness.")
t1=time.time()
print("On the way to the ship, you wonder why they have entrusted such an important mission to a beginner like you and why they hired you in the first place, you have not done any type of training or preparation with them before this.")
print("")
print("From what you could hear and find out, apparently they have been observing you for some time and have seen that you have a lot of potential to be a good spy, although you can't come up with an explanation as to why this is so.")
print("")
print("You have safely arrived at the position of the USS Horizon, now you have to find a way to infiltrate the ship...")
print("")
print("Before trying anything with your small infiltration ship, you decide to give a quick check on your most important skills as a spy:")
print("")
c=Combat(0)
print("1-",c[0])
print("")
p=Persuasion(0)
print("2-",p[0])
print("")
s=Stealth(0)
print("3-",s[0])
print("")
InfiltrationEncounter(c[1],p[1],s[1])
print("")
print("In the landing zone, you luckily find medicine, this will heal you up, in case you are physically injured.")
Meds=1
print("")
print("When you leave undetected the mostly lifeless landing zone, you find yourself in",LandingZone(0))
print("")
Object1=Objects1(0)
print("You search the room and find",Object1[0])
c[1]+=Object1[1]
p[1]+=Object1[2]
s[1]+=Object1[3]
print("")
Meds=Encounters(c[1],p[1],s[1],Meds)
print("")
print("Medicine left =",Meds)
print("")
print("You continue moving through the ship in search of the command room, where you assume the documents you are looking for will be and with a little luck, captain Red Fox himself.")
print("")
print("You open a door and enter",Places1(0))
print("")
Object2=Objects2(0)
print("You search the room and find",Object2[0])
c[1]+=Object2[1]
p[1]+=Object2[2]
s[1]+=Object2[3]
print("")
Meds=Encounters(c[1],p[1],s[1],Meds)
print("")
print("Medicine left =",Meds)
print("")
print("You feel that with every step you take, you get closer to the command room and your goals. You wonder why you don't meet more Tarkanian soldiers in the rooms you visit, it seems as if they are all in a specific place on the ship that you haven't found yet and that you hope not to find.")
print("")
print("You open the door in front of you and enter",Places2(0))
print("")
Object3=Objects3(0)
print("You search the room and find",Object3[0])
c[1]+=Object3[1]
p[1]+=Object3[2]
s[1]+=Object3[3]
print("")
Meds=Encounters(c[1],p[1],s[1],Meds)
print("")
print("Medicine left =",Meds)
print("")
print("You finally reach the command room, but there are three Tarkanian guards at the door, bickering between each other. You have to get rid of them somehow...")
print("")
Meds=DoorCommandRoomEncounter(c[1],p[1],s[1],Meds)
print("")
print("Medicine left =",Meds)
print("")
print("With the now unprotected door in front of you and a little grin on your face, you step into the command room.")
print("")
RedFoxEncounter(c[1],p[1],s[1],Meds)
print("")
print("With Red Fox lying dead, you search for the important documents, finding them shortly afterwards. Following that, you leave the USS Horizon alive and well.")
print("")
print("**********************SIMULATION COMPLETED***************************")
print("")
print("Well done, you have completed your training and initiation, welcome to ESA, Agent.")
print("")
print("Achievement unlocked: Simulated Universe Savior -> Simulation completed.")
print("")
lw+=1
print("----->CONGRATULATIONS! YOU HAVE COMPLETED SPACE COMMANDO!<-------")
print("")
t2=time.time()
print("Duration of the playthrough =",t2-t1,"seconds.")
print("")
if t2-t1<60:
    print("Achievement unlocked: Efficient Agent -> Finish the game in less than a minute.")
    print("")
    lw+=1
if t2-t1<47:
    print("Achievement unlocked: Agent 47 -> Finish the game in less than 47 seconds.")
    print("")
    lw+=1
if Meds!=0:
    print("Achievement unlocked: Untouchable -> Finish the game without using medicine.")
    print("")
    lw+=1
print("*******************ACHIEVEMENTS UNLOCKED************************")
print("")
Tot_Ach=7
print("Achievement unlocked: Simulated Universe Savior -> Simulation completed.")
print("")
if DumbMeter>4 and DumbMeter<10:
    print("Wake-up call -> You have changed your attitude.")
    print("")
if q==1:
    print("Positivity -> Always yes to everything.")
    print("")
if w==1:
    print("Eager to play -> Be complimented for your impatience at the beginning of the game.")
    print("")
if t2-t1<60:
    print("Efficient Agent -> Finish the game in less than a minute.")
    print("")
if t2-t1<47:
    print("Agent 47 -> Finish the game in less than 47 seconds.")
    print("")
if Meds!=0:
    print("Untouchable -> Finish the game without using medicine.")
    print("")
print("Achievements Unlocked =",lw)
print("")
print("Total Achievements =",Tot_Ach)
print("")
if lw==Tot_Ach:
    print("Congratulations! You have completed all the achievements in Space Commando! You should feel proud of yourself and maybe a little bit worried.")
else:
    x=["A. Unlocked","A. Locked"]
    y=[lw,Tot_Ach-lw]
    fig,Ej0=plt.subplots(figsize=(6,6))
    Ej0.bar(x,y,align="center",width=0.7,alpha=1)
    Ej0.set_title("Achievements")
    fig.tight_layout()
    plt.show()
print("**********************STATISTICS***************************")
print("")
y=[randrange(100),randrange(100),randrange(100),randrange(100),randrange(100),randrange(100),randrange(100)]
x=["Infiltration","Ship encounter 1","Ship encounter 2","Ship encounter 3","Command room door","Captain Red Fox","Playthrough finished"]
fig,h0=plt.subplots(figsize=(11,6))
h0.bar(x,y,align="center",width=0.7,alpha=0.7,color="g")
h0.set_title("Player Statistics (Deaths or playthroughs finished)")
fig.tight_layout()
plt.show()