"""
Created on Wed Sep 25 14:02:59 2019
Cheap Psychologist (The process of Happiness)
@author: Lucas Romero Fernández
"""
#main_program
print("“There is no path to happiness; happiness is the path.” Buddha")
print("")
print("Hello...")
print("")
print("Take a sit...")
print("")
print("Before we start with the questions, try to be clear, concise and upfront with your answers, this is the best way to make this process deliver the most positive results.")
print("")
print("Lets start with the the only and most important question of them all:")
print("")
Positives=["yes","Yes","YES","Yeah","yeah","YEAH"]
Negatives=["No","no","NO","Nope","nope","NOPE"]
PosibAnswers=[]
PosibAnswers.extend(Positives)
PosibAnswers.extend(Negatives)
PosibChanges=["Change something about your life, then.","Try changing even more things about your life.","Keep changing things about your life, then.","Then, you can try changing everything about your life, as a last resort."]
def changes(x):
    change=PosibChanges[0+x]
    return change
while True:
    question="Do you feel happy at this moment in your life? "
    answer=input(question)
    if answer in PosibAnswers:
        if answer in Positives:
            print("")
            print("Great, keep doing what you do then, you do not need this.")
            print("")
            print("Live well and enjoy the happy life you made for yourself.")
            print("")
            print("Goodbye.")
            break
        elif answer in Negatives:
            print("")
            print(changes(0))
            print("")
            question="Have you done it? "
            answer=input(question)
            if answer in Positives:
                print("")
                question="Do you feel happy now? "
                answer=input(question)
                if answer in Positives:
                    print("")
                    print("Great, keep doing what you do then, you do not need this anymore.")
                    print("")
                    print("Live well and enjoy the happy life you made for yourself.")
                    print("")
                    print("Goodbye.")
                    break
                elif answer in Negatives:
                    print("")
                    print(changes(1))
                    print("")
                    question="Do you feel happy now? "
                    answer=input(question)
                    if answer in Positives:
                        print("")
                        print("Great, keep doing what you do then, you do not need this anymore.")
                        print("")
                        print("Live well and enjoy the happy life you made for yourself.")
                        print("")
                        print("Goodbye.")
                        break
                    elif answer in Negatives:
                        print("")
                        print(changes(2))
                        print("")
                        question="Do you feel happy now? "
                        answer=input(question)
                        if answer in Positives:
                            print("")
                            print("Great, keep doing what you do then, you do not need this anymore.")
                            print("")
                            print("Live well and enjoy the happy life you made for yourself.")
                            print("")
                            print("Goodbye.")
                            break
                        elif answer in Negatives:
                            print("")
                            print(changes(3))
                            print("")
                            question="Do you feel happy now? "
                            answer=input(question)
                            if answer in Positives:
                                print("")
                                print("Great, keep doing what you do then, you do not need this anymore.")
                                print("")
                                print("Live well and enjoy the happy life you made for yourself.")
                                print("")
                                print("Goodbye.")
                                break
                            elif answer in Negatives:
                                print("")
                                print("Then, keep changing aspects in your life until you feel better than before the changes. When you feel like you have done it, be proud of yourself and if you want to, come back here and tell me about the experience, I would gladly appreciate it...")
                                print("")
                                print("Goodbye.")
                                break
                            else:
                               print("")
                               print("Write any of this answers ",PosibAnswers,", please.",sep="")
                               print("")
                               print("After this, it is better if we go back to the beginning...")
                               print("")
                               continue
                        else:
                           print("")
                           print("Write any of this answers ",PosibAnswers,", please.",sep="")
                           print("")
                           print("After this, it is better if we go back to the beginning...")
                           print("")
                           continue
                    else:
                       print("")
                       print("Write any of this answers ",PosibAnswers,", please.",sep="")
                       print("")
                       print("After this, it is better if we go back to the beginning...")
                       print("")
                       continue
            elif answer in Negatives:
                print("")
                print("What are you waiting for, then? Go and change something!")
                print("")
                question="Have you done it? "
                answer=input(question)
                if answer in Positives:
                    print("")
                    question="Do you feel happy now? "
                    answer=input(question)
                    if answer in Positives:
                        print("")
                        print("Great, keep doing what you do then, you do not need this anymore.")
                        print("")
                        print("Live well and enjoy the happy life you made for yourself.")
                        print("")
                        print("Goodbye.")
                        break
                    elif answer in Negatives:
                        print("")
                        print(changes(1))
                        print("")
                        question="Do you feel happy now? "
                        answer=input(question)
                        if answer in Positives:
                            print("")
                            print("Great, keep doing what you do then, you do not need this anymore.")
                            print("")
                            print("Live well and enjoy the happy life you made for yourself.")
                            print("")
                            print("Goodbye.")
                            break
                        elif answer in Negatives:
                            print("")
                            print(changes(2))
                            print("")
                            question="Do you feel happy now? "
                            answer=input(question)
                            if answer in Positives:
                                print("")
                                print("Great, keep doing what you do then, you do not need this anymore.")
                                print("")
                                print("Live well and enjoy the happy life you made for yourself.")
                                print("")
                                print("Goodbye.")
                                break
                            elif answer in Negatives:
                                print("")
                                print(changes(3))
                                print("")
                                question="Do you feel happy now? "
                                answer=input(question)
                                if answer in Positives:
                                    print("")
                                    print("Great, keep doing what you do then, you do not need this anymore.")
                                    print("")
                                    print("Live well and enjoy the happy life you made for yourself.")
                                    print("")
                                    print("Goodbye.")
                                    break
                                elif answer in Negatives:
                                    print("")
                                    print("Then, keep changing aspects in your life until you feel better than before the changes. When you feel like you have done it, be proud of yourself and if you want to, come back here and tell me about the experience, I would gladly appreciate it...")
                                    print("")
                                    print("Goodbye.")
                                    break
                                else:
                                   print("")
                                   print("Write any of this answers ",PosibAnswers,", please.",sep="")
                                   print("")
                                   print("After this, it is better if we go back to the beginning...")
                                   print("")
                                   continue
                            else:
                                print("")
                                print("Write any of this answers ",PosibAnswers,", please.",sep="")
                                print("")
                                print("After this, it is better if we go back to the beginning...")
                                print("")
                                continue
                elif answer in Negatives:
                    print("")
                    print("With that attitude, I cannot help you.")
                    print("")
                    print("Come back when you want to change.")
                    print("")
                    print("Goodbye.")
                    break
                else:
                   print("")
                   print("Write any of this answers ",PosibAnswers,", please.",sep="")
                   print("")
                   print("After this, it is better if we go back to the beginning...")
                   print("")
                   continue
            else:
               print("")
               print("Write any of this answers ",PosibAnswers,", please.",sep="")
               print("")
               print("After this, it is better if we go back to the beginning...")
               print("")
               continue 
    else:
        print("")
        print("Write any of this answers ",PosibAnswers,", please.",sep="")
        print("")
        continue