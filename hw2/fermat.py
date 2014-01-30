# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 20:26:37 2014

@author: cauerswald
"""

def check_fermat(a,b,c,n):
    if n<2:
        pass
    else:
        if a**n + b**n == c**n:
            print "Holy smokes! Fermat was wrong!"
        else:
            print "No, it doens't work."

def you_check():
    a=int(raw_input("What is a?"))
    b=int(raw_input("What is b?"))
    c=int(raw_input("What is c?"))
    n=int(raw_input("What is n?"))
    check_fermat(a,b,c,n)