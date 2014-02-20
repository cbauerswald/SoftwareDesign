# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 17:52:54 2014

@author: pruvolo
"""

def cumulative_sum(L):
    """ returns a list where each element in the returned list is the
    	cumulative sum of the elements up to the corresponding element in
	the original list.
	
	L: the original list
	returns: a new list where element i is equal to the sum of element
		 0 through i in the original list """
    sum_L=[L[0]]
    for i in range(len(L)-1):
        print i
        sum_L.append(L[i+1] + sum_L[i])
    return sum_L

if __name__ == '__main__':
    print cumulative_sum([4, 2, 3,1])
