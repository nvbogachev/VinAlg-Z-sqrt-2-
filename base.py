from __future__ import division
from itertools import product
import numpy as np
#import math
from math import *

class Q2(object):
    def __init__(self,a,b):
        self.a=a
        self.b=b
    def __getitem__(self, key):
        if key==0:
            return self.a
        else:
            return self.b
    def __repr__(self):
        return str((self.a,self.b))
    def __mul__(v1,v2):
        return Q2( v1[0]*v2[0] + 2*v1[1]*v2[1], v1[0]*v2[1] + v1[1]*v2[0] )    
    def __add__(v1,v2):
        return Q2( v1[0]+v2[0], v1[1]+v2[1])    
    def __neg__(self):
        return Q2(-self.a, -self.b)
    def __sub__(v1,v2):
        return Q2( v1[0]-v2[0], v1[1]-v2[1])
    def __truediv__(v1,v2):
        return (Q2( (v1[0]*v2[0] - 2*v1[1]*v2[1])/(v2[0]*v2[0] - 2*v2[1]*v2[1]), (- v1[0]*v2[1] + v1[1]*v2[0])/(v2[0]*v2[0] - 2*v2[1]*v2[1]) ))
    def square(v):
        return v*v
    def real(v):
        return v[0]+v[1]*sqrt(2)
    def antireal(v):
        return v[0]-v[1]*sqrt(2)

def div(v1,v2):
    return (Q2( (v1[0]*v2[0] - 2*v1[1]*v2[1])/(v2[0]*v2[0] - 2*v2[1]*v2[1]), (- v1[0]*v2[1] + v1[1]*v2[0])/(v2[0]*v2[0] - 2*v2[1]*v2[1]) ))

def divide(v1,v2):
    return((v1[0]*v2[0] - 2*v1[1]*v2[1])%(v2[0]*v2[0] - 2*v2[1]*v2[1]) == 0 and (- v1[0]*v2[1] + v1[1]*v2[0])%(v2[0]*v2[0] - 2*v2[1]*v2[1]) == 0 )




#print(minor(d,d))
    
def s(minor_,k):
    return floor(sqrt(abs(minor_-k)))




def Transpose(m):
    return map(list,zip(*m))

def Minor(m,i,j):
    return [row[:j] + row[j+1:] for row in (m[:i]+m[i+1:])]

def Det(m):
    #base case for 2x2 matrix
    if len(m) == 2:
        return m[0][0]*m[1][1]-m[0][1]*m[1][0]
    
    determinant = Q2(0,0)
    for c in range(len(m)):
        determinant = determinant + Q2((-1)**c,0)*m[0][c]*Det(Minor(m,0,c))
    return determinant

def Inv(m):
    determinant = Det(m)
    #print(Det(m))
    #assert(Det(m)!=Q2(0,0))
    #special case for 2x2 matrix:
    if len(m) == 2:
        return [[m[1][1]/determinant, -1*m[0][1]/determinant],
                [-1*m[1][0]/determinant, m[0][0]/determinant]]
    
    #find matrix of cofactors
    cofactors = []
    for r in range(len(m)):
        cofactorRow = []
        for c in range(len(m)):
            minor = Minor(m,r,c)
            cofactorRow.append(Q2((-1)**(r+c),0) * Det(minor))
        cofactors.append(cofactorRow)
    cofactors = Transpose(cofactors)
    for r in range(len(cofactors)):
        for c in range(len(cofactors)):
            #print(determinant)
            cofactors[r][c] = cofactors[r][c]/determinant
    return cofactors

def Mult(A, B):
    C = [[Q2(0,0) for row in range(len(A))] for col in range(len(B[0]))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                C[i][j] = C[i][j] + A[i][k]*B[k][j]
    return C


