from __future__ import division
import itertools
from itertools import product
import numpy as np
from math import *
import coxiter
#import numpy.matrix




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
    def __sub__(v1,v2):
        return Q2( v1[0]-v2[0], v1[1]-v2[1])
    def __truediv__(v1,v2):
        return (Q2( (v1[0]*v2[0] - 2*v1[1]*v2[1])/(v2[0]*v2[0] - 2*v2[1]*v2[1]), (- v1[0]*v2[1] + v1[1]*v2[0])/(v2[0]*v2[0] - 2*v2[1]*v2[1]) ))
    def __neg__(self):
        return Q2(-self.a, -self.b)
    def square(v):
        return v*v
    def real(v):
        return v[0]+v[1]*sqrt(2)
    def antireal(v):
        return v[0]-v[1]*sqrt(2)


def div(v1,v2):
    return (Q2( (v1[0]*v2[0] - 2*v1[1]*v2[1])/(v2[0]*v2[0] - 2*v2[1]*v2[1]), (-v1[0]*v2[1] + v1[1]*v2[0])/(v2[0]*v2[0] - 2*v2[1]*v2[1]) ))

def divide(v1,v2):
    return((v1[0]*v2[0] - 2*v1[1]*v2[1])%(v2[0]*v2[0] - 2*v2[1]*v2[1]) == 0 and (- v1[0]*v2[1] + v1[1]*v2[0])%(v2[0]*v2[0] - 2*v2[1]*v2[1]) == 0 )
#        return True


d = Q2(3,4)
lr = [(1,0), (2,0), (2,1), (5,-1), (10,-2), (8,3)]
lr = [Q2(*l) for l in lr]
m = 0
M = 1
#print(d.square())
		
def minor(*arg):
	sum = - arg[0].square()*d
    
	for v in arg[1:]:
		sum = sum + v.square()
	return sum

e0 = [Q2(0,0), Q2(0,0), Q2(0,0), Q2(-1,0)]
e1 = [Q2(0,0), Q2(0,0), Q2(-1,0), Q2(1,0)]
e2 = [Q2(0,0), Q2(-1,0), Q2(1,0), Q2(0,0)]



#e1 = [(-1, 1), (1, 2), (1, 0), (-1, 1)]
#e2 = [(-1, 1), (3, 0), (1, 1), (1, 0)]

#e1 = [Q2(*l) for l in e1]
#e2 = [Q2(*l) for l in e2]

def inner(arg1,arg2):
    sum = - (arg1[0]*arg2[0])*d
    n = len(arg1)
    for j in range(1,n):
        sum = sum + arg1[j]*arg2[j]
    return sum

#print(inner(e1,e1).real())

#print(minor(d,d))
	
def s(minor_,k):
	return int(floor(sqrt(abs(minor_/2-k))))
	
def Box(vectors, k):
	M = minor(*vectors)
	s_plus=s(M.real(),k.real())
	s_minus=s(M.antireal(),k.real())
	return (Q2(a,b) for a in range(-s_minus, s_plus+s_minus + 3) for b in range(int(floor(-a/sqrt(2) - sqrt(2))), int(floor(s_plus - a/sqrt(2) + 3))) )


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



#print(minor(Q2(0,0), Q2(0,0), Q2(0,0), Q2(0,-1)))

#([(-1, 1), (1, 2), (1, 0), (-1, 1)], (10, -2), (-1, 1), 0.01196201155856884)

#([(-1, 1), (3, 0), (1, 1), (1, 0)], (10, -2), (-1, 1), 0.01196201155856884)


#for e in eee:
#    e = [Q2(*l) for l in e]
#print(len(eee))


eee = [[(2, 1), (4, 2), (4, 2), (2, 1)], #(2, 0), (2, 1), 2.914213562373095)
[(2, 1), (4, 4), (2, 1), (0, 0)], #(2, 0), (2, 1), 2.914213562373095)
       [(2, 1), (6, 3), (0, 0), (0, 0)]] #(2, 0), (2, 1), 2.914213562373095)


e3 = [(2, 1), (4, 4), (2, 1), (0, 0)]
e3 = [div(Q2(*l),Q2(0,1)) for l in e3]

e4 = [(2, 1), (4, 2), (4, 2), (2, 1)]
e4 = [div(Q2(*l),Q2(0,1)) for l in e4]

e5 = [(2, 1), (6, 3), (0, 0), (0, 0)]
e5 = [div(Q2(*l),Q2(0,1)) for l in e5]

e6 = [(2, 1), (4, 3), (2, 2), (2, 1)]
e6 = [div(Q2(*l),Q2(0,1)) for l in e6]

e7 = [(2, 2), (6, 4), (4, 3), (0, 1)]
e7 = [(Q2(*l)/Q2(0,1)) for l in e7]


e8 = [(2, 2), (6, 4), (4, 3), (0, 1)]
e8 = [div(Q2(*l),Q2(0,1)) for l in e8]

e9 = [(2, 2), (6, 5), (2, 2), (2, 1)]
e9 = [div(Q2(*l),Q2(0,1)) for l in e9]

e10 = [(2, 2), (4, 3), (4, 3), (4, 3)]
e10 = [div(Q2(*l),Q2(0,1)) for l in e10]

e11 = [(2, 2), (6, 4), (4, 3), (0, 0)]
e11 = [div(Q2(*l),Q2(0,1)) for l in e11]

e12 = [(4, 2), (6, 4), (6, 4), (6, 4)]
e12 = [div(Q2(*l),Q2(0,1)) for l in e12]

e13 = [(4, 2), (8, 6), (4, 3), (4, 3)]
e13 = [div(Q2(*l),Q2(0,1)) for l in e13]

e14 = [(4, 2), (8, 6), (6, 4), (0, 0)]
e14 = [div(Q2(*l),Q2(0,1)) for l in e14]

e15 = [(4, 3), (8, 7), (6, 4), (6, 4)]
e15 = [div(Q2(*l),Q2(0,1)) for l in e15]

e16 = [(4, 3), (8, 7), (8, 6), (0, 0)]
e16 = [div(Q2(*l),Q2(0,1)) for l in e16]

e17 = [(4, 3), (12, 8), (2, 2), (2, 2)]
e17 = [div(Q2(*l),Q2(0,1)) for l in e17]

e18 = [(6, 4), (16, 12), (4, 2), (4, 2)]
e18 = [div(Q2(*l),Q2(0,1)) for l in e18]




#for e in eee:
#    e = [Q2(*l) for l in e]
#    if(inner(e,e3).real() <= 0):
#        print e


#print(inner(e12, e12).real())




g = [e0,e1,e2,e3,e4,e5,e6,e7,e9,e10,e15,e16,e17]

#print(cos(pi/8)*cos(pi/8))

#print([[int(floor(inner(u,v).real())) for u in g] for v in g])

#M = np.matrix([[inner(u,v).real() for u in g] for v in g])

#print(M)



N = [[inner(u,v).real() for u in g] for v in g]

#print(len(N))

#for p in N:
#    print(p)

#for k in range(len(N)):
#for j in range(k):
#        if N[k][j] > 0:
#            print((k,j))

#print(Q2(1,1)-Q2(1,0))

g1 = [e0,e1,e4,e5]

N1 = [[inner(u,v) for u in g1] for v in g1]

#for u in g1:
#    print(u)

#print(Det(N1))

for w in N1:
    print(w)

#print(len(N1))

#print(coxiter.gram_matrix2graph(N1, 4))

#print(Q2(1,1)*Q2(1,1) + Q2(1,1))

#print(Mult(g1,Transpose(g1)))

#print(e7)

g2 = [e1,e0,e4,e5]
N2 = [[inner(u,v) for u in g2] for v in g2]

A1 = Transpose(g1)

A2 = Transpose(g2)

#print(A1)

print(Det(A1))
#print(Det(A2))

#print(Minor(A1,0,0))

#print(Det(Minor(A1,0,0)))


#print(e0)

#print(e1)

#print(e6)

#print(e7[3].real())

#vv = [Q2(3,2),Q2(7,5),Q2(4,3),Q2(2,1)]
#print(inner(vv,vv))

C = Mult(A2,Inv(A1))

print(C)

print(Det(C))
#
#B1 = Transpose([e0,e1,e2,e6])

#B2 = Transpose([e0,e1,e2,e5])

#print(A1)

#print(Det(B1))
#print(Det(B2))



B = product(*[g for i in range(4)])

D = product(*[g for i in range(4)])

print('Equiv Matrices:')
'''
for b1 in B:
    X1 = [[inner(u,v).real() for u in b1] for v in b1]
    if (np.linalg.det(np.matrix(X1)) != 0):
        print(X1)
'''

for b1 in B:
    X1 = [[inner(u,v).real() for u in b1] for v in b1]
    if (np.linalg.det(np.matrix(X1)) != 0):
        #print(X1)
        for b2 in D:
            #X1 = [[inner(u,v).real() for u in b1] for v in b1]
            X2 = [[inner(u,v).real() for u in b2] for v in b2]
            #print(X1)
            #print(np.linalg.det(np.matrix(X1)))
            if(X1 == X2):
                print('000000000000000000')
                print('b1 =', b1)
                print('b2 =', b2)
                Z = np.matrix([[Mult(b2,Inv(b1))[i][j].real() for i in range(4)] for j in range(4)])
                print(np.matrix(Mult(b2,Inv(b1))))
                print(Det(Mult(b2,Inv(b1))))
                w,v = np.linalg.eig(Z)
                print('w =', w)
                print('v =', v)
                print('Mult2 =', Mult(Mult(b2,Inv(b1)),Mult(b2,Inv(b1))))



