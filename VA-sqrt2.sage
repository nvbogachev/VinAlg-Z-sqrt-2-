from itertools import product
import numpy as np
from math import *
from sage.geometry.polyhedron.constructor import Polyhedron
import coxiter
import time
import base
from base import Q2
from base import div
from base import divide
from base import s, Transpose, Minor, Det, Inv, Mult


t0 = time.time()

d = Q2(0,1)
d1 = Q2(3,1)

lr = [(1,0), (2,0), (2,1), (3,1), (4,2), (6,2), (8,5), (16,10)]
lr = [Q2(*l) for l in lr]

m=0
M=0.01

#m = 25.7
#M = 25.75

def minor(*arg):
    sum = 0
    if (len(arg) == 1):
        sum = - arg[0].square()*d
    if (len(arg) > 1):
        sum = - arg[0].square()*d + arg[1].square()*d1
        for v in arg[2:]:
            sum = sum + v.square()
    return sum

'''
def minor(*arg):
    sum = 0
    if (len(arg) == 1):
        sum = - arg[0].square()*d
    if (len(arg) > 1):
        sum = - arg[0].square()*d + arg[1].square()
        for v in arg[2:]:
            sum = sum + v.square()
    return sum
'''



def inner(arg1,arg2):
    sum = - (arg1[0]*arg2[0])*d + (arg1[1]*arg2[1])*d1
    n = len(arg1)
    for j in range(2,n):
        sum = sum + arg1[j]*arg2[j]
    return sum



def Box110(vectors, k): ## for y1/sqrt(2) and y0 = (0,0) and y1 ><= 0
    M = minor(*vectors)
    s_plus=floor(s(M.real(),k.real())/d1.real())
    s_minus=floor(s(M.antireal(),k.antireal())/d1.antireal())
    #print s_plus, s_minus
    return (Q2(a,b) for a in range(-int(floor(s_minus+s_plus))-2, int(floor(s_plus+s_minus)) + 3) for b in range(-int(floor(a/sqrt(2) + s_plus + sqrt(2))), int(floor(s_plus - a/sqrt(2) + 3))) )


def Box11(vectors, k): ## for y1/sqrt(2) and y1 > 0
    M = minor(*vectors)
    s_plus=floor(s(M.real(),k.real())/d1.real())
    s_minus=floor(s(M.antireal(),k.antireal())/d1.antireal())
    return (Q2(a,b) for a in range(-1-int(floor(s_minus)), int(floor(s_plus+s_minus)) + 3) for b in range(-int(floor(a/sqrt(2) + sqrt(2))), int(floor(s_plus - a/sqrt(2) + 3))) )






def Box10(vectors, k): ## for yj/sqrt(2)  and y0 = (0,0)
    M = minor(*vectors)
    s_plus=s(M.real(),k.real())
    s_minus=s(M.antireal(),k.antireal())
    return (Q2(a,b) for a in range(-int(floor(s_minus+s_plus))-2, int(floor(s_plus+s_minus)) + 3) for b in range(-int(floor(a/sqrt(2) + s_plus + sqrt(2))), int(floor(s_plus - a/sqrt(2) + 3))) )

def Box0(vectors, k): ## for yj and y0 = (0,0)
    M = minor(*vectors)
    s_plus=s(M.real(),k.real())
    s_minus=s(M.antireal(),k.antireal())
    return (Q2(a,b) for a in range(-int(floor((s_minus+s_plus)/sqrt(2)))-1, int(floor(s_plus/sqrt(2)+s_minus/sqrt(2))) + 3) for b in range(-int(floor((a+s_plus)/sqrt(2) + 1)), int(floor(s_plus/sqrt(2) - a/sqrt(2) + 3))) )


def Box1(vectors, k): ## for yj/sqrt(2) > 0
    M = minor(*vectors)
    s_plus=s(M.real(),k.real())
    s_minus=s(M.antireal(),k.antireal())
    return (Q2(a,b) for a in range(-1-int(floor(s_minus)), int(floor(s_plus+s_minus)) + 3) for b in range(-int(floor(a/sqrt(2) + sqrt(2))), int(floor(s_plus - a/sqrt(2) + 3))) )

def Box(vectors, k): ## for yj > 0
    M = minor(*vectors)
    s_plus=s(M.real(),k.real())
    s_minus=s(M.antireal(),k.antireal())
    return (Q2(a,b) for a in range(-1-int(floor((s_minus)/sqrt(2))), int(floor(s_plus/sqrt(2)+s_minus/sqrt(2))) + 3) for b in range(-int(floor((a)/sqrt(2) + 1)), int(floor(s_plus/sqrt(2) - a/sqrt(2) + 3))) )

#print(minor(Q2(0,0), Q2(0,0), Q2(0,0), Q2(0,-1)))

print inner([Q2(5,4), Q2(3,2), Q2(0,0), Q2(3,1)], [Q2(5,4), Q2(3,2), Q2(0,0), Q2(3,1)])


print (div(Q2(5,4)*Q2(5,4),Q2(6,2)).real())

def Dist(kor):
    return (kor[0].square().real())/minor(kor)


## Let's go!!!

'''
predfuncone = []



for k in lr:
    for y_1 in Box110([Q2(0,0)],k):
        if (divide((d1*Q2(0,1)*y_1), k)):
            #print 'y_1 =', y_1
            for y_2 in Box10([Q2(0,0), div(y_1,Q2(0,1))], k):
                if (y_1.a%2 == 0 and y_2.a%2 == 0 and (y_1.b + y_2.b)%2 == 0 and divide(Q2(0,1)*y_2, k)):
                    #print 'y_2 =', y_2
                    for y_3 in Box0([Q2(0,0), div(y_1,Q2(0,1)), div(y_2,Q2(0,1))],k):
                        if (divide(Q2(2,0)*y_3, k) and minor(Q2(0,0), div(y_1,Q2(0,1)), div(y_2,Q2(0,1)), y_3).real() == k.real() ):
                            #print [Q2(0,0), div(y_1,Q2(0,1)), div(y_2,Q2(0,1)), y_3]
                            predfuncone.append([Q2(0,0), div(y_1,Q2(0,1)), div(y_2,Q2(0,1)), y_3])


print 'predfuncone =', predfuncone
print len(predfuncone)





pr = [[j.real() for j in r] for r in predfuncone]
print pr



print time.time() - t0




K = QuadraticField(2)
OK = K.ring_of_integers()

hh = Polyhedron(ieqs = [tuple([0]+pr[0])], base_ring = RDF)
#print hh.Hrepresentation()


cone = Polyhedron(ieqs = [tuple([0]+[0]*4)], base_ring = RDF)
#print 'cone0 =', cone


conus = []
for v in pr:
    hh = Polyhedron(ieqs = [tuple([0]+ [-t for t in v])], base_ring = RDF)
    if cone.intersection(hh).dim() == 4:
        cone = cone.intersection(hh)
        conus.append(v)
        #    else:
#cone = cone.intersection(Cone([-root]).dual())
print cone.Hrepresentation()

print float(sqrt(2)/2)

print float(1 + sqrt(2)/2)
'''

#print (cone.intersection(halfplane).dim() )
#        cone = cone.intersection(halfplane)


e0 = [Q2(0,0), Q2(0,0), Q2(0,0), Q2(-1,0)]
e1 = [Q2(0,0), Q2(0,0), Q2(0,-1), Q2(0,0)]
e2 = [Q2(0,0), Q2(0,-1), Q2(0,0), Q2(0,0)]



e3 = [(1, 1), (0, 0), (0, 0), (2, 1)]
e3 = [Q2(*l) for l in e3]
e4 = [(1, 1), (0, 0), (2, 1), (0, 0)]
e4 = [Q2(*l) for l in e4]
e5 = [(2, 1), (0, 0), (2, 1), (2, 1)]
e5 = [Q2(*l) for l in e5]
e6 = [(1, 1), (1, 0), (1, 1), (0, 0)]
e6 = [Q2(*l) for l in e6]
e7 = [(2, 3/2), (1, 1/2), (1, 1/2), (2, 1)]
e7 = [Q2(*l) for l in e7]
e8 = [(5, 4), (2, 1), (0, 0), (5, 4)]
e8 = [Q2(*l) for l in e8]
#e6 = [(4.0, 3.0), (5.0, 4.0), (3, 1), (0, 0)]
#e6 = [Q2(*l) for l in e6]
#e7 = [(4.0, 3.0), (5.0, 3.0), (3, 2), (1, 1)]
#e7 = [Q2(*l) for l in e7]
#e8 = [(4.0, 3.0), (5.0, 4.0), (1, 1), (1, 1)]
#e8 = [Q2(*l) for l in e8]


roots = [e0, e1, e2, e3, e4, e5, e6, e7, e8]
    #N = [[inner(u,v).real() for u in roots] for v in roots]

#for r in N:
# print r

#print(coxiter.run(N, 4))

#print float( ( div(inner(e3,e7)*inner(e3,e7), (inner(e3,e3)*inner(e7,e7))) ).real() )


for k in lr:
    for a0 in range(int(floor((sqrt(m*k.real()) - sqrt(-k.antireal()/d.antireal()))/sqrt(2))), 2+int(floor((sqrt(M*k.real()) + sqrt(-k.antireal()/d.antireal()))/sqrt(2)))):
        for b0 in range(int(floor(sqrt(m*k.real()) - a0/sqrt(2))), int(floor(sqrt(M*k.real()) - a0/sqrt(2))) + 2):
            y_0 = Q2(a0,b0)
            if ((y_0).real() > 0 and y_0.square().real()/(2*k.real()) <= M and m <= y_0.square().real()/(2*k.real())):
                for y_1 in Box11([div(y_0,Q2(0,1))],k):
                    if (divide((d1*Q2(0,1)*y_1), k) and y_1.real() >= 0 and (y_0.a - y_1.a)%2 == 0):
                        for y_2 in Box1([div(y_0,Q2(0,1)), div(y_1,Q2(0,1))], k):
                            if (divide(Q2(0,1)*y_2, k) and y_2.real() >= 0 and (y_2.a + y_1.a)%2 == 0 and (y_2.b + y_1.b)%2 == 0):
                                for y_3 in Box([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), div(y_2,Q2(0,1))],k):
                                    if (y_3.real() >= 0 and divide(Q2(2,0)*y_3, k) and minor(div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), div(y_2,Q2(0,1)), y_3).real() == k.real() and inner([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), div(y_2,Q2(0,1)), y_3], e2).real() <= 0 and inner([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), div(y_2,Q2(0,1)), y_3], e3).real() <= 0 and inner([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), div(y_2,Q2(0,1)), y_3], e4).real() <= 0 and inner([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), div(y_2,Q2(0,1)), y_3], e5).real() <= 0 and inner([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), div(y_2,Q2(0,1)), y_3], e6).real() <= 0 and inner([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), div(y_2,Q2(0,1)), y_3], e7).real() <= 0 and inner([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), div(y_2,Q2(0,1)), y_3], e8).real() <= 0):
                                        print([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), div(y_2,Q2(0,1)), y_3], k, y_0, y_0.square().real()/(2*k.real()) )
#roots.append([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), y_2, y_3])
#prin


print time.time() - t0

t1 = time.time()

N = [[inner(u,v).real() for u in roots] for v in roots]

print 'N58 =', float((N[5][8]*N[5][8])/(N[5][5]*N[8][8]))


for r in N:
    print r

print(coxiter.run(N, 4))

print time.time() - t1

#print 'cos2 =', float(cos(pi/8)*cos(pi/8)), cos(pi/8)

#print float(2*(2*cos(pi/8)*cos(pi/8) - 1)*(2*cos(pi/8)*cos(pi/8) - 1))








'''

def FundCone(s):
    V1_roots = [v for k in s.root_lengths for v in s.Roots_decomposed_into(s.V([0]*s.n), k) if s.IsRoot(v)]
        print('roots in V1: {}'.format(V1_roots))
        cone = Cone([[0]*s.n]).dual()
        #print('cone', cone.rays())
        for root in V1_roots:
            halfplane = Cone([root]).dual()
            #print('halfplane', halfplane.rays())
            if cone.intersection(halfplane).dim() == s.n:
                cone = cone.intersection(halfplane)
            else:
                cone = cone.intersection(Cone([-root]).dual())
    #print('cone', cone.rays())
    print('FundCone returned',[s.V(r) for r in cone.dual().rays()])
        return [s.V(r) for r in cone.dual().rays()]

'''



'''
print('Box =', [i for i in Box([Q2(0,0), Q2(0,0)], Q2(2,0))])

#K = lr[-1]
#print(range(int(floor((sqrt(m*K.real()) - sqrt(-K.antireal()/d.antireal()))/sqrt(2))), 2+int(floor((sqrt(M*K.real()) + sqrt(-K.antireal()/d.antireal()))/sqrt(2)))))
for k in lr:
    #print k
    #for a0 in range(int(floor((sqrt(m*k.real()) - sqrt(-k.antireal()/d.antireal()))/sqrt(2))), 2+int(floor((sqrt(M*k.real()) + sqrt(-k.antireal()/d.antireal()))/sqrt(2)))):
    #for b0 in range(int(floor(sqrt(m*k.real()) - a0/sqrt(2))), int(floor(sqrt(M*k.real()) - a0/sqrt(2))) + 2):
                y_0 = Q2(0,0)
                #if ((div(y_0,Q2(0,1))).real() > 0):
                if (y_0.square().real()/(2*k.real()) <= M and m <= y_0.square().real()/(2*k.real())):
                    #print('k, y0:', k, y_0)
                    for y_1 in Box1([div(y_0,Q2(0,1))],k):
                        if ((y_0.a + y_1.a)%2 == 0 and (y_0.b + y_0.a + y_1.b)%2 == 0 and divide(Q2(0,1)*y_1, k)):
                            #print('00000')
                            for y_2 in Box([div(y_0,Q2(0,1)), div(y_1,Q2(0,1))], k):
                                if (divide(Q2(2,0)*y_2, k)):
                                    #print('000000000000000000')
                                    for y_3 in Box([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), y_2],k):
                                        if (divide(Q2(2,0)*y_3, k) and minor(div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), y_2, y_3).real() == k.real() ):
                                            #print('0000000')
                                            #if (inner([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), y_2, y_3], e3).real() <= 0 and inner([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), y_2, y_3], e4).real() <= 0): #and inner([div(y_0,Q2(0,1)), div(y_1,Q2(2,0)), y_2, y_3], e5).real() <= 0 and inner([div(y_0,Q2(0,1)), div(y_1,Q2(2,0)), y_2, y_3], e6).real() <= 0  and inner([div(y_0,Q2(0,1)), div(y_1,Q2(2,0)), y_2, y_3], e7).real() <= 0 and inner([div(y_0,Q2(0,1)), div(y_1,Q2(2,0)), y_2, y_3], e8).real() <= 0):
                                                    print([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), y_2, y_3], k, y_0, div(y_0, Q2(0,1)), y_0.square().real()/(2*k.real()) )
                                                    roots.append([div(y_0,Q2(0,1)), div(y_1,Q2(0,1)), y_2, y_3])
#print(inner([div(y_0,Q2(0,1)), div(y_1,Q2(2,0)), y_2, y_3],[div(y_0,Q2(0,1)), div(y_1,Q2(2,0)), y_2, y_3]))


Root = min(roots, key = Dist)



# ROOTS not v0
# (div(y_1,Q2(0,1))).real() >= y_2.real() and y_2.real() >= 0 and 
# y_2.real() >= y_3.real() and y_3.real() >= 0 and


g = [e0,e1,e2,e3,e4]

#print(cos(pi/8)*cos(pi/8))

#print([[int(floor(inner(u,v).real())) for u in g] for v in g])

#M = np.matrix([[inner(u,v).real() for u in g] for v in g])

#print(M)

print(div(Q2(1,1), Q2(0,1)))

#print(Q2(1,0)*Q2(0,1))

N1 = [[inner(u,v) for u in g] for v in g]

N = [[inner(u,v).real() for u in g] for v in g]

for r in N1:
    print(r)

print(coxiter.run(N, 4))

print time.time() - t0
t1 = time.time()

B = product(*[g for i in range(4)])

D = product(*[g for i in range(4)]) 

print('Equiv Matrices:')
'''
    
'''
    for b1 in B:
    X1 = [[inner(u,v).real() for u in b1] for v in b1]
    if (np.linalg.det(np.matrix(X1)) != 0):
    print(X1)
    

for b1 in B:
    X1 = [[inner(u,v).real() for u in b1] for v in b1]
    if (np.linalg.det(np.matrix(X1)) != 0):
        #print(X1)
        for b2 in D:
            #X1 = [[inner(u,v).real() for u in b1] for v in b1]
            X2 = [[inner(u,v).real() for u in b2] for v in b2]
            #print(X1)
            #print(np.linalg.det(np.matrix(X1)))
            if(X1 == X2 and b1 != b2):
                print('000000000000000000')
                print('b1 =', b1)
                print('b2 =', b2)
                #Z = np.matrix([[Mult(b2,Inv(b1))[i][j].real() for i in range(4)] for j in range(4)])
                #print(np.matrix(Mult(b2,Inv(b1))))
                print(Det(Mult(b2,Inv(b1))))
                #w,v = np.linalg.eig(Z)
                #print('w =', w)
                #print('v =', v)
                print('Mult2 =', Mult(Mult(b2,Inv(b1)),Mult(b2,Inv(b1))))

'''


'''

print time.time() - t1
print time.time() - t0

'''

