# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 10:28:29 2016

@author: coulaud
"""

import random
import numpy
import math


meshFile="Vega_Z09_RR1_sans_fils.geod"
fmmFile="Vega_Z09_RR1_sans_fils.fma"

ptfile = open(fmmFile,'w')
ptfile.write("8  4 \n")

Fichier = open(meshFile,'r')
line = Fichier.readline()
line = Fichier.readline()
# mear size triangles an points
line = Fichier.readline()
size = line.rstrip('\n\r').split()
print(line)
Npt = int(size[1])
NT = int(size[0])
x = numpy.zeros([Npt,3])
for i in range(Npt):
    line = Fichier.readline()
    size = line.rstrip('\n\r').split()
    x[i,0] = float(size[1])
    x[i,1] = float(size[2])
    x[i,2] = float(size[3])
    
a = numpy.amin(x,axis=0)
b = numpy.amax(x,axis=0)
print(a)
print(b)
length = math.ceil(max(b-a))
centre= (a+b)/2
print('Centre: ',centre)
print('length: ',length,max(b-a))
ptfile.write(str(Npt)+'  ' + str((length)/2) + '  '+ str(centre[0])
    + '  '+ str(centre[1])+ '  '+ str(centre[2])  +"\n"  )   

for i in range(Npt):
    rho = 2*random.random()-1
    str1 = str(x[i,0])+ '  ' + str(x[i,1])+'  ' +str(x[i,2])+ '  '  + str(rho) +"\n"
    ptfile.write(str1)
ptfile.close()