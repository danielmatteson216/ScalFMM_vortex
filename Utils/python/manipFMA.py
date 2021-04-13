# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 17:25:37 2014

@author: coulaud
"""
#Read file comming from Stamp 

import numpy
#import sys

#Filename='/Users/coulaud/Dev/src/ScalFMM/scalfmmT/Data/cea-200-per.xyzqf'

#
# Read Header
#
#global nbRecordPerline,dataType,NbParticles,BoxWidth,Centre
dataType        = 8
nbRecordPerline = 4
NbParticles     = 0
BoxWidth        = 1.0
Centre= (0.0,0.0,0.0)
#
def readFMAHeader(Filename):
    """
    This function read the header of the FMA file Filename
    """
    global nbRecordPerline,dataType,NbParticles,BoxWidth,Centre
    f = open(Filename, 'r')
    header1 = f.readline()
    header2 = f.readline()
    s1,s2=(item.strip() for item in header1.split())
    dataType = int(s1)
    nbRecordPerline = int(s2)
    s1,s2,s3,s4,s5 = (item.strip() for item in header2.split())
    NbParticles=int(s1)
    BoxWidth=float(s2)*2
    print "BoxWidth= ", BoxWidth 
    Centre = [float(s3),float(s4),float(s5)]
    f.close()
    return dataType,nbRecordPerline,NbParticles,BoxWidth,Centre
#
def readFMAfile(fileName):
    """
        Read data from an ascci file il the FMA format
        Return data (x,y,z,q) or (x,y,z,q,p,fx,fy,fz)
    """
    global nbRecordPerline,dataType,NbParticles,BoxWidth,Centre
    readFMAHeader(fileName)
    print NbParticles
    # Read data in a non compact mode
    #
    data = numpy.loadtxt(fileName,skiprows=2)

#    if nbRecordPerline == 4 or nbRecordPerline == 8 :
    data = numpy.loadtxt(fileName,skiprows=2)
#    else:
        #print " wrong data nbRecordPerline (",nbRecordPerline,"). Should be 4 or 8"
 #       os._exit()
    return data,BoxWidth,Centre

#
def writeFMAfilefmt(fileName,data,boxWidth, centre,myformat):
    """
        write data from an ascci file il the FMA format
        Return data (x,y,z,q) or (x,y,z,q,p,fx,fy,fz)
    """
  #  f = open(fileName, "w")
    s = ' 8  '+ str(numpy.shape(data)[1]) + '\n'
 #   f.write(s)
    s= s + str(numpy.shape(data)[0])+ '  '+ str(boxWidth/2)+ '  '+ str(Centre[0])+ '  '+ str(Centre[1])+ '  '+ str(Centre[2])
  #  f.write(s)
  #  f.close()
    numpy.savetxt(fileName,data,delimiter='  ',header=s,comments=' ',fmt=myformat)

def writeFMAfile(fileName,data,boxWidth, centre):
    """
        write data from an ascci file il the FMA format
        Return data (x,y,z,q) or (x,y,z,q,p,fx,fy,fz)
    """
  #  f = open(fileName, "w")
    s = ' 8  '+ str(numpy.shape(data)[1]) + '\n'
 #   f.write(s)
    s= s + str(numpy.shape(data)[0])+ '  '+ str(boxWidth/2)+ '  '+ str(Centre[0])+ '  '+ str(Centre[1])+ '  '+ str(Centre[2])
  #  f.write(s)
  #  f.close()
    numpy.savetxt(fileName,data,delimiter='  ',header=s,comments=' ')

