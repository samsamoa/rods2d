#!/usr/bin/python

"""
Measures the purity of clusters
"""
from __future__ import division
import struct
from optparse import OptionParser
import numpy
import matplotlib.pyplot

import readInfo

parser = OptionParser()
parser.add_option("--dir", dest="dir", default="./")

(options, args) = parser.parse_args()

dir=options.dir

info=readInfo.read(dir)

#if (info['absStopTime'] - info['absStartTime']) != 10000000:
#		print "-1"
#		exit()

equilibrationTime = (info['absStartTime'] + info['absStopTime'])/2
minimumClusterSize=3

"""
weightedClusterCount=0
weightedImpureCount=0

f=open(dir+"/aliveness.txt","r")

for line in f:
    q=line.split()
    time=int(q[0])
    if (time < equilibrationTime):
    	continue
    cSize=int(q[2])
    if cSize<minimumClusterSize:
		continue
    aliveC=int(q[3])
    weightedClusterCount = weightedClusterCount + cSize
    if (0.2 < (aliveC/cSize) < 0.8):
        weightedImpureCount = weightedImpureCount + aliveC


print 1-(weightedImpureCount/weightedClusterCount)
f.close()
"""

clumpInfoRecord=numpy.dtype([('time','i4'), ('sizes','i4',info['partNum']), ('aliveSizes','i4',info['partNum']), ('deadSizes','i4',info['partNum']), ('purities','i4',100), ('averageClusterSize','f8'), ('weightedAverageClusterSize','f8'), ('cutoffAverageClusterSize','f8')]);
clumpInfoData=numpy.fromfile(dir + "/otherClusterInfo.dat", dtype=clumpInfoRecord, count=-1)
