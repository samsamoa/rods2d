#!/usr/bin/python
from __future__ import division
import struct
import operator
import math
import os.path
from optparse import OptionParser
import numpy
from progressbar import ProgressBar

import readInfo

class struct():
	def __init__(self, *args, **kwargs):
		for k, v in kwargs.items():
			setattr(self, k, v)

parser = OptionParser()
parser.add_option("--dir", dest="dir", default="./")
parser.add_option("--startTime", dest="startTime", type="int", default=-1)
parser.add_option("--stopTime", dest="stopTime", type="int", default=-1)
parser.add_option("--speedup", dest="speedup", type="int", default=1)
parser.add_option("--singleFrame", action="store_true", dest="singleFrame", default=False)
parser.add_option("--debug", action="store_true", dest="debug", default=False)
parser.add_option("--allTime", action="store_true", dest="allTime", default=False)
parser.add_option("--endTimeFraction", dest="endTimeFraction", type="float", default=0.0)
parser.add_option("--track", action="store_true", dest="track", default=False)
parser.add_option("--minClumpSize", dest="minClumpSize", type="int", default=1)
parser.add_option("--dimensions", action="store_true", dest="storeDimensions", default=False)
parser.add_option("--noPBC", action="store_true", dest="noPBC", default=False)

(options, args) = parser.parse_args()

dir=options.dir
startTime=options.startTime
stopTime=options.stopTime
speedup=options.speedup
singleFrame=options.singleFrame
debug=options.debug
allTime=options.allTime
endTimeFraction=options.endTimeFraction
track=options.track
minClumpSize=options.minClumpSize
storeDimensions=options.storeDimensions
noPBC=options.noPBC

info=readInfo.read(dir)
startSl=info['startSl']
rodLength=info['rodLength']
frameSkip=info['frameSkip']
absStartTime=info['absStartTime']
absStopTime=info['absStopTime']
partNum=info['partNum']
timeStep=info['timeStep']



if allTime:
	#print "All time."
	startTime = absStartTime
	stopTime = absStopTime - 1
elif endTimeFraction != 0.0:
	#print "Fraction of time."
	stopTime = absStopTime - 1
	startTime = int(absStopTime - (absStopTime - absStartTime)*endTimeFraction)
elif (startTime == -1 or stopTime == -1):
	print str(absStartTime), " ", str(absStopTime), " ", str(frameSkip)
	startTime = int(raw_input("start: "))
	stopTime = int(raw_input("stop: "))
	
	

halfSideLength = startSl/2

cutoff = 4 # The cutoff for clumping
roiMaxX=startSl-rodLength-cutoff
roiMaxY=startSl-rodLength-cutoff
roiMinX=0+rodLength+cutoff
roiMinY=0+rodLength+cutoff


positionFile = open(os.path.join(dir,"coords-bin.dat"), "rb")
clumpFile = open(os.path.join(dir,"clumps2.dat"),"rb")

rodRecordType=numpy.dtype((numpy.double, (partNum,3)))
rodRecordSize = 24 * partNum
clumpRecordType = numpy.dtype(('i', (partNum)))
clumpRecordSize = 4 * partNum


if singleFrame:
	timeRange = [startTime]
else:
	timeRange = numpy.array(range(startTime, stopTime-1, frameSkip*speedup))


#clumpInfo=numpy.zeros(len(timeRange), dtype=object)
clumpInfo=range(len(timeRange))
type=numpy.dtype({'names': ['cm','rodIDs','tracked','timeIdx','rodAvgSpeed','hardDimensions','softDimensions','polarVector'],
				'formats': [object,object,object,object,object,object,object,object]})

allClumpData = numpy.fromfile(clumpFile, dtype=clumpRecordType, count=-1)

timeStepBetweenFrames = frameSkip*speedup*timeStep


progress = ProgressBar(maxval=len(timeRange))
for tI,t in enumerate(timeRange):
	progress.update(tI)
	#positionFile.seek(rodRecordSize*((t-absStartTime)//frameSkip))
	
	clumpData = allClumpData[(t-absStartTime)//frameSkip]
	cMax = numpy.max(clumpData)
	
	clumpList= sorted([numpy.flatnonzero(clumpData == selectedId) for selectedId in range(cMax+1)], key=len, reverse=True)
	if (minClumpSize != 1):
		clumpList = [ clump for clump in clumpList if clump.size >= minClumpSize ]
	
	clumpInfo[tI] = [[clumpList[i],-1] for i in range(len(clumpList))]
					
	
"""
	positionData = numpy.fromfile(positionFile, dtype=rodRecordType, count=1)[0]

	clumpCoords = [positionData[clump] for clump in clumpList]
	
	clumpHardDimensions = numpy.zeros((len(clumpList),2))
	clumpSoftDimensions = numpy.zeros((len(clumpList),2))
	clumpPolarVectors = numpy.zeros((len(clumpList),2))
	clumpMassCenters = numpy.zeros((len(clumpList),2))
	clumpPeriodic = numpy.zeros(len(clumpList), dtype=bool)
	
	for i in range(len(clumpList)): # Check for clumps that are "too big" to find a CM
		clumpPolarVectors[i][0] = numpy.mean(numpy.cos(clumpCoords[i][:,2]))
		clumpPolarVectors[i][1] = numpy.mean(numpy.sin(clumpCoords[i][:,2]))
		
		polarAxis = clumpPolarVectors[i]/numpy.linalg.norm(clumpPolarVectors[i])
		if not noPBC:
			if max(clumpCoords[i][:,0]) >= roiMaxX and min(clumpCoords[i][:,0]) <= roiMinX:
				clumpCoords[i][:,0] = ((clumpCoords[i][:,0]+halfSideLength) % startSl) - halfSideLength
				if max(clumpCoords[i][:,0]) >= (roiMaxX-startSl/2.0) and min(clumpCoords[i][:,0]) <= (roiMinX-startSl/2.0):
					clumpPeriodic[i] = True
					#if debug:
						#print "Periodic clump found"
			if max(clumpCoords[i][:,1]) >= roiMaxY and min(clumpCoords[i][:,1]) <= roiMinY:
				clumpCoords[i][:,1] = ((clumpCoords[i][:,1]+halfSideLength) % startSl) - halfSideLength
				if max(clumpCoords[i][:,1]) >= (roiMaxY-startSl/2.0) and min(clumpCoords[i][:,1]) <= (roiMinY-startSl/2.0):
					clumpPeriodic[i] = True
					#if debug:
					#	print "Periodic clump found"
		
		clumpMassCenters[i] = clumpCoords[i][:,0:2].mean(axis=0)
		
		if clumpPeriodic[i]:
			clumpHardDimensions[i] = clumpSoftDimensions[i] = (0,0)
		else:
			axisVector = clumpPolarVectors[i]/numpy.linalg.norm(clumpPolarVectors[i])
			offsets = clumpCoords[i][:,0:2]-clumpMassCenters[i]
			longOffsets = numpy.dot(offsets, axisVector)
			shortOffsets = numpy.cross(offsets, axisVector)
			
			clumpWidth = numpy.max(shortOffsets)-numpy.min(shortOffsets)
			clumpLength = numpy.max(longOffsets)-numpy.min(longOffsets)
			
			clumpWidth2 = numpy.std(shortOffsets)
			clumpLength2 = numpy.std(longOffsets)
			
			clumpHardDimensions[i] = (clumpWidth, clumpLength)
			clumpSoftDimensions[i] = (clumpWidth2*4, clumpLength2*4)

	if tI == 0:
		clumpRodSpeeds = numpy.zeros(len(clumpList))
	else:
		rodCoordDisp = numpy.abs(positionData[:,0:2] - lastPositions[:,0:2])
		rodCoordDispPBC = numpy.minimum(rodCoordDisp, startSl-rodCoordDisp)
		rodSpeeds = numpy.sqrt(numpy.sum(rodCoordDispPBC**2, axis=1))/timeStepBetweenFrames
		clumpRodSpeeds = [numpy.mean(rodSpeeds[clump]) for clump in clumpList]
	
	clumpInfo[tI] = [struct(
						cm = clumpMassCenters[i],
						rodIDs = clumpList[i],
						tracked = False,
						timeIdx = tI,
						rodAvgSpeed = clumpRodSpeeds[i],
						hardDimensions = clumpHardDimensions[i],
						softDimensions = clumpSoftDimensions[i],
						polarVector = clumpPolarVectors[i]
					) for i in range(len(clumpList))]
	
	lastPositions = positionData
"""	


progress.finish()

positionFile.close()
clumpFile.close()


if singleFrame:
	for clump in clumpInfo[0]:
		print len(clump.rodIDs), numpy.mod(clump.cm[0], startSl), numpy.mod(clump.cm[1], startSl), clump.polarVector[0], clump.polarVector[1], clump.hardDimensions[0], clump.hardDimensions[1]
	exit()

clumpSizeFile = open(dir+"/clumpSizes","w")
clumpSizes = [[clump[0].size for clump in clumps] for clumps in clumpInfo]	
for tI,tClumps in enumerate(clumpSizes):
	clumpSizeFile.write(str(timeRange[tI]) + " ")
	clumpSizeFile.write(' '.join([str(s) for s in tClumps]) + "\n")
clumpSizeFile.close()

if storeDimensions:
	a=numpy.array([(timeRange[tIdx],c[1].size,c.hardDimensions[0],c.hardDimensions[1]) for tIdx,p in enumerate(clumpInfo) for c in p if c[1].size>1])
	f=open (dir+"/clumpDimensions","w")
	for p in a:
		f.write(str(p[0])+" "+str (p[1])+" "+str (p[2])+" "+str(p[3])+"\n")

	f.close()

if debug:
	import matplotlib.pyplot as pyplot
	meanSizes = [numpy.mean(c) for c in clumpSizes]
	pyplot.figure(1)
	pyplot.plot(timeRange, meanSizes)
	
	pyplot.figure(2)
	a=numpy.array([(c[0].size, c.hardDimensions[1]) for p in  clumpInfo[(len(clumpInfo)//3)-1:] for c in p if c[0].size > 1])
	pyplot.scatter(a[:,0],a[:,1])
	
	pyplot.show()


idCounter = 0;



if track:
	progress = ProgressBar(maxval=len(timeRange))
	
	trackingData = list()
	
	idCounter = -1;
	
	for tI, t in enumerate(timeRange):
		progress.update(tI)
		for startClump in clumpInfo[tI]:
			if startClump[1] != -1:
				continue
			
			idCounter = idCounter + 1
			
			clumpTrack = list()
			prevClump = startClump
			for i in range(tI+1, len(timeRange)):
				prevClump[1] = idCounter
				clumpTrack.append(prevClump)
				newClumpCandidates=sorted(clumpInfo[i], key=lambda clump: numpy.intersect1d(prevClump[0],clump[0], assume_unique=True).size, reverse=True)
				newClump=newClumpCandidates[0]
				oldSize = prevClump[0].size
				commonSize = numpy.intersect1d(prevClump[0],newClump[0]).size
				newSize = newClump[0].size
				if ((commonSize <= 0.8*oldSize) or (newSize >= 1.2*oldSize)):
					break
				prevClump=newClump
			
			if len(clumpTrack)>1:
				trackingData.append(clumpTrack)
	
	"""
	allSpeedsFile = open(dir+"/allSpeeds","w")
	
	for track in trackingData:
		cmCoordDisp = numpy.abs(numpy.array([track[i].cm - track[i-1].cm for i in range(1,len(track))]))
		cmCoordDispPBC = numpy.minimum(cmCoordDisp, startSl-cmCoordDisp)
		cmDisp=numpy.array([numpy.linalg.norm(cmCoordDispPBC[i]) for i in range(len(track)-1)])
		
		sizes=[(track[i].rodIDs.size + track[i-1].rodIDs.size)/(2.0) for i in range(1,len(track))]
		speeds=[disp/timeStepBetweenFrames for disp in cmDisp]
		
		rodSpeeds = [track[i].rodAvgSpeed for i in range(1,len(track))]
		for p in zip(sizes, speeds, rodSpeeds):
			allSpeedsFile.write(str(p[0]) + "," + str(p[1])+"," + str(p[2]) + "\n")

	allSpeedsFile.close()
	"""
	
	progress.finish()
	
offset=(startTime-absStartTime)//frameSkip
newClumpData = allClumpData.copy()
newClumpData.fill(-1)
for tI, t in enumerate(timeRange):
	for clump in clumpInfo[tI]:
		newClumpData[offset+tI][clump[0]] = clump[1]