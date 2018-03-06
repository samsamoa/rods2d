#!/usr/bin/python

"""
Identifies swarms and makes images for g(r)
python -i ~/Work-Hagan/Rods2d/trunk/scripts/clumpCalc.py --dir ~/Work-Hagan/sim/sca-n100-f0.0-l20-s3/d0.0100-n1600/f1.0-a1.0/ --endTimeFraction 0.1 --speedup 10 --save --analysis --width 100 --height 100
"""
from __future__ import division
import struct
from optparse import OptionParser
import numpy
import subprocess
import readline
import os.path
import string

#import scipy.interpolate
#import scipy.optimize

import readInfo

parser = OptionParser()
parser.add_option("--dir", dest="baseDir", default="./")
parser.add_option("--startTime", dest="startTime", type="int", default=-1)
parser.add_option("--stopTime", dest="stopTime", type="int", default=-1)
parser.add_option("--speedup", dest="speedup", type="int", default=1)
parser.add_option("--analysisFrameSkip", dest="analysisFrameSkip", type="int", default=-1)
parser.add_option("--frameGroup", dest="frameGroup", type="int", default=0)
parser.add_option("--debug", action="store_true", dest="debug", default=False)
parser.add_option("--allTime", action="store_true", dest="allTime", default=False)
parser.add_option("--endTimeFraction", dest="endTimeFraction", type="float", default=0.0)
parser.add_option("--endFrameCount", dest="endFrameCount", type="int", default=0)
parser.add_option("--equilibrationTime", dest="equilibrationTime", type="int", default=-1)
parser.add_option("--calculateClumps", action="store_true", dest="calculateClumps", default=False)
parser.add_option("--storeClumps", action="store_true", dest="storeClumps", default=False)
parser.add_option("--analysis", action="store_true", dest="analysis", default=False)
parser.add_option("--save", action="store_true", dest="save", default=False)
parser.add_option("--width", dest="width", type="int", default=0)
parser.add_option("--height", dest="height", type="int", default=0)
parser.add_option("--resolution", dest="resolution", type="int", default=1)
parser.add_option("--noPBC", action="store_true", dest="noPBC", default=False)
parser.add_option("--rerun", action="store_true", dest="rerun", default=False)
parser.add_option("--onlyUniform", action="store_true", dest="onlyUniform", default=False)
parser.add_option("--fakeAlive", dest="fakeAlive", default="")
parser.add_option("--aliveFraction", dest="aliveFraction", type="float", default=0.5)


(options, args) = parser.parse_args()

baseDir=options.baseDir
startTime=options.startTime
stopTime=options.stopTime
speedup=options.speedup
analysisFrameSkip=options.analysisFrameSkip
debug=options.debug
allTime=options.allTime
storeClumps=options.storeClumps
calculateClumps=options.calculateClumps
analysis=options.analysis
endTimeFraction=options.endTimeFraction
endFrameCount=options.endFrameCount
equilibrationTime=options.equilibrationTime
save=options.save
imWidth=options.width
imHeight=options.height
noPBC=options.noPBC
resolution=options.resolution
rerun=options.rerun
fakeAlivePrefix=options.fakeAlive
aliveFraction=options.aliveFraction

baseDir=os.path.expanduser(baseDir)

import matplotlib
if save:
	matplotlib.use("PDF")
import matplotlib.pyplot as plt
matplotlib.rc('font', family='serif')

info=readInfo.read(baseDir,True)
startSl=info['startSl']
rodLength=info['rodLength']
frameSkip=info['frameSkip']
absStartTime=info['absStartTime']
absStopTime=info['absStopTime']
partNum=info['partNum']
timeStep=info['timeStep']

if (imWidth == 0):
	imWidth = numpy.floor(startSl/numpy.sqrt(8.0)-rodLength/2.0)

if (imHeight == 0):
	imHeight = numpy.floor(startSl/numpy.sqrt(8.0)-rodLength/2.0)





if allTime:
	#print "All time."
	startTime = absStartTime
	stopTime = absStopTime - 1

if endTimeFraction != 0.0:
	#print "Fraction of time."
	stopTime = absStopTime - 1
	startTime = int(absStopTime - (absStopTime - absStartTime)*endTimeFraction)

if endFrameCount != 0:
	stopTime = absStopTime - 1
	startTime = absStopTime - endFrameCount	

if equilibrationTime != -1:
	stopTime = absStopTime - 1
	startTime = absStartTime + equilibrationTime

if (startTime == -1 or stopTime == -1):
	print str(absStartTime), " ", str(absStopTime), " ", str(frameSkip)
	startTime = int(raw_input("start: "))
	stopTime = int(raw_input("stop: "))


minX = -imWidth
maxX = imWidth
minY = -imHeight
maxY = imHeight

pbcOption=""
if noPBC:
	pbcOption=" --noPBC "


allClumpInfo=[]

corrDataList=list()
corrDataOrientedList=list()
corrDataOrientedDirectionalList=list()
corrDataUniformList=list()


for subInfo in info['allInfo']:
	subAbsStartTime = subStartTime = subInfo['startTime']
	subStopTime = subInfo['stopTime']
	
	if (not allTime):
		if (startTime >= subStopTime):
			continue
		elif (subStartTime >= stopTime):
			continue
	
	frameSkip=info['frameSkip']
	dir = subInfo['dir']
	
	if (analysisFrameSkip != -1):
		if (analysisFrameSkip % frameSkip != 0):
			print "Frame skip and simulation frame skip not compatible."
			exit()
		speedup = int(analysisFrameSkip/frameSkip)
	
	purityCommand = "--purityFile " + dir + "/purity.txt"
	
	if options.frameGroup == 0:
		startTimeList = numpy.array([subStartTime,])
	else:
		startTimeList = numpy.arange(subAbsStartTime, subStopTime+1, frameSkip*options.frameGroup)
	
	uniformDists = []
	
	for runStartTime in startTimeList:
		if (options.frameGroup != 0):
			subStartTime = runStartTime
			subStopTime = runStartTime + (options.frameGroup-1)*frameSkip
		
		if (not os.path.exists(dir+"/coords-bin.dat")):
			continue
		
		if (fakeAlivePrefix != ""):
			aliveFile = fakeAlivePrefix+"/alive.dat"
		else:
			aliveFile = dir+"/alive.dat"
		
		command = ("~/Work-Hagan/Rods2d/trunk/build/clump -n " + str(partNum) + " -l " + str(rodLength) +
			 " -L " + str(startSl) + " -t " + str(subAbsStartTime) + " -T " + str(subStopTime) +
			 " -s " + str(frameSkip) + " -c " + str(subStartTime) + " -f " + dir + "/coords-bin.dat" +
			 " --speedup " + str(speedup) + pbcOption + " -a " + aliveFile + " " + purityCommand)
		
		#l=9, (5, 7, 15, 20)
		#l=20, (5, 7, 12, 12)
		#l=3, (10, 14, 24, 24)
		if (storeClumps or calculateClumps):
			command = command + (" -D 2.0 -F 0 --clumpInfoFile " + dir + "/otherClusterInfo.dat" + " --aliveInfoFile " + dir + "/aliveness.txt --doClumping")
			if (storeClumps):
				command = command + (" --clumpFile " + dir + "/clumps2.dat")
		elif (analysis):
			command = command + (" --analysis --GofRFile " + dir + "/GofR-cm.dat -R " + str(resolution) + " -w " +  str(imWidth) + " -h " + str(imHeight))
		else:
			print "Nothing to do."
			exit()
		
		print command
		
		if ((not rerun) and analysis and os.path.exists(dir + "/GofR-cm.dat")):
			print "using previous data"
		elif ((not rerun) and storeClumps and os.path.exists(dir + "/clumps2.dat")):
			print "using previous data"
		elif ((not rerun) and calculateClumps and os.path.exists(dir + "/otherClusterInfo.dat")):
			print "using previous data"
		else:
			output=subprocess.Popen(command,shell=True).communicate()
		
		imageExtent=[minX,maxX,minY,maxY]
				
		if analysis:			
			corrFile=open(dir + "/GofR-cm.dat")
			corrRecord=numpy.dtype("("+str(int(2*imWidth*resolution))+","+str(int(2*imHeight*resolution))+")i4");
			corrRecordOriented=numpy.dtype("("+str(int(2*imWidth*resolution))+","+str(int(2*imHeight*resolution))+")f8");
			corrRecordUniform=numpy.dtype(str(int(imWidth*resolution))+"i4");
			
			corrDataList.append(numpy.fromfile(corrFile, dtype=corrRecord,count=4))
			corrDataOrientedList.append(numpy.fromfile(corrFile, dtype=corrRecordOriented,count=4))
			corrDataOrientedDirectionalList.append(numpy.fromfile(corrFile, dtype=corrRecordOriented,count=4))
			corrDataUniformList.append(numpy.fromfile(corrFile, dtype=corrRecordUniform,count=4))
				
		if (storeClumps or calculateClumps):
			clumpInfoRecord=numpy.dtype([('time','i4'), ('sizes','i4',info['partNum']), ('aliveSizes','i4',info['partNum']), ('deadSizes','i4',info['partNum']), ('purities','i4',100), ('averageClusterSize','f8'), ('weightedAverageClusterSize','f8'), ('cutoffAverageClusterSize','f8')]);
			clumpInfoData=numpy.fromfile(dir + "/otherClusterInfo.dat", dtype=clumpInfoRecord, count=-1)
			allClumpInfo.append(clumpInfoData)



if analysis:
	density=partNum/(startSl**2)
	frameCount = int(numpy.floor((stopTime-startTime)/(speedup*frameSkip) + 1))
	
	uniformBinMinima = numpy.arange(0,imWidth,1.0/resolution)
	uniformBinMaxima = uniformBinMinima+(1.0/resolution)
	uniformNorm = numpy.pi * (uniformBinMaxima**2 - uniformBinMinima**2)
	
	twoDNorm = resolution*resolution/(frameCount*(density*aliveFraction)*(partNum*aliveFraction))#*(rodLength+1))
	oneDNorm = 1.0/(frameCount*(partNum*aliveFraction)*(density*aliveFraction)*uniformNorm)#*(rodLength+1))
	
	corrData = numpy.sum(numpy.array(corrDataList), axis=0) * twoDNorm
	corrDataOriented = numpy.sum(numpy.array(corrDataOrientedList), axis=0) * twoDNorm
	corrDataOrientedDirectional = numpy.sum(numpy.array(corrDataOrientedDirectionalList), axis=0) * twoDNorm
	corrDataUniform = numpy.sum(numpy.array(corrDataUniformList), axis=0) *oneDNorm
	
	corrData = (corrData[:,::-1]+corrData)/2.0
	corrDataOriented = (corrDataOriented[:,::-1]+corrDataOriented)/2.0
	corrDataOrientedDirectional = (corrDataOrientedDirectional[:,::-1]+corrDataOrientedDirectional)/2.0
	
	"""	
	horizontalSection=numpy.mean(corrData[3].transpose()[maxX-3-1:maxX+3-1],axis=0)
	horizontalSection=((horizontalSection+horizontalSection[::-1])/2.0)[(horizontalSection.shape[0]/2):]
	positions=numpy.arange(1.0/resolution,maxX+1.0/resolution,1.0/resolution)

	matplotlib.pyplot.figure(11)
	matplotlib.pyplot.loglog(positions, horizontalSection)
	matplotlib.pyplot.savefig(dir+"/horizontalSectionLogLog.pdf")

	matplotlib.pyplot.figure(12)
	matplotlib.pyplot.plot(positions, horizontalSection)
	matplotlib.pyplot.savefig(dir+"/horizontalSection.pdf")

	matplotlib.pyplot.figure(13)
	plt.plot(corrDataUniform[0])
	plt.plot(corrDataUniform[1])
	plt.plot(corrDataUniform[2])
	plt.plot(corrDataUniform[3])
	"""	

	matplotlib.pyplot.figure(13)	
	plt.plot(uniformBinMinima,corrDataUniform[3]/corrDataUniform[3][(uniformBinMinima==5.0).argmax()])
	#plt.plot(uniformBinMinima,corrDataUniform[0])
	plt.ylim((0,2))

	if (options.frameGroup != 0):
		matplotlib.pyplot.savefig(dir+"/G(r) uniform divided "+str(startTime)+".pdf")
	else:
		matplotlib.pyplot.savefig(dir+"/G(r) uniform divided.pdf")

	matplotlib.pyplot.figure(14)	
	plt.plot(uniformBinMinima,corrDataUniform[3])
	#plt.plot(uniformBinMinima,corrDataUniform[0])
	plt.ylim((0,2))

	if (options.frameGroup != 0):
		matplotlib.pyplot.savefig(dir+"/G(r) uniform "+str(startTime)+".pdf")
	else:
		matplotlib.pyplot.savefig(dir+"/G(r) uniform.pdf")

	uniformDists.append(corrDataUniform)

	"""
		plt.figure()
		#plt.xlabel("r/L(t)   (n(L(t), t) = 60)")
		#plt.ylabel("n(r/L(t), t)")
		plt.xlabel("r")
		plt.ylabel("n(r, t)")
		tmp=numpy.array(uniformDists)
		for corrDataUniform in tmp:
		#interpFunction = scipy.interpolate.interp1d(uniformBinMinima, corrDataUniform[3]*uniformNorm, kind="linear")
		#solveFunction = lambda x: interpFunction(x) - 60.0
		#lengthScale = scipy.optimize.fsolve(solveFunction, 8.0)
		#plt.plot(uniformBinMinima / lengthScale,corrDataUniform[3]*uniformNorm,"-^")
		plt.plot(uniformBinMinima ,corrDataUniform[3]*uniformNorm)
		
		matplotlib.pyplot.savefig(dir+"/G(r) uniform test.pdf")
		"""		
	#	plt.plot(uniformBinMinima, uniformNorm*a[3]/(a[3]+a[2]))

	#nearestNeighborData=numpy.fromfile(corrFile, dtype=numpy.int32, count=4)
	#deadPurity=(nearestNeighborData[0]/(nearestNeighborData[0]+nearestNeighborData[1]))
	#alivePurity=(nearestNeighborData[3]/(nearestNeighborData[2]+nearestNeighborData[3]))
	#corrFile.close()

	#print alivePurity, deadPurity

	fileNames=["dead-dead","dead-alive","alive-dead","alive-alive"]
	fileNamesGofR=["dead-dead","dead-alive","alive-dead","alive-alive"]
	fileNamesNematic=["dead-dead","dead-alive","alive-dead","alive-alive"]
	#			fileNamesGofR=["dead-dead","","","alive-alive"]
	#			fileNamesNematic=["dead-dead","","",""]
	fileNamesPolar=["","dead-alive","","alive-alive"]

	i=0

	for (data, dataName, dataIndex) in ((corrData, "G(r)", 1),
		#(corrDataOriented, "G(r) nematic orientation-density", 2), 
		#(corrDataOrientedDirectional, "G(r) oriented polar orientation-density", 3),
		(numpy.nan_to_num(corrDataOriented/corrData), "G(r) nematic orientation", 4),
		(numpy.nan_to_num(corrDataOrientedDirectional/corrData), "G(r) polar orientation", 5)):
		
		maxColor=numpy.max(corrData)
		minColor=numpy.min(corrData)
		
		#matplotlib.pyplot.figure(dataIndex*2)
		#matplotlib.pyplot.suptitle(dataName)
		
		if (dataName == "G(r) polar orientation"):
			fileNames=fileNamesPolar
		elif (dataName == "G(r) nematic orientation"):
			fileNames=fileNamesNematic
		elif (dataName == "G(r)"):
			fileNames=fileNamesGofR
		
		
		for (index, name) in enumerate(fileNames):
			if (name == ""):
				continue
			
			posImage=(data[index].transpose()[::-1])
			
			fig=matplotlib.pyplot.figure(i, figsize=(5,3.7))
			matplotlib.pyplot.xlabel("x")
			matplotlib.pyplot.ylabel("y")
			fig.subplots_adjust(left=0.2)
			#matplotlib.pyplot.margins(tight=True)
			#matplotlib.pyplot.title(dataName+" "+name)
			i += 1
			
			if (dataName == "G(r) polar orientation"):
				#imgMax = 1.0
				#if (numpy.min(posImage)>=0):
				#	imgMin = 0.0
				#else:
				#	imgMin = -1.0
				if (name == "alive-alive"):
					imgMax = 1.0
					imgMin = -0.1
				elif (name == "dead-alive"):
					imgMax = 0.3
					imgMin = -0.3
				imgNorm = matplotlib.colors.Normalize(vmin=imgMin,vmax=imgMax)
			elif (dataName == "G(r) nematic orientation"):
				imgMax = 1.0
				imgMin = 0.0
				imgNorm = matplotlib.colors.Normalize(vmin=imgMin,vmax=imgMax)
			elif (dataName == "G(r)"):
				imgMax = numpy.max(posImage) #maxColor #25 #maxColor
				imgMin = 0.0 #0.1
				imgNorm = matplotlib.colors.Normalize(vmin=imgMin,vmax=imgMax)
			#imgNorm = matplotlib.colors.LogNorm(vmin=imgMin,vmax=imgMax)
			
			
			#matplotlib.pyplot.figure(dataIndex*2)
			#matplotlib.pyplot.subplot(220+index+1,title=name)
			#matplotlib.pyplot.title(name)
			
			matplotlib.pyplot.imshow(posImage, interpolation='bilinear',
									 norm=imgNorm,
									 cmap=matplotlib.cm.gist_rainbow,#matplotlib.cm.gist_rainbow,
									 extent=imageExtent)
			matplotlib.pyplot.colorbar()
			
			if (dataIndex==1):
				xValues = numpy.arange(minX+1.0/(2*resolution),maxX-1.0/(2*resolution)+1.0/resolution,1.0/resolution)
				yValues = -numpy.arange(minY+1.0/(2*resolution),maxY-1.0/(2*resolution)+1.0/resolution,1.0/resolution)
			
			smoothData = (numpy.roll(numpy.roll(posImage,1,axis=0),1,axis=1) + numpy.roll(numpy.roll(posImage,1,axis=0),1,axis=-1) +
						  numpy.roll(numpy.roll(posImage,1,axis=0),-1,axis=1) + numpy.roll(numpy.roll(posImage,1,axis=0),-1,axis=-1) +
						  2*numpy.roll(posImage,1,axis=0) + 2*numpy.roll(posImage,-1,axis=0) + 
						  2*numpy.roll(posImage,1,axis=1) + 2*numpy.roll(posImage,-1,axis=1) + 
						  3*posImage)/15.0
	
			if (dataName == "G(r)"):
				if (name == "alive-alive"):
					contourLabels = (1,1.5,2,3,4)
				elif (name == "dead-dead"):
					contourLabels = (1,1.1,1.2,1.5)
				elif (name == "alive-dead"):
					contourLabels = (0.6,0.8,1)
				elif (name == "dead-alive"):
					contourLabels = (0.6,0.8,1,1.1)
			elif (dataName == "G(r) polar orientation"):
				if (name == "alive-alive"):
					contourLabels = (0.0,0.1,0.5,0.7)
				elif (name == "dead-alive"):
					contourLabels = (-0.1,0,0.1)
			elif (dataName == "G(r) nematic orientation"):
				contourLabels = (0.2,0.4,0.6,0.8,1.0)
			
			countours = matplotlib.pyplot.contour(xValues, yValues, smoothData, contourLabels, colors='k',
									  antialiased=True, linewidths=0.5)

			plt.clabel(countours, fontsize=12, inline=1, fmt="%1.1f")
			
			matplotlib.pyplot.savefig(dir+"/"+string.replace(string.replace(dataName,"G(r)","GofR")," ","-")+"-"+name+".pdf")



if (storeClumps):
	plt.figure()
	fullData=numpy.concatenate(allClumpInfo)
	goodData = (startTime <= fullData['time']) & (stopTime > fullData['time'])
	fullData = fullData[goodData]

	"""
	plt.scatter(fullData['time'],fullData['averageClusterSize'])
	if allTime:
		plt.savefig(baseDir+"/averageClusterSize.pdf")
		numpy.savetxt(baseDir+"/averageClusterSize.txt",fullData[['time','averageClusterSize','weightedAverageClusterSize','cutoffAverageClusterSize']],fmt=("%i","%f","%f","%f"))
	else:
		plt.savefig(baseDir+"/averageClusterSize-selected.pdf")
		numpy.savetxt(baseDir+"/averageClusterSize-selected.txt",fullData[['time','averageClusterSize','weightedAverageClusterSize','cutoffAverageClusterSize']],fmt=("%i","%f","%f","%f"))
	"""

	plt.figure()
	clumpSizes = numpy.arange(1,info['partNum']+1)


	allClumpSizeBins = fullData['sizes'].sum(axis=0)
	aliveClumpSizeBins = fullData['aliveSizes'].sum(axis=0)
	deadClumpSizeBins = fullData['deadSizes'].sum(axis=0)

	
	meanClumpSize = numpy.sum(clumpSizes*allClumpSizeBins/sum(allClumpSizeBins))
	massMeanClumpSize = numpy.sum(clumpSizes*clumpSizes*allClumpSizeBins/sum(allClumpSizeBins*clumpSizes))
	
	clumpWidths = numpy.sqrt(clumpSizes * (rodLength+1.0) / 2.0)
	boxSize = 10.0
	clumpAreas = clumpWidths*clumpWidths*2.0
	frameCount = len(fullData)
	totalArea = numpy.sum(allClumpSizeBins * clumpAreas / frameCount)
	
	interiorAreaFractions=list()
	interiorAreaFractions.append((meanClumpSize, massMeanClumpSize))
	for boxSize in (4.35, 5.02, 8.46, 13.05, 45.68, 9.62, 12.12, 14.84, 20.98, 31.48, 19.11, 22.07, 36.04, 40.78, 56.07):
		clumpInteriorAreas = (clumpWidths - boxSize)*(clumpWidths*2.0 - boxSize) * (clumpWidths > boxSize)
		interiorArea = numpy.sum(allClumpSizeBins * clumpInteriorAreas / frameCount)
		interiorAreaFraction = interiorArea/totalArea
		interiorAreaFractions.append((boxSize, interiorAreaFraction))
	
	numpy.savetxt(baseDir+"/fullAverageClusterSize-selected.txt", interiorAreaFractions, fmt=("%f","%f"))

	maxClumpSize = numpy.max(numpy.nonzero(allClumpSizeBins))

	numpy.savetxt(baseDir+"/clusterSizeDistribution.txt", numpy.transpose(numpy.array([clumpSizes,allClumpSizeBins])), fmt=("%d","%d"))

	#if (numpy.max(deadClumpSizeBins)>0):
	#	plt.loglog(clumpSizes,deadClumpSizeBins)
	#
	#if (numpy.max(aliveClumpSizeBins)>0):
	#	plt.loglog(clumpSizes,aliveClumpSizeBins)
	
	if (numpy.max(allClumpSizeBins)>0):
		plt.loglog(clumpSizes[0:maxClumpSize],allClumpSizeBins[0:maxClumpSize]/allClumpSizeBins.sum())

	if allTime:
		plt.savefig(baseDir+"/clumpDist.pdf")
	else:
		plt.savefig(baseDir+"/clumpDist-selected.pdf")
	
	
	plt.figure()
	
	weightedClumpSizeBins = allClumpSizeBins[0:maxClumpSize] * clumpSizes[0:maxClumpSize]
	if (numpy.max(allClumpSizeBins)>0):
		plt.loglog(clumpSizes[0:maxClumpSize],weightedClumpSizeBins[0:maxClumpSize]/weightedClumpSizeBins.sum())

	if allTime:
		plt.savefig(baseDir+"/weightedClumpDist.pdf")
	else:
		plt.savefig(baseDir+"/weightedClumpDist-selected.pdf")
