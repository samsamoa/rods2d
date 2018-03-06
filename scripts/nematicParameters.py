#!/usr/bin/python

"""
Calculates the polar and nematic order parameters over boxes in each 
frame, as well as over the whole frame.
"""
from __future__ import division
from optparse import OptionParser
import numpy
import readInfo
from progressbar import ProgressBar

parser = OptionParser()
parser.add_option("--dir", dest="dir", default="./")
parser.add_option("--startTime", dest="startTime", type="int", default=-1)
parser.add_option("--noBoxes", dest="noBoxes", action="store_true", default=False)
parser.add_option("--sideNumBoxes", dest="sideNumBoxes", type="int", default=10)
parser.add_option("--stopTime", dest="stopTime", type="int", default=-1)
parser.add_option("--analysisFrameSkip", dest="analysisFrameSkip", type="int", default=10000)
parser.add_option("--singleFrame", action="store_true", dest="singleFrame", default=False)
parser.add_option("--debug", action="store_true", dest="debug", default=False)
parser.add_option("--save", action="store_true", dest="save", default=False)

(options, args) = parser.parse_args()

dir=options.dir
startTime=options.startTime
stopTime=options.stopTime
analysisFrameSkip=options.analysisFrameSkip
singleFrame=options.singleFrame
debug=options.debug
sideNumBoxes=options.sideNumBoxes
save=options.save
doBoxes=(not options.noBoxes)

info=readInfo.read(dir,False)
startSl=info['startSl']
rodLength=info['rodLength']
frameSkip=info['frameSkip']
absStartTime=info['absStartTime']
partNum=info['partNum']
timeStep=info['timeStep']


import matplotlib
if (save):
	matplotlib.use("pdf")
import matplotlib.pyplot as pyplot


if (analysisFrameSkip < frameSkip):
	print "Frame skip too low, sorry. " , analysisFrameSkip, frameSkip
	exit()

if startTime==-1:
	startTime=info['absStartTime']
if stopTime==-1:
	stopTime=info['absStopTime']-1



file = open(dir+"/coords-bin.dat", "rb")
rodRecordType=numpy.dtype((numpy.double, (partNum,3)))
rodRecordSize = 24 # bytes

try:
	aliveRecord=numpy.dtype((bool,partNum))
	aliveData=numpy.fromfile(dir+"/alive.dat",dtype=aliveRecord,count=-1)[-1]
except IOError:
    aliveData=numpy.zeros(partNum, bool)

if startTime==stopTime:
	times = [startTime]
else:
	times = range(startTime, stopTime, analysisFrameSkip)

if (doBoxes):
	boxCounts = numpy.zeros((len(times), sideNumBoxes, sideNumBoxes))
	
	boxNematicTensors = numpy.zeros((len(times), sideNumBoxes, sideNumBoxes, 2, 2))
	boxNematicOrderParameters = numpy.zeros((len(times), sideNumBoxes, sideNumBoxes))
	boxNematicDirectors = numpy.zeros((len(times), sideNumBoxes, sideNumBoxes, 2))
	
	boxPolarOrderParameters = numpy.zeros((len(times), sideNumBoxes, sideNumBoxes))
	boxPolarDirectors = numpy.zeros((len(times), sideNumBoxes, sideNumBoxes, 2))

frameNematicTensors = numpy.zeros((len(times), 2, 2))
frameNematicOrderParameters = numpy.zeros((len(times)))

framePolarDirectors = numpy.zeros((len(times), 2))
framePolarOrderParameters = numpy.zeros((len(times)))


tmp = numpy.zeros(partNum, dtype=bool)

boxWidth = startSl/sideNumBoxes
centerWidth = startSl - boxWidth
sideInc = centerWidth/(sideNumBoxes-1)

progress = ProgressBar(maxval=len(times))
for tIdx, t in enumerate(times):
	progress.update(tIdx)
	file.seek(rodRecordSize*partNum*((t-absStartTime)//frameSkip))
	data = numpy.fromfile(file, dtype=rodRecordType, count=1)[0]
	pos = data[:,(0,1)]
	theta = data[:,2]
	sin = numpy.ma.array(numpy.sin(theta));
	cos = numpy.ma.array(numpy.cos(theta));

	frameNematicTensors[tIdx] = numpy.array([[(cos*cos).mean() - 0.5, (sin*cos).mean()], [(cos*sin).mean(), (sin*sin).mean() - 0.5]])*2
	frameNematicOrderParameters[tIdx] = numpy.max(numpy.linalg.eig(frameNematicTensors[tIdx])[0])

	framePolarDirectors[tIdx] = numpy.array([cos.mean(), sin.mean()])
	framePolarOrderParameters[tIdx] = numpy.linalg.norm(framePolarDirectors[tIdx])
	
	if (doBoxes):
		for i in range(sideNumBoxes):
			for j in range(sideNumBoxes):
				sin.mask = cos.mask = numpy.logical_not(numpy.all((pos > [sideInc*i,sideInc*j]) & (pos < [sideInc*i+boxWidth,sideInc*j+boxWidth]),1))
				if numpy.ma.count_masked(sin) == partNum:
					boxNematicTensors[tIdx, i, j] = numpy.array([[0,0],[0,0]])
					boxPolarDirectors[tIdx, i, j] = numpy.array([0,0])
				else:
					boxNematicTensors[tIdx, i, j] = numpy.array([[(cos*cos).mean() - 0.5, (sin*cos).mean()], [(cos*sin).mean(), (sin*sin).mean() - 0.5]])*2
					boxPolarDirectors[tIdx, i, j] = numpy.array([cos.mean(), sin.mean()])
				
				eig = numpy.linalg.eig(boxNematicTensors[tIdx,i,j])

			boxCounts[tIdx, i, j] = sin.count() # How many rods in each box
				
				boxPolarOrderParameters[tIdx, i, j] = numpy.linalg.norm(boxPolarDirectors[tIdx, i, j])
				boxNematicOrderParameters[tIdx, i, j] = numpy.max(eig[0])
				if (eig[0][0] > eig[0][1]):
					boxNematicDirectors[tIdx,i,j] = eig[1][:,0]
				else:
					boxNematicDirectors[tIdx,i,j] = eig[1][:,1]
			

file.close()
progress.finish()

if (doBoxes):
	#pyplot.figure(1)
	
	frame=len(times)-1
	X,Y = numpy.meshgrid( numpy.arange(0,startSl,sideInc)+sideInc/2.0,numpy.arange(0,startSl,sideInc)+sideInc/2.0 )
	pyplot.subplot(1,2,1)
	polar = matplotlib.pyplot.quiver( X, Y, boxPolarDirectors[frame,:,:,0], boxPolarDirectors[frame,:,:,1], pivot="middle", color="r")
	pyplot.subplot(1,2,2)
	nematic = matplotlib.pyplot.quiver( X, Y, boxNematicDirectors[frame,:,:,0], boxNematicDirectors[frame,:,:,1], pivot="middle", color="b", headlength=0, headwidth=0)
	

pyplot.figure(2)
pyplot.ylim(0,1)
pyplot.plot(times, frameNematicOrderParameters, 'r-', times, framePolarOrderParameters, 'b-')


if save:
	pyplot.savefig(dir+"/orderParam.pdf")
	numpy.savetxt(dir+"/orientationOrderParam.dat", numpy.dstack((times, frameNematicOrderParameters, framePolarOrderParameters))[0], fmt="%f")
else:
	pyplot.show()

#pickleFile = open(filename+"boxes.dat","wb")
#pickle.dump(dpts, pickleFile)
#pickleFile.close()
#for t in range(int(len(times)/10)):
#	dptsqs = dpts[:,t:t+10]**2
#	dptms = numpy.mean(dpts[:,t:t+10], axis=1)
#	dptsqms = numpy.mean(dptsqs, axis=1)
#	stdevs = numpy.sqrt(dptsqms - dptms**2)
#	
#	naverages = numpy.mean(numpy.mean(dptms, axis=1),axis=1)
#	nstdevavgs = numpy.mean(numpy.mean(stdevs, axis=1),axis=1)
#
#	for x in zip(naverages, nstdevavgs):
#		print t, x[0], x[1]
#

