"""

python -i ~/Work-Hagan/Rods2d/trunk/scripts/clumpCalc.py --dir ~/Work-Hagan/sim/sca-n100-f0.0-l20-s3/d0.0100-n1600/f1.0-a1.0/ --endTimeFraction 0.1 --speedup 10 --save --analysis --width 100 --height 100
"""
from __future__ import division
import struct
from optparse import OptionParser
import math
import subprocess
import readline

import readInfo

parser = OptionParser()
parser.add_option("--dir", dest="dir", default="./")
parser.add_option("--startTime", dest="startTime", type="int", default=-1)
parser.add_option("--stopTime", dest="stopTime", type="int", default=-1)
parser.add_option("--analysisFrameSkip", dest="analysisFrameSkip", type="int", default=10000)
parser.add_option("--boxes", dest="boxes", type="int", default=50)
parser.add_option("--debug", action="store_true", dest="debug", default=False)
parser.add_option("--allTime", action="store_true", dest="allTime", default=False)
parser.add_option("--endTimeFraction", dest="endTimeFraction", type="float", default=0.0)
parser.add_option("--fakeAlive", dest="fakeAlive", default="")
parser.add_option("--motileFraction", dest="motileFraction", default="0.5")

(options, args) = parser.parse_args()

dir=options.dir
startTime=options.startTime
stopTime=options.stopTime
analysisFrameSkip=options.analysisFrameSkip
debug=options.debug
allTime=options.allTime
endTimeFraction=options.endTimeFraction
fakeAlivePrefix=options.fakeAlive
boxes=options.boxes
motileFraction=options.motileFraction

info=readInfo.read(dir, False)
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
	stopTime = absStopTime
elif endTimeFraction != 0.0:
	#print "Fraction of time."
	stopTime = absStopTime
	startTime = absStopTime - (absStopTime - absStartTime)*endTimeFraction
elif (startTime == -1 or stopTime == -1):
	print str(absStartTime), " ", str(absStopTime), " ", str(frameSkip)
	startTime = int(raw_input("start: "))
	stopTime = int(raw_input("stop: "))

if (analysisFrameSkip % frameSkip != 0):
	print "Error: frame skips don't match"

# This is used so that we can pretend some of the rods are alive for passive simulations
if (fakeAlivePrefix != ""):
	fluctuationFile = "mixed-fluctuations-"+motileFraction+".dat"
	aliveFile = fakeAlivePrefix+"/alive.dat"
else:
	fluctuationFile = "fluctuations.dat"
	aliveFile = dir+"/alive.dat"

command = ("~/Work-Hagan/Rods2d/trunk/build/fluctuate -n " + str(partNum) + " -l " + str(rodLength) + " -L " + str(startSl) + " -t " + str(absStartTime) + " -T " + str(stopTime) +
		" -s " + str(frameSkip) + " -c " + str(startTime) + " -f " + dir + "/coords-bin.dat -b " + dir +
		"/boxes.dat -a " + aliveFile + " --speedup " + str(analysisFrameSkip/frameSkip) +
		" -F " + dir + "/" + fluctuationFile + " -G " + dir + "/GofR.dat" + " -E " + dir + "/fluctuationExponent.dat" + 
		" -A " + dir + "/meanSqDisp.dat -C " + str(boxes) + " -B 40")

print command

frameCount = math.floor((stopTime-startTime)/(analysisFrameSkip) + 1)
output=subprocess.Popen(command,shell=True).communicate()
