#!/usr/bin/python

import subprocess
import os.path

def read(dir, findMore):
	if (findMore):
		simulationDirectories = subprocess.Popen("ls -d "+dir+"*/info.dat", shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')
		previousSimulations = [x.replace("info.dat","").split(dir.split("/")[-1])[-1].replace("/","").replace("-","") for x in simulationDirectories if len(x)>0]
		
		if (len(previousSimulations)==0):
			lastSimulationIndex = 1
		else:
			previousSimulations = [int(p) for p in previousSimulations if len(p)>0]
			if len(previousSimulations)>0:
				lastSimulationIndex = max(previousSimulations)
			else:
				lastSimulationIndex = 1
	else:
		lastSimulationIndex = 1
	
	info = []
	for i in (range(lastSimulationIndex)):
		if i == 0:
			postfix = ""
		else:
			postfix = "-" + str(i+1)
		
		infoFile=open(os.path.expanduser(dir+postfix+"/info.dat"))

		for line in infoFile:
			parsed = line.split()
			info.append(dict(
				startTime=int(parsed[0]),
				stopTime=int(parsed[1]),
				particleCount=int(parsed[2]),
				endSl=float(parsed[3]),
				rodLength=float(parsed[4]),
				outputSkip=int(parsed[5]), 
				startSl=float(parsed[6]),
				timeStep=float(parsed[7]),
				circleRadius=float(parsed[8]),
				dir=dir+postfix,
		 		index=i+1
			))
		infoFile.close()
	return dict(
		startSl=info[0]['startSl'],
		rodLength=info[0]['rodLength'],
		frameSkip=info[0]['outputSkip'],
		absStartTime=info[0]['startTime'],
		partNum=info[0]['particleCount'],
		absStopTime=info[-1]['stopTime'],
		timeStep=info[0]['timeStep'],
		allInfo=info
	)


def whichDir(step, info):
	possibilities = [sub.index for sub in info if (step >= sub['absStartTime'] and step < sub['absStopTime'])]
	if (len(possibilities) != 1):
		return -1
	else:
		return possibilities[0]


