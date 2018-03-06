/*
 *  RodSimulation.cpp
 *  bd
 *
 *  Created by Sam McCandlish on 12/22/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

/* Non-dimensionalization:
(14.8nm)^2 * 2Pi * (10^(-3) pascal*second)(0.88 um) / (boltzmann constant 300K ln(60)) 
 
 tau = sigma^2 2π (viscosity)(length) / (kT ln(length/diameter))

Fd = 10pN * 14.8nm/kT
D=14.8nm

sigma/tau=D/sigma=kt/(sigma gamma)=kT/(2 Pi viscosity D^2 (L/sigma)/Log[L/sigma])

 x = distance/sigma
 sigma = molecular diameter
 t = time/tau
 tau = sigma^2 / D
 D = kt/squiggle
 squiggle = m gamma
 (2Dt) = sqrt(2) sigma sqrt(time/tau)
 f = force sigma/kT
 E = energy/kT
 Isotropic-Nematic Transition ** Von Otter, reflecting, overlapping ** 
 */

#include "RodSimulation.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <ctime>

# ifdef __INTEL_COMPILER 
# include <mathimf.h> 
# else 
# include <math.h> 
# endif 

//#include "Rand.cpp"

using namespace std;


RodSimulation::RodSimulation(int _particleCount, double _rodLength, double startDensity, double _boundaryRadius, // Initial parameters
					int randomSeed, double endDensity, double endBoundaryRadius, double _timeStep, int _stopTime, int _outputSkip, double _selfPushForce, double motorCutoff, double friction, double aliveFraction) {	// Simulation parameters
	particleCount = _particleCount;
	rodLength = _rodLength;
	timeStep = _timeStep;
	stopTime = _stopTime;
	outputSkip = _outputSkip;
	selfPushForce = _selfPushForce;
	motorForceCutoffSq = motorCutoff*motorCutoff;
	boundaryRadius = _boundaryRadius;
	alignmentParameter = friction;
	
	randomGen = new Random(randomSeed);
	int sideNumber = (int)sqrt(particleCount);
	double randomSideLength = sideNumber*(rodLength+2); // If we have enough space, orient the rods randomly
	
	particleList.reserve(particleCount);
	
	double initialBoxSideLength; // The size of the box we put particles in initially
	
	if (boundaryRadius == -1) {
		resizeBox = (startDensity != endDensity);	
		sideLength = sqrt((double)particleCount/startDensity);

		if (resizeBox)
			targetSideLength = sqrt((double)particleCount/endDensity);
		else
			targetSideLength = -1;
		
		initialBoxSideLength = sideLength;
	}
	else {
		resizeBox = (boundaryRadius != endBoundaryRadius);
		sideLength = boundaryRadius*2;
		
		if (resizeBox)
			targetSideLength = endBoundaryRadius*2;
		else
			targetSideLength = -1;

		initialBoxSideLength = sqrt(2.0)*boundaryRadius;  // We want to try to put the rods in the inscribed square		
	}


	initValues();

	double spacing = (initialBoxSideLength)/(double)sideNumber;
	bool useRandomOrientation = (initialBoxSideLength >= randomSideLength);

	Rod* tempP;
	for (int p=0; p<particleCount; p++) {
		double theta, x, y;
		if (useRandomOrientation)
			theta = randomGen->RealIE()*2*M_PI;
		else
			theta = 0.6;
		
		x = (p % sideNumber) * spacing;
		y = (p / sideNumber) * spacing;
		

		if (boundaryRadius != -1) {
			x -= (initialBoxSideLength*0.5);
			y -= (initialBoxSideLength*0.5);
			
			if (sideNumber % 2 == 0) {
				x += spacing/2.0;
				y += spacing/2.0;
			}
		}		
		
		tempP = new Rod(x, y, theta, rodLength2, p);
		tempP->calculateU();
		particleList.push_back(tempP);
	}
	
	if (aliveFraction != -1)
		killRods(aliveFraction);

	
	currentTime = 0;
	
	openOut(stepOut, "steps.dat");
	openOut(coordBinOut, "coords-bin.dat");
	openOut(forceOut, "forces.dat");
}

RodSimulation::RodSimulation(int randomSeed, double _endDensity, double endBoundaryRadius, double _timeStep, int _stopTime, int _outputSkip, double _selfPushForce, double motorCutoff,	double friction, double aliveFraction, // Simulation parameters
									istream &savedState, bool resetClock, int multiplySideLength) {	// Resume parameters
	randomGen = new Random(randomSeed);
	//randomGen->LoadState(savedState);
	alignmentParameter = friction;
		
	/* Read the saved state */
	loadState(savedState);
	if (multiplySideLength != 1) {
		cout << sideLength << "," << particleCount << endl;
		multiplyPBC(multiplySideLength);
		cout << sideLength << "," << particleCount << endl;
	}
	
	openOut(stepOut, "steps.dat");
	openOut(coordBinOut, "coords-bin.dat");
	openOut(forceOut, "forces.dat");
	
	if (_stopTime != -1) {
		if (resetClock)
			stopTime = currentTime + _stopTime;
		else
			stopTime = _stopTime;
	}
	if (_selfPushForce != -1)
		selfPushForce = _selfPushForce;
	if (_timeStep != -1)
		timeStep = _timeStep;
	if (_outputSkip != -1)
		outputSkip = _outputSkip;
	if (friction != -1)
		alignmentParameter = friction;
	
	motorForceCutoffSq = motorCutoff*motorCutoff;

	if (boundaryRadius == -1) {
		resizeBox = (_endDensity != -1);
		if (resizeBox)
			targetSideLength = pow((double)particleCount/_endDensity, 0.5);
		else
			targetSideLength = -1;
	}
	else {
		sideLength = boundaryRadius*2;
		resizeBox = (endBoundaryRadius != -1);	

		if (resizeBox)
			targetSideLength = endBoundaryRadius * 2;
		else
			targetSideLength = -1;
	}
	
	if (aliveFraction != -1)
		killRods(aliveFraction);
	
	initValues();
}

RodSimulation::~RodSimulation() {
	stepOut.close();
	coordBinOut.close();
}

void RodSimulation::initValues() {
	if (selfPushForce == -1)
		selfPushForce = 0;
	if (alignmentParameter == -1)
		alignmentParameter = 0;
	frictionCutoffSq = 1.5*1.5; // This will need to be tuned!!	

	halfSideLength = sideLength*0.5;
	rodLength2 = rodLength/2.0;
	cutoff = pow(2.0,1.0/6.0);
	cutoffSq = pow(2.0,1.0/3.0);
	quickCutoffSq = pow(rodLength + cutoff, 2.0);
	
	periodicBoundary = (boundaryRadius == -1); // If there is no boundary radius, use PBC
	initSideLength = sideLength;
	
	
	boundaryRadiusSq = boundaryRadius*boundaryRadius;
	
	XParrFrictionCoeff = 1;
	XPerpFrictionCoeff = 2; // = 2 * XParrFrictionCoeff;
	TFrictionCoeff = rodLength*rodLength/6.0; // = rodlength^2/6 * XParrFrictionCoeff
	
	//randomXParrStdev = sqrt(2.0 * timeStep / XParrFrictionCoeff);
	//randomXPerpStdev = sqrt(2.0 * timeStep / XPerpFrictionCoeff);
	//randomTStdev = sqrt(2.0 * timeStep / TFrictionCoeff);
	
	double neighborListBuffer = 2;
	neighborListBufferSq = neighborListBuffer*neighborListBuffer;
	neighborListBufferSq4 = neighborListBufferSq*0.25;
	
	double neighborListCutoff = neighborListBuffer + cutoff;
	neighborListCutoffSq = neighborListCutoff*neighborListCutoff;
	
	double neighborListQuickCutoff = rodLength + neighborListCutoff;
	neighborListQuickCutoffSq = neighborListQuickCutoff*neighborListQuickCutoff;
	
	motorGrid = new MotorGrid(this, 100, selfPushForce, 20);
}


int RodSimulation::runSimulation () {
	int status = 0;
	//debugOut.open("debug.dat", ios::out | ios::binary);
	ofstream descOut;
	openOut(descOut, "desc.txt");
	descOut << "Starting simulation:" << endl
		<< "Particle Count: " << particleCount << endl
		<< "Rod Length: " << rodLength << endl
		<< "Side Length: " << sideLength << endl
		<< "Time: " << currentTime << endl
		<< "Stop Time: " << stopTime << endl
		<< "Time Step: " << timeStep << endl
		<< "Circular Boundary: " << boundaryRadius << endl
		<< "Output Skip: " << outputSkip << endl
		<< "Self-Propulsion Force: " << selfPushForce << endl
		<< "Initial Side Length: " << initSideLength << endl
		<< "Target Side Length: " <<targetSideLength << endl
		<< "Resize Box: " << resizeBox << endl
		<< "Max Rotation: " << thetaDisp_max << endl
		<< "Max Displacement (Squared): " << posSqDisp_max << endl
		<< "Motor Force Cutoff: " << motorForceCutoffSq << endl
		<< "Parallel Friction Coefficient: " << XParrFrictionCoeff << endl
		<< "Random Parallel Displacement Stdev: " << randomXParrStdev << endl
		<< "Perpendicular Friction Coefficient: " << XPerpFrictionCoeff << endl
		<< "Random Perpendicular Displacement Stdev: " << randomXPerpStdev << endl
		<< "Theta Friction Coefficient: " << TFrictionCoeff << endl
		<< "Random Theta Stdev: " << randomTStdev << endl
		<< "Friction: " << alignmentParameter << endl
		<< endl << "...starting..." << endl;
			
	
	cout << "Starting simulation with non-dimensionalized parameters:\n";
	cout << "Particle Count: " << particleCount 
		<< ", Rod Length: " << rodLength 
		<< ", Side Length: " << sideLength
		<< ", Time: " << currentTime
		<< ", Stop Time: " << stopTime << "\n";
	
	/* This uses the following formula:
	 * tau = sigma^3 2π (viscosity)(length) / (kT ln(length))
	 * with sigma=14.8nm, viscosity=10^-3Pa, T=300K
	 */
	double realTimeStep = timeStep*4.917678*rodLength/log(rodLength);
	cout << "Actual time step: " << realTimeStep << " microseconds.\n";
	cout << "Actual stop time: " << stopTime*realTimeStep*0.000001 << " seconds.\n";
	
	ofstream aliveOut;
	openOut(aliveOut, "alive.dat");
	for (int i=0; i<particleCount; i++) {
		aliveOut.write((char *)&particleList[i]->alive, sizeof(bool));
	}
	aliveOut.close();

	
	//openOut(clumpOut, "clumps.dat");
	
	// Main loop
	int startTime = currentTime;
	stop = false;

	if (!checkConsistency()) {
		cerr << "Initial consistency check failed. Exiting...\n";
		return 1;
	}
	
	failCount = 0;
	failFactor = 1;
	
	saveState(true);
	
	for (; currentTime < stopTime; currentTime++) {
		if (currentTime % 10000 == 0 && failCount == 0) {
			// We've been doing well.
			saveState(true);
			if (failFactor > 1)
				failFactor -= 1;
		}
		else if (failCount > 0) {
			// We're not doing so well. Do the failsafe.
			//cerr << "Initiating failsafe...\n";
			revertState();
			if (failFactor <= 5)
				failFactor = failFactor + 1;
			else if (failFactor <= 10)
				failFactor = failFactor + 2;
			else
				failFactor = failFactor + 3;
			failCount = 0;
		}
		
		if (failFactor > 20) {
			cerr << "Recovery using failsafe has stalled. Bailing...\n";
			status=1;
			break;
		}
		
		
		if (currentTime % outputSkip == 0)
			output = true;
		else
			output = false;
		
		
		if (stop) {
			cout << "Stopping at end of step " << currentTime-1 << "\n";
			break;
		}


		for (int f=0; f<failFactor; f++) {
			if (resizeBox) {
				// Go to the target linearly (0=initial, stopTime)
				// Also, eat 2 slices of cake
				// Note that the last scaling occurs at stopTime-1
				//  (which is why this is so messy).
				double newSideLength = (double)(currentTime + (double)f/failFactor - startTime)*targetSideLength/(stopTime - startTime - 1) + (double)(stopTime - 1 - (currentTime + (double)f/failFactor))*initSideLength/(stopTime - startTime - 1);
				scaleBoxTo(newSideLength);
			}
			
			/*** Second Order Predictor Corrector ***/
			calculateForces(true);
			applyForces(true); // Make first guess
			if (failCount > 0)
				break;
			calculateForces(false);
			applyForces(false); // Correct guess
			if (failCount > 0)
				break;
			/****************************************/
		}
		
		if (failCount > 0)
			continue;
		
		if (output) {
			/** Status Update **/
			cout << "Step " << currentTime << " tsReduction " << failFactor << " %done " << (float)(currentTime-startTime)/(stopTime-startTime) << endl;
			
			/** Some calculations **/
			stepOut	<< currentTime*timeStep << " "
				<< totalDisp << " " 
				<< orientAllCorr << " "
				<< orientInitCorr << "\n";
			
			//updateClumps();
			
			/** Output the current position **/
			Rod *r;
			for (int p=0; p<particleCount; p++) {
				r = particleList[p];
				coordBinOut.write((char*)(r->pos), sizeof(double)*2);
				coordBinOut.write((char*)&(r->theta), sizeof(double));
			    forceOut.write((char*)r->force, sizeof(double)*2);
			    forceOut.write((char*)&(r->torque), sizeof(double));
			    //clumpOut.write((char*)&(r->clumpId), sizeof(int));
				
			}
			
			if (!checkConsistency()) {
				failCount += 1;
				continue;
			}
		}
		if (currentTime % 2000000 == 0) {
			// Back up, just in case
			saveState(false);
		}
	}
	
	
	coordBinOut.close();
	forceOut.close();
	//debugOut.close();
	stepOut.close();
	
	saveState(false);
	
	descOut << endl
		<< "Actual End Time: " << currentTime << endl
		<< "Actual End Side Length: " << sideLength << endl
		<< "-----------------------------------------------" << endl;
	descOut.close();
	
	ofstream infoOut;
	openOut(infoOut, "info.dat");	
	infoOut << startTime << " " << currentTime << " " << particleCount << " " << sideLength << " " << rodLength << " " << outputSkip << " " << initSideLength << " " << timeStep << " " << boundaryRadius << "\n";
	infoOut.close();

	cerr << "Done" << endl;
    return status;
}

void RodSimulation::saveState(bool saveTemp) {
	stringstream fileName;
	if (saveTemp) {
		fileName << "temp-state.dat";
		coordBinPos = coordBinOut.tellp();
		forceOutPos = forceOut.tellp();
		//clumpOutPos = clumpOut.tellp();
		stepOutPos = stepOut.tellp();
	}
	else
		fileName << "state-" << currentTime << ".dat";	

	ofstream savedState;
	savedState.open(fileName.str().c_str(), ios::out | ios::binary | ios::trunc);
	
	//randomGen->SaveState(savedState);
	
	// Changes:
	//   Version 2: - Add initial information and aliveness to rod info
	//				- Store the force and friction (for easy resuming)
	int version = 2;
	
	savedState.write((char*)&version, sizeof(int));
	
	savedState.write((char*)&particleCount, sizeof(int));
	savedState.write((char*)&rodLength, sizeof(double));
	savedState.write((char*)&sideLength, sizeof(double));
	
	for (int i=0; i<particleCount; i++) {
		particleList[i]->saveState(savedState, version);
	}
	
	savedState.write((char*)&currentTime, sizeof(int));
	savedState.write((char*)&stopTime, sizeof(int));
	savedState.write((char*)&timeStep, sizeof(double));
	savedState.write((char*)&boundaryRadius, sizeof(double));
	savedState.write((char*)&outputSkip, sizeof(int));
	savedState.write((char*)&selfPushForce, sizeof(double));
	savedState.write((char*)&alignmentParameter, sizeof(double));
	
	savedState.close();
}

void RodSimulation::revertState() {
	ifstream resumeFile;
	resumeFile.open("temp-state.dat", ios::in | ios::binary);
	if (resumeFile.fail()) {
		cerr << "Unable to open temporary resume file." << endl;
		abort();
	}
	loadState(resumeFile);
	resumeFile.close();
	
	coordBinOut.seekp(coordBinPos);
	forceOut.seekp(forceOutPos);
	//clumpOut.seekp(clumpOutPos);
	stepOut.seekp(stepOutPos);
	
	updateNeighbors();
}

void RodSimulation::loadState(istream &savedState) {
	int version;
	savedState.read((char*)&version, sizeof(int));
	if (version < 1 || version > 2) {
		cerr << "Could not read resume file (bad version)" << endl;
		abort();
	}
	
	savedState.read((char*)&particleCount, sizeof(int));
	
	for (int i=0; i<particleList.size(); i++) {
		delete particleList[i];
	}
	
	particleList.clear();
	particleList.reserve(particleCount);
	
	savedState.read((char*)&rodLength, sizeof(double));
	savedState.read((char*)&sideLength, sizeof(double));
	
	rodLength2 = rodLength*0.5;
	for (int i=0; i<particleCount; i++) {
		particleList.push_back(new Rod(savedState, rodLength2, i, version));
	}
	
	savedState.read((char*)&currentTime, sizeof(int));
	savedState.read((char*)&stopTime, sizeof(int));

	savedState.read((char*)&timeStep, sizeof(double));
	
	savedState.read((char*)&boundaryRadius, sizeof(double));
	savedState.read((char*)&outputSkip, sizeof(int));
	
	if (version == 1) {
		selfPushForce = 0;
		alignmentParameter = 0;
	}
	else if (version == 2) {
		savedState.read((char*)&selfPushForce, sizeof(double));
		savedState.read((char*)&alignmentParameter, sizeof(double));
	}
	
	if (savedState.fail()) {
		cerr << "Failed to load resume file." << endl;
		abort();
	}
}

// Make some copies of the current state and glue them together, thus
// creating a new state with side length sideFactor times longer.
void RodSimulation::multiplyPBC(int sideFactor) {
	vector<Rod*> newParticleList;

	Rod *tempP;
	
	for (int i=0; i<sideFactor; i++) {
		for (int j=0; j<sideFactor; j++) {
			for (int p=0; p<particleCount; p++) {
				tempP = new Rod(particleList[p]->pos[0] + i*sideLength,
								particleList[p]->pos[1] + j*sideLength,
								particleList[p]->theta, rodLength2,
								particleList[p]->idNum + i*particleCount + sideFactor*j*particleCount);
				
				tempP->alive = particleList[p]->alive;
				
				tempP->calculateU();
				
				newParticleList.push_back(tempP);
			}
		}
	}
	
	particleList = newParticleList;
	sideLength = sideLength * sideFactor;
	particleCount = particleCount * sideFactor * sideFactor;
}

void RodSimulation::killRods(double aliveFraction) {
	if (aliveFraction >= 1) {
		// Jesus rods!
		for (int i=0; i<particleCount; i++)
			particleList[i]->alive = true;
		return;
	}
	
	bool alive[particleCount];	
	int numberAlive = particleCount*aliveFraction;
	int aliveCount = 0;
	
	// This is not the right way to pick randomly; however, it's good enough.
	for (int i=0; i<particleCount; i++) {
		if (aliveCount == numberAlive) {
			alive[i] = false;
		}
		else {
			bool tmpAlive = (randomGen->RealIE() < aliveFraction);
			if (tmpAlive) {
				aliveCount += 1;
			}
			alive[i] = tmpAlive;
		}
	}
	
	for (int i=0; i<particleCount; i++) {
		if (aliveCount == numberAlive) {
			break;
		}
		else {
			if (alive[i] == (aliveCount > numberAlive)) {
				continue;
			}
			else {
				if (aliveCount > numberAlive) {
					aliveCount -= 1;
					alive[i] = false;
				}
				else {
					aliveCount += 1;
					alive[i] = true;
				}
			}
		}
	}
	
	for (int i=0; i<particleCount; i++)
		particleList[i]->alive = alive[i];
}

void RodSimulation::stopSimulation() {
	stop = true;
}

/* Scale positions to get a higher density equilibrium */
void RodSimulation::scaleBoxTo(double newSideLength) {
	double scale = newSideLength / sideLength;
	sideLength = newSideLength;
	halfSideLength = sideLength/2;
	
	if (boundaryRadius != -1) {
		boundaryRadius = halfSideLength;
		boundaryRadiusSq = boundaryRadius*boundaryRadius;
	}
	
	for (int i=0; i<particleCount; i++) {
		Rod *p = particleList[i];
		for (int dim=0; dim<2; dim++) {
			p->pos[dim] *= scale;
		}
	}
}


/* Check whether the specified rods are crossing */
bool RodSimulation::crossedRods(Rod *p, Rod *q) {
	double dir = cross(p->u(),q->u());
	double dx[2] = {q->pos[0] - p->pos[0], q->pos[1] - p->pos[1]};
	
	// Make sure we have the right pair in the boundary conditions
	// (closest pair)
	if (periodicBoundary) {
		if (dx[0] < -halfSideLength)
			dx[0]+=sideLength;
		if (dx[1] < -halfSideLength)
			dx[1]+=sideLength;
		if (dx[0] > halfSideLength)
			dx[0]-=sideLength;
		if (dx[1] > halfSideLength)
			dx[1]-=sideLength;
	}
		
	double u1 = cross(q->u(), dx)/dir;
	double u2 = cross(p->u(), dx)/dir;
	
	return (u1 < 1 && u1 > -1 && u2 < 1 && u2 > -1);
}

/* Check whether the simulation is in a good state */
bool RodSimulation::checkConsistency() {
	for (int i=0; i<particleCount; i++) {
		Rod *p = particleList[i];
		for (int j=i+1; j<particleCount; j++) {
			Rod *q = particleList[j];
			if (crossedRods(p,q)) {
				cerr << "Rods " << p->idNum << " and " << q->idNum << " crossed\n";
				return false; // Uh oh...
			}
		}
	}
	
	return true; // All is well!
}

/* Apply Forces+Torques using Forward Euler */
void RodSimulation::applyForces (bool temp) {
	totalDisp = 0;
	orientAllCorr = 0;
	orientInitCorr = 0;
	
	thetaDisp_max = (0.4/rodLength)/failFactor; // (Max end displacement of 0.2, or 1/10 of a diameter)
	posSqDisp_max = 0.04/(failFactor*failFactor); // (Max particle displacement of 0.2, or 1/10 of a diameter))
	randomXParrStdev = sqrt(2.0 * timeStep / (XParrFrictionCoeff * failFactor));
	randomXPerpStdev = sqrt(2.0 * timeStep / (XPerpFrictionCoeff * failFactor));
	randomTStdev = sqrt(2.0 * timeStep / (TFrictionCoeff * failFactor));
	
	
	for (int i=0; i<particleCount; i++) {
		Rod *p = particleList[i];
		/** Calculate displacement by using the rod's coordinate system **/
		double *parrvec = p->uhat();
		double perpvec[2] = {-parrvec[1], parrvec[0]};
		
		double ldx;
		double dx[2];
		
		ldx = dot(p->force,parrvec)*timeStep/(XParrFrictionCoeff*failFactor);
		if (!temp) {
			ldx /= 2; // We added both forces in here, so divide by 2 to average them
		}
		dx[0] = ldx*parrvec[0];
		dx[1] = ldx*parrvec[1];
				
		ldx = dot(p->force,perpvec)*timeStep/(XPerpFrictionCoeff*failFactor);
		if (!temp) {
			ldx /= 2;
		}
		dx[0] += ldx*perpvec[0];
		dx[1] += ldx*perpvec[1];
		
		
		if (temp) {
			// Calculate the random force in the predictor step, then save it
			// for the corrector step.
			double randomDispMag;
			
			// Parallel component
			randomDispMag = randomGen->RandomNormal()*randomXParrStdev;
			p->randomDisp[0] = randomDispMag*parrvec[0];
			p->randomDisp[1] = randomDispMag*parrvec[1];
			
			// Perpendicular component
			randomDispMag = randomGen->RandomNormal()*randomXPerpStdev;
			p->randomDisp[0] += randomDispMag*perpvec[0];
			p->randomDisp[1] += randomDispMag*perpvec[1];
		}

		dx[0] += p->randomDisp[0];
		dx[1] += p->randomDisp[1];		
		
		double disp = dx[0]*dx[0]+dx[1]*dx[1];
		
		
		double dtheta = p->torque*timeStep/(TFrictionCoeff*failFactor);
		if (!temp) {
			dtheta /= 2;
		}
		
		if (temp) {
			p->randomThetaDisp = randomGen->RandomNormal()*randomTStdev;
		}
		dtheta += p->randomThetaDisp;

		
		/** Cut off the force at the given maximum displacement value **/
		//if (disp > posSqDisp_max) {
			//cout << "Step " << currentTime << ": Displacement too large on particle " << p->idNum << ", " << disp << "\n";
			if (disp > posSqDisp_max*10) {
				failCount += 1;
				//cout << "Too big - returning\n";
				return;
			}
			//continue;
			//disp = sqrt(disp/posSqDisp_max);
			//dx[0] /= disp;
			//dx[1] /= disp;
		//}
		
		//if (dtheta > thetaDisp_max) {
			if (dtheta > thetaDisp_max*10) {
				failCount += 1;
				//cout << "Too big - returning\n";
				return;
			}
			//cout << "Step " << currentTime << ": Torque too large on particle " << p->idNum << ", " << dtheta << "\n";
			//dtheta = thetaDisp_max;
		//}
		//else if (dtheta < -thetaDisp_max) {
			if (dtheta < -thetaDisp_max*10) {
				failCount += 1;
				//cout << "Too big - returning\n";
				return;
			}
			//cout << "Step " << currentTime << ": Torque too large on particle " << p->idNum << ", " << dtheta << "\n";
			//dtheta = -thetaDisp_max;
		//}
		
				
		/* If this is temporary (for the predictor corrector), make a note of 
		 * the changes so they can be undone.
		 */
		if (temp) {
			p->tempPosDisp[0] = dx[0];
			p->tempPosDisp[1] = dx[1];
			p->tempThetaDisp = dtheta;
		}
		else {
			p->revertTemp();

			if (output)
				totalDisp += p->positionSqDisp();
		}

		p->theta += dtheta;		
		
		for (int dim=0; dim<2; dim++) {			
			p->pos[dim] += dx[dim];
			p->absPos[dim] += dx[dim];
			
			if (periodicBoundary) {
				if (p->pos[dim] < 0)
					p->pos[dim] += sideLength;
				else if (p->pos[dim] >= sideLength)
					p->pos[dim] -= sideLength;
			}
		}
		
		
		if (p->theta < 0)
			p->theta += 2*M_PI;
		else if (p->theta > 2*M_PI)
			p->theta -= 2*M_PI;
		

		if (temp) {
			p->calculateTempU();
		}
		else {
			p->calculateU();
			
			if (output) {
				orientInitCorr += p->orientChange();

				for (int j=i+1; j<particleCount; j++) {
					Rod *q = particleList[j];
					
					orientAllCorr = cos(2*(p->theta - q->theta));
				}
			}
		}
	}
	
	if (output) {
		totalDisp /= particleCount;
		orientInitCorr /= particleCount;
		if (particleCount != 1) { // Can't do this if there's only one rod...not sure why you'd have that
			// There were n(n-1)/2 comparisons  (1 + 2 + ... + n-1)
			orientAllCorr *= 2;
			orientAllCorr /= (particleCount)*(particleCount-1);
		}
	}
}

/* Calculate Forces+Torques */
void RodSimulation::calculateForces (bool reset) {
	/** Zero the forces, torques, energy **/
	for (int i=0; i<particleCount; i++) {
		if (reset) {
			particleList[i]->resetForces();
		}
		particleList[i]->currentForce[0] = 0;
		particleList[i]->currentForce[1] = 0;

		particleList[i]->frictionForce[0] = 0;
		particleList[i]->frictionForce[1] = 0;
	}
	
	/** Calculate the distance to each rod. P is rod 1, Q is rod 2.  **/
	bool needUpdate = false;
	for (int i=0; i<particleCount; i++) {
		if (particleList[i]->maxDisplacement()>neighborListBufferSq4) {
			needUpdate = true;
			break;
		}
	}
		
	
	if (needUpdate) {
		updateNeighbors();
	}

	InteractionPair *ip;
	
	for (int i=0; i<neighborList.size(); i++) {
		// Calculate initial forces - mark which rods interact
		ip = &(neighborList[i]);
		calculateForcesForRods(ip);
	}
	
	for (int i=0; i<particleCount; i++) {
        // Add in the self propulsion force, if the force is not too big already.
		Rod *r = particleList[i];
		//motorGrid->addForceToRod(r);
		if (r->alive) {
			r->currentForce[0] += selfPushForce*r->uhat()[0];
			r->currentForce[1] += selfPushForce*r->uhat()[1];
		}
		
		if (boundaryRadius != -1) {
			double *u = r->u();
			double endA[2] = {r->pos[0]+u[0], r->pos[1]+u[1]};
			double endB[2] = {r->pos[0]-u[0], r->pos[1]-u[1]};
			double radSqA = endA[0]*endA[0]+endA[1]*endA[1];
			double radSqB = endB[0]*endB[0]+endB[1]*endB[1];
			if ( radSqA > boundaryRadiusSq ) {
				double rad = sqrt(radSqA);
				double dist = cutoff-(rad - boundaryRadius);
				double dx[2] = {-(endA[0])*dist/rad,-(endA[1])*dist/rad};
			
				double df = 48 * (pow(dist, -14) - 0.5 *pow(dist, -8));
				double dfx[2];

				for (int dim=0; dim<2; dim++) {
					dfx[dim] = dx[dim] * df;
					r->currentForce[dim] += dfx[dim];
				}

				r->torque -= cross(dfx, u);
			}
			if ( radSqB > boundaryRadiusSq ) {
				double rad = sqrt(radSqB);
				double dist = cutoff-(rad - boundaryRadius);
				double dx[2] = {-(endB[0])*dist/rad,-(endB[1])*dist/rad};
			
				double df = 48 * (pow(dist, -14) - 0.5 *pow(dist, -8));
				double dfx[2];

				for (int dim=0; dim<2; dim++) {
					dfx[dim] = dx[dim] * df;
					r->currentForce[dim] += dfx[dim];
				}
				
				r->torque += cross(dfx, u);	
			}
		}
	}
	
	if (alignmentParameter != 0) {
		cout << "Adding friction" << endl;
		for (int i=0; i<neighborList.size(); i++) {
			// Calculate "friction" forces
			// This is a bit weird, because it is velocity dependent.  I think it still works with the
			// "predictor-corrector" though.
			ip = &(neighborList[i]);
			Rod *rod1 = ip->rod1;
			Rod *rod2 = ip->rod2;
			if (ip->interactedThisTime == true) {
				double parallelness = fabs(dot(rod1->uhat(),rod2->uhat()));
				double forceDifference[2] = {rod1->currentForce[0]-rod2->currentForce[0], rod1->currentForce[1]-rod2->currentForce[1]};
				double interactionStrength = -sqrt(ip->currentDistSq/frictionCutoffSq) + 1.0;
				//double forceDifferenceSize = dot(forceDifference, forceDifference);
				double force1 = -dot(forceDifference, rod1->uhat())*parallelness*alignmentParameter*interactionStrength;
				double force2 =  dot(forceDifference, rod2->uhat())*parallelness*alignmentParameter*interactionStrength;

				// f1 = -(f1-f2).(u1) * u1
				// f2 = -(f2-f1).(u2) * u2
				
				rod1->frictionForce[0] = force1*rod1->uhat()[0];
				rod1->frictionForce[1] = force1*rod1->uhat()[1];
				
				rod2->frictionForce[0] = force2*rod2->uhat()[0];
				rod2->frictionForce[1] = force2*rod2->uhat()[1];
			}
			else {
				rod1->frictionForce[0] = 0;
				rod1->frictionForce[1] = 0;
				
				rod2->frictionForce[0] = 0;
				rod2->frictionForce[1] = 0;
			}
		}
	}
	
	for (int i=0; i<particleCount; i++) {
		Rod *r = particleList[i];
		
		// Add the current interaction and "friction" forces into the total force
		r->force[0] += r->frictionForce[0];
		r->force[1] += r->frictionForce[1];

		r->force[0] += r->currentForce[0];
		r->force[1] += r->currentForce[1];
	}
	
}

void RodSimulation::calculateForcesForRods(InteractionPair *ip) {
	if (ip->isZero) {
		return;
		ip->interactedThisTime = false;
	}

	Rod *p = ip->rod1;
	Rod *q = ip->rod2;

	double *endvec1 = p->u();
	double *endvec2 = q->u();

	double dx[2];
	double torqueR1[2], torqueR2[2];
	double distSq = 0;
	
	/*** Calculate the shortest distance between two rods ***/
	/* Try the ends. A is the + end, B is the - end. */
	
	double qloc[2] = {q->pos[0], q->pos[1]};
	
	if (periodicBoundary) {
		if (qloc[0]-p->pos[0] < -halfSideLength)
			qloc[0]+=sideLength;
		else if (qloc[0]-p->pos[0] > halfSideLength)
			qloc[0]-=sideLength;
	
		if (qloc[1]-p->pos[1] < -halfSideLength)
			qloc[1]+=sideLength;
		else if (qloc[1]-p->pos[1] > halfSideLength)
			qloc[1]-=sideLength;
	}
	
	double end1a[2] = {p->pos[0] + endvec1[0], p->pos[1] + endvec1[1]};
	double end1b[2] = {p->pos[0] - endvec1[0], p->pos[1] - endvec1[1]};
	double end2a[2] = {qloc[0] + endvec2[0], qloc[1] + endvec2[1]};
	double end2b[2] = {qloc[0] - endvec2[0], qloc[1] - endvec2[1]};


	double end_ab[2], end_ba[2], end_aa[2], end_bb[2];
	double x_12a[2], x_12b[2], x_21a[2], x_21b[2];
	
	double dend_ab=DINFINITY, dend_ba=DINFINITY, dend_aa=DINFINITY, dend_bb=DINFINITY;
	double dmid_12a=DINFINITY, dmid_12b=DINFINITY, dmid_21a=DINFINITY, dmid_21b=DINFINITY;

	
	//dend_xy = distance between end x of rod 1 and end y of rod 2
	if (ip->inSubList(endA_t, endA_t)) {
		end_aa[0] = end1a[0] - end2a[0]; end_aa[1] = end1a[1] - end2a[1];
		dend_aa = end_aa[0]*end_aa[0] + end_aa[1]*end_aa[1];
	}
	if (ip->inSubList(endB_t, endB_t)) {
		end_bb[0] = end1b[0] - end2b[0]; end_bb[1] = end1b[1] - end2b[1];
		dend_bb = end_bb[0]*end_bb[0] + end_bb[1]*end_bb[1];
	}
	if (ip->inSubList(endA_t, endB_t)) {
		end_ab[0] = end1a[0] - end2b[0]; end_ab[1] = end1a[1] - end2b[1];
		dend_ab = end_ab[0]*end_ab[0] + end_ab[1]*end_ab[1];
	}
	if (ip->inSubList(endB_t, endA_t)) {
		end_ba[0] = end1b[0] - end2a[0]; end_ba[1] = end1b[1] - end2a[1];
		dend_ba = end_ba[0]*end_ba[0] + end_ba[1]*end_ba[1];
	}
	
	//dmid_1a2 = distance between rod 1 and end a of rod 2
	double u; // paramaterize the rod as (Xcm + u*endvec) u = -1...1

	double cmDisp[2];
	
	if (ip->inSubList(mid_t, endA_t)) {
	cmDisp[0] = end2a[0]-p->pos[0]; cmDisp[1] = end2a[1]-p->pos[1];
		u = (cmDisp[0]*endvec1[0] + cmDisp[1]*endvec1[1])/(rodLength2*rodLength2);
		if (u < -1 || u > 1)
			dmid_12a = DINFINITY;
		else {
			x_12a[0] = p->pos[0] + u*endvec1[0];
			x_12a[1] = p->pos[1] + u*endvec1[1];
			dmid_12a = (end2a[0] - x_12a[0])*(end2a[0] - x_12a[0]) + (end2a[1] - x_12a[1])*(end2a[1] - x_12a[1]);
		}
	}
	
	if (ip->inSubList(mid_t, endB_t)) {
		cmDisp[0] = end2b[0]-p->pos[0]; cmDisp[1] = end2b[1]-p->pos[1];
		u = (cmDisp[0]*endvec1[0] + cmDisp[1]*endvec1[1])/(rodLength2*rodLength2);
		if (u < -1 || u > 1)
			dmid_12b = DINFINITY;
		else {
			x_12b[0] = p->pos[0] + u*endvec1[0];
			x_12b[1] = p->pos[1] + u*endvec1[1];
			dmid_12b = (end2b[0] - x_12b[0])*(end2b[0] - x_12b[0]) + (end2b[1] - x_12b[1])*(end2b[1] - x_12b[1]);
		}
	}
	
	if (ip->inSubList(endA_t, mid_t)) {
		cmDisp[0] = end1a[0]-qloc[0]; cmDisp[1] = end1a[1]-qloc[1];
		u = (cmDisp[0]*endvec2[0] + cmDisp[1]*endvec2[1])/(rodLength2*rodLength2);
		if (u < -1 || u > 1)
			dmid_21a = DINFINITY;
		else {
			x_21a[0] = qloc[0] + u*endvec2[0];
			x_21a[1] = qloc[1] + u*endvec2[1];
			dmid_21a = (end1a[0] - x_21a[0])*(end1a[0] - x_21a[0]) + (end1a[1] - x_21a[1])*(end1a[1] - x_21a[1]);
		}
	}

	if (ip->inSubList(endB_t, mid_t)) {
		cmDisp[0] = end1b[0]-qloc[0]; cmDisp[1] = end1b[1]-qloc[1];
		u = (cmDisp[0]*endvec2[0] + cmDisp[1]*endvec2[1])/(rodLength2*rodLength2);
		if (u < -1 || u > 1)
			dmid_21b = DINFINITY;
		else {
			x_21b[0] = qloc[0] + u*endvec2[0];
			x_21b[1] = qloc[1] + u*endvec2[1];
			dmid_21b = (end1b[0] - x_21b[0])*(end1b[0] - x_21b[0]) + (end1b[1] - x_21b[1])*(end1b[1] - x_21b[1]);
		}
	}
	
	
	if (!isnan(dend_aa) && dend_aa < dend_bb && dend_aa < dend_ab && dend_aa < dend_ba
		&& dend_aa < dmid_12a && dend_aa < dmid_12b && dend_aa < dmid_21a && dend_aa < dmid_21b) {
		torqueR1[0] = end1a[0] - p->pos[0]; torqueR1[1] = end1a[1] - p->pos[1];
		torqueR2[0] = end2a[0] - qloc[0]; torqueR2[1] = end2a[1] - qloc[1];
		dx[0] = end_aa[0]; dx[1] = end_aa[1];
	}
	else if (!isnan(dend_bb) && dend_bb < dend_ab && dend_bb < dend_ba
			 && dend_bb < dmid_12a && dend_bb < dmid_12b && dend_bb < dmid_21a && dend_bb < dmid_21b) {
		torqueR1[0] = end1b[0] - p->pos[0]; torqueR1[1] = end1b[1] - p->pos[1];
		torqueR2[0] = end2b[0] - qloc[0]; torqueR2[1] = end2b[1] - qloc[1];
		dx[0] = end_bb[0]; dx[1] = end_bb[1];
	}
	else if (!isnan(dend_ab) && dend_ab < dend_ba
			 && dend_ab < dmid_12a && dend_ab < dmid_12b && dend_ab < dmid_21a && dend_ab < dmid_21b) {
		torqueR1[0] = end1a[0] - p->pos[0]; torqueR1[1] = end1a[1] - p->pos[1];
		torqueR2[0] = end2b[0] - qloc[0]; torqueR2[1] = end2b[1] - qloc[1];
		dx[0] = end_ab[0]; dx[1] = end_ab[1];
	}
	else if (!isnan(dend_ba) && dend_ba < dmid_12a && dend_ba < dmid_12b && dend_ba < dmid_21a && dend_ba < dmid_21b) {
		torqueR1[0] = end1b[0] - p->pos[0]; torqueR1[1] = end1b[1] - p->pos[1];
		torqueR2[0] = end2a[0] - qloc[0]; torqueR2[1] = end2a[1] - qloc[1];
		dx[0] = end_ba[0]; dx[1] = end_ba[1];
	}
	else if (!isnan(dmid_12a) && dmid_12a < dmid_12b && dmid_12a < dmid_21a && dmid_12a < dmid_21b) {
		torqueR1[0] = x_12a[0] - p->pos[0]; torqueR1[1] = x_12a[1] - p->pos[1];
		torqueR2[0] = end2a[0] - qloc[0]; torqueR2[1] = end2a[1] - qloc[1];
		dx[0] = x_12a[0] - end2a[0]; dx[1] = x_12a[1] - end2a[1];
	}
	else if (!isnan(dmid_12b) && dmid_12b < dmid_21a && dmid_12b < dmid_21b) {
		torqueR1[0] = x_12b[0] - p->pos[0]; torqueR1[1] = x_12b[1] - p->pos[1];
		torqueR2[0] = end2b[0] - qloc[0]; torqueR2[1] = end2b[1] - qloc[1];
		dx[0] = x_12b[0] - end2b[0]; dx[1] = x_12b[1] - end2b[1];
	}
	else if (!isnan(dmid_21a) && dmid_21a < dmid_21b) {
		torqueR1[0] = end1a[0] - p->pos[0]; torqueR1[1] = end1a[1] - p->pos[1];
		torqueR2[0] = x_21a[0] - qloc[0]; torqueR2[1] = x_21a[1] - qloc[1];
		dx[0] = end1a[0] - x_21a[0]; dx[1] = end1a[1] - x_21a[1];
	}
	else if (!isnan(dmid_21b)){
		torqueR1[0] = end1b[0] - p->pos[0]; torqueR1[1] = end1b[1] - p->pos[1];
		torqueR2[0] = x_21b[0] - qloc[0]; torqueR2[1] = x_21b[1] - qloc[1];
		dx[0] = end1b[0] - x_21b[0]; dx[1] = end1b[1] - x_21b[1];
	} 
	else {
		// This shouldn't happen.
		torqueR1[0] = torqueR1[1] = torqueR2[0] = torqueR2[1] = 0;
		dx[0] = dx[1] = DINFINITY;
	}
	
	/********************************************************/
	
	for (int dim=0; dim<2; dim++) {					
		distSq += dx[dim]*dx[dim];
	}
	
	/**** Use this to see what the forces look like ****
	if (debugOn) {
		double debugInfo[4] = {0,0,0,0};

		cout << "Step " << currentTime << ": " << distSq << "\n";
		
		debugInfo[0] = p->pos[0]+torqueR1[0];
		debugInfo[1] = p->pos[1]+torqueR1[1];
		debugInfo[2] = debugInfo[0]-dx[0];
		debugInfo[3] = debugInfo[1]-dx[1];
		
		debugOut.write((char*)debugInfo, sizeof(double)*4);
	}
	***************************************************/
	
	/** Cutoff distance is 2^(1/6), the end of repulsion **/
	if (distSq > cutoffSq) {
		ip->interactedThisTime = (distSq < frictionCutoffSq);
		return;
	}
	else {
		ip->interactedThisTime = true;
		ip->currentDistSq = distSq;
	}

	
	if (distSq < 0.6) {
//		cerr << "Step " << currentTime << ": Small displacement particle " << p->idNum << "->" << q->idNum << ": (" << dx[0] << ", " << dx[1] << ") =" << distSq <<"\n";
	}
	
	
	/** Calculate the forces, using LJ
	 ** F_x = 48 (r^-14 - r^-8 /2)r_x
	 **/
	//if (potentialEnergy != NULL)
	//	*potentialEnergy += pow(distSq,-6.0) - pow(distSq,-3.0) + 0.25;
	
	double df = 48 * (pow(distSq, -7) - 0.5 *pow(distSq, -4));
	double dfx[2];
	
	for (int dim=0; dim<2; dim++) {
		dfx[dim] = dx[dim] * df;
		p->currentForce[dim] += dfx[dim];
		q->currentForce[dim] -= dfx[dim];
	}
	
	p->torque -= cross(dfx, torqueR1);
	q->torque += cross(dfx, torqueR2);
}


void RodSimulation::updateNeighbors() {
	neighborList.clear();
	
	for (int i=0; i<particleCount; i++) {
		particleList[i]->neighbors.clear();
		particleList[i]->clumpId=-1;
	}
	
	Rod *p, *q;
	for (int i=0; i<particleCount; i++) {
		p = particleList[i];
		p->storePosition();
		for (int j=i+1; j<particleCount; j++) {		
			q = particleList[j];
			
			InteractionPair ip = InteractionPair(p, q);
			
			double *endvec1 = p->u();
			double *endvec2 = q->u();
			
			/*** Calculate the shortest distance between two rods ***/
			/* Try the ends. A is the + end, B is the - end. */

			double qloc[2] = {q->pos[0], q->pos[1]};

			if (periodicBoundary) {
				if (qloc[0]-p->pos[0] < -halfSideLength)
					qloc[0]+=sideLength;
				if (qloc[1]-p->pos[1] < -halfSideLength)
					qloc[1]+=sideLength;
				if (qloc[0]-p->pos[0] > halfSideLength)
					qloc[0]-=sideLength;
				if (qloc[1]-p->pos[1] > halfSideLength)
					qloc[1]-=sideLength;
			}
			
			double cmDistSq = (qloc[0]-p->pos[0])*(qloc[0]-p->pos[0]) + (qloc[1]-p->pos[1])*(qloc[1]-p->pos[1]);
			if (cmDistSq > neighborListQuickCutoffSq)
				continue;

			double end1a[2] = {p->pos[0] + endvec1[0], p->pos[1] + endvec1[1]};
			double end1b[2] = {p->pos[0] - endvec1[0], p->pos[1] - endvec1[1]};
			double end2a[2] = {qloc[0] + endvec2[0], qloc[1] + endvec2[1]};
			double end2b[2] = {qloc[0] - endvec2[0], qloc[1] - endvec2[1]};

			double end_ab[2], end_ba[2], end_aa[2], end_bb[2];
			double dend_ab, dend_ba, dend_aa, dend_bb;

			//dend_xy = distance between end x of rod 1 and end y of rod 2
			end_aa[0] = end1a[0] - end2a[0]; end_aa[1] = end1a[1] - end2a[1];
			end_bb[0] = end1b[0] - end2b[0]; end_bb[1] = end1b[1] - end2b[1];
			end_ab[0] = end1a[0] - end2b[0]; end_ab[1] = end1a[1] - end2b[1];
			end_ba[0] = end1b[0] - end2a[0]; end_ba[1] = end1b[1] - end2a[1];

			dend_aa = (end_aa[0])*(end_aa[0]) + (end_aa[1])*(end_aa[1]);
			dend_bb = (end_bb[0])*(end_bb[0]) + (end_bb[1])*(end_bb[1]);
			dend_ab = (end_ab[0])*(end_ab[0]) + (end_ab[1])*(end_ab[1]);
			dend_ba = (end_ba[0])*(end_ba[0]) + (end_ba[1])*(end_ba[1]);


			double x_12a[2]={0,0}, x_12b[2]={0,0}, x_21a[2]={0,0}, x_21b[2]={0,0};
			double dmid_12a, dmid_12b, dmid_21a, dmid_21b;
			//dmid_21a = distance between rod 1 and end a of rod 2
			double u; // paramaterize the rod as (Xcm + u*endvec) u = -1...1
			double cmDisp[2];

			cmDisp[0] = end2a[0]-p->pos[0]; cmDisp[1] = end2a[1]-p->pos[1];
			u = (cmDisp[0]*endvec1[0] + cmDisp[1]*endvec1[1])/(rodLength2*rodLength2);
			if (u < -1 || u > 1)
				dmid_12a = DINFINITY;
			else {
				x_12a[0] = p->pos[0] + u*endvec1[0];
				x_12a[1] = p->pos[1] + u*endvec1[1];
				dmid_12a = ((end2a[0] - x_12a[0])*(end2a[0] - x_12a[0]) + (end2a[1] - x_12a[1])*(end2a[1] - x_12a[1]));
			}

			cmDisp[0] = end2b[0]-p->pos[0]; cmDisp[1] = end2b[1]-p->pos[1];
			u = (cmDisp[0]*endvec1[0] + cmDisp[1]*endvec1[1])/(rodLength2*rodLength2);
			if (u < -1 || u > 1)
				dmid_12b = DINFINITY;
			else {
				x_12b[0] = p->pos[0] + u*endvec1[0];
				x_12b[1] = p->pos[1] + u*endvec1[1];
				dmid_12b = ((end2b[0] - x_12b[0])*(end2b[0] - x_12b[0]) + (end2b[1] - x_12b[1])*(end2b[1] - x_12b[1]));
			}

			cmDisp[0] = end1a[0]-qloc[0]; cmDisp[1] = end1a[1]-qloc[1];
			u = (cmDisp[0]*endvec2[0] + cmDisp[1]*endvec2[1])/(rodLength2*rodLength2);
			if (u < -1 || u > 1)
				dmid_21a = DINFINITY;
			else {
				x_21a[0] = qloc[0] + u*endvec2[0];
				x_21a[1] = qloc[1] + u*endvec2[1];
				dmid_21a = ((end1a[0] - x_21a[0])*(end1a[0] - x_21a[0]) + (end1a[1] - x_21a[1])*(end1a[1] - x_21a[1]));
			}

			cmDisp[0] = end1b[0]-qloc[0]; cmDisp[1] = end1b[1]-qloc[1];
			u = (cmDisp[0]*endvec2[0] + cmDisp[1]*endvec2[1])/(rodLength2*rodLength2);
			if (u < -1 || u > 1)
				dmid_21b = DINFINITY;
			else {
				x_21b[0] = qloc[0] + u*endvec2[0];
				x_21b[1] = qloc[1] + u*endvec2[1];
				dmid_21b = ((end1b[0] - x_21b[0])*(end1b[0] - x_21b[0]) + (end1b[1] - x_21b[1])*(end1b[1] - x_21b[1]));
			}


			// There should be some way of using the some distances to avoid calculating later ones.

			// Now let's see which pairs of rod parts are within the cutoff distance
			if (dend_aa < neighborListCutoffSq)
				ip.setSubList(endA_t, endA_t);
			if (dend_ab < neighborListCutoffSq)
				ip.setSubList(endA_t, endB_t);
			if (dend_ba < neighborListCutoffSq)
				ip.setSubList(endB_t, endA_t);
			if (dend_bb < neighborListCutoffSq)
				ip.setSubList(endB_t, endB_t);
				
			if (dmid_21a < neighborListCutoffSq)
				ip.setSubList(endA_t, mid_t);
			if (dmid_21b < neighborListCutoffSq)
				ip.setSubList(endB_t, mid_t);
			if (dmid_12a < neighborListCutoffSq)
				ip.setSubList(mid_t, endA_t);
			if (dmid_12b < neighborListCutoffSq)
				ip.setSubList(mid_t, endB_t);
			/***************************************************/
			//			cout << "cutoff is " << neighborListCutoffSq << endl;
			if (!ip.isZero) {
				neighborList.push_back(ip);
				p->neighbors.push_back(q);
				q->neighbors.push_back(p);
				
				/*cout << "1: " << p->idNum << ", 2: " << q->idNum << endl;
				cout << "\tdend_aa " << dend_aa << " " << (dend_aa<neighborListCutoffSq) << endl
					<< "\tdend_ab " << dend_ab << " " << (dend_ab<neighborListCutoffSq)<< endl
					<< "\tdend_ba " << dend_ba << " " << (dend_ba<neighborListCutoffSq) << endl
					<< "\tdend_bb " << dend_bb << " " << (dend_bb<neighborListCutoffSq) << endl
					<< "\tdmid_12a " << dmid_12a << " " << (dmid_12a<neighborListCutoffSq) << endl
					<< "\tdmid_21a " << dmid_21a << " " << (dmid_21a<neighborListCutoffSq) << endl
					<< "\tdmid_12b " << dmid_12b << " " << (dmid_12b<neighborListCutoffSq) << endl
					<< "\tdmid_21b " << dmid_21b << " " << (dmid_21b<neighborListCutoffSq) << endl;*/
			}
		}		
	}
}

// Calculate clumps by seeing which rods are interacting.
void RodSimulation::updateClumps() {
	int newId = 0;
	for (int i=0; i<particleCount; i++) {
		if (drillDownClump(newId, particleList[i]))
			newId++;
	}
}

bool RodSimulation::drillDownClump(int newId, Rod *r) {
	if (r->clumpId != -1)
		return false;
	
	r->clumpId = newId;
	
	for (int i=0; i < r->neighbors.size(); i++) {
		drillDownClump(newId, r->neighbors[i]);
	}
	
	return true;
}

double RodSimulation::dot(double v1[2], double v2[2]) {
	return (v1[0]*v2[0] + v1[1]*v2[1]);
}

double RodSimulation::cross(double v1[2], double v2[2]) {
	return (v1[0]*v2[1]-v1[1]*v2[0]);
}

double RodSimulation::wrap(double x) {
	if (periodicBoundary) {
		double newx = x;
		if (newx < -halfSideLength)
			newx += sideLength;
		else if (newx > halfSideLength)
			newx -= sideLength;
		return newx;
	}
	else {
		return x;
	}
}


/*
 * This is really annoying.  I basically guessed until I found the combination of
 * flags that worked.  The effect of the first open/close is to create the file if
 * it doesn't exist, without accidentally truncating it.
 */
void RodSimulation::openOut(ofstream &file, const char *name) {
	file.open(name, ios::out | ios::app);
	file.close();
	
	file.open(name, ios::out | ios::ate | ios::binary | ios::in);
}
