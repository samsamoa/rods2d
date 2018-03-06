/*
 *  Rod.cpp
 *  bd
 *
 *  Created by Sam McCandlish on 12/22/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "Rod.h"

# ifdef __INTEL_COMPILER 
# include <mathimf.h> 
# else 
# include <math.h> 
# endif 

Rod::Rod(istream &stateFile, double _halfLength, int _idNum, int version) {
	halfLength = _halfLength;
	idNum = _idNum;
	
	stateFile.read((char*)pos, sizeof(double)*2);
	stateFile.read((char*)&theta, sizeof(double));
	
	calculateU();
	
	if (version == 1) {
		initPos[0] = absPos[0] = pos[0];
		initPos[1] = absPos[1] = pos[1];
		
		initTheta = theta;
		initUhat[0] = _uhat[0];
		initUhat[1] = _uhat[1];
		
		alive = true;
	}
	else if (version == 2) {
		stateFile.read((char*)absPos, sizeof(double)*2);
		stateFile.read((char*)initPos, sizeof(double)*2);
		stateFile.read((char*)&initTheta, sizeof(double));
		stateFile.read((char*)&alive, sizeof(bool));
		
		initUhat[0] = cos(initTheta);
		initUhat[1] = sin(initTheta);
	}
}

void Rod::reload(istream &stateFile) {
	stateFile.read((char*)pos, sizeof(double)*2);
	stateFile.read((char*)&theta, sizeof(double));
	absPos[0]=pos[0];
	absPos[0]=pos[0];
	calculateU();
}

Rod::Rod(double _x, double _y, double _theta, double _halfLength, int _idNum) {
	pos[0] = initPos[0] = absPos[0] = _x;
	pos[1] = initPos[1] = absPos[1] = _y;
	theta = initTheta = _theta;
	
	halfLength = _halfLength;
	idNum = _idNum;
	
	calculateU();
	initUhat[0] = _uhat[0];
	initUhat[1] = _uhat[1];
	
	alive = true;
}

void Rod::resetForces() {
	torque = 0;
	force[0] = 0;
	force[1] = 0;
	frictionForce[0] = 0;
	frictionForce[1] = 0;
}

void Rod::calculateU() {
	_uhat[0] = cos(theta);
	_uhat[1] = sin(theta);
	
	_u[0] = halfLength*_uhat[0];
	_u[1] = halfLength*_uhat[1];
}

void Rod::calculateTempU() {
	_uhatprev[0] = _uhat[0];
	_uhatprev[1] = _uhat[1];
	_uprev[0] = _u[0];
	_uprev[1] = _u[1];
	
	calculateU();
}

void Rod::revertTemp() {
	pos[0] -= tempPosDisp[0];
	pos[1] -= tempPosDisp[1];
	absPos[0] -= tempPosDisp[0];
	absPos[1] -= tempPosDisp[1];
	theta -= tempThetaDisp;
		
	_uhat[0] = _uhatprev[0];
	_uhat[1] = _uhatprev[1];
	_u[0] = _uprev[0];
	_u[1] = _uprev[1];
}

double *Rod::u() {
	return _u;
}

double *Rod::uhat() {
	return _uhat;
}

void Rod::storePosition() {
	neighborListPos[0] = absPos[0];
	neighborListPos[1] = absPos[1];
	neighborListU[0] = _u[0];
	neighborListU[1] = _u[1];
}

double Rod::maxDisplacement() {
	double endADisp[2] =  { absPos[0]+_u[0]-neighborListPos[0]-neighborListU[0],
				absPos[1]+_u[1]-neighborListPos[1]-neighborListU[1]};
	double endADispSq = endADisp[0]*endADisp[0] + endADisp[1]*endADisp[1];
	double endBDisp[2] =  {	absPos[0]-_u[0]-neighborListPos[0]+neighborListU[0],
				absPos[1]-_u[1]-neighborListPos[1]+neighborListU[1]};
	double endBDispSq = endBDisp[0]*endBDisp[0] + endBDisp[1]*endBDisp[1];
	
	if (endADispSq>endBDispSq)
		return endADispSq;
	else
		return endBDispSq;
}

void Rod::saveState(ostream &stateFile, int version) {
	stateFile.write((char*)pos, sizeof(double)*2);
	stateFile.write((char*)&theta, sizeof(double));
	
	if (version == 1) {
		// Do nothing
	}
	else if (version == 2) {
		stateFile.write((char*)absPos, sizeof(double)*2);
		stateFile.write((char*)initPos, sizeof(double)*2);
		stateFile.write((char*)&initTheta, sizeof(double));
		stateFile.write((char*)&alive, sizeof(bool));
	}
}

double Rod::orientChange() {
	return (_uhat[0]*initUhat[0] + _uhat[1]*initUhat[1]);
}

double Rod::positionSqDisp() {
	return pow(absPos[0] - initPos[0],2) + pow(absPos[1] - initPos[1],2);
}
