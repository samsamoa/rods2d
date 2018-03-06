#pragma once
/*
 *  MotorGrid.h
 *  bd
 *
 *  Created by Sam McCandlish on 12/22/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "Rod.h"

class RodSimulation;

using namespace std;

class MotorGrid {
public:
	MotorGrid(RodSimulation *sim, double motorCount, double motorForce, double captureRadius);
	
	void addForceToRod(Rod *p);

private:
	int sideCount;
	double motorSpacing;
	double motorForce;
	double captureRadiusSq;
	double quickCutoff;
	double quickCutoffSq;
	RodSimulation *sim;
	
	double distanceToMotor(Rod *p, double motorPos[2]);
};
