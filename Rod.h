#pragma once
/*
 *  Rod.h
 *  bd
 *
 *  Created by Sam McCandlish on 12/22/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

class Rod {
public:
	/** State Values **/
	double pos[2];			// X,Y coordinates
	double theta;			// Rotational position
		
	double tempPosDisp[2];	// Displacement on pos for a temporary step
	double tempThetaDisp;	// Displacement on theta for a temporary step
	
	double absPos[2];		// pos, without periodic boundary conditions applied
	
	double initPos[2];		// Initial position, for total displacement calculation
	double initTheta;		// Initial theta value, for total rotation calculation
	double initUhat[2];
	
	/** Forces **/
	double torque;			// Current total torque on this rod
	double force[2];		// Current total force on this rod
	
	double currentForce[2];	// The force in the current predictor-corrector state
							// (This is necessary to calculate the friction correctly)
	double frictionForce[2];// Friction force - calculated separately
	
	double randomDisp[2];  // Noise force - in the simulation frame
	double randomThetaDisp;
	
	/** Constant Values **/
	double halfLength;		// Half the rod length
	int idNum;				// A number to identify this rod, not used in sim
	bool alive;				// Does this rod get propelled?
	
	vector<Rod *> neighbors;
	vector<Rod *> clusterNeighbors;
	int clumpId;
	

	Rod(istream &stateFile, double _halfLength, int _idNum, int version);
		// Read the position from the given input stream
	
	Rod(double _x, double _y, double _theta, double _halfLength, int _idNum);
		// New rod
	
	
	void resetForces();		// Set forces to 0
	double* u();			// Vector from center to + end of rod
	double* uhat();			// Unit vector from center to + end of rod
	void calculateU();		// Calculate the u and uhat vectors
	void calculateTempU();	// Do as above, but store the previous value (for speed)
	void revertTemp();		// Revert the position and u,uhat vectors 
	void storePosition();	// Store the position for calculating max displacement for neighborlist
	double maxDisplacement();	// Which end has moved the farthest, and what is the distsq?
		
	double orientChange();	// Dot between uHat and initUHat
	double positionSqDisp();// Position displacement, squared
	
	void saveState(ostream &stateFile, int version);
		// Save this rod to the given stream
	void reload(istream &stateFile);

	double _u[2];
	double _uhat[2];
	
	double _uprev[2];
	double _uhatprev[2];
	
	double neighborListPos[2];
	double neighborListU[2];

};
