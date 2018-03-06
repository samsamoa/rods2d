#pragma once
/*
 *  RodSimulation.h
 *  bd
 *
 *  Created by Sam McCandlish on 12/22/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "Rand.h"
#include "Rod.h"
#include "InteractionPair.h"
#include "MotorGrid.h"

#include <vector>
#include <fstream>

#include <limits>
#define DINFINITY std::numeric_limits<double>::infinity()


class RodSimulation {
public:
	
	RodSimulation(int particleCount, double rodLength, double startDensity, double circularBoundary, // Initial parameters
				int randomSeed, double endDensity, double endCircularBoundary, double timeStep, int stopTime, int outputSkip, double selfPushForce, double motorCutoff, double friction, double aliveFraction); // Simulation parameters

		// Use these values for a new simulation
	

	RodSimulation(int randomSeed, double endDensity, double endCircularBoundary, double timeStep, int stopTime, int outputSkip, double selfPushForce, double motorCutoff, double friction, double aliveFraction,	// Simulation parameters
									istream &savedState, bool resetClock, int multiplySideLength);	// Resume parameters

		// Use these values for a saved simulation
	
	~RodSimulation();
	
	int runSimulation();
	void stopSimulation();
	
private:
	/*
	 * Note: everything is in non-dimensionalized units.
	 */
	
	/** User Set Constant Values **/
	int particleCount;
	double rodLength;
	double timeStep;
	int stopTime;
	int outputSkip;			// Only output positions every x frames
	
	bool resizeBox;			// Are we shrinking down to the specified density?
	
	bool periodicBoundary;		// Do we apply periodic boundary conditions?
	double boundaryRadius;		// The radius of the circular boundary, if we are using it
					// instead of periodic boundary conditions (0 otherwise)
	double boundaryRadiusSq;

//	bool motorGrid;			// Are we using a motor grid?
//	double motorDensity;		// If we are using a motor grid, the density of those motors
	double selfPushForce;		// If we are not using a motor grid, the constant force
					// added in the + direction of each rod
	
	

	
	
	/** Calculated Constant Values **/
public:
	double rodLength2;		// Half rod length
	double cutoff;			// Actual cutoff for force (2^(1/6))
	double cutoffSq;
	double quickCutoffSq;		// Maximum distance for interacting rods
	double initSideLength;		// The starting side length (for a shrinking box)
	double targetSideLength;
	double sideLength;		// Side length of the square box
	double halfSideLength;		// Half the side length
	
	double thetaDisp_max;		// Cutoff for rotation in one timestep
	double posSqDisp_max;		// Cutoff for displacement in one timestep
	double motorForceCutoffSq;	// Don't apply the motor force if the other forces are more than this

	double XParrFrictionCoeff;	// Friction coefficient, moving parallel to rod
	double randomXParrStdev;	// Stdev of the random parallel displacement
	double XPerpFrictionCoeff;	// Friction coefficient, moving perpendicular
	double randomXPerpStdev;	// Stdev of the random perpendicular displacement
	double TFrictionCoeff;		// Friction coefficient, rotating
	double randomTStdev;		// Stdev of the random rotational displacement
	
	double neighborListBufferSq;	// The additional buffer distance for neighbor lists 
	double neighborListBufferSq4;	// The buffer divided by 4 (the max distSq a particle can move before
					// recalculating)
	double neighborListCutoffSq;	// The buffer added to the force cutoff
	double neighborListQuickCutoffSq;// If two rods are farther than this from each other,
					// they can't possibly be in each other's neighbor list.
		
	double alignmentParameter;	// Strength of the friction interaction
	double frictionCutoffSq;	// Max distance for friction interaction
	
	/** State Values **/
	int currentTime;
	bool stop;			// Set by the signal callback, save state and quit
	vector<Rod*> particleList;	// Contains pointers to all the rods
	vector<InteractionPair> neighborList;
	Random *randomGen;
	MotorGrid *motorGrid;
	bool output;			// Are we outputing data in this round?
	
	double totalDisp;		// Displacement from the initial positions
	double orientInitCorr;		// Correlation between current and initial rod orientations
	double orientAllCorr;		// Correlation between current rod orientations
	
	ofstream coordBinOut;		// Binary output of particle coordinates
	ofstream stepOut;		// ASCII output of interesting information per step
	ofstream debugOut;		// Used for whatever you want
	ofstream forceOut;
	ofstream clumpOut;
	
	bool failSafeMode;
	int failCount;
	int failFactor;
	streampos coordBinPos;
	streampos forceOutPos;
	streampos clumpOutPos;
	streampos stepOutPos;

private:	
	/** Private Functions **/
	void initValues();			// Calculate the constant values, etc.
	
	void calculateForces(bool reset);
		// Add on the current forces, optionally resetting them first

	void calculateForcesForRods(InteractionPair *ip);
		// Add on the forces between the two specified rods
	
	void applyForces (bool temp);
		// Apply the forces, optionally only temporarily.
		// Divides by two if not temporary, for the 2nd order predictor corrector
	
	void scaleBoxTo(double newSideLength);
	
	bool crossedRods(Rod *p, Rod *q);
	bool checkConsistency();
	void updateNeighbors();
	
	void updateClumps();
	bool drillDownClump(int newId, Rod *r);

	void saveState(bool saveTemp);
		// Save the simulation state to a file
	void loadState(istream &savedState);
		// Reload a state from the given file
	void revertState();
	void multiplyPBC(int sideFactor);
	void killRods(double aliveFraction);
	
	void openOut(ofstream &file, const char *name);
	
public:
	/* Convenience Functions */
	double dot(double v1[2], double v2[2]);	// Dot Product
	double cross(double v1[2], double v2[2]); // 2D Cross Products
	double wrap(double x);
};
