#include "MotorGrid.h"
#include "RodSimulation.h"

# ifdef __INTEL_COMPILER 
# include <mathimf.h> 
# else 
# include <math.h> 
# endif 


MotorGrid::MotorGrid(RodSimulation *_sim, double motorCount, double _motorForce, double captureRadius) {
	motorForce = _motorForce;
	sim = _sim;
	
	sideCount = floor(sqrt(motorCount));
	motorSpacing = sim->sideLength / sideCount;
	captureRadiusSq = captureRadius*captureRadius;
	quickCutoff = captureRadius + sim->rodLength2;
	quickCutoffSq = pow(captureRadius + sim->rodLength2,2);
}
	
void MotorGrid::addForceToRod(Rod *p) {
	double motorPos[2];
	
	for (int i=0; i<sideCount; i++) {
		for (int j=0; j<sideCount; j++) {
			/*if (sim->dot(p->uhat(), p->force) < -motorForce) {
				// If the rod is being pushed in the wrong direction
				// too hard, detach the motor.
				continue;
			}*/
			
			motorPos[0] = i*motorSpacing;
			motorPos[1] = j*motorSpacing;
			
			if (distanceToMotor(p, motorPos) <= captureRadiusSq) {
				p->force[0] += motorForce*p->uhat()[0];
				p->force[1] += motorForce*p->uhat()[1];
			}
		}
	}
}

double MotorGrid::distanceToMotor(Rod *p, double motorPos[2]) {
	double distSq;
	double cmDisp[2];
	
	cmDisp[0] = sim->wrap(motorPos[0]-p->pos[0]); cmDisp[1] = sim->wrap(motorPos[1]-p->pos[1]);
	
	// Don't do the hard part if we're farther than the quick cutoff
	/*if (cmDisp[0] > quickCutoff || cmDisp[0] < -quickCutoff 
	  || cmDisp[1] > quickCutoff || cmDisp[1] < -quickCutoff
	  || cmDisp[0]*cmDisp[0] + cmDisp[1]*cmDisp[1] > quickCutoffSq) {
		return DINFINITY;
	}*/
	
	
	double u = (cmDisp[0]*p->u()[0] + cmDisp[1]*p->u()[1])/(p->halfLength*p->halfLength);
	if (u < -1 || u > 1) {
		double distSqEndA = pow(cmDisp[0]-p->u()[0],2)+pow(cmDisp[1]-p->u()[1],2);
		double distSqEndB = pow(cmDisp[0]+p->u()[0],2)+pow(cmDisp[1]+p->u()[1],2);

		distSq = (distSqEndA < distSqEndB) ? distSqEndA : distSqEndB; // Use the smaller one
	}
	else {
		distSq = pow(cmDisp[0]-u*p->u()[0], 2) + pow(cmDisp[1]-u*p->u()[1], 2);
	}
	
	return distSq;
}
