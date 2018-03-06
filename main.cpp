#include <getopt.h>
#include <csignal>
#include <fstream>

#include "RodSimulation.h"

RodSimulation *sim;

void interrupt (int sig) {
	cout << "\nWaiting for current step to finish...\n";
	sim->stopSimulation();
}

void doNothing (int sig) {
	cout << "Caught signal " << sig << ", doing nothing...\n";
}

void usage() {
	cout << "Usage: " << 
		"bd --number=... --length=... --startDensity=... [--endDensity=...] --timestep=... --ticks=... --skipFrames=... --seed=... --push=...\n"
		<< "or " << "bd --resume=... [--endDensity=...] --timestep=... --ticks=... [--reset-clock] --skipFrames=... --seed=... --push=...\n";
}


int main (int argc, char * const argv[]) {
	int particleCount;
	double rodLength; // In units of sigma
	double startDensity;
	double endDensity;
	double timeStep;
	int stopTime;
	int outputSkip;
	int randomSeed;
	double selfPushForce;
	bool resetClock;
	double cutoff;
	double circularBoundary;
	double endCircularBoundary;
	ifstream resumeFile;
	int multiplySideLength;
	double friction;
	bool sameFile;
	double aliveFraction;
	
	/*** Default Values ***/
	particleCount = -1;		//
	rodLength = -1;			// In units of sigma
	startDensity = -1;		// Number density
	endDensity = -1;		// this will be reset to the start density if it remains at -1
	timeStep = -1;			//
	stopTime = -1;			// Time to stop at
	outputSkip = -1;		// Output every x frames
	randomSeed = 14;		// 
	selfPushForce = -1;		// Self-propelling force
	cutoff = 0;
	resetClock = true;		// Reset the clock to 0, if resuming
	circularBoundary = -1;	// Don't put this in a petri dish, by default
	endCircularBoundary = -1;
	multiplySideLength = 1; // Don't multiply the system, by default
	friction = -1;
	sameFile = false;		// Use the old coords-bin file (this option does nothing yet)
	aliveFraction = -1;
	/**********************/
	
	if (argc < 2) {
		usage();
		exit(1);
	}
	
	while (1)
	{
		static struct option long_options[] =
		{
			{"help",		no_argument	 , 0, '?'},
			{"number",		required_argument, 0, 'n'},
			{"length",		required_argument, 0, 'l'},
			{"startDensity",	required_argument, 0, 'd'},
			{"endDensity",		required_argument, 0, 'D'},
			{"timestep",		required_argument, 0, 'h'},
			{"ticks",		required_argument, 0, 't'},
			{"skipFrames",		required_argument, 0, 'f'},
			{"friction",		required_argument, 0, 'F'},
			{"seed",		required_argument, 0, 's'},
			{"push",		required_argument, 0, 'p'},
			{"resume",		required_argument, 0, 'r'},
			{"reset-clock",		no_argument	 , 0, 'R'},
			{"cutoff",		required_argument, 0, 'c'},
			{"circularBoundary",	required_argument, 0, 'C'},
			{"endCircularBoundary",	required_argument, 0, 'E'},
			{"multiplySideLength", required_argument, 0, 'm'},
			{"sameFile",			no_argument, 0, 'S'},
			{"aliveFraction",			required_argument, 0, 'a'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		int c = getopt_long (argc, argv, "?n:l:d:D:h:t:f:F:s:p:r:c:C:E:m:RSa:",
							 long_options, &option_index);
		
		/* Detect the end of the options. */
		if (c == -1)
			break;
		
		switch (c)
		{
			case '?':
				usage();
				exit(1);

			case 'n':
				particleCount = atoi(optarg);
				break;
				
			case 'l':
				rodLength = atof(optarg);
				break;
				
			case 'd':
				startDensity = atof(optarg);
				break;

			case 's':
				randomSeed = atoi(optarg);
				break;
				
			case 'D':
				endDensity = atof(optarg);
				break;

			case 'b':
				endDensity = atof(optarg);
				break;
				
			case 'h':
				timeStep = atof(optarg);
				break;
				
			case 't':
				stopTime = atoi(optarg);
				break;
				
			case 'f':
				outputSkip = atoi(optarg);
				break;
			
				
			case 'p':
				selfPushForce = atof(optarg);
				break;
			
			case 'r':
				resumeFile.open(optarg, ios::in | ios::binary);
				if (resumeFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
					return 1;
				}
				break;
				
			case 'c':
				cutoff = atof(optarg);
				break;
			
			case 'C':
				circularBoundary = atof(optarg);
				break;
			
			case 'E':
				endCircularBoundary = atof(optarg);
				break;
				
			case 'm':
				multiplySideLength = atoi(optarg);
				break;
			
			case 'R':
				resetClock = true;
				break;
				
			case 'F':
				friction = atof(optarg);
				break;
			
			case 'S':
				sameFile = true;
				break;
				
			case 'a':
				aliveFraction = atof(optarg);
				break;
			
			default:
				abort();
		}
	}
	
	if (circularBoundary != -1 && (endDensity != -1 || startDensity != -1)) {
		cerr << "Cannot set density with a circular boundary.\n";
		abort();
	}
	if (resumeFile.is_open() && (startDensity != -1 || circularBoundary != -1 || rodLength != -1 || particleCount != -1)) {
		cerr << "Cannot set initial conditions for resume.\n";
		abort();
	} 
	if (!resumeFile.is_open() && ( particleCount == -1 || rodLength == -1 || timeStep == -1 ||
								  stopTime == -1 || outputSkip == -1 || selfPushForce == -1 || friction == -1)) {
		cerr << "Did not set all the initial conditions.\n";
		abort();
	}
	if (!resumeFile.is_open() && sameFile) {
		cerr << "Can't use an old coordinate file for new simulation.\n";
		abort();
	}
	
	if (endDensity == -1)
		endDensity = startDensity;
	if (endCircularBoundary == -1)
		endCircularBoundary = circularBoundary;
	
	if (resumeFile.is_open()) {
		sim = new RodSimulation(randomSeed, endDensity, endCircularBoundary, timeStep, stopTime, outputSkip, selfPushForce, cutoff, friction, aliveFraction, // Simulation parameters
					resumeFile, resetClock, multiplySideLength); // Resume parameters
	}
	else {
		sim = new RodSimulation(particleCount, rodLength, startDensity, circularBoundary, // Initial parameters
					randomSeed, endDensity, endCircularBoundary, timeStep, stopTime, outputSkip, selfPushForce, cutoff, friction, aliveFraction); // Simlation parameters
	}
	
	signal(SIGINT, interrupt);
	//signal(SIGUSR2, interrupt);
	//signal(SIGUSR1, doNothing);
	
	return sim->runSimulation();
	
}
