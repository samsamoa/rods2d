#include <getopt.h>
#include <csignal>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <stdlib.h>

#include "Rod.h"

void usage() {
	cerr << "Usage: " << 
	"bd --number=... --length=... --startDensity=... [--endDensity=...] --timestep=... --ticks=... --skipFrames=... --seed=... --push=...\n"
	<< "or " << "bd --resume=... [--endDensity=...] --timestep=... --ticks=... [--reset-clock] --skipFrames=... --seed=... --push=...\n";
}


// Start Time, End Time, output skip, rod length, particle count, side length
int particleCount;
double rodLength, rodLength2; // In units of sigma
double startSl;
int absStartTime;
int absStopTime;
int currentTime;
int outputSkip;
int startTime;
double sideLength;
double halfSideLength;
double quickCutoffSq;
double windowQuickCutoffSq;
double neighborQuickCutoffSq;
double neighborCutoffSq;
bool periodicBoundary = true;
double orientationCutoff;
double displacementDifferenceSq;
vector<Rod> *particleList;
int currentIndex, comparisonIndex;
vector<Rod> lastParticleList;
vector< vector<Rod *> > clumpList;
bool analysis = false;
bool doClumping = false;
bool doDisplacementComparison = false;
bool doOrientationCutoff = false;
bool onlyAliveClumping ;
ofstream GofRFile;
int windowHeight, windowWidth;
bool *alive;
bool aliveKnown;
int resolution;
int speedup;
int *GofR;
int *GofRUniform;
double *GofROriented, *GofROrientedDirectional;
int *clusterSizes;
int *aliveClusterSizes;
int *deadClusterSizes;
int neighborCounts[4];
int clusterPurities[100];
double pixelSize;

void updateNeighbors();
void updateClumps();
bool drillDownClump(int newId, Rod *r, vector<Rod *> &clump);
double dot(double v1[2], double v2[2]);
double cross(double v1[2], double v2[2]);

using namespace std;

int main (int argc, char * const argv[]) {
	ifstream coordsFile;
	ofstream clumpFile;
	ofstream aliveInfoFile;
	ifstream aliveFile;
	ofstream purityFile;
	ofstream clumpInfoFile;
	
	/***** Default Values ***/
	particleCount = -1;		//
	rodLength = -1;			// In units of sigma
	startSl = -1;
	absStartTime = -1;
	absStopTime = -1;			// Time to stop at
	currentTime = -1;
	outputSkip = -1;		// Output every x frames
	windowHeight = 10;
	windowWidth = 5;
	speedup = 1;
	resolution = 3;
	/** Important Parameters **/
	//double neighborCutoff = 2; // use for purity
	double neighborCutoff = 3.0; // use for cluster identification
	double displacementDifference = 2.0;
	orientationCutoff = cos(M_PI/6.0); // Gompper value
	int frameDifference = 0;
	doOrientationCutoff = true;
	onlyAliveClumping = false;
	/**************************/
	
	
	if (argc < 6) {
		usage();
		return 1;
	}
	
	while (1)
	{
		static struct option long_options[] =
		{
			{"help",			no_argument	 , 0, '?'},
			{"particleCount",	required_argument, 0, 'n'},
			{"rodLength",		required_argument, 0, 'l'},
			{"startSl",			required_argument, 0, 'L'},
			{"absStartTime",	required_argument, 0, 't'},
			{"absStopTime",		required_argument, 0, 'T'},
			{"outputSkip",		required_argument, 0, 's'},
			{"time",			required_argument, 0, 'c'},
			{"coordsFile",		required_argument, 0, 'f'},
			{"GofRFile",		required_argument, 0, 'C'},
			{"clumpFile",		required_argument, 0, 'p'},
			{"doClumping",		no_argument, 0, 'm'},
			{"clumpInfoFile",	required_argument, 0, 'i'},
			{"aliveInfoFile",	required_argument, 0, 'd'},
			{"analysis",		no_argument		 , 0, 'A'},
			{"noPBC",			no_argument		, 0, 'P'},
			{"windowWidth",		required_argument, 0, 'w'},
			{"windowHeight",	required_argument, 0, 'h'},
			{"aliveness",		required_argument, 0, 'a'},
			{"speedup",			required_argument, 0, 'S'},
			{"resolution",		required_argument, 0, 'R'},
			{"purityFile",		required_argument, 0, 'u'},
			{"frameDifference",	required_argument, 0, 'F'},
			{"displacementDifference",required_argument, 0, 'D'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		int c = getopt_long (argc, argv, "?n:l:L:t:T:c:s:f:p:AC:w:h:a:S:Pd:R:u:i:F:D:m",
							 long_options, &option_index);
		
		/* Detect the end of the options. */
		if (c == -1)
			break;
		
		switch (c)
		{
			case '?':
				usage();
				return 1;
				
			case 'n':
				particleCount = atoi(optarg);
				break;
				
			case 'l':
				rodLength = atof(optarg);
				break;
			
			case 'L':
				startSl = atof(optarg);
				break;
				
			case 't':
				absStartTime = atoi(optarg);
				break;

			case 'T':
				absStopTime = atoi(optarg);
				break;
				
			case 'c':
				startTime = atoi(optarg);
				break;
				
			case 's':
				outputSkip = atoi(optarg);
				break;
			
			case 'f':
				coordsFile.open(optarg, ios::in | ios::binary);
				if (coordsFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
					return 1;
				}
				break;
			
			case 'p':
				clumpFile.open(optarg, ios::out | ios::binary);
				if (clumpFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
					return 1;
				}
				doClumping = true;
				break;
			
			case 'C':
				GofRFile.open(optarg, ios::out | ios::binary);
				if (GofRFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
					return 1;
				}
				break;

			case 'm':
				doClumping = true;	
				break;

			case 'a':
				aliveFile.open(optarg, ios::in | ios::binary);
				if (aliveFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
					return 1;
				}
				break;

			case 'd':
				aliveInfoFile.open(optarg, ios::out | ios::binary);
				if (aliveInfoFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
					return 1;
				}
				break;
				
			case 'i':
				clumpInfoFile.open(optarg, ios::out | ios::binary);
				if (clumpInfoFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
					return 1;
				}
				break;
				
			case 'u':
				purityFile.open(optarg, ios::out | ios::binary);
				if (purityFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
					return 1;
				}
				break;
				
			case 'A':
				analysis = true;
				break;

			case 'P':
				periodicBoundary = false;
				break;
			
			case 'w':
				windowWidth = atoi(optarg);
				break;
			
			case 'h':
				windowHeight = atoi(optarg);
				break;
			
			case 'S':
				speedup = atoi(optarg);
				break;
			
			case 'R':
				resolution = atoi(optarg);
				break;

			case 'D':
				displacementDifference = atof(optarg);
				break;
				
			case 'F':
				frameDifference = atoi(optarg);
				break;
			
			default:
				cerr << "Invalid argument: " << (char)c << "\n";
				abort();
		}
	}
	
	
	alive = new bool[particleCount];
	if (!aliveFile.is_open() || !aliveFile.good()) {
		cerr << "No aliveness information..." << endl;
		aliveKnown = false;
	}
	else {
		cerr << "Reading aliveness..." << endl;
		aliveFile.seekg(-sizeof(bool)*particleCount, std::ios::end);
		for (int i=0; i<particleCount; i++) {
			aliveFile.read((char *)(&alive[i]), sizeof(bool));
		}
		aliveKnown = true;
	}
	aliveFile.close();

	
	sideLength = startSl;
	halfSideLength = sideLength*0.5;
		
	neighborCutoffSq=neighborCutoff*neighborCutoff;
	neighborQuickCutoffSq = pow(rodLength + neighborCutoff, 2.0);

	if (analysis) {
		windowQuickCutoffSq = (windowHeight+rodLength2+2)*(windowHeight+rodLength2+2)+(windowWidth+rodLength2+2)*(windowWidth+rodLength2+2);
		
		pixelSize = 1.0/resolution;
		
		GofR = new int[2*2*2*2*windowWidth*resolution*windowHeight*resolution];
		GofROriented = new double[2*2*2*2*windowWidth*resolution*windowHeight*resolution];
		GofROrientedDirectional = new double[2*2*2*2*windowWidth*resolution*windowHeight*resolution];
		GofRUniform = new int[2*2*windowWidth*resolution];
		for (int i=0; i<2*2*2*2*windowWidth*resolution*windowHeight*resolution; i++) {
			GofR[i]=0;
			GofROriented[i]=0;
			GofROrientedDirectional[i]=0;
		}
		for (int i=0; i<2*2*windowWidth*resolution; i++) {
			GofRUniform[i]=0;
		}		
	}
	
	clusterSizes = new int[particleCount];
	deadClusterSizes = new int[particleCount];
	aliveClusterSizes = new int[particleCount];
	
	displacementDifferenceSq = displacementDifference*displacementDifference;

	rodLength2 = rodLength*0.5;
	
	if (frameDifference == 0)
		doDisplacementComparison = false;
	else
		doDisplacementComparison = true;
	
	particleList = new vector<Rod>[frameDifference+1];
	
	for (int i=0; i<=frameDifference; i++) {
		streampos s = ((streampos)(startTime-absStartTime))*24*particleCount/outputSkip;
		coordsFile.seekg(s);
		particleList[i].clear();
		particleList[i].reserve(particleCount);
		for (int j=0; j<particleCount; j++) {
			particleList[i].push_back(Rod(coordsFile, rodLength2, j, 1));
		}
	}
	
	
	for (currentTime=absStartTime; currentTime < absStopTime; currentTime += outputSkip) {
		if (currentTime < startTime) {
			if (doClumping && clumpFile.is_open()) {
				int cID=0;
				for (int i=0; i<particleCount; i++) {
					clumpFile.write((char*)&cID, sizeof(int));
				}
			}
			continue;
		}
		if ((currentTime - absStartTime) % (outputSkip*speedup) != 0) {
			if (doClumping && clumpFile.is_open()) {
				for (int i=0; i<particleCount; i++) {
					clumpFile.write((char*)&particleList[currentIndex][i].clumpId, sizeof(int));
				}
			}
			continue;			
		}
		
		for (int i=0; i<particleCount; i++) {
			clusterSizes[i] = 0;
			deadClusterSizes[i] = 0;
			aliveClusterSizes[i] = 0;
		}
		for (int i=0; i<100; i++) {
			clusterPurities[i] = 0;
		}

		
		currentIndex = (currentTime/(outputSkip*speedup)) % (frameDifference+1);
		comparisonIndex = (currentIndex+1)%(frameDifference+1);	
		
		streampos s = ((streampos)(currentTime-absStartTime))*24*particleCount/outputSkip;
		coordsFile.seekg(s);
		particleList[currentIndex].clear();
		particleList[currentIndex].reserve(particleCount);
		for (int i=0; i<particleCount; i++) {
			particleList[currentIndex].push_back(Rod(coordsFile, rodLength2, i, 1));
		}
		
		neighborCounts[0]=neighborCounts[1]=neighborCounts[2]=neighborCounts[3]=0;
		
		updateNeighbors();
		if (doClumping) {
			updateClumps();
		}
		
		if (aliveInfoFile.is_open()) {
			for (int clumpId=0; clumpId < clumpList.size(); clumpId++) {
				/*int aliveCount=0;
				for (int i=0; i<clumpList[clumpId].size(); i++) {
					if (alive[clumpList[clumpId][i]->idNum])
						aliveCount += 1;
				}
				aliveInfoFile << currentTime << " " << clumpId << " " << clumpList[clumpId].size() << " " << aliveCount << "\n";
				*/
				clusterSizes[clumpList[clumpId].size()-1] += 1;
				/*if (clumpList[clumpId].size() > 2) {
					int purityBin = (int)floor(100*((double)aliveCount/clumpList[clumpId].size()));
					if (purityBin == 100)
						purityBin = 99;
					clusterPurities[purityBin] += clumpList[clumpId].size();
					
					if (purityBin >= 80)
						aliveClusterSizes[clumpList[clumpId].size()-1] += 1;
					else if (purityBin <= 19)
						deadClusterSizes[clumpList[clumpId].size()-1] += 1;
				}*/
				
				/*** This code can be used to compule how polar the clumps are ***/
				 
				/*double polarVector[2] = {0,0};
				for (int i=0; i<clumpList[clumpId].size(); i++) {
					Rod *r = clumpList[clumpId][i];
					polarVector[0] += r->_uhat[0];
					polarVector[1] += r->_uhat[1];
				}

				double norm = sqrt(dot(polarVector, polarVector));
							
				polarVector[0] /= norm;
				polarVector[1] /= norm;
				
				//cout << "c  " << clumpId << ": " << clumpList[clumpId].size() << " " << norm/clumpList[clumpId].size() << endl;
				//cout << "   " << longAxis[1]-longAxis[0] << " x " << shortAxis[1]-shortAxis[0] << endl;*/
				
				/*******************************************************************/
			}
		}
		
		if (clumpInfoFile.is_open()) {
			int clusterCount = 0;
			int cutoffClusterCount = 0;
			int averageClusterCount = 0;
			int weightedAverageClusterCount = 0;
			int cutoffAverageClusterCount = 0;
			
			for (int i=0; i<particleCount; i++) {
				clusterCount += clusterSizes[i];
				averageClusterCount += clusterSizes[i]*(i+1);
				weightedAverageClusterCount += clusterSizes[i]*(i+1)*(i+1);
				if (i>1) {
					cutoffClusterCount += clusterSizes[i];
					cutoffAverageClusterCount += clusterSizes[i]*(i+1);
				}
			}
			
			double averageClusterSize = (double)averageClusterCount/clusterCount;
			double weightedAverageClusterSize = (double)weightedAverageClusterCount/(clusterCount*particleCount);
			double cutoffAverageClusterSize = (double)cutoffAverageClusterCount/cutoffClusterCount;
			
			clumpInfoFile.write((char *)&currentTime, sizeof(int));
			clumpInfoFile.write((char *)clusterSizes, sizeof(int)*particleCount);
			clumpInfoFile.write((char *)aliveClusterSizes, sizeof(int)*particleCount);
			clumpInfoFile.write((char *)deadClusterSizes, sizeof(int)*particleCount);
			clumpInfoFile.write((char *)clusterPurities, sizeof(int)*100);
			clumpInfoFile.write((char *)&averageClusterSize, sizeof(double));
			clumpInfoFile.write((char *)&weightedAverageClusterSize, sizeof(double));
			clumpInfoFile.write((char *)&cutoffAverageClusterSize, sizeof(double));
		}
		
		if (purityFile.is_open()) {
			// Time, dead purity, alive purity.
			purityFile << currentTime << " " << (double)neighborCounts[0]/(neighborCounts[0]+neighborCounts[1]) << " " << (double)neighborCounts[3]/(neighborCounts[2]+neighborCounts[3]) << "\n";
		}
		
		if (doClumping && clumpFile.is_open()) {
			for (int i=0; i<particleCount; i++) {
				clumpFile.write((char*)&particleList[currentIndex][i].clumpId, sizeof(int));
			}
		}
//		if ((currentTime - startTime)/(outputSkip*speedup) % 20 == 0)
			cerr << currentTime << endl;
	}
	
	if (analysis) {
		GofRFile.write((char*)GofR, sizeof(int)*2*2*2*2*windowHeight*windowWidth*resolution*resolution);
		GofRFile.write((char*)GofROriented, sizeof(double)*2*2*2*2*windowHeight*windowWidth*resolution*resolution);
		GofRFile.write((char*)GofROrientedDirectional, sizeof(double)*2*2*2*2*windowHeight*windowWidth*resolution*resolution);
		GofRFile.write((char*)GofRUniform, sizeof(int)*2*2*windowWidth*resolution);
	}

	
	GofRFile.close();
	clumpFile.close();
	aliveInfoFile.close();
	coordsFile.close();
	purityFile.close();
	clumpInfoFile.close();
}

inline void fixPeriodicBoundary(double *a) {
	if (a[0] < -halfSideLength)
		a[0] +=sideLength;
	if (a[1] < -halfSideLength)
		a[1] +=sideLength;
	if (a[0] > halfSideLength)
		a[0] -=sideLength;
	if (a[1] > halfSideLength)
		a[1] -=sideLength;
}

inline void fixPeriodicBoundary(double *a, double *b) {
	if (a[0]-b[0] < -halfSideLength)
		a[0] +=sideLength;
	if (a[1]-b[1] < -halfSideLength)
		a[1] +=sideLength;
	if (a[0]-b[0] > halfSideLength)
		a[0] -=sideLength;
	if (a[1]-b[1] > halfSideLength)
		a[1] -=sideLength;	
}


inline void addToGofR(double *displacement, bool alive1, bool alive2, int *grid) {
	if (-windowWidth <= displacement[0] && displacement[0] <= windowWidth
		&& -windowHeight <= displacement[1] && displacement[1] <= windowHeight) {
		
	grid[ alive1*(windowWidth*windowHeight*resolution*resolution*2*2*2)
	+ alive2*(windowWidth*windowHeight*resolution*resolution*2*2)
	+ (windowWidth*resolution+(int)floor(displacement[0]/pixelSize))*(windowHeight*2*resolution)
	+ (windowHeight*resolution+(int)floor(displacement[1]/pixelSize))] += 1;
	}
}

void updateNeighbors() {
	
	for (int m=0; m<particleCount; m++) {
		particleList[currentIndex][m].neighbors.clear();
		particleList[currentIndex][m].clusterNeighbors.clear();
		particleList[currentIndex][m].clumpId=-1;
	}
	
	Rod *p, *q, *pPrev, *qPrev;
	double pDisp[2], qDisp[2], dispDifference[2];
	for (int i=0; i<particleCount; i++) {
		if (onlyAliveClumping && !alive[i])
			continue;
		
		p = &(particleList[currentIndex][i]);
		double *uhat1 = p->_uhat;
		
		/** Uncomment to enable "old neighbors" method **/
		/* bool oldNeighbors[particleCount];
		for (int j=0; j<particleCount; j++)
			oldNeighbors[j]=false;
		
		for (vector<Rod*>::iterator it = pPrev->neighbors.begin(); it!=pPrev->neighbors.end(); ++it) {
			oldNeighbors[(*it)->idNum] = true;
		}*/
		
		if (doDisplacementComparison) {
			pPrev = &particleList[comparisonIndex][i];
			pDisp[0] = p->pos[0] - pPrev->pos[0];
			pDisp[1] = p->pos[1] - pPrev->pos[1];
			fixPeriodicBoundary(pDisp);
		}

		for (int j=i+1; j<particleCount; j++) {
			if (onlyAliveClumping && !alive[j])
				continue;
			
			q = &(particleList[currentIndex][j]);
						
			double *uhat2 = q->_uhat;			
			
			double orientationCorrelation = dot(uhat1, uhat2);
			//if (aliveKnown && (!alive[i] || !alive[j]))
			//	orientationCorrelation = fabs(orientationCorrelation);
			
			/** Uncomment to enable orientation cutoff **/
			if (!analysis && doOrientationCutoff && (orientationCorrelation < orientationCutoff)) // Difference in direction is <20º
				continue;
			

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
			
			double cmDist[2] = {qloc[0]-p->pos[0], qloc[1]-p->pos[1]};
			double cmDistSq = dot(cmDist, cmDist);
						
			if (analysis && cmDistSq > windowQuickCutoffSq)
				continue;
			
			if (doDisplacementComparison) {
				qPrev = &particleList[comparisonIndex][j];
				qDisp[0] = q->pos[0] - qPrev->pos[0];
				qDisp[1] = q->pos[1] - qPrev->pos[1];
				fixPeriodicBoundary(qDisp);
			}
			
			if (analysis) {
				double projectedDistance1[2] = {(cross(cmDist, uhat1)), dot(cmDist, uhat1)}; // Rod p at center
				double projectedDistance2[2] = {(-cross(cmDist, uhat2)), -dot(cmDist, uhat2)}; // Rod q at center
				
				//to cmDist, we are adding either uhat1 or uhat2, in order to separate the rods
				// into "balls"


				double tempDistance1[2], tempDistance2[2];
				
				int *thisGridLocation1 = &GofR[alive[i]*(windowWidth*windowHeight*resolution*resolution*2*2*2)
											   + alive[j]*(windowWidth*windowHeight*resolution*resolution*2*2)];
				int *thisGridLocation2 = &GofR[alive[j]*(windowWidth*windowHeight*resolution*resolution*2*2*2)
											   + alive[i]*(windowWidth*windowHeight*resolution*resolution*2*2)];

				double *thisGridLocationOriented1 = &GofROriented[alive[i]*(windowWidth*windowHeight*resolution*resolution*2*2*2)
											   + alive[j]*(windowWidth*windowHeight*resolution*resolution*2*2)];
				double *thisGridLocationOriented2 = &GofROriented[alive[j]*(windowWidth*windowHeight*resolution*resolution*2*2*2)
											   + alive[i]*(windowWidth*windowHeight*resolution*resolution*2*2)];

				double *thisGridLocationOrientedDirectional1 = &GofROrientedDirectional[alive[i]*(windowWidth*windowHeight*resolution*resolution*2*2*2)
																  + alive[j]*(windowWidth*windowHeight*resolution*resolution*2*2)];
				double *thisGridLocationOrientedDirectional2 = &GofROrientedDirectional[alive[j]*(windowWidth*windowHeight*resolution*resolution*2*2*2)
																  + alive[i]*(windowWidth*windowHeight*resolution*resolution*2*2)];
				
				int *thisGridLocationUniform1 = &GofRUniform[alive[i]*(windowWidth*resolution*2)
													  + alive[j]*(windowWidth*resolution)];
				int *thisGridLocationUniform2 = &GofRUniform[alive[j]*(windowWidth*resolution*2)
													  + alive[i]*(windowWidth*resolution)];
				

				double increment1[2] = {(cross(uhat2, uhat1)), dot(uhat2, uhat1)};
				double increment2[2] = {(cross(uhat1, uhat2)), dot(uhat1, uhat2)};
				/*projectedDistance1[0] -= increment1[0]*rodLength2;
				projectedDistance1[1] -= increment1[1]*rodLength2;
				projectedDistance2[0] -= increment2[0]*rodLength2;
				projectedDistance2[1] -= increment2[1]*rodLength2;				
				for (int ix=0; ix<=rodLength; ix++) {*/
				for (int ix=0; ix<=0; ix++) {
					tempDistance1[0] = projectedDistance1[0] + increment1[0]*ix;
					tempDistance1[1] = projectedDistance1[1] + increment1[1]*ix;

					tempDistance2[0] = projectedDistance2[0] + increment2[0]*ix;
					tempDistance2[1] = projectedDistance2[1] + increment2[1]*ix;
					
					double tempAbsDistance1 = sqrt(dot(tempDistance1, tempDistance1));
					double tempAbsDistance2 = sqrt(dot(tempDistance2, tempDistance2));
						
					if (-windowWidth <= tempDistance1[0] && tempDistance1[0] <= windowWidth
						&& -windowHeight <= tempDistance1[1] && tempDistance1[1] <= windowHeight) {
						thisGridLocation1[(windowWidth*resolution+(int)floor(tempDistance1[0]/pixelSize))*(windowHeight*2*resolution)
										  + (windowHeight*resolution+(int)floor(tempDistance1[1]/pixelSize))] += 1;
						thisGridLocationOriented1[(windowWidth*resolution+(int)floor(tempDistance1[0]/pixelSize))*(windowHeight*2*resolution)
										  + (windowHeight*resolution+(int)floor(tempDistance1[1]/pixelSize))] += orientationCorrelation*orientationCorrelation*2-1.0;
						thisGridLocationOrientedDirectional1[(windowWidth*resolution+(int)floor(tempDistance1[0]/pixelSize))*(windowHeight*2*resolution)
										  + (windowHeight*resolution+(int)floor(tempDistance1[1]/pixelSize))] += orientationCorrelation;
					}
					if (-windowWidth <= tempDistance2[0] && tempDistance2[0] <= windowWidth
						&& -windowHeight <= tempDistance2[1] && tempDistance2[1] <= windowHeight) {
						thisGridLocation2[(windowWidth*resolution+(int)floor(tempDistance2[0]/pixelSize))*(windowHeight*2*resolution)
										  + (windowHeight*resolution+(int)floor(tempDistance2[1]/pixelSize))] += 1;
						thisGridLocationOriented2[(windowWidth*resolution+(int)floor(tempDistance2[0]/pixelSize))*(windowHeight*2*resolution)
										  + (windowHeight*resolution+(int)floor(tempDistance2[1]/pixelSize))] += orientationCorrelation*orientationCorrelation*2-1.0;
						thisGridLocationOrientedDirectional2[(windowWidth*resolution+(int)floor(tempDistance2[0]/pixelSize))*(windowHeight*2*resolution)
										  + (windowHeight*resolution+(int)floor(tempDistance2[1]/pixelSize))] += orientationCorrelation;
					}
					
					if (tempAbsDistance1 < windowWidth) {
						thisGridLocationUniform1[(int)floor(tempAbsDistance1/pixelSize)] += 1;
					}
					if (tempAbsDistance2 < windowWidth) {
						thisGridLocationUniform2[(int)floor(tempAbsDistance2/pixelSize)] += 1;
					}
				}
				
				continue;
			}
			
						
			if (cmDistSq > neighborQuickCutoffSq)
				continue;
			
			if (doDisplacementComparison) {
				dispDifference[0] = qDisp[0]-pDisp[0];
				dispDifference[1] = qDisp[1]-pDisp[1];
				
				if (doClumping && (dot(dispDifference, dispDifference) > displacementDifferenceSq))
					continue;
			}

			double *endvec1 = p->_u;
			double *endvec2 = q->_u;
			
			
			//double neighborDisp[2];
						
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
			
			bool success = false;
			
			if (!isnan(dend_aa) && dend_aa < neighborCutoffSq) {
				success = true;
				//neighborDisp[0]=end_aa[0];
				//neighborDisp[1]=end_aa[1];
			}
			else if (!isnan(dend_bb) && dend_bb < neighborCutoffSq) {
				success = true;
				//neighborDisp[0]=end_bb[0];
				//neighborDisp[1]=end_bb[1];
			}
			else if (!isnan(dend_ab) && dend_ab < neighborCutoffSq) {
				success = true;
				//neighborDisp[0]=end_ab[0];
				//neighborDisp[1]=end_ab[1];
			}
			else if (!isnan(dend_ba) && dend_ba < neighborCutoffSq) {
				success = true;
				//neighborDisp[0]=end_ba[0];
				//neighborDisp[1]=end_ba[1];
			}
			
			if (!success) {			
				double x_12a[2]={0,0}, x_12b[2]={0,0}, x_21a[2]={0,0}, x_21b[2]={0,0};
				double dmid_12a, dmid_12b, dmid_21a, dmid_21b;
				double u; // paramaterize the rod as (Xcm + u*endvec) u = -1...1
				double cmDisp[2];
				
				cmDisp[0] = end2a[0]-p->pos[0]; cmDisp[1] = end2a[1]-p->pos[1];
				u = (cmDisp[0]*endvec1[0] + cmDisp[1]*endvec1[1])/(rodLength2*rodLength2);
				if (u < -1 || u > 1)
					dmid_12a = NAN;
				else {
					x_12a[0] = p->pos[0] + u*endvec1[0];
					x_12a[1] = p->pos[1] + u*endvec1[1];
					dmid_12a = ((end2a[0] - x_12a[0])*(end2a[0] - x_12a[0]) + (end2a[1] - x_12a[1])*(end2a[1] - x_12a[1]));
				}
				
				cmDisp[0] = end2b[0]-p->pos[0]; cmDisp[1] = end2b[1]-p->pos[1];
				u = (cmDisp[0]*endvec1[0] + cmDisp[1]*endvec1[1])/(rodLength2*rodLength2);
				if (u < -1 || u > 1)
					dmid_12b = NAN;
				else {
					x_12b[0] = p->pos[0] + u*endvec1[0];
					x_12b[1] = p->pos[1] + u*endvec1[1];
					dmid_12b = ((end2b[0] - x_12b[0])*(end2b[0] - x_12b[0]) + (end2b[1] - x_12b[1])*(end2b[1] - x_12b[1]));
				}
				
				cmDisp[0] = end1a[0]-qloc[0]; cmDisp[1] = end1a[1]-qloc[1];
				u = (cmDisp[0]*endvec2[0] + cmDisp[1]*endvec2[1])/(rodLength2*rodLength2);
				if (u < -1 || u > 1)
					dmid_21a = NAN;
				else {
					x_21a[0] = qloc[0] + u*endvec2[0];
					x_21a[1] = qloc[1] + u*endvec2[1];
					dmid_21a = ((end1a[0] - x_21a[0])*(end1a[0] - x_21a[0]) + (end1a[1] - x_21a[1])*(end1a[1] - x_21a[1]));
				}
				
				cmDisp[0] = end1b[0]-qloc[0]; cmDisp[1] = end1b[1]-qloc[1];
				u = (cmDisp[0]*endvec2[0] + cmDisp[1]*endvec2[1])/(rodLength2*rodLength2);
				if (u < -1 || u > 1)
					dmid_21b = NAN;
				else {
					x_21b[0] = qloc[0] + u*endvec2[0];
					x_21b[1] = qloc[1] + u*endvec2[1];
					dmid_21b = ((end1b[0] - x_21b[0])*(end1b[0] - x_21b[0]) + (end1b[1] - x_21b[1])*(end1b[1] - x_21b[1]));
				}
				
				if (!isnan(dmid_12a) && dmid_12a < neighborCutoffSq) {
					success = true;
				}
				else if (!isnan(dmid_12b) && dmid_12b < neighborCutoffSq)
					success = true;
				else if (!isnan(dmid_21a) && dmid_21a < neighborCutoffSq)
					success = true;
				else if (!isnan(dmid_21b) && dmid_21b < neighborCutoffSq)
					success = true;
			}
			
			if (success) {
				neighborCounts[alive[i]*2+alive[j]] += 1;
				neighborCounts[alive[j]*2+alive[i]] += 1;
				
				// Only look at pairs that did not move too much relative to each other.
				if (dot(dispDifference, dispDifference) < displacementDifferenceSq) {
					/** Remove comments here for "old neighbors" method **/
					//p->neighbors.push_back(q);
					//q->neighbors.push_back(p);
					
					//if (oldNeighbors[q->idNum]) {
						p->clusterNeighbors.push_back(q);
						q->clusterNeighbors.push_back(p);
					//}
				}
			}
		}	
	}
}

// Calculate clumps by seeing which rods are interacting.
void updateClumps() {
	int newId = 0;
	clumpList.clear();
	clumpList.push_back(vector<Rod *>());
	
	for (int i=0; i<particleCount; i++) {
		if (onlyAliveClumping && !alive[i])
			continue;
		
		if (drillDownClump(newId, &particleList[currentIndex][i], clumpList[newId])) {
			newId++;
			clumpList.push_back(vector<Rod *>());
		}
	}
	
	//cout << "id " << newId << endl;
	
	if (clumpList.back().size() == 0) {
		clumpList.pop_back();
	}
}

bool drillDownClump(int newId, Rod *r, vector<Rod *> &clump) {
	if (r->clumpId != -1)
		return false;
	
	r->clumpId = newId;
	clump.push_back(r);
	
	for (int i=0; i < r->clusterNeighbors.size(); i++) {
		drillDownClump(newId, r->clusterNeighbors[i], clump);
	}
	
	return true;
}

inline double dot(double v1[2], double v2[2]) {
	return (v1[0]*v2[0] + v1[1]*v2[1]);
}

inline double cross(double v1[2], double v2[2]) {
	return (v1[0]*v2[1]-v1[1]*v2[0]);
}
