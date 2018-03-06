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
int rodLengthInt;
double startSl;
int absStartTime;
int absStopTime;
int currentTime;
int outputSkip;
int startTime;
double sideLength;
double halfSideLength;
bool periodicBoundary = true;
vector<Rod> particleList;
bool *alive;
bool aliveKnown;
double aliveFraction;
int speedup;
int sideBoxCount = 50;
int boxesPerN = 5;

double dot(double v1[2], double v2[2]);
double cross(double v1[2], double v2[2]);

using namespace std;

inline void fixPeriodicBoundary(double *a) {
	if (a[0] < 0)
		a[0] +=sideLength;
	if (a[1] < 0)
		a[1] +=sideLength;
	if (a[0] >= sideLength)
		a[0] -= sideLength;
	if (a[1] >= sideLength)
		a[1] -=sideLength;
}

inline void minimalImageBox(double *a) {
	if (a[0] < -sideBoxCount/2)
		a[0] += sideBoxCount;
	if (a[1] < -sideBoxCount/2)
		a[1] += sideBoxCount;
	if (a[0] >= sideBoxCount/2)
		a[0] -= sideBoxCount;
	if (a[1] >= sideBoxCount/2)
		a[1] -= sideBoxCount;
}

int main (int argc, char * const argv[]) {
	ifstream coordsFile;
	ifstream aliveFile;
	ofstream boxFile;
	ofstream fluctuationFile;
	ofstream fluctuationExpFile;
	ofstream GofRFile;
	ofstream meanSqDispFile;
	
	/*** Default Values ***/
	particleCount = -1;		//
	rodLength = -1;			// In units of sigma
	startSl = -1;
	absStartTime = -1;
	absStopTime = -1;			// Time to stop at
	currentTime = -1;
	outputSkip = -1;		// Output every x frames
	speedup = 1;
	/**********************/
	
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
			{"aliveness",		required_argument, 0, 'a'},
			{"speedup",			required_argument, 0, 'S'},
			{"boxFile",			required_argument, 0, 'b'},
			{"sideBoxCount",	required_argument, 0, 'C'},
			{"boxesPerN",		required_argument, 0, 'B'},
			{"fluctuationFile",	required_argument, 0, 'F'},
			{"GofRFile",		required_argument, 0, 'G'},
			{"fluctuationExpFile",required_argument, 0, 'E'},
			{"meanSqDispFile",required_argument, 0, 'A'},
			{0, 0, 0, 0}
		};
		/* getopt_long stores the option index here. */
		int option_index = 0;
		
		int c = getopt_long (argc, argv, "?n:l:L:t:T:s:c:f:a:S:b:C:B:F:G:E:A:",
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

				
			case 'b':
				boxFile.open(optarg, ios::out | ios::binary);
				if (boxFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
					return 1;
				}
				break;
				
			case 'a':
				aliveFile.open(optarg, ios::in | ios::binary);
				if (aliveFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
				}
				break;

			case 'F':
				fluctuationFile.open(optarg, ios::out | ios::binary);
				if (fluctuationFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
				}
				break;
			
			case 'G':
				GofRFile.open(optarg, ios::out | ios::binary);
				if (GofRFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
				}
				break;

			case 'E':
				fluctuationExpFile.open(optarg, ios::out | ios::binary);
				if (fluctuationExpFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
				}
				break;
				
			case 'A':
				meanSqDispFile.open(optarg, ios::out | ios::binary);
				if (meanSqDispFile.fail()) {
					cerr << "Unable to open " << optarg << ".\n";
				}
				break;

			case 'S':
				speedup = atoi(optarg);
				break;
				
			case 'C':
				sideBoxCount = atoi(optarg);
				break;

			case 'B':
				boxesPerN = atoi(optarg);
				break;
				
			default:
				cerr << "Invalid argument: " << (char)c << "\n";
				abort();
		}
	}
	
	
	alive = new bool[particleCount];
	int numberAlive = 0;
	if (!aliveFile.is_open() || !aliveFile.good()) {
		cerr << "No aliveness information..." << endl;
		aliveKnown = false;
	}
	else {
		cerr << "Reading aliveness..." << endl;
		aliveFile.seekg(-sizeof(bool)*particleCount, std::ios::end);
		for (int i=0; i<particleCount; i++) {
			aliveFile.read((char *)(&alive[i]), sizeof(bool));
			if (alive[i])
				numberAlive++;
		}
		aliveFraction = (double)numberAlive/particleCount;
		aliveKnown = true;
	}
	aliveFile.close();

	
	sideLength = startSl;
	halfSideLength = sideLength*0.5;
	
	rodLength2 = rodLength*0.5;
	rodLengthInt = (int)rodLength;
	
	double boxSize = sideLength/sideBoxCount;
	double boxArea = boxSize*boxSize;
	double normalization = boxArea * (rodLength + 1);
	
	int *aliveBoxCounts, *deadBoxCounts;
	
	aliveBoxCounts = new int[sideBoxCount*sideBoxCount*((absStopTime-startTime)/(outputSkip*speedup))];
	deadBoxCounts = new int[sideBoxCount*sideBoxCount*((absStopTime-startTime)/(outputSkip*speedup))];
	
	int frameCount = ((absStopTime-startTime)/(outputSkip*speedup));

	double *polarization = new double[frameCount*sideBoxCount*sideBoxCount*2];
	
	cerr << ((absStopTime-startTime)/(outputSkip*speedup)) << " frames" << endl;
	
	/**************************
	 *   Histogram the Rods   * -> boxFile
	 * ************************/
	vector<Rod> initParticleList;
	vector<Rod> previousParticleList;
	int XwrapAround[particleCount];
	int YwrapAround[particleCount];
	
	for (int i=0; i<particleCount; i++) {
		XwrapAround[i]=0;
		YwrapAround[i]=0;
	}
	
	for (currentTime=absStartTime; currentTime < absStopTime; currentTime += outputSkip*speedup) {
		if (currentTime < startTime) {
			continue;
		}
		
		int frame = ((currentTime-startTime)/(outputSkip*speedup));
		int frameOffset = frame*sideBoxCount*sideBoxCount;
		
		for (int i=0; i<sideBoxCount*sideBoxCount; i++) {
			aliveBoxCounts[frameOffset+i]=0;	
			deadBoxCounts[frameOffset+i]=0;	
		}
		
		streampos s = ((streampos)(currentTime-absStartTime))*24*particleCount/outputSkip;
		coordsFile.seekg(s);

		previousParticleList.clear();
		previousParticleList=particleList;
		particleList.clear();
		particleList.reserve(particleCount);
		for (int i=0; i<particleCount; i++) {
			particleList.push_back(Rod(coordsFile, rodLength2, i, 1));
		}
		if (currentTime == startTime) {
			initParticleList=particleList;
			previousParticleList=particleList;
		}
		
		
		double *direction;
		double divisionPos[2];
		
		double meanSqDisp = 0, aliveMeanSqDisp=0, deadMeanSqDisp=0;
		
		for (int i=0; i<particleCount; i++) {
			direction = particleList[i].uhat();
			divisionPos[0] = particleList[i].pos[0] - direction[0]*rodLength2;
			divisionPos[1] = particleList[i].pos[1] - direction[1]*rodLength2;
			
			for (int j=0; j <= rodLengthInt; j++) {
				fixPeriodicBoundary(divisionPos);				
				
				if ((!aliveKnown) || alive[i]) {
					aliveBoxCounts[frameOffset + (int)floor(divisionPos[0]/boxSize)*sideBoxCount + (int)floor(divisionPos[1]/boxSize)] += 1;
					double *u = &polarization[frameOffset*2+(int)floor(divisionPos[0]/boxSize)*sideBoxCount*2+(int)floor(divisionPos[1]/boxSize)*2];
					u[0] += direction[0];
					u[1] += direction[1];
				}
				else
					deadBoxCounts[frameOffset + (int)floor(divisionPos[0]/boxSize)*sideBoxCount + (int)floor(divisionPos[1]/boxSize)] += 1;
				
				divisionPos[0] += direction[0];
				divisionPos[1] += direction[1];
			}
			
			// Calculate the absolute position
			double dX = particleList[i].pos[0]-previousParticleList[i].pos[0];
			double dY = particleList[i].pos[1]-previousParticleList[i].pos[1];			
			if (dX > halfSideLength)
				XwrapAround[i] -= 1;
			else if (dX < -halfSideLength)
				XwrapAround[i] += 1;
			
			if (dY > halfSideLength)
				YwrapAround[i] -= 1;
			else if (dY < -halfSideLength)
				YwrapAround[i] += 1;
			
			double absdX = particleList[i].pos[0] - initParticleList[i].pos[0] + sideLength*XwrapAround[i];
			double absdY = particleList[i].pos[1] - initParticleList[i].pos[1] + sideLength*YwrapAround[i];
			meanSqDisp += absdX*absdX+absdY*absdY;
			if (alive[i]) {
				aliveMeanSqDisp += absdX*absdX+absdY*absdY;
			}
			else {
				deadMeanSqDisp += absdX*absdX+absdY*absdY;
			}

		}
		
		meanSqDisp /= particleCount;
		aliveMeanSqDisp /= numberAlive;
		deadMeanSqDisp /= (particleCount - numberAlive);
		
		for (int x=0; x<sideBoxCount; x++) {
			for (int y=0; y<sideBoxCount; y++) {
				double *u = &polarization[frameOffset*2+x*sideBoxCount*2+y*2];
				double norm = sqrt(dot(u,u));
				if (norm == 0)
					continue;
				u[0] /= norm;
				u[1] /= norm;
			}
		}
		
		boxFile.write((char *)(aliveBoxCounts+frameOffset), sizeof(int)*sideBoxCount*sideBoxCount);
		boxFile.write((char *)(deadBoxCounts+frameOffset), sizeof(int)*sideBoxCount*sideBoxCount);
		
		meanSqDispFile << currentTime << " " << meanSqDisp << " " << aliveMeanSqDisp << " " << deadMeanSqDisp << "\n";
		
		if (currentTime % 500000 == 0)
			cerr << currentTime << endl;
	}
	
	meanSqDispFile.close();
	boxFile.write((char *)polarization, sizeof(double)*frameCount*sideBoxCount*sideBoxCount*2);
	
	/********************************
	 * Density Fluctuation Exponent * -> fluctuationExpFile
	 * ******************************/
	int boxIncr = sideBoxCount/boxesPerN;
	for (int boxNumber = 1; boxNumber <= sideBoxCount/3; boxNumber++) {
		//cerr << n << endl;
		double meanMean = 0, aliveMeanMean=0, deadMeanMean=0;
		double meanStdev = 0, aliveMeanStdev=0, deadMeanStdev=0;
		
		long int totalSumSq=0, aliveTotalSumSq=0, deadTotalSumSq=0;
		long int totalSum=0, aliveTotalSum=0, deadTotalSum=0;
		
		int distribution[particleCount];
		for (int i=0; i<particleCount; i++) {
			distribution[i] = 0;
		}
		
		for (int sx = 0; sx < boxesPerN; sx++)  {
			for (int sy = 0; sy < boxesPerN; sy++)  {
				// We are looking at groups of n^2 boxes, formed in a square, with the lower left corner at 
				// the coordinates (sx, sy) (in # of boxes).  The boxes also wrap around the PBCs.
				long int sum=0, aliveSum=0, deadSum=0;
				long int sumSq=0, aliveSumSq=0, deadSumSq=0;
				
				for (int frame=0; frame<frameCount; frame++) {
					long int frameSum = 0, aliveFrameSum=0, deadFrameSum=0;
					// Loop over all sub-boxes in this box
					for (int bx = sx*boxIncr; bx<(sx*boxIncr+boxNumber); bx++) {
						for (int by = sy*boxIncr; by<(sy*boxIncr+boxNumber); by++) {
							frameSum += aliveBoxCounts[frame*sideBoxCount*sideBoxCount + (bx % sideBoxCount)*sideBoxCount + (by%sideBoxCount)];
							frameSum += deadBoxCounts[frame*sideBoxCount*sideBoxCount + (bx % sideBoxCount)*sideBoxCount + (by%sideBoxCount)];

							aliveFrameSum += aliveBoxCounts[frame*sideBoxCount*sideBoxCount + (bx % sideBoxCount)*sideBoxCount + (by%sideBoxCount)];
							deadFrameSum += deadBoxCounts[frame*sideBoxCount*sideBoxCount + (bx % sideBoxCount)*sideBoxCount + (by%sideBoxCount)];
						}
					}
					//if (sx == 0 && sy == 0)
					//	cout << frameSum << " ";
					distribution[(int)(frameSum/(rodLength+1))] += 1;
					sum += frameSum;
					aliveSum += aliveFrameSum;
					deadSum += deadFrameSum;
					sumSq += frameSum*frameSum;
					aliveSumSq += aliveFrameSum*aliveFrameSum;
					deadSumSq += deadFrameSum*deadFrameSum;
				}
				/*
				double mean = (double)sum/(frameCount*(rodLength+1));
				double meanSquare = (double)sumSq/(frameCount*(rodLength+1)*(rodLength+1));
				double squareMean = mean*mean;
				double stdev = sqrt(meanSquare-squareMean);

				double aliveMean = (double)aliveSum/(frameCount*(rodLength+1));
				double aliveMeanSquare = (double)aliveSumSq/(frameCount*(rodLength+1)*(rodLength+1));
				double aliveSquareMean = aliveMean*aliveMean;
				double aliveStdev = sqrt(aliveMeanSquare-aliveSquareMean);

				double deadMean = (double)deadSum/(frameCount*(rodLength+1));
				double deadMeanSquare = (double)deadSumSq/(frameCount*(rodLength+1)*(rodLength+1));
				double deadSquareMean = deadMean*deadMean;
				double deadStdev = sqrt(deadMeanSquare-deadSquareMean);
				*/
				totalSumSq += sumSq;
				totalSum += sum;

				aliveTotalSumSq += aliveSumSq;
				aliveTotalSum += aliveSum;

				deadTotalSumSq += deadSumSq;
				deadTotalSum += deadSum;
				
				/*
				meanMean += mean;
				meanStdev += stdev;

				aliveMeanMean += aliveMean;
				aliveMeanStdev += aliveStdev;

				deadMeanMean += deadMean;
				deadMeanStdev += deadStdev;
				 */

				//if (sx == 0 && sy == 0)
				//	cout << "\n";
				
			}
		}
		
		for (int i=0; i<particleCount; i++) {
			fluctuationExpFile << distribution[i] << " ";
		}
		fluctuationExpFile << "\n";
		double totalMean = (double)totalSum/(frameCount*boxesPerN*boxesPerN*(rodLength+1));
		double totalMeanSq = (double)totalSumSq/(frameCount*boxesPerN*boxesPerN*(rodLength+1)*(rodLength+1));
		double totalStdev = sqrt(totalMeanSq - totalMean*totalMean);

		double aliveTotalMean = (double)aliveTotalSum/(frameCount*boxesPerN*boxesPerN*(rodLength+1));
		double aliveTotalMeanSq = (double)aliveTotalSumSq/(frameCount*boxesPerN*boxesPerN*(rodLength+1)*(rodLength+1));
		double aliveTotalStdev = sqrt(aliveTotalMeanSq - aliveTotalMean*aliveTotalMean);
		
		double deadTotalMean = (double)deadTotalSum/(frameCount*boxesPerN*boxesPerN*(rodLength+1));
		double deadTotalMeanSq = (double)deadTotalSumSq/(frameCount*boxesPerN*boxesPerN*(rodLength+1)*(rodLength+1));
		double deadTotalStdev = sqrt(deadTotalMeanSq - deadTotalMean*deadTotalMean);
		
		meanStdev /= (boxesPerN*boxesPerN);
		meanMean /= (boxesPerN*boxesPerN);
		//cout << meanMean << " " << meanStdev << endl;
		fluctuationExpFile << totalMean << " " << totalStdev << " " << totalMeanSq << " "
			<< aliveTotalMean << " " << aliveTotalStdev << " " << aliveTotalMeanSq << " "
			<< deadTotalMean << " " << deadTotalStdev << " " << deadTotalMeanSq << endl;
	}
	
	fluctuationExpFile.close();
	
	
	/************************************
	 * Single Time Density Fluctuations * -> fluctuationFile
	 * **********************************/
	int aliveCount, deadCount;
	for (int frame=0; frame<frameCount; frame++) {
		for (int boxNumber = 1; boxNumber <= sideBoxCount/3; boxNumber++) {
			long int aliveTotal=0, aliveTotalSq=0;
			long int deadTotal=0, deadTotalSq=0;
			long int onlyDeadTotal=0, onlyDeadTotalSq=0, onlyDeadBoxCount=0;
			long int total=0, totalSq=0;
			
			double aliveFractionDev;
			double aliveFractionTotalSqDev=0, aliveFractionWeightedTotalSqDev=0, aliveFractionBoxCount=0;
			double aliveFractionTotalAbsDev=0, aliveFractionWeightedTotalAbsDev=0;
			
			for (int sx = 0; sx < boxesPerN; sx++)  {
				for (int sy = 0; sy < boxesPerN; sy++)  {
					// We are looking at groups of n^2 boxes, formed in a square, with the lower left corner at 
					// the coordinates (sx*boxIncr, sy*boxIncr) (in # of boxes).  The boxes also wrap around the PBCs.

					// Loop over all sub-boxes in this box
					aliveCount=0;
					deadCount=0;
					
					for (int bx = sx*boxIncr; bx<(sx*boxIncr+boxNumber); bx++) {
						for (int by = sy*boxIncr; by<(sy*boxIncr+boxNumber); by++) {
							aliveCount += aliveBoxCounts[frame*sideBoxCount*sideBoxCount + (bx % sideBoxCount)*sideBoxCount + (by%sideBoxCount)];
							deadCount += deadBoxCounts[frame*sideBoxCount*sideBoxCount + (bx % sideBoxCount)*sideBoxCount + (by%sideBoxCount)];
						}
					}
					
					aliveTotal += aliveCount;
					aliveTotalSq += aliveCount*aliveCount;
					
					deadTotal += deadCount;
					deadTotalSq += deadCount*deadCount;
					
					if (aliveCount == 0) {
						onlyDeadTotal += deadCount;
						onlyDeadTotalSq += deadCount*deadCount;
						onlyDeadBoxCount += 1;
					}
					
					total += aliveCount;
					total += deadCount;
					totalSq += (aliveCount+deadCount)*(aliveCount+deadCount);
					
					if (!(aliveCount==0 && deadCount==0)) {
						aliveFractionDev = (double)(aliveCount)/(aliveCount+deadCount) - aliveFraction ;
						
						aliveFractionTotalSqDev += pow(aliveFractionDev,2);
						aliveFractionWeightedTotalSqDev += pow(aliveFractionDev,2)*(aliveCount+deadCount);
						
						aliveFractionTotalAbsDev += fabs(aliveFractionDev);
						aliveFractionWeightedTotalAbsDev += fabs(aliveFractionDev)*(aliveCount+deadCount);

						aliveFractionBoxCount += 1;
					}
					// else: No rods in the box.  Adds nothing to weighted.
					// Ignore the box for unweighted (think about this, it makes sense)
				}
			}
			
			double aliveStdev = sqrt((double)aliveTotalSq/(boxesPerN*boxesPerN) - (double)(aliveTotal*aliveTotal)/(boxesPerN*boxesPerN*boxesPerN*boxesPerN))/(normalization*boxNumber*boxNumber);
			double deadStdev = sqrt((double)deadTotalSq/(boxesPerN*boxesPerN) - (double)(deadTotal*deadTotal)/(boxesPerN*boxesPerN*boxesPerN*boxesPerN))/(normalization*boxNumber*boxNumber);
			double onlyDeadStdev = sqrt((double)onlyDeadTotalSq/(onlyDeadBoxCount) - (double)(onlyDeadTotal*onlyDeadTotal)/(onlyDeadBoxCount*onlyDeadBoxCount))/(normalization*boxNumber*boxNumber);
			
			double aliveFractionStdev = sqrt((double)aliveFractionTotalSqDev/(aliveFractionBoxCount));
			double aliveFractionWeightedStdev = sqrt((double)aliveFractionWeightedTotalSqDev/(total));
			double aliveFractionAbsDev = ((double)aliveFractionTotalAbsDev/(aliveFractionBoxCount));
			double aliveFractionWeightedAbsDev = ((double)aliveFractionWeightedTotalAbsDev/(total));

			
			double totalStdev = sqrt((double)totalSq/(boxesPerN*boxesPerN) - (double)(total*total)/(boxesPerN*boxesPerN*boxesPerN*boxesPerN))/(normalization*boxNumber*boxNumber);
			
			fluctuationFile << startTime+outputSkip*speedup*frame << " " << boxNumber << " " << totalStdev << " " << aliveStdev << " " << deadStdev << " " << onlyDeadStdev << " " << aliveFractionStdev << " " << aliveFractionWeightedStdev << " " << aliveFractionAbsDev << " " << aliveFractionWeightedAbsDev << "\n";			
		}
		
		if (frame % 50 == 0)
			cerr << frame << endl;
	}
	
	fluctuationFile.close();
	
	/*
	
	long long int AliveGofR[sideBoxCount];
	long long int AliveGofRNorm[sideBoxCount];
	long long int DeadGofR[sideBoxCount];
	long long int DeadGofRNorm[sideBoxCount];
	
	long long int AliveGofXY[sideBoxCount][sideBoxCount];
	
	for (int i=0; i<sideBoxCount; i++) {
		AliveGofR[i]=0;
		AliveGofRNorm[i]=0;
		DeadGofR[i]=0;
		DeadGofRNorm[i]=0;
		
		for (int j=0; j<sideBoxCount; j++) {
			AliveGofXY[i][j]=0;
		}
	}
	
	int maxD = (int)floor(sideBoxCount/2.0);
	int maxDSq = maxD*maxD;
	long int XYnorm = 0;
	
	for (int frame=frameCount*4/5; frame<frameCount; frame+= 50) {
		int frameOffset = frame*sideBoxCount*sideBoxCount;
		
		for (int sx = 0; sx<sideBoxCount; sx++)  {
			for (int sy = 0; sy<sideBoxCount; sy++)  {

				XYnorm += 2;
				int centerAliveCount = aliveBoxCounts[frameOffset + sx*sideBoxCount + sy];
				int centerDeadCount = deadBoxCounts[frameOffset + sx*sideBoxCount + sy];
				double *uCenter = &polarization[frameOffset*2+sx*sideBoxCount*2+sy*2];
				
				for (int dx = 0; dx<(sideBoxCount-sx); dx++) {
					for (int dy = 0; dy<sideBoxCount; dy++) {
						int dSq = (min(sideBoxCount-dx,dx)*min(sideBoxCount-dx,dx)+min(sideBoxCount-dy,dy)*min(sideBoxCount-dy,dy));
						if (dSq > maxDSq)
							continue;
						int d = (int)floor(sqrt(dSq));						
						
						int ox = (dx + sx) % sideBoxCount;
						int oy = (dy + sy) % sideBoxCount;						
						double cmDist[2] = {ox-sx, oy-sy};
						minimalImageBox(cmDist);
						
						int totalOffset = frameOffset + (ox)*sideBoxCount + (oy);
						int outAliveCount    = aliveBoxCounts[totalOffset];
						int outDeadCount    = deadBoxCounts[totalOffset];
						double *uOut = &polarization[totalOffset*2];
						
						double projectedDistance1[2] = {fabs(cross(cmDist, uCenter)), fabs(dot(cmDist, uCenter))};
						double projectedDistance2[2] = {fabs(-cross(cmDist, uOut)), fabs(-dot(cmDist, uOut))};
						
						AliveGofXY[(int)floor(projectedDistance1[0])][(int)floor(projectedDistance1[1])] += centerAliveCount*outAliveCount;
						AliveGofXY[(int)floor(projectedDistance2[0])][(int)floor(projectedDistance2[1])] += centerAliveCount*outAliveCount;
						
						AliveGofR[d] += 2*centerAliveCount*outAliveCount;
						AliveGofRNorm[d] += (outAliveCount+centerAliveCount);
						
						DeadGofR[d] += 2*centerDeadCount*outDeadCount;
						DeadGofRNorm[d] += (outDeadCount+centerDeadCount);							
					}
				}
			}
		}
		if ((frame) % 10 == 0) {
			cerr << frame << endl;
		}			
	}
	
	GofRFile << "{{";
	for (int i=0; i<(sideBoxCount); i++) {
		if (AliveGofRNorm[i]==0) break;
		if (i != 0)
			GofRFile << ",";
		GofRFile << ((double)AliveGofR[i]/AliveGofRNorm[i])/normalization;
	}
	GofRFile << "},\n\n{";
	for (int i=0; i<(sideBoxCount); i++) {
		if (DeadGofRNorm[i]==0) break;
		if (i != 0)
			GofRFile << ",";
		GofRFile << ((double)DeadGofR[i]/DeadGofRNorm[i])/normalization;
	}
	GofRFile << "},\n\n" << endl;
	int boxLimit = (sideBoxCount/(2.0*sqrt(2.0)));
	for (int i=0; i<boxLimit; i++) {
		GofRFile << "{";
		for (int j=0; j<boxLimit; j++) {
			GofRFile << ((double)AliveGofXY[i][j]/XYnorm)/normalization;
			if (j != boxLimit-1)
				GofRFile << ",";
		}
		GofRFile << "}\n";
		if (i != boxLimit-1)
			GofRFile << ",";
	}
	GofRFile << "}";
	
	cerr << "done!" << endl;
	
	GofRFile.close();
	 
	*/
	
	boxFile.close();
	fluctuationFile.close();
	coordsFile.close();
}

double dot(double v1[2], double v2[2]) {
	return (v1[0]*v2[0] + v1[1]*v2[1]);
}

double cross(double v1[2], double v2[2]) {
	return (v1[0]*v2[1]-v1[1]*v2[0]);
}
