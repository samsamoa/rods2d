#include "InteractionPair.h"

InteractionPair::InteractionPair(Rod *_rod1, Rod *_rod2) {
	rod1 = _rod1;
	rod2 = _rod2;
	isZero = true;
	subNeighbors[0] = subNeighbors[1] = subNeighbors[2] = 
		subNeighbors[3] = subNeighbors[4] = subNeighbors[5] =
		subNeighbors[6] = subNeighbors[7] = subNeighbors[8] = false;
}

bool InteractionPair::inSubList(rodLoc_t aLoc, rodLoc_t bLoc) {
	return subNeighbors[aLoc + 3*bLoc];
}

void InteractionPair::setSubList(rodLoc_t aLoc, rodLoc_t bLoc) {
	subNeighbors[aLoc + 3*bLoc]=true;
	if (aLoc == endA_t || aLoc == endB_t)
		subNeighbors[mid_t + 3*bLoc]=true;
	if (bLoc == endA_t || bLoc == endB_t)
		subNeighbors[aLoc + 3*mid_t]=true;
	
	isZero=false;
}

