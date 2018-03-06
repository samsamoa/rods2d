#pragma once
/*
 *  InteractionPair.h
 *  bd
 *
 *  Created by Sam McCandlish on 12/22/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "Rod.h"

using namespace std;

class InteractionPair;

typedef enum {
	endA_t = 0,
	endB_t = 1,
	mid_t = 2,
} rodLoc_t;


class InteractionPair {
public:
	Rod *rod1;
	Rod *rod2;
	
	bool subNeighbors[9];
	bool isZero;
	
	bool interactedThisTime;
	double currentDistSq;
	
	InteractionPair(Rod *_rod1, Rod *_rod2);

	bool inSubList(rodLoc_t aLoc, rodLoc_t bLoc);
	void setSubList(rodLoc_t aLoc, rodLoc_t bLoc);
};
