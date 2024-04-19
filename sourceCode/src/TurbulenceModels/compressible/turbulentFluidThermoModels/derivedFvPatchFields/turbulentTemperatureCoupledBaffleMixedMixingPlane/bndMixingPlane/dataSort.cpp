#include <cmath>
#include <algorithm>
#include <vector>
#include <iterator>
#include "bndStructs.h"
using namespace std;

bool dataSwap(int l, int k, vector<dataOutput>* dataOutPtr)
{
	if (l == k)
	{
		return true;
	}
	else
	{
		tempStr temp = {};

		//Save data in the temporary matrix - Tnew is unpopulated so ignore
		temp.x = (*dataOutPtr)[l].xOut;
		temp.y = (*dataOutPtr)[l].yOut;
		temp.z = (*dataOutPtr)[l].zOut;
		temp.T = (*dataOutPtr)[l].TOut;
		temp.origIndex = (*dataOutPtr)[l].origIndexOut;
		temp.r = (*dataOutPtr)[l].rOut;
		temp.theta = (*dataOutPtr)[l].thetaOut;

		//Overwrite lower matrix i with upper matrix j
		(*dataOutPtr)[l].xOut = (*dataOutPtr)[k].xOut;
		(*dataOutPtr)[l].yOut = (*dataOutPtr)[k].yOut;
		(*dataOutPtr)[l].zOut = (*dataOutPtr)[k].zOut;
		(*dataOutPtr)[l].TOut = (*dataOutPtr)[k].TOut;
		(*dataOutPtr)[l].origIndexOut = (*dataOutPtr)[k].origIndexOut;
		(*dataOutPtr)[l].rOut = (*dataOutPtr)[k].rOut;
		(*dataOutPtr)[l].thetaOut = (*dataOutPtr)[k].thetaOut;

		//Overwrite upper matrix j with saved lower matrix i
		(*dataOutPtr)[k].xOut = temp.x;
		(*dataOutPtr)[k].yOut = temp.y;
		(*dataOutPtr)[k].zOut = temp.z;
		(*dataOutPtr)[k].TOut = temp.T;
		(*dataOutPtr)[k].origIndexOut = temp.origIndex;
		(*dataOutPtr)[k].rOut = temp.r;
		(*dataOutPtr)[k].thetaOut = temp.theta;

		//Clear structure
		temp = {};

		return true;
	}
}

bool dataSort(int numOfRows, vector<dataOutput>* dataOutPtr)
{
	//Create temporary vector
	bool swapped;

	for (int i = 0; i < (numOfRows - 1); i++)
	{
		swapped = false;
		for (int j = 0; j < (numOfRows - i - 1); j++)
		{
			int m = j + 1;
			if ((*dataOutPtr)[j].rOut > (*dataOutPtr)[m].rOut)
			{
				swapped = dataSwap(j, m, dataOutPtr);
			}
			else if ((*dataOutPtr)[j].rOut == (*dataOutPtr)[m].rOut)
			{
				if ((*dataOutPtr)[j].xOut > (*dataOutPtr)[m].xOut)
				{
					swapped = dataSwap(j, m, dataOutPtr);
				}
				else if ((*dataOutPtr)[j].xOut == (*dataOutPtr)[m].xOut)
				{
					if ((*dataOutPtr)[j].thetaOut > (*dataOutPtr)[m].thetaOut)
					{
						swapped = dataSwap(j, m, dataOutPtr);
					}
				}
			}
		}

		if (swapped == false)
		{
			break;
		}
	}

	//Return true to mixingPlane
	return true;
}

