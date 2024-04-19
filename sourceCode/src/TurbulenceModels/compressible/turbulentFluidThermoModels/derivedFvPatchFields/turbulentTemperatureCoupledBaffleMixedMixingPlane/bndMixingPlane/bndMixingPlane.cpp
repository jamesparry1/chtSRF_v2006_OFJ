// bndMixingPlane.cpp
//This contains the main function.

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "mixingPlane.h"
#include "bndStructs.h"
#include "scalarField.H"

using namespace std;

scalarField bndMixingPlane(scalarField iList, scalarField xList, scalarField yList, scalarField zList, scalarField TList)
{
	//Initialise Vectors
	vector<dataInput> dataIn;
	vector<dataOutput> dataOut;
        bool processDataResult;

        //Convert Data to Structure
        for(int i = 0; i < iList.size(); i++)
        {
            dataIn[i].xIn = xList[i];
            dataIn[i].yIn = yList[i];
            dataIn[i].zIn = zList[i];
            dataIn[i].TIn = TList[i];
            dataIn[i].indexIn = iList[i];
        }

	//Run mixingPlane
	processDataResult = mixingPlane(dataIn, &dataOut);
	//Export
        float TOut = dataOut.TNewOut;

	//Tidy variables
	dataIn.clear();
	dataOut.clear();

	return TOut;
}
