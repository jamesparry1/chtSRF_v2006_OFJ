//bndStructs.h
#ifndef MY_STRUCTURES
#define MY_STRUCTURES

#include <vector>
using namespace std;

//Data input
// Input data structure
typedef struct dataInput {
	// Create fields
        double xIn;
        double yIn;
        double zIn;
        double TIn;
        double procIndexIn;
        double indexIn;

	// Constructor
	dataInput(
                double x,
                double y,
                double z,
                double T,
                double procIndex,
                double index
        ) : xIn(x),
                yIn(y),
                zIn(z),
                TIn(T),
                procIndexIn(procIndex),
                indexIn(index)
	{}
};

//Data output
// Output data structure
typedef struct dataOutput {
	// Create fields
        double xOut;
        double yOut;
        double zOut;
        double TOut;
        double procOrigIndexOut;
        double origIndexOut;
        double rOut;
        double thetaOut;
        double TnewOut;

	// Constructor
	dataOutput(
                double x,
                double y,
                double z,
                double T,
                double procOrigIndex,
                double origIndex,
                double r,
                double theta,
                double Tnew
	) : xOut(x),
		yOut(y),
		zOut(z),
		TOut(T),
                procOrigIndexOut(procOrigIndex),
		origIndexOut(origIndex),
		rOut(r),
		thetaOut(theta),
		TnewOut(Tnew)
	{}
};

//Temporary matrix
// Temporary
typedef struct tempStr {
        double x;
        double y;
        double z;
        double T;
        double procOrigIndex;
        double origIndex;
        double r;
        double theta;
        double Tnew;
};

#endif
