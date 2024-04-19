//mixingPlane.H
#include <cmath>
#include <algorithm>
#include <vector>
#include <iterator>
#include "bndStructs.H"
#include "IOstream.H"
#include "turbulentTemperatureCoupledBaffleMixedMixingPlaneFvPatchScalarField.H"
namespace Foam
{
//using namespace std;

//Data swapping function for use within the mixingPlane function
bool dataSwap(int l, int k, std::vector<dataOutput>* dataOutPtr);

//Data sort function for use within mixingPlane function
bool dataSort(int numOfRows, std::vector<dataOutput>* dataOutPtr);

//mixing plane function for use on boudnary data
bool mixingPlane(std::vector<dataInput> dataIn, std::vector<dataOutput>* dataOutPtr);
}

