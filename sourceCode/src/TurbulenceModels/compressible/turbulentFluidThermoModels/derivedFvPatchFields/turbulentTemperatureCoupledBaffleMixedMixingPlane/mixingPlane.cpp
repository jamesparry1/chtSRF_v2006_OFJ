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

bool dataSwap(int l, int k, std::vector<dataOutput>* dataOutPtr)
{
    if (l == k)
    {
        return true;
    }
    else
    {
        //Initialise temporary vector
        tempStr temp = {};
        //Save data in the temporary matrix - Tnew is unpopulated so ignore
        temp.x = (*dataOutPtr)[l].xOut;
        temp.y = (*dataOutPtr)[l].yOut;
        temp.z = (*dataOutPtr)[l].zOut;
        temp.T = (*dataOutPtr)[l].TOut;
        temp.procOrigIndex = (*dataOutPtr)[l].procOrigIndexOut;
        temp.origIndex = (*dataOutPtr)[l].origIndexOut;
        temp.r = (*dataOutPtr)[l].rOut;
        temp.theta = (*dataOutPtr)[l].thetaOut;

        //Overwrite lower matrix i with upper matrix j
        (*dataOutPtr)[l].xOut = (*dataOutPtr)[k].xOut;
        (*dataOutPtr)[l].yOut = (*dataOutPtr)[k].yOut;
        (*dataOutPtr)[l].zOut = (*dataOutPtr)[k].zOut;
        (*dataOutPtr)[l].TOut = (*dataOutPtr)[k].TOut;
        (*dataOutPtr)[l].procOrigIndexOut = (*dataOutPtr)[k].procOrigIndexOut;
        (*dataOutPtr)[l].origIndexOut = (*dataOutPtr)[k].origIndexOut;
        (*dataOutPtr)[l].rOut = (*dataOutPtr)[k].rOut;
        (*dataOutPtr)[l].thetaOut = (*dataOutPtr)[k].thetaOut;

        //Overwrite upper matrix j with saved lower matrix i
        (*dataOutPtr)[k].xOut = temp.x;
        (*dataOutPtr)[k].yOut = temp.y;
        (*dataOutPtr)[k].zOut = temp.z;
        (*dataOutPtr)[k].TOut = temp.T;
        (*dataOutPtr)[k].procOrigIndexOut = temp.procOrigIndex;
        (*dataOutPtr)[k].origIndexOut = temp.origIndex;
        (*dataOutPtr)[k].rOut = temp.r;
        (*dataOutPtr)[k].thetaOut = temp.theta;

        //Clear structure
        temp = {};

        return true;
    }
}

bool dataSort(int numOfRows, std::vector<dataOutput>* dataOutPtr)
{
    //Initialise swap check
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

bool mixingPlane(std::vector<dataInput> dataIn, std::vector<dataOutput>* dataOutPtr)
{
    //Begin averaging process
    try {
        // Populate output matrix
        float k = 0.0;
        float dp = 10000.0; //set number of decimal places

//        if(debug)
//        {
//            //DEBUG FLAG UPDATE
//            Info<< "Entering mixingPlane Function" << endl;
//            //Output dataIn
//            Info<< "Data In =" << dataIn.size() << endl;
//            Info<< "Data Out =" << dataOutPtr << endl;
//        }

        for (int i = 0; i < dataIn.size(); i++)
        {
            // Fill dataOut with dataIn information
            dataOutPtr->push_back(
                dataOutput(
                    dataIn[i].xIn,	//Fill x
                    dataIn[i].yIn,	//Fill y
                    dataIn[i].zIn,	//Fill z
                    dataIn[i].TIn,	//Fill T
                    dataIn[i].procIndexIn, //Fill processor original index
                    i,		//Fill origIndex
                    round((sqrt(pow(dataIn[i].zIn, 2.0) + pow(dataIn[i].yIn, 2.0)))*dp)/dp,	//Fill r
                    round((atan2(dataIn[i].yIn, dataIn[i].zIn))*dp)/dp,	//Fill theta
                    k		//Fill blank for TNew
                ));
        }

        //DEBUG FLAG UPDATE
        //Foam::Info<< "dataOut Allocation Successful" << endl;

        // Create temporary vector
        tempStr temp = {};
        const int numOfRows = dataIn.size();

        //Sort data
        bool sortOut = dataSort(numOfRows, dataOutPtr);

        // Create average for TNew
        const int numOfDiv = 360; //NOTE: CHANGE FOR NUMBER OF DIVISIONS 
        const int numOfAvg = (numOfRows / numOfDiv);
        float sum, avg;
        for (int w = 0; w < numOfAvg; w++)
        {
            //Re-initialise Sum
            sum = 0.0;
            for (int q = 0; q < numOfDiv; q++)
            {
                //Create sum of temperatures
                sum += (*dataOutPtr)[q + (numOfDiv * w)].TOut;
            }
            //Create Average
            avg = sum / numOfDiv;
            //Repopulate TNew in dataOut
            //std::fill((*dataOut).[(numOfDiv * l)], (*dataOut).[(numOfDiv * l + 1)], avg);

            //Alternate repopulation
            for (int e = 0; e < numOfDiv; e++)
            {
                (*dataOutPtr)[e + (numOfDiv * w)].TnewOut = avg;
            }
        }

        //Re-sort by Original Index
        for (int i = 0; i < (numOfRows - 1); i++)
        {
            for (int j = 0; j < (numOfRows - i - 1); j++)
            {
                int m = j + 1;
                if ((*dataOutPtr)[j].origIndexOut > (*dataOutPtr)[m].origIndexOut)
                {
                    //Save data in the temporary matrix - Tnew is unpopulated so ignore
                    temp.x = (*dataOutPtr)[j].xOut;
                    temp.y = (*dataOutPtr)[j].yOut;
                    temp.z = (*dataOutPtr)[j].zOut;
                    temp.T = (*dataOutPtr)[j].TOut;
                    temp.procOrigIndex = (*dataOutPtr)[j].procOrigIndexOut;
                    temp.origIndex = (*dataOutPtr)[j].origIndexOut;
                    temp.r = (*dataOutPtr)[j].rOut;
                    temp.theta = (*dataOutPtr)[j].thetaOut;
                    temp.Tnew = (*dataOutPtr)[j].TnewOut;

                    //Overwrite lower matrix i with upper matrix j
                    (*dataOutPtr)[j].xOut = (*dataOutPtr)[m].xOut;
                    (*dataOutPtr)[j].yOut = (*dataOutPtr)[m].yOut;
                    (*dataOutPtr)[j].zOut = (*dataOutPtr)[m].zOut;
                    (*dataOutPtr)[j].procOrigIndexOut = (*dataOutPtr)[m].procOrigIndexOut;
                    (*dataOutPtr)[j].origIndexOut = (*dataOutPtr)[m].origIndexOut;
                    (*dataOutPtr)[j].TOut = (*dataOutPtr)[m].TOut;
                    (*dataOutPtr)[j].rOut = (*dataOutPtr)[m].rOut;
                    (*dataOutPtr)[j].thetaOut = (*dataOutPtr)[m].thetaOut;
                    (*dataOutPtr)[j].TnewOut = (*dataOutPtr)[m].TnewOut;

                    //Overwrite upper matrix j with saved lower matrix i
                    (*dataOutPtr)[m].xOut = temp.x;
                    (*dataOutPtr)[m].yOut = temp.y;
                    (*dataOutPtr)[m].zOut = temp.z;
                    (*dataOutPtr)[m].procOrigIndexOut = temp.procOrigIndex;
                    (*dataOutPtr)[m].origIndexOut = temp.origIndex;
                    (*dataOutPtr)[m].TOut = temp.T;
                    (*dataOutPtr)[m].rOut = temp.r;
                    (*dataOutPtr)[m].thetaOut = temp.theta;
                    (*dataOutPtr)[m].TnewOut = temp.Tnew;

                    //clear temp structure
                    temp = {};
                }
            }
        }

        return true; //true if successful
    }
    catch (exception ex)
    {
        return false; //false if failed
    }
}

}

