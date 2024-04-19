//dataSort.h
#include <cmath>
#include <algorithm>
#include <vector>
#include <iterator>
#include "bndStructs.h"
using namespace std;

void swap(dataOutput* iData, dataOutput* jData);

bool dataSort(int numOfRows, vector<dataOutput>* dataOut);
