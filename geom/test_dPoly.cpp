// MIT License Terms (http://en.wikipedia.org/wiki/MIT_License)
// 
// Copyright (C) 2011 by Oleg Alexandrov
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <limits>
#include <algorithm>
#include "dPoly.h"

using namespace std;
using namespace utils;

int main(int argc, char** argv){

  dPoly poly;
  if (argc <= 1){
    cout << "Usage: " << argv[0] << " file.xg" << endl;
    exit(1);
  }
  
  char * filename = argv[1];
  cout << "Reading " << filename << endl;
  
  bool isPointCloud = false;

  if (! poly.readPoly(filename, isPointCloud) ) exit(1);

  double * xv   = (double*)poly.get_xv(); 
  double * yv   = (double*)poly.get_yv();
  int      numV = poly.get_totalNumVerts();

  double xmin = *min_element(xv, xv + numV) - 10;
  double xmax = *max_element(xv, xv + numV) + 10;
  double ymin = *min_element(yv, yv + numV) - 10;
  double ymax = *max_element(yv, yv + numV) + 10;

  char outFile[] = "pointsInside.xg";
  cout << "Writing " << outFile << endl;
  ofstream of(outFile);
  for (int i = (int)xmin; i <= (int)xmax; i++){
    for (int j = (int)ymin; j <= (int)ymax ; j++){
      if (isPointInPolyOrOnEdges(i, j, numV, xv, yv)){
        of << i  << ' ' << j << endl;
      }
    }
  }
  
#if 0
  poly.sortFromLargestToSmallest();
  double * xv   = (double*)poly.get_xv(); // to do: fix this hack
  double * yv   = (double*)poly.get_yv();
  int      numV = poly.get_totalNumVerts();

  bool isClosedPolyLine = true;
  utils::snapPolyLineTo45DegAngles(isClosedPolyLine, numV, xv, yv);
  
  const char * outFile = "out.xg";
  cout << "Writing to " << outFile << endl;
  poly.writePoly(outFile);

#endif
  return 0;
  
}
