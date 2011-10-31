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
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <limits>
#include <fstream>
#include "dPoly.h"
#include "cutPoly.h"

using namespace std;
using namespace utils;

int main(int argc, char** argv){

  if (argc <= 1){
    cout << "Usage: " << argv[0] << " file.xg" << endl;
    exit(1);
  }
  
  char * filename = argv[1];
  bool isPointCloud = false;
  cout << "Reading " << filename << endl;

  dPoly poly;
  if (! poly.readPoly(filename, isPointCloud) ) exit(1);

  const double * xv       = poly.get_xv();
  const double * yv       = poly.get_yv();
  const int    * numVerts = poly.get_numVerts();
  int numPolys            = poly.get_numPolys();
  //int totalNumVerts     = poly.get_totalNumVerts();

//   double xll = 100, xur = 426;
//   double yll = 100, yur = 501.214;

  double w[] = {373.255, -54.6423, 1086, 4837};
  
  vector<double> cutX, cutY;
  vector<int> cutNumPolys;
  
  cutPoly(numPolys, numVerts, xv, yv, w[0], w[1], w[2], w[3], // inputs
          cutX, cutY, cutNumPolys                             // outputs
          );


  const char * outFile = "clipped.xg";
  cout << "Writing to " << outFile << endl;
  //double scale = 1.0;
  vector<string> colorsOut, layersOut; 
  std::vector<anno>  annotations;
  // writePoly(outFile, isPointCloud,
  // cutX, cutY, cutNumPolys, cutNumPolys.size(), cutX.size(),
  // colorsOut, "yellow", scale, layersOut, annotations);

  ofstream win("window.xg");
  cout << "Writing to window.xg" << endl;
  win << "color = cyan" << endl;
  win << w[0] << ' ' << w[1] << endl;
  win << w[2] << ' ' << w[1] << endl;
  win << w[2] << ' ' << w[3] << endl;
  win << w[0] << ' ' << w[3] << endl;
  win.close();
  return 0;
  
}
