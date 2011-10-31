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
#include "dPoly.h"
#include "dTree.h"
#include "geomUtils.h"

using namespace std;
using namespace utils;

int main(int argc, char** argv){

  if (argc <= 2){
    cout << "Usage: " << argv[0] << " poly.xg box.xg" << endl;
    exit(1);
  }
  
  char * polyFile   = argv[1];
  char * boxFile    = argv[2];
  bool isPointCloud = false;

  dPoly poly;
  if (! poly.readPoly(polyFile, isPointCloud) ) exit(1);

  dPoly box;
  if (! box.readPoly(boxFile, isPointCloud) ) exit(1);
  double xl, yl, xh, yh;
  box.bdBox(xl, yl, xh, yh);

  edgeTree T;
  T.putPolyEdgesInTree(poly);

  vector<seg> edgesInBox;
  T.findPolyEdgesInBox(xl, yl, xh, yh, // inputs
                       edgesInBox      // output
                       );
  
  const char * outFile = "edgesInBox.xg";
  cout << "Writing to " << outFile << endl;
  ofstream fs(outFile);
  fs << "color = red" << endl;
  for (int t = 0; t < (int)edgesInBox.size(); t++){
    const seg & S = edgesInBox[t];
    fs << S.begx  << ' ' << S.begy << endl;
    fs << S.endx  << ' ' << S.endy << endl;
    fs << "NEXT" << endl;
  }
  fs.close();

  return 0;
  
}
