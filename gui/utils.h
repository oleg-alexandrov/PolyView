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
#ifndef UTILS_H
#define UTILS_H
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <geom/polyUtils.h>

enum closedPolyInfo{
  // If an array of points as read from file has the first vertex equal to the last
  // one, we treat it is a closed polygon (last vertex joins the first vertex).
  // If the user wants to override this behavior, the first two fields below
  // become necessary.
  forceClosedPoly, forceNonClosedPoly, readClosedPolyInfoFromFile
};


struct polyOptions{
  // Each polygon file has these options
  bool            plotAsPoints;
  bool            isPolyFilled;
  bool            clockwisePoly;
  closedPolyInfo  isPolyClosed;
  bool            useCmdLineColor;
  int             fontSize;
  int             lineWidth;
  bool            isGridOn;
  double          gridSize;
  int             gridWidth;
  bool            readPolyFromDisk;
  double          panelRatio;
  std::string     bgColor;
  std::string     fgColor;
  std::string     cmdLineColor;
  std::string     gridColor;
  std::string     polyFileName;

  polyOptions(){
    plotAsPoints     = false;
    isPolyFilled     = false;
    clockwisePoly    = false;
    isPolyClosed     = readClosedPolyInfoFromFile;
    fontSize         = 10;
    lineWidth        = 1;
    useCmdLineColor  = false;
    isGridOn         = false;
    gridWidth        = 1;
    gridSize         = -1;
    readPolyFromDisk = true;
    panelRatio       = 0.2;
    bgColor          = "black";
    fgColor          = "white";
    cmdLineColor     = "green";
    gridColor        = "white";
    polyFileName     = "unnamed.xg";
  }

};

struct cmdLineOptions{
  std::vector<polyOptions> polyOptionsVec;
};

namespace utils{

  std::string getDocText();

  void extractWindowDims(// inputs
                         int numArgs, char ** args,
                         // outputs
                         int & windowWidX, int & windowWidY
                         );

  void parseCmdOptions(//inputs
                       int argc, char** argv, std::string exeName,
                       // outputs
                       int & windowWidX, int & windowWidY, cmdLineOptions & options
                       );

  std::string inFileToOutFile(const std::string & inFile);

  void printUsage(std::string progName);

  std::string getFilenameExtension(std::string filename);
  std::string replaceAll(std::string result,
                         const std::string & replaceWhat,
                         const std::string & replaceWithWhat);

}

#endif
