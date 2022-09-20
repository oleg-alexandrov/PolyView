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
#include <limits>
#include <cstring>
#include <algorithm>
#include <gui/utils.h>

using namespace std;

void utils::printUsage(std::string progName){
  cout << "Usage: " << progName << " "
       << "[ -geo[metry] 1000x800 ] [-bg | -backgroundColor black ] "
       << "[ -c | -color yellow ] "
       << "[ -nc | -noColorOverride ] "
       << "[ -fs | -fontSize 10 ] "
       << "[ -lw | -lineWidth 2 ] "
       << "[ -p | -points ] [ -cp | -closedPoly ] [ -ncp | -nonClosedPoly ] "
       << "[ -f | -filledPoly ] [ -nf | -nonFilledPoly ] [ -cw | clockwisePoly ]"
       << "[ -grid on | off ] [ -gridSize 10 ] [ -gridWidth 1 ] "
       << "[ -gridColor green ] "
       << "file_1.xg ... file_N.xg " << endl;
}

void utils::extractWindowDims(// inputs
                              int numArgs, char ** args,
                              // outputs
                              int & windowWidX, int & windowWidY
                              ){

  // Parse the command line arguments '-geo[metry] 500x600'

  windowWidX = 900; windowWidY = 700; // defaults

  for (int s = 1; s < numArgs; s++){

    if ( !strstr(args[s-1], "-geo") ) continue;

    string lineStr = args[s];
    char * line    = (char*) lineStr.c_str();

    // Blank the geometry settings once located
    // to not confuse other parsers.
    args[s-1][0] = '\0';
    args[s  ][0] = '\0';

    char * pch;
    char delimiter[] = "x";

    pch = strtok (line, delimiter);
    if (pch == NULL) continue;
    int windowWidX_tmp = atoi(pch);

    pch = strtok (NULL, delimiter);
    if (pch == NULL) continue;
    int windowWidY_tmp = atoi(pch);

    if (windowWidX_tmp > 0 && windowWidY_tmp > 0){
      windowWidX = windowWidX_tmp;
      windowWidY = windowWidY_tmp;
    }

  }

}

void utils::parseCmdOptions(//inputs
                            int argc, char** argv, std::string exeName,
                            // outputs
                            int & windowWidX, int & windowWidY, cmdLineOptions & options){

  options.polyOptionsVec.clear();

  polyOptions opt; // Each polygon file will have one such entry

  // Skip argv[0] as that's the program name
  extractWindowDims(argc - 1, argv + 1, windowWidX, windowWidY);

  for (int argIter = 1; argIter < argc; argIter++){

    char * currArg = argv[argIter];

    if (currArg == NULL || strlen(currArg) == 0) continue;

    if (currArg[0] == '-'){
      // Transform -P into -p, etc.
      transform(currArg, currArg + strlen(currArg), currArg, ::tolower);
    }

    if (strcmp( currArg, "-h"     ) == 0 || strcmp( currArg, "--h"    ) == 0 ||
        strcmp( currArg, "-help"  ) == 0 || strcmp( currArg, "--help" ) == 0 ||
        strcmp( currArg, "-?"     ) == 0 || strcmp( currArg, "--?"    ) == 0 ){
      printUsage(exeName);
      exit(0);
    }

    if ( strcmp(currArg, "-p") == 0 || strcmp(currArg, "-points") == 0 ){
      opt.plotAsPoints = !opt.plotAsPoints;
      continue;
    }

    if ( strcmp(currArg, "-f") == 0 || strcmp(currArg, "-filledpoly") == 0 ){
      opt.isPolyFilled = true;
      continue;
    }

    if ( strcmp(currArg, "-nf") == 0 || strcmp(currArg, "-nonfilledpoly") == 0 ){
      opt.isPolyFilled = false;
      continue;
    }

    if ( strcmp(currArg, "-cw") == 0 || strcmp(currArg, "-clockwisepoly") == 0 ){
      opt.clockwisePoly = true;
      continue;
    }

    if ( strcmp(currArg, "-cp") == 0 || strcmp(currArg, "-closedpoly") == 0 ){
      // Plot as closed polygons
      opt.isPolyClosed = forceClosedPoly;
      continue;
    }

    if ( strcmp(currArg, "-ncp") == 0 || strcmp(currArg, "-nonclosedpoly") == 0 ){
      // Plot as polygonal lines
      opt.isPolyClosed = forceNonClosedPoly;
      continue;
    }

    if ( (strcmp(currArg, "-bg") == 0 || strcmp(currArg, "-backgroundcolor") == 0 )
         &&
         argIter < argc - 1
         ){
      opt.bgColor = argv[argIter + 1];
      argIter++;
      continue;
    }

    if ( (strcmp(currArg, "-fs"      ) == 0 ||
          strcmp(currArg, "-fontsize") == 0 )
         &&
         argIter < argc - 1
         ){
      int fs = (int)round(atof(argv[argIter + 1]));
      if (fs > 0) opt.fontSize = fs;
      argIter++;
      continue;
    }

    if ( (strcmp(currArg, "-lw"       ) == 0 ||
          strcmp(currArg, "-linewidth") == 0 )
         &&
         argIter < argc - 1
         ){
      int lw = (int)round(atof(argv[argIter + 1]));
      if (lw > 0) opt.lineWidth = lw;
      argIter++;
      continue;
    }

    if ( strcmp(currArg, "-gridsize") == 0 &&
         argIter < argc - 1
         ){
      double gs = atof(argv[argIter + 1]);
      if (gs > 0) opt.gridSize = gs;
      argIter++;
      continue;
    }

    if ( strcmp(currArg, "-gridwidth") == 0  &&
         argIter < argc - 1
         ){
      int gw = (int)round(atof(argv[argIter + 1]));
      if (gw > 0) opt.gridWidth = gw;
      argIter++;
      continue;
    }

    if ( strcmp(currArg, "-gridcolor") == 0  &&
         argIter < argc - 1
         ){
      opt.gridColor = argv[argIter + 1];
      argIter++;
      continue;
    }

    if ( strcmp(currArg, "-grid") == 0  &&
         argIter < argc - 1             &&
         strcmp(argv[argIter + 1], "on") == 0
         ){
      opt.isGridOn = true;
      argIter++;
      continue;
    }

    if ( (strcmp(currArg, "-c"    ) == 0 ||
          strcmp(currArg, "-color") == 0 )
         && argIter < argc - 1){
      opt.useCmdLineColor = true;
      opt.cmdLineColor    = argv[argIter + 1];
      argIter++;
      continue;
    }

    if ((strcmp(currArg, "-nc") == 0 ||
         strcmp(currArg, "-nocoloroverride") == 0)){
      opt.useCmdLineColor = false;
      continue;
    }

    // Other command line options are ignored
    if (currArg[0] == '-') continue;

    opt.polyFileName = currArg;

    options.polyOptionsVec.push_back(opt);
  }

  // Push one more time, to guarantee that the options vector is
  // non-empty even if no polygons were provided as input, and to make
  // sure we also parsed the options after the last polygon filename.
  options.polyOptionsVec.push_back(opt);

  return;
}

std::string utils::inFileToOutFile(const std::string & inFile){

  string outFile = "";

  bool lastDot = true;
  for (int s = (int)inFile.length() - 1; s >= 0; s--){

    string currChar = inFile.substr(s, 1);
    if (currChar == "." && lastDot){
      outFile = string("_out") + currChar + outFile;
      lastDot = false;
    }else{
      outFile = currChar + outFile;
    }

  }

  if (outFile.length() == 0){
    cerr << "Invalid filename" << endl;
  }

  return outFile;

}

std::string utils::getFilenameExtension(std::string filename){

  std::string::size_type idx;
  idx = filename.rfind('.');

  if(idx != std::string::npos) return filename.substr(idx+1);
  else                         return "";
}

std::string utils::replaceAll(std::string result,
                              const std::string & replaceWhat,
                              const std::string & replaceWithWhat){

  while(1){
    const int pos = result.find(replaceWhat);
    if (pos == -1) break;
    result.replace(pos,replaceWhat.size(),replaceWithWhat);
  }
  return result;
}
