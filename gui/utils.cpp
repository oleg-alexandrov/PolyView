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

#ifdef POLYVIEW_USE_OPENMP
#include <omp.h>
#endif
// For Visual studio
#ifdef _MSC_VER 
#define strcasecmp _stricmp
#endif

using namespace std;

void utils::printUsage(std::string progName){
  cout <<endl;
  cout << "USAGE: " <<endl<< progName
      <<" [app_options] [file_options] file1.xg [file_options] file2.xg"<<  endl<<endl;

  cout <<"app_options:"<<endl;
  cout <<"     -geo | -geometry 1200x800 " <<endl;
  cout <<"     -bg  | -backgroundColor black "<<endl;
  cout <<"     -grid (on or off) "<<endl;
  cout <<"     -gridSize 10 "<<endl;
  cout <<"     -gridWidth 1 "<<endl;
  cout <<"     -gridColor green "<<endl;
  cout <<"     -panelRatio 0.2  ([0-1] defines the ratio of the menu size to the display size)"<<endl;
#ifdef POLYVIEW_USE_OPENMP
  cout <<"     -nt  | -numThreads    number of threads to use for openmp loops"<<endl;
#endif
  cout <<endl;

  cout <<"file_options:"<<endl;
  cout <<"     -c   | -color            yellow "<<endl;
  cout <<"     -nc  | -noColorOverride     (use color in file) "<<endl;
  cout <<"     -lw  | -lineWidth        1  (thickness of the drawn lines)"<<endl;


  cout <<endl<<"   polygon options:"<<endl;

  cout <<"     -cp  | -closedPoly " <<endl;
  cout <<"     -ncp | -nonClosedPoly "<<endl;
  cout <<"     -f   | -filledPoly "<<endl;
  cout <<"     -nf  | -nonFilledPoly "<<endl;
  cout <<"     -cw  | -clockwisePoly    (if polygon orientation is clockwise)"<<endl;
  cout <<"     -tr  | -transparency [0-1] (transparency for filled polygons)"<<endl;

  cout <<endl<<"   point options:"<<endl;
  cout <<"     -p   | -points     (read as point cloud, not polygons)" <<endl;
  cout <<"     -sh  | -shape   o  ([x, +, o, s, t] defines shape of the points in point mode display"<<endl;
  cout <<"     -si  | -size    3  ([1-8]  defines size of the points in point mode display"<<endl;


  cout <<endl<<"   annotation options:"<<endl;
  cout <<"     -fs  | -fontSize         10        (for annotations)"<<endl;
  cout <<"     -ha  | -hideAnno                   (do not show annotations in file) "<<endl;
  cout <<"     -sa  | -scatterAnno                (plot annotation values as scattered points)"<<endl;
  cout <<"     -cs  | -colorScale min_val max_val (fixed color scale for scattered plot)"<<endl;

  cout <<endl<<"   image options:"<<endl;
  cout <<"     -cm  | -colorMap                    (Display grayscale image with built in colored map)"<<endl;
  cout <<"     -cs  | -colorScale  min_val max_val (color scale for the grayscale image, default 0->black 1->white )"<<endl;

  cout<<endl;

}

namespace {

double clamp(double v)
{
  const double t = v < 0 ? 0 : v;
  return t > 1.0 ? 1.0 : t;
}

}

void utils::getRGBColor(double v, double vmin, double vmax, double &r, double &g, double &b){
  if (vmax <= vmin){
    r = 0;
    b = 0;
    g = 1;
  } else {
    v = std::max(v, vmin);
    v = std::min(v, vmax);
    double t = -1 + 2*(v-vmin)/(vmax-vmin);
    r = clamp(1.5 - std::abs(2.0 * t - 1.0));
    g = clamp(1.5 - std::abs(2.0 * t));
    b = clamp(1.5 - std::abs(2.0 * t + 1.0));
  }
}


void utils::extractWindowDims(// inputs
                              int numArgs, char ** args,
                              // outputs
                              int & windowWidX, int & windowWidY){

  // Parse the command line arguments '-geo[metry] 500x600'

  windowWidX = 1200; windowWidY = 800; // defaults

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

    if (strcasecmp( currArg, "-h"     ) == 0 || strcasecmp( currArg, "--h"    ) == 0 ||
        strcasecmp( currArg, "-help"  ) == 0 || strcasecmp( currArg, "--help" ) == 0 ||
        strcasecmp( currArg, "-?"     ) == 0 || strcasecmp( currArg, "--?"    ) == 0 ){
      printUsage(exeName);
      exit(0);
    }

    if ( strcasecmp(currArg, "-p") == 0 || strcasecmp(currArg, "-points") == 0 ){
      opt.plotAsPoints = !opt.plotAsPoints;
      continue;
    }

    if ( strcasecmp(currArg, "-f") == 0 || strcasecmp(currArg, "-filledpoly") == 0 ){
      opt.isPolyFilled = true;
      continue;
    }

    if ( strcasecmp(currArg, "-nf") == 0 || strcasecmp(currArg, "-nonfilledpoly") == 0 ){
      opt.isPolyFilled = false;
      continue;
    }

    if ( strcasecmp(currArg, "-cw") == 0 || strcasecmp(currArg, "-clockwisepoly") == 0 ){
      opt.clockwisePoly = true;
      continue;
    }

    if ( strcasecmp(currArg, "-cp") == 0 || strcasecmp(currArg, "-closedpoly") == 0 ){
      // Plot as closed polygons
      opt.isPolyClosed = forceClosedPoly;
      continue;
    }

    if ( strcasecmp(currArg, "-ncp") == 0 || strcasecmp(currArg, "-nonclosedpoly") == 0 ){
      // Plot as polygonal lines
      opt.isPolyClosed = forceNonClosedPoly;
      continue;
    }

    if ((strcasecmp(currArg, "-bg") == 0 || strcasecmp(currArg, "-backgroundcolor") == 0) &&
        argIter < argc - 1) {
      opt.bgColor = argv[argIter + 1];
      argIter++;
      continue;
    }

    if ((strcasecmp(currArg, "-fs") == 0 || strcasecmp(currArg, "-fontsize") == 0 ) &&
        argIter < argc - 1) {
      int fs = (int)round(atof(argv[argIter + 1]));
      if (fs > 0) opt.fontSize = fs;
      argIter++;
      continue;
    }

    if ((strcasecmp(currArg, "-lw") == 0 || strcasecmp(currArg, "-linewidth") == 0) &&
         argIter < argc - 1) {
      int lw = (int)round(atof(argv[argIter + 1]));
      if (lw > 0) opt.lineWidth = lw;
      argIter++;
      continue;
    }

    if (strcasecmp(currArg, "-gridsize") == 0 && argIter < argc - 1){
      double gs = atof(argv[argIter + 1]);
      if (gs > 0) opt.gridSize = gs;
      argIter++;
      continue;
    }

    if ( strcasecmp(currArg, "-gridwidth") == 0  &&
         argIter < argc - 1
         ){
      int gw = (int)round(atof(argv[argIter + 1]));
      if (gw > 0) opt.gridWidth = gw;
      argIter++;
      continue;
    }

    if ( strcasecmp(currArg, "-gridcolor") == 0  &&
         argIter < argc - 1
         ){
      opt.gridColor = argv[argIter + 1];
      argIter++;
      continue;
    }

    if (strcasecmp(currArg, "-grid") == 0  &&
        argIter < argc - 1             &&
        strcasecmp(argv[argIter + 1], "on") == 0){
      opt.isGridOn = true;
      argIter++;
      continue;
    }

    if ((strcasecmp(currArg, "-c"    ) == 0 ||
         strcasecmp(currArg, "-color") == 0 )
        && argIter < argc - 1){
      opt.useCmdLineColor = true;
      opt.cmdLineColor    = argv[argIter + 1];
      argIter++;
      continue;
    }

    if ((strcasecmp(currArg, "-sh"    ) == 0 ||
         strcasecmp(currArg, "-shape") == 0 )
        && argIter < argc - 1){
       char shape = tolower(argv[argIter + 1][0]);
       if (shape == 'x'){
         opt.pointShape = 0;
       } else if (shape == '+'){
         opt.pointShape = 1;
       }else if (shape == 's'){
         opt.pointShape = 2;
       }else if (shape == 'o'){
         opt.pointShape = 3;
       }else if (shape == 't'){
         opt.pointShape = 4;
       } else {
         cout <<"ERROR un-supported shape specified, using triangles"<<endl;
         opt.pointShape = 5;
       }
      argIter++;
      continue;
    }

    if ((strcasecmp(currArg, "-si"    ) == 0 ||
        strcasecmp(currArg, "-size") == 0 )
        && argIter < argc - 1){
      opt.pointSize = (int)round(atof(argv[argIter + 1]));
      argIter++;
      continue;
    }

    if ((strcasecmp(currArg, "-tr"    ) == 0 ||
        strcasecmp(currArg, "-transparency") == 0 )
        && argIter < argc - 1){
      opt.transparency = atof(argv[argIter + 1]);
      argIter++;
      continue;
    }

    if ((strcasecmp(currArg, "-panelRatio") == 0) &&
        argIter < argc - 1) {
      opt.panelRatio = std::max(std::min(atof(argv[argIter + 1]), 0.9), 0.0);
      argIter++;
      continue;
    }
    
    if ((strcasecmp(currArg, "-nc") == 0 ||
         strcasecmp(currArg, "-nocoloroverride") == 0)){
      opt.useCmdLineColor = false;
      continue;
    }

    if ((strcasecmp(currArg, "-ha") == 0 ||
         strcasecmp(currArg, "-hideAnno") == 0)){
      opt.hideAnnotation = !opt.hideAnnotation;
      continue;
    }

    if ((strcasecmp(currArg, "-sa") == 0 ||
         strcasecmp(currArg, "-scatterAnno") == 0)){
      opt.scatter_annotations = !opt.scatter_annotations;
      continue;
    }

    if (strcasecmp(currArg, "-cm") == 0 || strcasecmp(currArg, "-colorMap") == 0){
      opt.useColorMap = !opt.useColorMap;
      continue;
    }

    if ((strcasecmp(currArg, "-cs") == 0 ||
        strcasecmp(currArg, "-colorScale") == 0) && argIter < argc - 2){

      opt.colorScale.resize(2);
      argIter++;
      opt.colorScale[0] = atof(argv[argIter]);
      argIter++;
      opt.colorScale[1] = atof(argv[argIter]);
      if (opt.colorScale[0] >= opt.colorScale[1]){
        cout <<"ERROR: min_val in color scale must be less than max_val"<<endl;
        opt.colorScale.clear();
      }

      continue;
    }

#ifdef POLYVIEW_USE_OPENMP

    if ((strcasecmp(currArg, "-nt") == 0 || strcasecmp(currArg, "-numThreads") == 0) &&
            argIter < argc - 1) {
         int nt = atoi(argv[argIter + 1]);
         if (nt > 0) {
           omp_set_num_threads(nt);
         }
         argIter++;
         continue;
       }
#endif

    // Other command line options are ignored
    if (currArg[0] == '-') continue;

    opt.polyFileName = currArg;

    options.polyOptionsVec.push_back(opt);
  }

  // Push one more time, to guarantee that the options vector is
  // non-empty even if no polygons were provided as input, and to make
  // sure we also parsed the options after the last polygon filename.
  // TODO(oalexan1): This is awkward, it better be stored separately
  // to start with.
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

// Get filename extension and make it lowercase
std::string utils::getFilenameExtension(std::string filename){

  std::string ext;
  std::string::size_type idx;
  idx = filename.rfind('.');

  if (idx != std::string::npos)
    ext = filename.substr(idx+1);

  std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
  
  return ext;
}

// Remove everything after the last dot and the dot itself
std::string utils::removeExtension(std::string filename){

  std::string ext;
  std::string::size_type idx;
  idx = filename.rfind('.');

  if (idx == std::string::npos)
    return filename;

  return filename.substr(0, idx);
}

bool utils::isImage(std::string const& filename) {
  string type = utils::getFilenameExtension(filename);

  return (type == "jpg" || type == "jpeg" || type == "png" || type == "tif" || type == "gif" ||
          type == "bmp" || type == "xpm");

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


bool utils::readImagePosition(std::string const& filename, std::vector<double> & pos) {

  std::ifstream fh(filename.c_str());
  if (!fh) {
    cerr << "Error: Could not open metadata file: " << filename << endl;
    return false;
  }

  pos.clear();
  double val = -1;
  while (fh >> val) {
    pos.push_back(val);
    if (pos.size() >= 4) {
      break;
    }
  }

  if (pos.size() < 4) {
    std::cerr << "Could not read four values from metadata file " << filename << ".\n";
    return false;
  }

  if (pos[2] == 0 || pos[3] == 0) {
    std::cerr << "Expecting last two values to be non-zero in " << filename << ".\n";
    return false;
  }
  
  return true;
}

// Convert from world coordinates to this image's pixel coordinates.
void utils::worldToImage(double wx, double wy, utils::PositionedImage const& img, // inputs
                  double & ix, double & iy) { // outputs
   // half grid shift for pixel center
   ix = (wx - img.pos[0] + img.pos[2]/2) / img.pos[2];
   iy = (wy - img.pos[1] + img.pos[3]/2) / img.pos[3];

   // Flip in y
   iy = img.qimg.height() - 1 - iy;
}

// The inverse of worldToImage()
void utils::imageToWorld(double ix, double iy, utils::PositionedImage const& img,
                         double & wx, double & wy) { // outputs
  
   // Flip in y
  iy = img.qimg.height() - 1 - iy;
  
  // half grid shift for pixel center
  wx = ix * img.pos[2] + img.pos[0] - img.pos[2]/2;
  wy = iy * img.pos[3] + img.pos[1] - img.pos[3]/2;
}

// Find the box containing all polygons and images
void utils::setUpViewBox(// inputs
                         const std::vector<dPoly> & polyVec,
                         // Ensure the produced box contains this box
                         double in_xll, double in_yll, double in_xur, double in_yur,
                         // outputs
                         double & xll,  double & yll,
                         double & widx, double & widy) {

  double xur, yur; // local variables

  bdBox(polyVec,             // inputs
        xll, yll, xur, yur); // outputs

  xll = std::min(xll, in_xll); yll = std::min(yll, in_yll);
  xur = std::max(xur, in_xur); yur = std::max(yur, in_yur);
  
  // Treat the case of empty polygons
  if (xur < xll || yur < yll) {
    xll = 0.0; yll = 0.0; xur = 1000.0; yur = 1000.0;
  }

  // Treat the case when the polygons are degenerate
  if (xur == xll) { xll -= 0.5; xur += 0.5; }
  if (yur == yll) { yll -= 0.5; yur += 0.5; }

  widx = xur - xll; assert(widx > 0.0);
  widy = yur - yll; assert(widy > 0.0);

  // Expand the box slightly for plotting purposes
  double factor = 0.05;
  xll -= widx*factor; xur += widx*factor; widx *= 1.0 + 2*factor;
  yll -= widy*factor; yur += widy*factor; widy *= 1.0 + 2*factor;

  return;
}

// The bounding box of an image in world coordinates
void utils::imageToWorldBdBox(// inputs
                              utils::PositionedImage const& positioned_img,
                              // outputs
                              double & xll, double & yll, double & xur, double & yur) {
  
  std::vector<int> x = {0, positioned_img.qimg.width(), positioned_img.qimg.width(), 0};
  std::vector<int> y = {0, 0, positioned_img.qimg.height(), positioned_img.qimg.height()};
  
  double big = std::numeric_limits<double>::max();
  xll = big; yll = big; xur = -big; yur = -big;
  for (size_t c = 0; c < x.size(); c++) {
    double wx, wy;
    utils::imageToWorld(x[c], y[c], positioned_img, wx, wy);
    
    xll = std::min(xll, wx); yll = std::min(yll, wy);
    xur = std::max(xur, wx); yur = std::max(yur, wy);
  }

  return;
}

// The bounding box of a sequence of images in world coordinates
void utils::imageToWorldBdBox(// inputs
                              const std::vector<dPoly> & polyVec,
                              // outputs
                              double & image_xll, double & image_yll,
                              double & image_xur, double & image_yur) {

  double big = DBL_MAX;
  image_xll = big; image_yll = big; image_xur = -big; image_yur = -big;
  
  for (int vi = 0; vi < (int)polyVec.size(); vi++) {
    
    if (polyVec[vi].img == NULL)
      continue;
    
    // Recover the image. This may crash if it was not populated correctly.
    utils::PositionedImage const & positioned_img
      = *(utils::PositionedImage*)(polyVec[vi].img);
    
    double xll, yll, xur, yur;
    utils::imageToWorldBdBox(// inputs
                             positioned_img,  
                             // outputs
                             xll, yll, xur, yur);
    
    image_xll = std::min(image_xll, xll); image_yll = std::min(image_yll, yll);
    image_xur = std::max(image_xur, xur); image_yur = std::max(image_yur, yur);
  }

  return;
}

