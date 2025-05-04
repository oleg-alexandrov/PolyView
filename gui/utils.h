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
#include <QImage>

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
  bool            hideAnnotation;
  closedPolyInfo  isPolyClosed;
  bool            useCmdLineColor;
  int             fontSize;
  int             pointShape;
  int             pointSize;
  int             lineWidth;
  bool            isGridOn;
  double          gridSize;
  int             gridWidth;
  bool            readPolyFromDisk;
  double          panelRatio;
  bool            scatter_annotations;
  bool            useColorMap;
  double          transparency;
  std::string     bgColor;
  std::string     fgColor;
  std::string     cmdLineColor;
  std::string     gridColor;
  std::string     markColor;
  std::string     polyFileName;
  std::vector<double> colorScale;

  polyOptions(){
    plotAsPoints     = false;
    isPolyFilled     = false;
    clockwisePoly    = false;
    hideAnnotation   = false;
    isPolyClosed     = readClosedPolyInfoFromFile;
    fontSize         = 10;
    pointShape       = -1;
    pointSize        = 0;
    lineWidth        = 1;
    useCmdLineColor  = false;
    isGridOn         = false;
    gridWidth        = 1;
    gridSize         = -1;
    readPolyFromDisk = true;
    panelRatio       = 0.2;
    scatter_annotations = false;
    useColorMap        = false;
    transparency      = 0;
    bgColor          = "black";
    fgColor          = "white";
    cmdLineColor     = "green";
    gridColor        = "white";
    markColor        = "magenta";
    polyFileName     = "unnamed.xg";
  }

};

struct cmdLineOptions{
  std::vector<polyOptions> polyOptionsVec;
};

namespace utils{

  // An image buffer and its world position. The vector 'pos'
  // has world coordinates of the center of the image lower-left
  // pixel, and the vector from the lower-left image pixel to the next
  // pixel along the 45 degree diagonal, so that pixel's dimensions.
  struct PositionedImage {
    QImage qimg;
    std::vector<double> pos;
  };
  
  std::string getDocText();
  void getRGBColor(double v, double vmin, double vmax, double &r, double &g, double &b);

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

  // Get filename extension and make it lowercase
  std::string getFilenameExtension(std::string filename);

  // Remove everything after the last dot and the dot itself
  std::string removeExtension(std::string filename);

  bool isImage(std::string const& filename);

  std::string replaceAll(std::string result,
                         const std::string & replaceWhat,
                         const std::string & replaceWithWhat);

  // Read a file with four numbers determining how the image will be plotted.
  bool readImagePosition(std::string const& filename, std::vector<double> & pos);

  // Convert from world coordinates to this image's pixel coordinates.
  void worldToImage(double wx, double wy, PositionedImage const& img, // inputs
                    double & ix, double & iy); // outputs

  // The inverse of worldToImage()
  void imageToWorld(double ix, double iy, PositionedImage const& img,
                    double & wx, double & wy);  // outputs

  // Find the box containing all polygons and images
  void setUpViewBox(// inputs
                    const std::vector<dPoly> & polyVec,
                    // Ensure the produced box contains this box
                    double in_xll, double in_yll, double in_xur, double in_yur,
                    // outputs
                    double & xll,  double & yll,
                    double & widx, double & widy);

  // The bounding box of an image in world coordinates
  void imageToWorldBdBox(// inputs
                         PositionedImage const& positioned_img,
                         // outputs
                         double & xll, double & yll, double & xur, double & yur);

  // The bounding box of a sequence of images in world coordinates
  void imageToWorldBdBox(// inputs
                         const std::vector<dPoly> & polyVec,
                         // outputs
                         double & xll, double & yll, double & xur, double & yur);
  
}

#endif
