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
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <cassert>
#include <algorithm>
#include "edgeUtils.h"
#include "dPoly.h"
#include "polyUtils.h"
#include "geomUtils.h"
#include "kdTree.h"
#include "dTree.h"

using namespace std;
using namespace utils;

void utils::findClosestPolyVertex(// inputs
                                  double x0, double y0,
                                  const std::vector<dPoly> & polyVec,
                                  // outputs
                                  int & polyVecIndex,
                                  int & polyIndexInCurrPoly,
                                  int & vertIndexInCurrPoly,
                                  double & minX, double & minY,
                                  double & minDist
                                  ){

  // Find the closest point in a given vector of polygons to a given point.
  
  polyVecIndex = -1; polyIndexInCurrPoly = -1; vertIndexInCurrPoly = -1;
  minX = x0; minY = y0; minDist = DBL_MAX;
  
  for (int s = 0; s < (int)polyVec.size(); s++){

    double minX0, minY0, minDist0;
    int polyIndex, vertIndex;
    polyVec[s].findClosestPolyVertex(// inputs
                                     x0, y0,    
                                     // outputs
                                     polyIndex, vertIndex, minX0, minY0, minDist0
                                     );

    if (minDist0 <= minDist){
      polyVecIndex  = s;
      polyIndexInCurrPoly = polyIndex;
      vertIndexInCurrPoly = vertIndex;
      minDist       = minDist0;
      minX          = minX0;
      minY          = minY0;
    }
    
  }

  return;
}

void utils::findClosestPolyEdge(// inputs
                                double x0, double y0,
                                const std::vector<dPoly> & polyVec,
                                // outputs
                                int & vecIndex, int & polyIndex,
                                int & vertIndex,
                                double & minX, double & minY, double & minDist
                                ){

  // Find the closest edge in a given vector of polygons to a given point.
  
  vecIndex  = -1;
  polyIndex = -1;
  vertIndex = -1;
  minX      = DBL_MAX;
  minY      = DBL_MAX;
  minDist   = DBL_MAX;
  
  for (int vecIter = 0; vecIter < (int)polyVec.size(); vecIter++){

    double lx, ly, ldist;
    int pIndex, vIndex;
    polyVec[vecIter].findClosestPolyEdge(x0, y0,                        // in
                                         pIndex, vIndex, lx, ly, ldist  // out
                                         );

    if (ldist <= minDist){
      vecIndex  = vecIter;
      polyIndex = pIndex;
      vertIndex = vIndex;      
      minX      = lx;
      minY      = ly;
      minDist   = ldist;
    }
    
  }

  return;
}

void utils::alignPoly1ToPoly2(dPoly       & poly1,
                              const dPoly & poly2,
                              utils::linTrans & T // save the applied transform
                              ){

  // Find the closest pair of vertices from poly1 to poly2. Shift poly1
  // so that the first vertex is on top of the second one.
  vector<segDist> distVec;

  findDistanceFromVertsOfPoly1ToVertsPoly2(// inputs
                                           poly1, poly2,  
                                           // outputs
                                           distVec
                                           );
  int len = distVec.size();
  if (len == 0) return; // one of the polygons is empty
    
  segDist S = distVec[len -1]; // this corresponds to shortest distance
  poly1.applyTransform(1, 0, 0, 1, S.endx - S.begx, S.endy - S.begy, T);

  return;
}

void utils::findDistanceFromVertsOfPoly1ToVertsPoly2(// inputs
                                                     const dPoly & poly1,
                                                     const dPoly & poly2,
                                                     // outputs
                                                     std::vector<segDist> & distVec
                                                     ){
  
  // Given two sets of polygons, for each vertex in the first set of
  // polygons find the distance to the closest vertex in the second
  // set of polygons, and the segment with the smallest distance. Sort
  // these segments in decreasing value of their lengths.

  // The complexity of this algorithm is roughly
  // size(poly1)*log(size(poly2)).
  
  distVec.clear();

  const double * x1 = poly1.get_xv();
  const double * y1 = poly1.get_yv();
  int numVerts1     = poly1.get_totalNumVerts();
  const double * x2 = poly2.get_xv();
  const double * y2 = poly2.get_yv();
  int numVerts2     = poly2.get_totalNumVerts();

  if (numVerts1 == 0 || numVerts2 == 0) return; // no vertices

  // Put the edges of the second polygon in a tree for fast access
  kdTree T;
  T.formTreeOfPoints(numVerts2, x2, y2); 

  for (int t = 0; t < numVerts1; t++){

    double x = x1[t], y = y1[t];
    double closestDist;
    PointWithId closestPt;
    T.findClosestVertexToPoint(// inputs
                               x, y,  
                               // outputs
                               closestPt, closestDist
                               );
    distVec.push_back(segDist(x, y, closestPt.x, closestPt.y, closestDist));
    
  }

  sort(distVec.begin(), distVec.end(), segDistGreaterThan);

  return;
}

void utils::findDistanceBwPolys(// inputs
                                const dPoly & poly1,
                                const dPoly & poly2,
                                // outputs
                                std::vector<segDist> & distVec
                                ){

  // Find the distances from poly1 to poly2, then from poly2 to poly1.
  // Sort them in decreasing order of their lengths. See
  // findDistanceFromPoly1ToPoly2 for more info.
  
  findDistanceFromPoly1ToPoly2(poly1, poly2, // inputs  
                               distVec       // outputs
                               );
  
  vector<segDist> l_distVec;
  findDistanceFromPoly1ToPoly2(poly2, poly1, // inputs  
                               l_distVec     // outputs
                               );
  
  for (int s = 0; s < (int)l_distVec.size(); s++) distVec.push_back(l_distVec[s]);
  
  sort(distVec.begin(), distVec.end(), segDistGreaterThan);

  return;
}

void utils::findDistanceFromPoly1ToPoly2(// inputs
                                         const dPoly & poly1,
                                         const dPoly & poly2,
                                         // outputs
                                         std::vector<segDist> & distVec
                                         ){
  
  // Given two sets of polygons, for each vertex in the first set of
  // polygons find the distance to the closest point (may be on edge)
  // in the second set of polygons, and the segment with the smallest
  // distance. Sort these segments in decreasing value of their
  // lengths.

  // The complexity of this algorithm is roughly
  // size(poly1)*log(size(poly2)).
  
  distVec.clear();

  const double * x1 = poly1.get_xv();
  const double * y1 = poly1.get_yv();
  int numVerts1     = poly1.get_totalNumVerts();
  int numVerts2     = poly2.get_totalNumVerts();

  if (numVerts1 == 0 || numVerts2 == 0) return; // no vertices

  // Put the edges of the second polygon in a tree for fast access
  edgeTree T;
  T.putPolyEdgesInTree(poly2);

  for (int t = 0; t < numVerts1; t++){

    double x = x1[t], y = y1[t];
    double closestX, closestY, closestDist;
    seg closestEdge;
    T.findClosestEdgeToPoint(x, y,                                        // inputs 
                             closestEdge, closestDist, closestX, closestY // outputs
                             );
    distVec.push_back(segDist(x, y, closestX, closestY, closestDist));
    
  }

  sort(distVec.begin(), distVec.end(), segDistGreaterThan);

  return;
}

void utils::findDistanceBwPolysBruteForce(// inputs
                                          const dPoly & poly1,
                                          const dPoly & poly2,
                                          // outputs
                                          std::vector<segDist> & distVec
                                          ){

  // A naive (but simple) implementation of findDistanceFromPoly1ToPoly2.
  
  distVec.clear();

  const double * x = poly1.get_xv();
  const double * y = poly1.get_yv();
  int numVerts1    = poly1.get_totalNumVerts();
  int numVerts2    = poly2.get_totalNumVerts();
  
  if (numVerts1 == 0 || numVerts2 == 0) return; // no vertices

  for (int t  = 0; t < numVerts1; t++){
    
    int minPolyIndex, minVertIndex;
    double minX, minY, minDist = DBL_MAX;
    poly2.findClosestPolyEdge(x[t], y[t],                                      // inputs
                              minPolyIndex, minVertIndex, minX, minY,  minDist // outputs
                              );

    
    if (minDist != DBL_MAX)
      distVec.push_back(segDist(x[t], y[t], minX, minY, minDist));
  }

  sort(distVec.begin(), distVec.end(), segDistGreaterThan);
  
  return;
}

void utils::putPolyInMultiSet(const dPoly & P, std::multiset<dPoint> & mP){

  const double * x  = P.get_xv();
  const double * y  = P.get_yv();
  int totalNumVerts = P.get_totalNumVerts();

  mP.clear();
  for (int v = 0; v < totalNumVerts; v++){
    dPoint P;
    P.x = x[v];
    P.y = y[v];
    mP.insert(P);
  }

  return;
}

void utils::findPolyDiff(const dPoly & P, const dPoly & Q, // inputs
                         std::vector<dPoint> & vP, std::vector<dPoint> & vQ // outputs
                         ){
    
  // Compare two polygons point-by-point. We assume that the polygons
  // may have collinear points. If one polygon has a point repeated
  // twice, but the second polygon has it repeated just once, this
  // will be flagged as a difference as well.

  // This utility will not be able to detect when two polygons are
  // different but contain exactly the same points.
  
  multiset<dPoint> mP; putPolyInMultiSet(P, mP);
  multiset<dPoint> mQ; putPolyInMultiSet(Q, mQ);

  // If a point is in mP, and also in mQ, mark it as being in mP and wipe it from mQ
  vector<dPoint> shared;
  shared.clear();
  multiset<dPoint>::iterator ip, iq;
  for (ip = mP.begin(); ip != mP.end(); ip++){
    iq = mQ.find(*ip);
    if ( iq != mQ.end() ){
      shared.push_back(*ip);
      mQ.erase(iq); // Erase just the current instance of the given value
    }
  }
  
  // Wipe it from mP as well
  for (int s = 0; s < (int)shared.size(); s++){
    ip = mP.find(shared[s]);
    if ( ip != mP.end() ){
      mP.erase(ip); // Erase just the current instance of the given value
    }
  }

  vP.clear(); vQ.clear();
  for (ip = mP.begin(); ip != mP.end(); ip++){
    dPoint p;
    p.x = ip->x;
    p.y = ip->y;
    vP.push_back(p);
  }
  for (iq = mQ.begin(); iq != mQ.end(); iq++){
    dPoint q;
    q.x = iq->x;
    q.y = iq->y;
    vQ.push_back(q);
  }
  
  return;
}

void utils::bdBox(const std::vector<dPoly> & polyVec,
                  // outputs
                  double & xll, double & yll,
                  double & xur, double & yur
                  ){

  double big = DBL_MAX;
  xll = big; yll = big; xur = -big; yur = -big;
  for (int p = 0; p < (int)polyVec.size(); p++){
    if (polyVec[p].get_totalNumVerts() == 0) continue;
    double xll0, yll0, xur0, yur0;
    polyVec[p].bdBox(xll0, yll0, xur0, yur0);
    xll = min(xll, xll0); xur = max(xur, xur0);
    yll = min(yll, yll0); yur = max(yur, yur0);
  }
  
  return;
}
  

void utils::setUpViewBox(// inputs
                         const std::vector<dPoly> & polyVec,
                         // outputs
                         double & xll,  double & yll,
                         double & widx, double & widy
                         ){
  
  // Given a set of polygons, set up a box containing these polygons.

  double xur, yur; // local variables
  
  bdBox(polyVec,           // inputs 
        xll, yll, xur, yur // outputs
        );
  
  // Treat the case of empty polygons
  if (xur < xll || yur < yll){
    xll = 0.0; yll = 0.0; xur = 1000.0; yur = 1000.0;
  }

  // Treat the case when the polygons are degenerate
  if (xur == xll){ xll -= 0.5; xur += 0.5; }
  if (yur == yll){ yll -= 0.5; yur += 0.5; }
    
  widx = xur - xll; assert(widx > 0.0);
  widy = yur - yll; assert(widy > 0.0);

  // Expand the box slightly for plotting purposes
  double factor = 0.05;
  xll -= widx*factor; xur += widx*factor; widx *= 1.0 + 2*factor;
  yll -= widy*factor; yur += widy*factor; widy *= 1.0 + 2*factor;
  
  return;
  
}


void utils::markPolysInHlts(// Inputs
                            const std::vector<dPoly> & polyVec,
                            const std::vector<dPoly> & highlights,
                            // Outputs
                            std::map< int, std::map<int, int> > & markedPolyIndices
                            ){
  
  markedPolyIndices.clear();

  for (int s = 0; s < (int)highlights.size(); s++){

    double xll, yll, xur, yur;
    assert(highlights[s].get_totalNumVerts() == 4);
    highlights[s].bdBox(xll, yll, xur, yur);

    map<int, int> mark;
    for (int t = 0; t < (int)polyVec.size(); t++){
      polyVec[t].markPolysIntersectingBox(xll, yll, xur, yur, // Inputs 
                                          mark                // Outputs
                                          );
      for (map<int, int>::iterator it = mark.begin(); it != mark.end(); it++){
        markedPolyIndices[t][it->first] = it->second;
      }
    }
  }
  
  return;
}

void utils::shiftMarkedPolys(// Inputs
                             std::map< int, std::map<int, int> > & markedPolyIndices,
                             double shift_x, double shift_y,
                             // Inputs-outputs
                             std::vector<dPoly> & polyVec
                             ){

  for (int pIter = 0; pIter < (int)polyVec.size(); pIter++){
    polyVec[pIter].shiftMarkedPolys(markedPolyIndices[pIter], shift_x, shift_y);
  }

  return;
}

void utils::scaleMarkedPolysAroundCtr(// Inputs
                                       std::map< int, std::map<int, int> > & markedPolyIndices,
                                       double scale,
                                       // Inputs-outputs
                                       std::vector<dPoly> & polyVec
                                       ){

  matrix2 M;
  M.a11 = scale; M.a12 = 0.0; M.a21 = 0.0; M.a22 = scale;
  transformMarkedPolysAroundCtr(// Inputs
                                markedPolyIndices, M,  
                                // Inputs-outputs
                                polyVec
                                );
  
  return;  
}

void utils::rotateMarkedPolysAroundCtr(// Inputs
                                       std::map< int, std::map<int, int> > & markedPolyIndices,
                                       double angle,
                                       // Inputs-outputs
                                       std::vector<dPoly> & polyVec
                                       ){

  double a = angle*M_PI/180.0, c = cos(a), s= sin(a);
  
  if (angle == round(angle) && int(angle)%90 == 0 ){
    // The special case of angle multiple of 90 degrees
    c = round(c), s = round(s);
  }

  matrix2 M;
  M.a11 = c; M.a12 = -s; M.a21 = s; M.a22 = c;
  transformMarkedPolysAroundCtr(// Inputs
                                markedPolyIndices, M,  
                                // Inputs-outputs
                                polyVec
                                );
  return;  
}

void utils::transformMarkedPolysAroundCtr(// Inputs
                                         std::map< int, std::map<int, int> > & markedPolyIndices,
                                         const utils::matrix2 & M,
                                         // Inputs-outputs
                                         std::vector<dPoly> & polyVec
                                         ){

  if (getNumElements(markedPolyIndices) == 0) return;

  vector<dPoly> extractedPolyVec;
  extractMarkedPolys(// Inputs
                     polyVec, markedPolyIndices,  
                     // Outputs
                     extractedPolyVec
                     );
  
  // Find the center of the bounding box of the marked polygons
  double xll, yll, xur, yur;
  bdBox(extractedPolyVec,   // inputs 
        xll, yll, xur, yur  // outputs
        );
  dPoint P;
  P.x = (xll + xur)/2.0;
  P.y = (yll + yur)/2.0;
    
  for (int pIter = 0; pIter < (int)polyVec.size(); pIter++){
    polyVec[pIter].transformMarkedPolysAroundPt(markedPolyIndices[pIter], M, P);
  }

  return;
}

void utils::eraseMarkedPolys(// Inputs
                               std::map< int, std::map<int, int> > & markedPolyIndices,
                               // Inputs-outputs
                               std::vector<dPoly> & polyVec
                               ){
  for (int pIter = 0; pIter < (int)polyVec.size(); pIter++){
    polyVec[pIter].eraseMarkedPolys(markedPolyIndices[pIter]);
  }

  return;
}

void utils::extractMarkedPolys(// Inputs
                               const std::vector<dPoly> & polyVec,
                               std::map< int, std::map<int, int> > & markedPolyIndices,
                               // Outputs
                               std::vector<dPoly> & extractedPolyVec
                               ){

  int pSize = polyVec.size();
  extractedPolyVec.resize(pSize);
  
  for (int pIter = 0; pIter < pSize; pIter++){
    polyVec[pIter].extractMarkedPolys(markedPolyIndices[pIter], // input
                                      extractedPolyVec[pIter]   // output
                                      );
  }
  
  return;
}

int utils::getNumElements(std::map< int, std::map<int, int> > & Indices){

  int num = 0;
  map< int, map<int, int> >::iterator it; 
  for (it = Indices.begin(); it != Indices.end(); it++){
    num += (it->second).size();
  }
  
  return num;
}
