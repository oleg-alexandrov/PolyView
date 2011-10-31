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
#include <fstream>
#include <algorithm>
#include <vector>
#include <cassert>
#include <cstring>
#include "edgeUtils.h"
#include "dTree.h"
#include "geomUtils.h"
#include "dPoly.h"
using namespace std;
using namespace utils;

// See the .h file for documentation.

void edgeTree::putPolyEdgesInTree(const dPoly & poly){

  const double * xv       = poly.get_xv();
  const double * yv       = poly.get_yv();
  const int    * numVerts = poly.get_numVerts();
  int numPolys            = poly.get_numPolys();
  int totalNumVerts       = poly. get_totalNumVerts();
  
  m_allEdges.resize(totalNumVerts);
  
  int start = 0;
  for (int pIter = 0; pIter < numPolys; pIter++){
      
    if (pIter > 0) start += numVerts[pIter - 1];

    for (int vIter = 0; vIter < numVerts[pIter]; vIter++){

      int vIter2 = (vIter + 1) % numVerts[pIter];
      double bx = xv[start + vIter ], by = yv[start + vIter ];
      double ex = xv[start + vIter2], ey = yv[start + vIter2];

      // Transform an edge into a box, with the id storing
      // the information necessary to reverse this later.
      dRectWithId R; 
      edgeToBox(//inputs
                bx, by, ex, ey, // inputs  
                R               // output
                );
      m_allEdges[start + vIter] = R;
      
    }
  }
  
  // Form the tree. Boxes will be reordered but otherwise unchanged
  // inside of this function. Do not modify the vector m_allEdges
  // afterward.
  m_boxTree.formTreeOfBoxes(m_allEdges);

  return;
}

void edgeTree::findPolyEdgesInBox(// inputs
                                  double xl, double yl,
                                  double xh, double yh,
                                  // outputs
                                  std::vector<utils::seg> & edgesInBox
                                  ){
  
  // Search the tree
  m_boxTree.getBoxesInRegion(xl, yl, xh, yh, m_boxesInRegion);

  // Save the edges in the box
  edgesInBox.clear();
  for (int s = 0; s < (int)m_boxesInRegion.size(); s++){

    const dRectWithId & R = m_boxesInRegion[s]; // alias
    double bx, by, ex, ey;
    boxToEdge(R,             // input
              bx, by, ex, ey // outputs
              );

    bool res = edgeIntersectsBox(bx, by, ex, ey,  // arbitrary edge (input)
                                 xl, yl, xh, yh   // box to intersect (input)
                                 );
    if (res) edgesInBox.push_back(seg(bx, by, ex, ey));
  }

  return;
}

void edgeTree::findClosestEdgeToPoint(// inputs
                                      double x0, double y0,
                                      // outputs
                                      utils::seg & closestEdge,
                                      double     & closestDist,
                                      // The location on the closest
                                      // edge where closestDist is achieved
                                      double & closestX, double & closestY
                                      ){

  // Fast searching for the closest edge to a given point. We assume
  // that the tree of edges is formed by now. The idea is the
  // following: putting the edges in a tree creates a hierarchical
  // partition of the plane. When searching for the closest edge
  // through the regions in this partition, search first the regions
  // closest to the point, and skip altogether the regions which are
  // too far.

  // This function returns DBL_MAX for the closest distance if there
  // are no edges to search.
  
  closestDist = DBL_MAX;

  boxNode<dRectWithId> * root = m_boxTree.getTreeRoot();
  if (root == NULL) return;
    
  findClosestEdgeToPointInternal(x0, y0, root,            // inputs 
                                 closestEdge, closestDist // outputs
                                 );

  // Find the point on the closest edge at which the closest distance is achieved
  minDistFromPtToSeg(// inputs
                     x0, y0,
                     closestEdge.begx, closestEdge.begy,
                     closestEdge.endx, closestEdge.endy,
                     // outputs
                     closestX, closestY, closestDist
                     );

  return;
}

void edgeTree::findClosestEdgeToPointInternal(// inputs
                                              double x0, double y0,
                                              boxNode<utils::dRectWithId> * root,
                                              // outputs
                                              utils::seg & closestEdge,
                                              double     & closestDist
                                              ){

  // Find the distance from the input point to the edge at the root
  // node. Decide if first to visit the left or right subtree from the
  // root depending on which looks more promising.
  
  assert (root != NULL);
  
  dRectWithId R = root->Rect;
  double bx, by, ex, ey;
  boxToEdge(R,             // input
            bx, by, ex, ey // outputs
            );

  double dist = DBL_MAX, xval, yval;
  minDistFromPtToSeg(// inputs
                     x0, y0, bx, by, ex, ey,
                     // outputs
                     xval, yval, dist
                     );
  
  if (dist < closestDist){
    closestEdge = seg(bx, by, ex, ey);
    closestDist = dist;
  }

  // Midpoint of the root edge 
  double midx = (bx + ex)/2.0, midy = (by + ey)/2.0;

  // Unify the left-right and down-up cases to avoid duplicating code.
  double lx0 = x0, ly0 = y0, lmidx = midx, lmidy = midy;
  if (!root->isLeftRightSplit){ // down-up split
    swap(lx0,   ly0);
    swap(lmidx, lmidy);
  }
  
  if (lx0 <= lmidx){
    // Search the entire left subtree first
    if (root->left != NULL)
      findClosestEdgeToPointInternal(x0, y0, root->left,      // inputs 
                                     closestEdge, closestDist // outputs
                                     );
    // Don't go right unless there's any chance on improving what we already found
    bool bad = (root->right == NULL || lx0 + closestDist <= root->minInRightChild);
    if (!bad) findClosestEdgeToPointInternal(x0, y0, root->right,     // inputs 
                                             closestEdge, closestDist // outputs
                                             );
  }else{
    // Search the entire right subtree first
    if (root->right != NULL)
      findClosestEdgeToPointInternal(x0, y0, root->right,     // inputs 
                                     closestEdge, closestDist // outputs
                                     );
    // Don't go left unless there's any chance on improving what we already found
    bool bad = (root->left == NULL || root->maxInLeftChild + closestDist <= lx0);
    if (!bad) findClosestEdgeToPointInternal(x0, y0, root->left,      // inputs 
                                             closestEdge, closestDist // outputs
                                             );
    
  }
  
  return;
}
