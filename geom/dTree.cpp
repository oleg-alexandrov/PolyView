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
#include <edgeUtils.h>
#include <dTree.h>
#include <geomUtils.h>
#include <dPoly.h>
using namespace std;
using namespace utils;

// See the .h file for documentation.

void edgeTree::putPolyEdgesInTree(const dPoly & poly){

  const double * xv       = poly.get_xv();
  const double * yv       = poly.get_yv();
  const int    * numVerts = poly.get_numVerts();
  int numPolys            = poly.get_numPolys();
  int totalNumVerts       = poly. get_totalNumVerts();
  const auto  &isclosed   = poly.get_isPolyClosed();

  std::vector<utils::dRectWithId>  allBoxes;
  allBoxes.resize(totalNumVerts);
  m_allEdges.resize(totalNumVerts);
  //m_polyEdgeIds.resize(totalNumVerts);
  
  int start = 0;
  for (int pIter = 0; pIter < numPolys; pIter++){
      
    if (pIter > 0) start += numVerts[pIter - 1];
    bool skip_last_edge = !isclosed[pIter] && (numVerts[pIter] > 1);

    for (int vIter = 0; vIter < numVerts[pIter]; vIter++){

      if ( skip_last_edge && vIter == (numVerts[pIter]-1)) continue;

      int vIter2 = (vIter + 1) % numVerts[pIter];
      double bx = xv[start + vIter ], by = yv[start + vIter ];
      double ex = xv[start + vIter2], ey = yv[start + vIter2];

      //m_polyEdgeIds[start + vIter] = make_pair(pIter, vIter);

      dRectWithId R; 
      edgeToBox(bx, by, ex, ey, R );
      R.id = start + vIter;
      allBoxes[start + vIter] = R;

      m_allEdges[start + vIter] = segWidthId(bx, by, ex, ey, start + vIter);

    }
  }
  
  // Form the tree. Boxes will be reordered but otherwise unchanged
  // inside of this function. Do not modify the vector m_allEdges
  // afterward.
  m_boxTree.formTreeOfBoxes(allBoxes);

  return;
}

void edgeTree::findPolyEdgesInBox(// inputs
                                  double xl, double yl,
                                  double xh, double yh,
                                  // outputs
                                  std::vector<utils::segWidthId> & edgesInBox
                                  ){
  
  // Search the tree
  m_boxTree.getBoxesInRegion(xl, yl, xh, yh, m_boxesInRegion);

  // Save the edges in the box
  edgesInBox.clear();
  for (int s = 0; s < (int)m_boxesInRegion.size(); s++){

    const dRectWithId & R = m_boxesInRegion[s]; // alias
    auto edge = m_allEdges[R.id];

    bool res = edgeIntersectsBox(edge.begx, edge.begy, edge.endx, edge.endy,  // arbitrary edge (input)
                                 xl, yl, xh, yh   // box to intersect (input)
                                 );
    if (res) edgesInBox.push_back(edge);
  }

  return;
}

int edgeTree::findClosestEdge( double x0, double y0, utils::seg &closestEdge, double &closestDist) const{

  closestDist = DBL_MAX;

  int root = m_boxTree.getTreeRoot();
  if (root == -1) return -1;
  int edge_id;
  findClosestEdgeToPointInternal(x0, y0, root,            // inputs
                                 edge_id, closestEdge, closestDist);

  closestDist = sqrt(closestDist);
  return edge_id;

}

void edgeTree::findClosestEdgeToPoint(// inputs
                                      double x0, double y0,
                                      // outputs
                                      utils::seg & closestEdge,
                                      double     & closestDist,
                                      // The location on the closest
                                      // edge where closestDist is achieved
                                      double & closestX, double & closestY
                                      ) const{

  // Fast searching for the closest edge to a given point. We assume
  // that the tree of edges is formed by now. The idea is the
  // following: putting the edges in a tree creates a hierarchical
  // partition of the plane. When searching for the closest edge
  // through the regions in this partition, search first the regions
  // closest to the point, and skip altogether the regions which are
  // too far.

  // This function returns DBL_MAX for the closest distance if there
  // are no edges to search.
  
  findClosestEdge(x0, y0, closestEdge, closestDist);


  // Find the point on the closest edge at which the closest distance is achieved
  double distSq;
  minDistSqFromPtToSeg(// inputs
                     x0, y0,
                     closestEdge.begx, closestEdge.begy,
                     closestEdge.endx, closestEdge.endy,
                     // outputs
                     closestX, closestY, distSq
                     );

  return;
}

void edgeTree::findClosestEdgeToPointInternal(// inputs
                                              double x0, double y0,
                                              int root,
                                              // outputs
                                              int &edge_id,
                                              utils::seg & closestEdge,
                                              double     & closestDistSq
                                              ) const{

  // Find the distance from the input point to the edge at the root
  // node. Decide if first to visit the left or right subtree from the
  // root depending on which looks more promising.

  assert (root != -1);
  const auto &bnode = m_boxTree.getBoxNode(root);
  
  const dRectWithId &R = bnode.Rect;

  double dist = DBL_MAX, xval, yval;
  auto edge = m_allEdges[R.id];
  utils::minDistSqFromPtToSeg(x0, y0, edge, xval, yval, dist);
  
  if (dist < closestDistSq){
    closestEdge =  edge;
    closestDistSq = dist;
    edge_id = R.id;
  }

  // Midpoint of the root edge 
  double midx = (edge.begx + edge.endx)/2.0, midy = (edge.begy + edge.endy)/2.0;

  // Unify the left-right and down-up cases to avoid duplicating code.
  double lx0 = x0, ly0 = y0, lmidx = midx, lmidy = midy;
  if (!bnode.isLeftRightSplit){ // down-up split
    swap(lx0,   ly0);
    swap(lmidx, lmidy);
  }
  
  if (lx0 <= lmidx){
    // Search the entire left subtree first
    if (bnode.left != -1)
      findClosestEdgeToPointInternal(x0, y0, bnode.left,      // inputs
                                     edge_id, closestEdge, closestDistSq // outputs
                                     );
    // Don't go right unless there's any chance on improving what we already found
    double dd = bnode.minInRightChild - lx0;
    double dd2 = dd*dd;
    bool bad = (bnode.right == -1 || (dd >= 0 && closestDistSq <= dd2));
    if (!bad) findClosestEdgeToPointInternal(x0, y0, bnode.right,     // inputs
                                             edge_id, closestEdge, closestDistSq // outputs
                                             );
  }else{
    // Search the entire right subtree first
    if (bnode.right != -1)
      findClosestEdgeToPointInternal(x0, y0, bnode.right,     // inputs
                                     edge_id, closestEdge, closestDistSq // outputs
                                     );
    // Don't go left unless there's any chance on improving what we already found
    double dd = lx0 - bnode.maxInLeftChild;
    double dd2 = dd*dd;
    bool bad = (bnode.left == -1 || (dd >= 0 && closestDistSq <= dd2));
    if (!bad) findClosestEdgeToPointInternal(x0, y0, bnode.left,      // inputs
                                             edge_id, closestEdge, closestDistSq // outputs
                                             );
    
  }
  
  return;
}
