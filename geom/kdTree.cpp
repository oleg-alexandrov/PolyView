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
#include <cfloat> // defines DBL_MAX
#include <kdTree.h>
#include <baseUtils.h>
#include <geomUtils.h>

using namespace std;
using namespace utils;

kdTree::kdTree(){
  reset();
  return;
}

void kdTree::reset(){
  m_freeNodeIndex = 0;
  m_root          = -1;
  m_nodePool.clear();
  return;
}


int kdTree::getNewNode(){
  // Get a node from the pool
  assert( m_freeNodeIndex < (int)m_nodePool.size() );
  int new_ind = m_freeNodeIndex;
  m_freeNodeIndex++;
  return new_ind;
}

void kdTree::formTreeOfPoints(int numPts, const double * xv, const double * yv){
  
  vector<PointWithId> Pts;
  Pts.resize(numPts);
  for (int s = 0; s < numPts; s++){
    Pts[s] = PointWithId(xv[s], yv[s], s);
  }
  
  formTreeOfPoints(// Pts will be reordered but otherwise
                   // unchanged inside this function
                   Pts
                   );
  return;
}

void kdTree::formTreeOfPoints(// Pts will be reordered but otherwise
                              // unchanged inside this function
                              std::vector<utils::PointWithId> & Pts 
                              ){

  reset();
  int numPts = Pts.size();
  m_nodePool.resize(numPts);

  bool isLeftRightSplit = true;
  formTreeOfPointsInternal(vecPtr(Pts), numPts, isLeftRightSplit, m_root);
  return;
}


void kdTree::formTreeOfPointsInternal(utils::PointWithId * Pts, int numPts, bool isLeftRightSplit,
                                      int &root){

  // To do: Implement this without recursion.
  // To do: No need to store the point P in the tree. Store just a pointer to P
  // as sorting the left and right halves does not change the midpoint P.
  // To do: No need even for the tree, it can be stored in-place in Pts,
  // in the same way as an in-place heap is stored.
  
  assert(numPts >= 0);
  
  if (numPts == 0){
    root = -1;
    return; 
  }
  
  root = getNewNode();
  utils::Node &node = getNode(root);
  node.isLeftRightSplit = isLeftRightSplit;
  
  if (isLeftRightSplit){ // Split points into left and right halves
    sort(Pts, Pts + numPts, leftLessThan);
  }else{               // Split points into bottom and top halves 
    sort(Pts, Pts + numPts, botLessThan);
  }

  int mid = numPts/2;
  assert( 0 <= mid && mid < numPts);
  
  node.P = Pts[mid]; // Must happen after sorting

  // At the next split we will split perpendicularly to the direction
  // of the current split.
  formTreeOfPointsInternal( Pts,            mid,
                            !isLeftRightSplit, node.left
                            );
  formTreeOfPointsInternal( Pts + mid + 1,  numPts - mid - 1,
                            !isLeftRightSplit, node.right
                            );

}

void kdTree::getPointsInBox(double xl, double yl, double xh, double yh, // input box
                            std::vector<utils::PointWithId> & outPts) const{

  outPts.clear();
  if (xl > xh || yl > yh) return;
  getPointsInBoxInternal(xl, yl, xh, yh, m_root, // inputs
                         outPts                  // outputs
                         );
  return;
}

void kdTree::getPointsInBoxInternal(//Inputs
                                    double xl, double yl, double xh, double yh,
                                    int root,
                                    // Outputs
                                    std::vector<utils::PointWithId> & outPts) const{

  // To do: Can this be done without recursion?
  
  if (root == -1) return;

  const utils::Node &node = getNode(root);

  const PointWithId & P = node.P; // alias
  
  if (xl <= P.x && P.x <= xh && yl <= P.y && P.y <= yh) outPts.push_back(P);
  
  if (node.isLeftRightSplit){
    if (xl <= P.x) getPointsInBoxInternal(xl, yl, xh, yh, node.left,  outPts);
    if (xh >= P.x) getPointsInBoxInternal(xl, yl, xh, yh, node.right, outPts);
  }else{
    if (yl <= P.y) getPointsInBoxInternal(xl, yl, xh, yh, node.left,  outPts);
    if (yh >= P.y) getPointsInBoxInternal(xl, yl, xh, yh, node.right, outPts);
  }
    
  return;
}

void kdTree::findClosestVertexToPoint(// inputs
                                      double x0, double y0,
                                      // outputs
                                      utils::PointWithId & closestVertex,
                                      double & closestDist
                                      ) const{

  // Fast searching for the closest vertex in the tree to a given
  // point. We assume that the tree of vertices is formed by now. The
  // idea is the following: putting the vertices in a tree creates a
  // hierarchical partition of the plane. When searching for the
  // closest vertex through the regions in this partition, search
  // first the regions closest to the point, and skip altogether the
  // regions which are too far.

  // This function returns DBL_MAX for the closest distance if there
  // are no vertices to search.
  
  closestDist = DBL_MAX;

  if (m_root == -1) return;
    
  findClosestVertexToPointInternal(x0, y0, m_root,            // inputs 
                                   closestVertex, closestDist // outputs
                                   );


  return;
}

void kdTree::findClosestVertexToPointInternal(// inputs
                                              double x0, double y0,
                                              int root,
                                              // outputs
                                              utils::PointWithId & closestVertex,
                                              double & closestDist
                                              ) const{

  // Find the distance from the input point to the vertex at the root
  // node. Decide if first to visit the left or right subtree from the
  // root depending on which looks more promising.
  
  assert (root != -1);
  
  const utils::Node &node = getNode(root);
  const PointWithId & P = node.P; // alias

  double dist = distance(x0, y0, P.x, P.y);
  
  if (dist < closestDist){
    closestVertex = P;
    closestDist   = dist;
  }

  // Unify the left-right and down-up cases to avoid duplicating code.
  double lx0 = x0, ly0 = y0, rootX = P.x, rootY = P.y;
  if (!node.isLeftRightSplit){ // down-up split
    swap(lx0,   ly0);
    swap(rootX, rootY);
  }
  
  if (lx0 <= rootX){
    // Search the entire left subtree first
    if (node.left != -1)
      findClosestVertexToPointInternal(x0, y0, node.left,        // inputs
                                       closestVertex, closestDist // outputs
                                       );
    // Don't go right unless there's any chance on improving what we already found
    bool bad = (node.right == -1 || lx0 + closestDist <= rootX);
    if (!bad) findClosestVertexToPointInternal(x0, y0, node.right,       // inputs
                                               closestVertex, closestDist // outputs
                                               );
  }else{
    // Search the entire right subtree first
    if (node.right != -1)
      findClosestVertexToPointInternal(x0, y0, node.right,       // inputs
                                       closestVertex, closestDist // outputs
                                       );
    // Don't go left unless there's any chance on improving what we already found
    bool bad = (node.left == -1 || rootX + closestDist <= lx0);
    if (!bad) findClosestVertexToPointInternal(x0, y0, node.left,        // inputs
                                               closestVertex, closestDist // outputs
                                               );
    
  }
  
  return;
}
