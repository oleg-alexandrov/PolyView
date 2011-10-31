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
#ifndef KDTREE_H
#define KDTREE_H

// A naive (but hopefully correct) implementation of kd-trees.

#include <vector>

namespace utils{ 

  struct PointWithId{
    double x, y;
    int id;
    PointWithId(): x(0), y(0), id(0){}
    PointWithId(double x_in, double y_in, int id_in): x(x_in), y(y_in), id(id_in){}
  };

  inline bool leftLessThan (PointWithId P, PointWithId Q){ return P.x < Q.x; }
  inline bool botLessThan  (PointWithId P, PointWithId Q){ return P.y < Q.y; }
  inline bool lexLessThan  (PointWithId P, PointWithId Q){
    return (P.x < Q.x) || ( (P.x == Q.x) && (P.y < Q.y) );
  }

  struct Node{
    Node * left;
    Node * right;
    utils::PointWithId  P;
    bool   isLeftRightSplit;
    Node(): left(NULL), right(NULL), isLeftRightSplit(false){}
  };
  
}

class kdTree{

public:
  kdTree();

  void formTreeOfPoints(int numPts, const double * xv, const double * yv);
  
  void formTreeOfPoints(// Pts will be reordered but otherwise
                        // unchanged inside this function
                        std::vector<utils::PointWithId> & Pts 
                        );
  void getPointsInBox(double xl, double yl, double xh, double yh, // input box
                      std::vector<utils::PointWithId> & outPts);

  void findClosestVertexToPoint(// inputs
                                double x0, double y0,
                                // outputs
                                utils::PointWithId & closestVertex,
                                double & closestDist
                                );

private:
  void findClosestVertexToPointInternal(// inputs
                                        double x0, double y0,
                                        utils::Node * root,
                                        // outputs
                                        utils::PointWithId & closestVertex,
                                        double & closestDist
                                        );

  void formTreeOfPointsInternal(utils::PointWithId * Pts, int numPts, bool isLeftRightSplit,
                                utils::Node *&  root);

  void getPointsInBoxInternal(//Inputs
                              double xl, double yl, double xh, double yh,
                              utils::Node * root,
                              // Outputs
                              std::vector<utils::PointWithId> & outPts);
  void reset();
  utils::Node * getNewNode();
  
  // Will get nodes from this pool (for performance reasons)
  std::vector<utils::Node> m_nodePool; 

  int m_freeNodeIndex;
  utils::Node * m_root;
    
};
  
#endif
