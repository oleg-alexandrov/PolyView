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
#include <ctime>
#include "dTree.h"
#include "dPoly.h"

using namespace std;
using namespace utils;

double rand_ab(double a, double b){
  assert(a <= b);
  return a + rand()%max(int(b - a), 1);
}

void saveBoxes(std::vector<dRect> & Boxes, std::string file, std::string color){
  cout << "Writing " << file << endl;
  ofstream of(file.c_str());
  of << "color = " << color << endl;
  for (int s = 0; s < (int)Boxes.size(); s++){
    const dRect & B = Boxes[s]; // alias
    of << B.xl << ' ' << B.yl << endl;
    of << B.xh << ' ' << B.yl << endl;
    of << B.xh << ' ' << B.yh << endl;
    of << B.xl << ' ' << B.yh << endl;
    of << B.xl << ' ' << B.yl << endl;
    of << "NEXT" << endl;
  }
  return;
}


int main(int argc, char** argv){

 srand(time(NULL));

 int numBoxes   = 1000;
 int L          = 100, w = 10; // To determine region size and box width 
 int numRepeats = 1;           // How many times to repeat each experiment
 int numRuns    = 10000;       // How many experiments
 bool doTiming  = false; 

 if (doTiming){
   numRuns = 1;
   numRepeats = 100000;
 }
 
 vector<dRect> Boxes;

 for (int q = 0; q < numRuns; q++){
   
    Boxes.clear();

   for (int s = 0; s < numBoxes; s++){
     double xl = rand_ab(-L, L), xh = xl + rand_ab(0, w);
     double yl = rand_ab(-L, L), yh = yl + rand_ab(0, w);
     Boxes.push_back(dRect(xl, yl, xh, yh));
   }
 
   boxTree<dRect> T;
   // Boxes will be reordered but otherwise unchanged inside of this function.
   // Do not modify this vector afterward.
   T.formTreeOfBoxes(Boxes); 

   double L2 = 2*L;
   double xl = rand_ab(-L2, L2), xh = xl + rand_ab(0, L2);
   double yl = rand_ab(-L2, L2), yh = yl + rand_ab(0, L2);

   time_t Start_t = 0, End_t = 0;
   double time_task;
   
   vector<dRect> outBoxes;
   if(doTiming) Start_t = time(NULL);    //record time that task 1 begins
   for (int v = 0; v < numRepeats; v++){
     T.getBoxesInRegion(xl, yl, xh, yh, outBoxes);
   }
   if(doTiming){
     End_t = time(NULL);  
     time_task = difftime(End_t, Start_t);
     cout << "Tree search time: " << time_task << endl;
   }

   sort(outBoxes.begin(), outBoxes.end(), lexLessThan<dRect>);
 
   vector<dRect> outBoxes2; // Must be different than Boxes
   outBoxes2.resize(numBoxes);
   
   if(doTiming) Start_t = time(NULL);
   for (int v = 0; v < numRepeats; v++){
     outBoxes2.clear();
     for (int s = 0; s < (int)Boxes.size(); s++){
       const dRect & B = Boxes[s];
       if (boxesIntersect(B.xl, B.yl, B.xh, B.yh, xl, yl, xh, yh)){
         outBoxes2.push_back(B);
       }
     }
   }
   if(doTiming){
     End_t = time(NULL);
     time_task = difftime(End_t, Start_t);
     cout << "Brute force search time: " << time_task << endl;
   }

   sort(outBoxes2.begin(), outBoxes2.end(), lexLessThan<dRect>);

   bool haveProblem = false;
 
   int len = outBoxes.size();
   if ( len != (int)outBoxes2.size() ){
     haveProblem = true;
   }else{
     for (int s = 0; s < len; s++){
       if ( outBoxes[s] != outBoxes2[s] ){
         haveProblem = true;
       }
     }
   }

   if (haveProblem){
     cerr << "Have a problem!" << endl;

     cout << "In:   " << Boxes.size()     << endl;
     cout << "Out:  " << outBoxes.size()  << endl;
     cout << "Out2: " << outBoxes2.size() << endl;
     
     vector<dRect> W;
     W.push_back(dRect(xl, yl, xh, yh));
   
     saveBoxes( Boxes,    "all.xg",       "blue"   );
     saveBoxes( W,        "window.xg",    "green"  );
     saveBoxes( outBoxes, "outBoxes.xg",  "white"  );
     saveBoxes( outBoxes, "outBoxes2.xg", "yellow" );
     exit(1);
   }

   if (q%1000 == 0){
     cout << "In:   " << Boxes.size()     << endl;
     cout << "Out:  " << outBoxes.size()  << endl;
     cout << "Out2: " << outBoxes2.size() << endl;
   }
   
 }

 return 0;
}
