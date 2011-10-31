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
#include <qstatusbar.h>
#include <qmessagebox.h>
#include <qmenubar.h>
#include <qapplication.h>
#include <qimage.h>
#include <iostream>
#include <stdlib.h>
#include "appWindow.h"
#include "polyView.h"
#include "utils.h"

using namespace std;
using namespace utils;

int main(int argc, char** argv){

  char * exeName = argv[0];

  int windowWidX, windowWidY;
  cmdLineOptions options;
  
  parseCmdOptions(// inputs
                  argc, argv, exeName,
                  // outputs
                  windowWidX, windowWidY, options
                  );
  
  QApplication app(argc, argv);
  string progName = "PolyView";
  
  appWindow m(NULL, progName, options, windowWidX, windowWidY);
  m.setCaption(progName.c_str());
  m.show();
  
  QObject::connect( qApp, SIGNAL(lastWindowClosed()), qApp, SLOT(quit()) );
  
  return app.exec();
}
