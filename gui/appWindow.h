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
#ifndef EXAMPLE_H
#define EXAMPLE_H

#include <qmainwindow.h>
#include <qlineedit.h>
#include <QEvent>
#include <string>
#include <vector>
class polyView;
class QCloseEvent;
class cmdLineOptions;
class QTextEdit;

class cmdLine : public QLineEdit {
  Q_OBJECT
  
public:
  cmdLine(QWidget* parent);
  virtual ~cmdLine();
};

class docWindow: public QMainWindow{
  Q_OBJECT

public:
  docWindow(QWidget * parent = NULL);
  virtual ~docWindow();
  void setText(const std::string & docText);
private:
  QTextEdit * m_textArea;
};

class appWindow : public QMainWindow {
  Q_OBJECT

public:
  appWindow(QWidget* parent, std::string progName,
            const cmdLineOptions & options, 
            int windowWidX, int windowWidY
            );
  ~appWindow();
  
protected:
  bool eventFilter(QObject *obj, QEvent *event);

private slots:
  void createMenusAndMainWidget(const cmdLineOptions & opt);
  void showDoc();
  void about();
  void procCmdLine();
  void shiftUp();
  void shiftDown();
  void forceQuit();
  
private:
  void closeEvent(QCloseEvent *);
  void insertCmdFromHist();
  
  polyView       * m_poly;
  cmdLine        * m_cmdLine;
  std::string      m_progName;
  std::vector<std::string> m_cmdHist;
  int m_histPos;
  docWindow m_docWindow;
};


#endif
