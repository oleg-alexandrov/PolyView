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
#ifndef CHOOSE_FILES_DLG_H
#define CHOOSE_FILES_DLG_H
#include <set>
#include <string>
#include <vector>
#include <QDialog>

class QWidget;
class QTableWidget;
struct polyOptions;

class chooseFilesDlg: public QDialog{
  Q_OBJECT
  
public:
  chooseFilesDlg(QWidget * parent = NULL);
  ~chooseFilesDlg();
  void chooseFiles(const std::vector<polyOptions> & optionsVec);
  QTableWidget * getFilesTable(){ return m_filesTable; }
  static QString selectFilesTag(){ return "Select files to hide/show"; }
private:
  QTableWidget * m_filesTable;

private slots:
};

#endif // CHOOSE_FILES_DLG_H
