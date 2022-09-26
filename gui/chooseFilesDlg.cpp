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
#include <cassert>
#include <iostream>
#include <QDialogButtonBox>
#include <QLabel>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QHeaderView>
#include <QWidget>
#include <chooseFilesDlg.h>
#include <utils.h>
using namespace std;

// Allow the user to choose which files to hide/show in the GUI.
// User's choice will be processed by polyView::showFilesChosenByUser().
chooseFilesDlg::chooseFilesDlg(QWidget * parent): QWidget(parent) {

  setWindowModality(Qt::ApplicationModal); 

  int spacing = 0;
  
  QVBoxLayout * vBoxLayout = new QVBoxLayout(this);
  vBoxLayout->setSpacing(spacing);
  vBoxLayout->setAlignment(Qt::AlignLeft);

  // The layout having the file names. It will be filled in dynamically later.
  m_filesTable = new QTableWidget();

  m_filesTable->verticalHeader()->hide();

  vBoxLayout->addWidget(m_filesTable);

  return;
}

chooseFilesDlg::~chooseFilesDlg(){}

void chooseFilesDlg::chooseFiles(const std::vector<polyOptions> & optionsVec){

  // See the top of this file for documentation.

  int numFiles = optionsVec.size();
  int numCols = 1;
  m_filesTable->setRowCount(numFiles);
  m_filesTable->setColumnCount(numCols);

  for (int fileIter = 0; fileIter < numFiles; fileIter++) {
    string fileName = optionsVec[fileIter].polyFileName;
    QTableWidgetItem *item = new QTableWidgetItem(fileName.c_str());
    item->data(Qt::CheckStateRole);
    item->setCheckState(Qt::Checked);
    item->setForeground(QColor::fromRgb(0, 0, 0));
    
    m_filesTable->setItem(fileIter, numCols - 1, item);
  }

#if 1
  // Sort this out later, for now hide the header caption
    m_filesTable->horizontalHeader()->hide();
#else
  // Horizontal header caption
  QTableWidgetItem *item = new QTableWidgetItem("Hide/show all");
  item->setFlags(Qt::NoItemFlags);
  item->setForeground(QColor::fromRgb(0, 0, 0));
  m_filesTable->setHorizontalHeaderItem(1, item);
#endif
  
  m_filesTable->setSelectionMode(QTableWidget::ExtendedSelection);
  std::string style = std::string("QTableWidget::indicator:unchecked ")
    + "{background-color:white; border: 1px solid black;}; " +
    "selection-background-color: rgba(128, 128, 128, 40);";
  m_filesTable->setStyleSheet(style.c_str());
  
  m_filesTable->resizeColumnsToContents();
  m_filesTable->resizeRowsToContents();
  
  // The processing of user's choice happens in polyView::showFilesChosenByUser()

  bool isVisible = (numFiles > 1);
  this->setVisible(isVisible);

  return;
}

void chooseFilesDlg::keyPressEvent(QKeyEvent * /*event*/) {
  // std::cout << "Key was pressed " << event->key() << std::endl;
}
