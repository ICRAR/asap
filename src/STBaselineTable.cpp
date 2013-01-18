//
// C++ Implementation: STBaselineTable
//
// Description:
//
//
// Author: Wataru Kawasaki <wataru.kawasaki@nao.ac.jp> (C) 2013
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <assert.h>

#include <casa/Exceptions/Error.h>
#include <casa/OS/Path.h>
#include <tables/Tables/TableDesc.h>
#include <tables/Tables/SetupNewTab.h>
#include <tables/Tables/ArrColDesc.h>
#include <tables/Tables/ScaColDesc.h>
#include <tables/Tables/TableRecord.h>
#include <measures/TableMeasures/TableMeasDesc.h>
#include <measures/TableMeasures/TableMeasRefDesc.h>
#include <measures/TableMeasures/TableMeasValueDesc.h>

#include "Scantable.h"
#include "STBaselineTable.h"


using namespace casa;

namespace asap {

const String STBaselineTable::name_ = "APPLY_BASELINE";

STBaselineTable::STBaselineTable(const Scantable& parent)
  : STApplyTable(parent, name_)
{
  setup();
}

STBaselineTable::STBaselineTable(const String &name)
  : STApplyTable(name)
{
  attachOptionalColumns();
}

STBaselineTable::~STBaselineTable()
{
}

void STBaselineTable::setup()
{
  table_.addColumn(ScalarColumnDesc<uInt>("NCHAN"));
  table_.addColumn(ScalarColumnDesc<uInt>("FUNCTION"));
  table_.addColumn(ScalarColumnDesc<uInt>("ORDER"));
  table_.addColumn(ScalarColumnDesc<uInt>("CLIP_ITERATION"));
  table_.addColumn(ScalarColumnDesc<Float>("CLIP_THRESHOLD"));
  table_.addColumn(ArrayColumnDesc<Float>("SECTION"));
  table_.addColumn(ArrayColumnDesc<Float>("PARAMETER"));
  table_.addColumn(ArrayColumnDesc<Float>("MASKLIST"));
  table_.addColumn(ScalarColumnDesc<Float>("RMS"));

  table_.rwKeywordSet().define("ApplyType", "BASELINE");

  attachOptionalColumns();
}

void STBaselineTable::attachOptionalColumns()
{
  nchanCol_.attach(table_, "NCHAN");
  funcCol_.attach(table_, "FUNCTION");
  orderCol_.attach(table_, "ORDER");
  clipiterCol_.attach(table_, "CLIP_ITERATION");
  clipthresCol_.attach(table_, "CLIP_THRESHOLD");
  sectCol_.attach(table_, "SECTION");
  paramCol_.attach(table_, "PARAMETER");
  maskCol_.attach(table_, "MASKLIST");
  rmsCol_.attach(table_, "RMS");
}

void STBaselineTable::save(const std::string &filename)
{
  String inname(filename);
  Path path(inname);
  inname = path.expandedName();
  table_.deepCopy(inname, Table::New);
}

void STBaselineTable::setdata(uInt irow, uInt scanno, uInt cycleno, 
                             uInt beamno, uInt ifno, uInt polno, uInt freqid,  
				   Double time, uInt nchan, STBaselineFunc::FuncName func, uInt order, 
                             uInt clipiter, Float clipthres,
                             Vector<Float> sect, Vector<Float> param,
                             Vector<Float> mask, Float rms)
{
  if (irow >= (uInt)nrow()) {
    throw AipsError("row index out of range");
  }

  if (!sel_.empty()) {
    os_.origin(LogOrigin("STBaselineTable","setdata",WHERE));
    os_ << LogIO::WARN << "Data selection is effective. Specified row index may be wrong." << LogIO::POST;
  }  

  setbasedata(irow, scanno, cycleno, beamno, ifno, polno, freqid, time);
  nchanCol_.put(irow, nchan);
  funcCol_.put(irow, uInt(func));
  orderCol_.put(irow, order);
  clipiterCol_.put(irow, clipiter);
  clipthresCol_.put(irow, clipthres);
  sectCol_.put(irow, sect);
  paramCol_.put(irow, param);
  maskCol_.put(irow, mask);
  rmsCol_.put(irow, rms);
}

void STBaselineTable::appenddata(uInt scanno, uInt cycleno, 
				      uInt beamno, uInt ifno, uInt polno, uInt freqid,
				      Double time, uInt nchan, STBaselineFunc::FuncName func, uInt order, 
				      uInt clipiter, Float clipthres, 
				      Vector<Float> sect, Vector<Float> param,
				      Vector<Float> mask, Float rms)
{
  uInt irow = nrow();
  table_.addRow(1, True);
  setdata(irow, scanno, cycleno, beamno, ifno, polno, freqid, time, 
	  nchan, func, order, clipiter, clipthres, sect, param, mask, rms);
}

Vector<STBaselineFunc::FuncName> STBaselineTable::getFunctionAsString()
{
  Vector<uInt> rawBlfuncColumn = funcCol_.getColumn();
  uInt n = rawBlfuncColumn.nelements();
  Vector<STBaselineFunc::FuncName> blfuncColumn(n);
  for (uInt i = 0; i < n; ++i) {
    blfuncColumn[i] = STBaselineFunc::FuncName(rawBlfuncColumn(i));
  }
  return blfuncColumn;
}

uInt STBaselineTable::nchan(uInt ifno)
{
  STSelector org = sel_;
  STSelector sel;
  sel.setIFs(vector<int>(1,(int)ifno));
  setSelection(sel);
  uInt n = nchanCol_(0);
  unsetSelection();
  if (!org.empty())
    setSelection(org);
  return n;
}
}
