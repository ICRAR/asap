#include "Scantable.h"
#include "STMath.h"
#include "STFillerWrapper.h"
#include <casa/Utilities/CountedPtr.h>
#include <vector>
#include <iostream>

using namespace asap;
using namespace std;

int main() {
  const string fname = "/Users/mar637/data/test.asap";
  //const string fname = "test/data/tid-t002.rpf";
  //STFillerWrapper stf(fname);
  //cout << "reading ..." << endl;
  //stf.read();
  
  casa::CountedPtr<Scantable>  st(new Scantable(fname, casa::Table::Plain));
  STMath stm;
  vector<casa::CountedPtr<Scantable> > on;
  on.push_back(st);
  cout << "averaging..." << endl;
  casa::CountedPtr<Scantable> avst = stm.average(on);
  /*
  cout << "quotient" << endl;
  casa::CountedPtr<Scantable> quot = stm.autoQuotient(avst);
  cout << quot->summary() << endl;
  */
  return 0;

}
