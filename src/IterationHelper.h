#include <vector>

#include <casa/Arrays/Vector.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/CountedPtr.h>

#include "Scantable.h"
#include "STTcal.h"
#include "STIdxIter.h"
#include "STSelector.h"

using namespace casa;
using namespace asap;

namespace {
vector<string> split(const string &str, char delim)
{
  vector<string> result;
  size_t current = 0;
  size_t found;
  while ((found = str.find_first_of(delim, current)) != string::npos) {
    result.push_back(string(str, current, found - current));
    current = found + 1;
  }
  result.push_back(string(str, current, str.size() - current));
  return result;
}
  
// Iteration Helper
template<class T>
class IterationHelper
{
public:
  static void Iterate(T &processor, const string cols_list)
  {
    vector<string> cols = split(cols_list, ',');
    // for (vector<string>::iterator i = cols.begin(); i != cols.end(); ++i)
    //   cout << *i << endl;
    STIdxIter2 iter(processor.target(), cols);
    STSelector sel ;
    while ( !iter.pastEnd() ) {
      const Record current = iter.currentValue() ;
      Vector<uInt> rows = iter.getRows( SHARE ) ;
      // any process
      processor.Process(cols, current, rows);
      // go next
      iter.next() ;
    }    
  }
};

  
} // anonymous namespace
