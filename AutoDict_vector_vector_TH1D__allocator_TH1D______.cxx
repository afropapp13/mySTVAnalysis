#include "vector"
class TH1D;
#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class vector<vector<TH1D*,allocator<TH1D*> > >+;
#pragma link C++ class vector<vector<TH1D*,allocator<TH1D*> > >::*;
#ifdef G__VECTOR_HAS_CLASS_ITERATOR
#pragma link C++ operators vector<vector<TH1D*,allocator<TH1D*> > >::iterator;
#pragma link C++ operators vector<vector<TH1D*,allocator<TH1D*> > >::const_iterator;
#pragma link C++ operators vector<vector<TH1D*,allocator<TH1D*> > >::reverse_iterator;
#endif
#endif
