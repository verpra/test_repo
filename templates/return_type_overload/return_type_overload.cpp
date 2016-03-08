#include <iostream>
#include "return_type_overload.hpp"

using namespace mynamespace;

int
main(int argc, char** argv) {
  myclass m;
  bar* ba;
  foo* fo;
  int *pi;
  float *fl;
  int* pint = m.foo(pi);
  foo* pfoo = m.foo(fo);
  bar* pbar = m.foo(ba);
  float *flo = m.foo(fl);
  
}
