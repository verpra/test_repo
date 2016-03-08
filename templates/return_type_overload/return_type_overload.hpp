#ifndef __return_type_overload_hpp
#define __return_type_overload_hpp
#include <typeinfo>

namespace mynamespace {
  class foo { };
  class bar { };

  template <typename T>
  struct return_type {typedef T type; };
  template<>
  struct return_type<foo> {typedef foo type; };
  template<>
  struct return_type<bar> {typedef bar type; };
  
  class myclass {
    public:

    template <typename T>
    typename return_type<T>::type* foo(T* f) {
      std::cout << "Template " << typeid(T).name() << " Called.." << std::endl;
      return new T;
    }
  
  };

  template <>
  bar* myclass::foo<bar>(bar* b) {
    std::cout << "Specialization " << typeid(bar).name() << " Called.." << std::endl;
    return new bar;
  }
  
}

#endif
