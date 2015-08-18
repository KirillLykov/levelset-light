//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <stdexcept>
#include <string>
#include <assert.h>

namespace ls
{
  class array_size_error : public std::logic_error
  {
  public:
    explicit array_size_error(const std::string& msg) throw()
      : std::logic_error(msg)
    { }
  };

#ifndef THROW_ON_INDEX_ERRORS
  #define INDEX_ASSERT(_condition) { assert(_condition); }
#else
  #define INDEX_ASSERT(_condition) { if (!(_condition)) throw std::out_of_range(""); }
#endif
}

#endif /* EXCEPTIONS_H_ */
