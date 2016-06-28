//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef ALLOCATION_UTILS_COMMON_H_
#define ALLOCATION_UTILS_COMMON_H_

#include <algorithm>
#include <stdexcept>

#ifdef __GNUG__
#define GCC_VERSION (__GNUC__ * 10000 \
                   + __GNUC_MINOR__ * 100 \
                   + __GNUC_PATCHLEVEL__)
#if GCC_VERSION < 40500
#define NO_CXX11_ADDRESSOF
#endif
#else //for other compilers assume there are no addressof
#define NO_CXX11_ADDRESSOF
#endif

#ifdef __clang__
#undef NO_CXX11_ADDRESSOF
#endif

#if defined(NO_CXX11_ADDRESSOF)
namespace std
{
  template<typename T>
  T* addressof(T& r) throw()
  {
    return reinterpret_cast<T*>(&const_cast<char&>(reinterpret_cast<const volatile char&>(r)));
  }
}
#endif


/* TODO Add to build script: for gcc 4.5-4.6, use option -std=c++0x. In case if gcc 4.7 or higher is used, use option -std=c++11 */

namespace ls
{
namespace allocation_utils
{
  /**
   * These routines are similar to stl routines but an allocator was added (postfix _a to denote that)
   */

  template<typename _ForwardIterator, typename _Allocator>
  void destroy_a(_ForwardIterator first, _ForwardIterator last, _Allocator& alloc)
  {
    for (; first != last; ++first)
      alloc.destroy( std::addressof(*first) );
  }

  template<typename _ForwardIterator, typename _Size, typename _T, typename _Allocator>
  void uninitialized_fill_n_a(_ForwardIterator first, _Size n, const _T& t, _Allocator& alloc)
  {
    _ForwardIterator cur = first;
    try
    {
      for (; n > 0; --n, ++cur)
        alloc.construct(std::addressof(*cur), t);
    }
    catch (...)
    {
      allocation_utils::destroy_a(first, cur, alloc);
      throw;
    }
  }

  template<typename _InputIterator, typename _ForwardIterator, typename _Allocator>
  _ForwardIterator uninitialized_copy_a(_InputIterator first, _InputIterator last,
         _ForwardIterator result, _Allocator& alloc)
  {
      _ForwardIterator cur = result;
      try
      {
        for (; first != last; ++first, ++cur)
          alloc.construct(std::addressof(*cur), *first);
        return cur;
      }
      catch (...)
      {
        allocation_utils::destroy_a(result, cur, alloc);
        throw;
      }
    }
}
}
#endif /* ALLOCATION_UTILS_COMMON_H_ */
