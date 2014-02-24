//  (C) Copyright Kirill Lykov 2012.
//
// Distributed under the FreeBSD Software License (See accompanying file license.txt)

#ifndef GRID_H_
#define GRID_H_

#include "allocation_utils_common.h"
#include <stdexcept>
#include <assert.h>
#include "exceptions.h"
#include "box.h"

namespace ls
{
//TODO use macro for handling exceptional behavior:
//it may throw an exception or an assertion. Exceptions are handy for testing

  template< class T, class allocator >
  class Grid_impl
  {
  protected:

    typedef typename allocator::pointer _pointer;
    typedef typename allocator::size_type _size_type;
    typedef Grid_impl<T, allocator> _TGridImpl;

    mutable allocator m_allocator;
    _size_type m_linearSz;
    typename allocator::pointer m_data;

    Grid_impl()
      : m_linearSz(0), m_data(0)
    {
    }

    Grid_impl(_size_type linesize)
       : m_linearSz(linesize), m_data(0)
     {
       if (linesize == 0)
         throw array_size_error("grid must not have 0 size");

       m_data = m_allocator.allocate(m_linearSz);
       allocation_utils::uninitialized_fill_n_a(m_data, m_linearSz, T(), m_allocator);
     }

    Grid_impl(const _TGridImpl& another)
      : m_linearSz(another.m_linearSz)
    {
      m_data = m_allocator.allocate(m_linearSz);
      std::uninitialized_copy(another.m_data, another.m_data + m_linearSz, this->m_data);
    }

    virtual ~Grid_impl()
    {
      allocation_utils::destroy_a(m_data, m_data + m_linearSz, m_allocator);
      m_allocator.deallocate(m_data, m_linearSz);
    }

    _TGridImpl& operator= (const _TGridImpl& another)
    {
      if (this == &another)
       return *this;

      std::copy(another.m_data, another.m_data + m_linearSz, this->m_data);

      return *this;
    }

    bool operator== (const _TGridImpl& another) const
    {
      return std::equal(m_data, m_data + m_linearSz, another.m_data);
    }

    bool operator!= (const _TGridImpl& another) const
    {
      return !this->operator== (another);
    }

    void swap(_TGridImpl& another)
    {
      std::swap(m_allocator, another.m_allocator);
      std::swap(m_data, another.m_data);
      std::swap(m_linearSz, another.m_linearSz);
    }

    void resize(_size_type newSize)
    {
      // limited implementation - supposed to be used only
      // when grid was constructed using default constructor
      // so is is called from grid with some data, data will be lost

      // allocate new memory block
      _pointer new_data(m_allocator.allocate(newSize));
      try
      {
        // fill new objects by default values
        allocation_utils::uninitialized_fill_n_a(new_data, newSize, T(), m_allocator);
        //TODO fix copying - code under comment is not correct for cube or square
        // copy old values to the first part of the buffer
        //if (m_data != 0)
        //  allocation_utils::uninitialized_copy_a(m_data, m_data + m_linearSz, __new_start, m_allocator);
      }
      catch(...)
      {
        m_allocator.deallocate(new_data, newSize);
        throw;
      }
      // remove data from the old buffer
      allocation_utils::destroy_a(m_data, m_data + m_linearSz, m_allocator);
      m_allocator.deallocate(m_data, m_linearSz);

      m_data = new_data;
      m_linearSz = newSize;
    }
  };

  template< class T, class allocator = std::allocator<T> >
  class Grid2D : public Grid_impl<T, allocator>
  {
    typedef Grid_impl<T, allocator> _TGridImpl;
    typedef typename _TGridImpl::_size_type _size_type;
    typedef typename allocator::reference _reference;
    typedef typename allocator::const_reference _const_reference;
    _size_type m_n1, m_n2;
    geometry_utils::Box2D m_bbox;
  public:

    typedef Grid2D<T, allocator> _TGrid;

    Grid2D()
      : m_n1(0), m_n2(0)
    {
    }

    Grid2D(_size_type n1, _size_type n2, geometry_utils::Box2D box = geometry_utils::Box2D())
      : _TGridImpl(n1 * n2), m_n1(n1), m_n2(n2), m_bbox(box)
    {
    }

    Grid2D(const _TGrid& another)
      : _TGridImpl(another), m_n1(another.m_n1), m_n2(another.m_n2), m_bbox(another.m_bbox)
    {
    }

    _TGrid& operator= (const _TGrid& another)
    {
      if (this == &another)
        return *this;

      if (m_n1 != another.m_n1 || m_n2 != another.m_n2)
        throw array_size_error("grids in the assignment must have the same sizes");

      m_bbox = another.m_bbox;
      _TGridImpl::operator= (another);

      return *this;
    }

    bool operator== (const _TGrid& another) const
    {
      if (m_n1 != another.m_n1 || m_n2 != another.m_n2)
        throw array_size_error("logical operators on grids may be applied only to grids having the same sizes");

      if (m_bbox != another.m_bbox) {
        return false; // if domains bounding boxes are different
      }

      return _TGridImpl::operator== (another);
    }

    bool operator!= (const _TGrid& another) const
    {
      return !this->operator== (another);
    }

    _const_reference operator() (_size_type i, _size_type j) const
    {
      assert(i < m_n1 && j < m_n2);
      return _TGridImpl::m_data[i + m_n1 * j];
    }

    _reference operator() (_size_type i, _size_type j)
    {
      assert(i < m_n1 && j < m_n2);
      return _TGridImpl::m_data[i + m_n1 * j];
    }

    void swap(_TGrid& another)
    {
      _TGridImpl::swap(another);
      std::swap(m_n1, another.m_n1);
      std::swap(m_n2, another.m_n2);
      std::swap(m_bbox, another.m_bbox);
    }

    _size_type size(_size_type index) const
    {
      INDEX_ASSERT(index <= 2);
      return index == 0 ? m_n1 : m_n2;
    }

    void resize(_size_type n1, _size_type n2)
    {
      m_n1 = n1;
      m_n2 = n2;
      _TGridImpl::resize(n1 * n2);
    }

    geometry_utils::Box2D getBoundingBox() const
    {
      return m_bbox;
    }

    void setBoundingBox(const geometry_utils::Box2D& box)
    {
      m_bbox = box;
    }
  };

  template< class T, class allocator = std::allocator<T> >
  class Grid3D : public Grid_impl<T, allocator>
  {
    typedef Grid_impl<T, allocator> _TGridImpl;
    typedef typename _TGridImpl::_size_type _size_type;
    typedef typename allocator::reference _reference;
    typedef typename allocator::const_reference _const_reference;
    _size_type m_n1, m_n2, m_n3;
    geometry_utils::Box3D m_bbox;
  public:

    typedef Grid3D<T, allocator> _TGrid;
    typedef T _DataType;

    Grid3D()
      : m_n1(0), m_n2(0), m_n3(0)
    {
    }

    Grid3D(_size_type n1, _size_type n2, _size_type n3, geometry_utils::Box3D box = geometry_utils::Box3D())
      : _TGridImpl(n1 * n2 * n3), m_n1(n1), m_n2(n2), m_n3(n3), m_bbox(box)
    {
    }

    Grid3D(const _TGrid& another)
      :  _TGridImpl(another), m_n1(another.m_n1), m_n2(another.m_n2),
         m_n3(another.m_n3), m_bbox(another.m_bbox)
    {
    }

    _TGrid& operator= (const _TGrid& another)
    {
      if (this == &another)
        return *this;

      if (m_n1 != another.m_n1 || m_n2 != another.m_n2 || m_n3 != another.m_n3)
        throw array_size_error("grids in the assignment must have the same sizes");

      m_bbox = another.m_bbox;
      _TGridImpl::operator= (another);

      return *this;
    }

    bool operator== (const _TGrid& another) const
    {
      if (m_n1 != another.m_n1 || m_n2 != another.m_n2 || m_n3 != another.m_n3)
        throw array_size_error("logical operators on grids may be applied only to grids having the same sizes");

      if (m_bbox != another.m_bbox) {
        return false; // if domains bounding boxes are different
      }

      return _TGridImpl::operator== (another);
    }

    bool operator!= (const _TGrid& another) const
    {
      return !this->operator== (another);
    }

    _const_reference operator() (_size_type i, _size_type j, _size_type k) const
    {
      assert(i < m_n1 && j < m_n2 && k < m_n3);
      return _TGridImpl::m_data[i + m_n1 * (j + m_n2 * k)];
    }

    _reference operator() (_size_type i, _size_type j, _size_type k)
    {
      assert(i < m_n1 && j < m_n2 && k < m_n3);
      return _TGridImpl::m_data[i + m_n1 * (j + m_n2 * k)];
    }

    void swap(_TGrid& another)
    {
      _TGridImpl::swap(another);
      std::swap(m_n1, another.m_n1);
      std::swap(m_n2, another.m_n2);
      std::swap(m_n3, another.m_n3);
      std::swap(m_bbox, another.m_bbox);
    }

    _size_type size(_size_type index) const
    {
      INDEX_ASSERT(index <= 3);
      return index == 0 ? m_n1 : index == 1 ? m_n2 : m_n3;
    }

    void resize(_size_type n1, _size_type n2, _size_type n3)
    {
      m_n1 = n1;
      m_n2 = n2;
      m_n3 = n3;
      _TGridImpl::resize(n1 * n2 * n3);
    }

    geometry_utils::Box3D getBoundingBox() const
    {
      return m_bbox;
    }

    void setBoundingBox(const geometry_utils::Box3D& box)
    {
      m_bbox = box;
    }
  };
}

namespace std {
  template<typename _Tp, typename _Alloc>
  inline void
  swap(ls::Grid2D<_Tp, _Alloc>& __x, ls::Grid2D<_Tp, _Alloc>& __y)
  { __x.swap(__y); }

  template<typename _Tp, typename _Alloc>
  inline void
  swap(ls::Grid3D<_Tp, _Alloc>& __x, ls::Grid3D<_Tp, _Alloc>& __y)
  { __x.swap(__y); }
}

#endif /* GRID_H_ */
