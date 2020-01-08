// fieldvector.h
//

#ifndef fieldvector_h
#define fieldvector_h

#include <complex>
#include <mpi.h>
#include <sstream>
#include <vector>

namespace Darius
{
  typedef double Double;

  typedef std::complex <double> Complex;

  extern MPI_Datatype MPI_TYPE;

  template <typename ElementType>
  class FieldVector
  {
  public:
    ElementType (a) [3];
    FieldVector ();
    template <typename T>
    FieldVector (T x);
    ElementType & operator [] (unsigned int n);
    ElementType norm ();
    template <typename T1, typename T2>
    void mv (T1 y, FieldVector <T2> & x);
    template <typename T1, typename T2>
    void pmv (T1 y, FieldVector <T2> & x);
    template <typename T1, typename T2>
    void mmv (T1 y, FieldVector <T2> & x);
    template <typename T1, typename T2>
    void dv (T1 y, FieldVector <T2> & x);
    template <typename T1, typename T2>
    void pdv (T1 y, FieldVector <T2> & x);
    template <typename T1, typename T2>
    void mdv (T1 y, FieldVector <T2> & x);
    void operator = (ElementType & y);
    template <typename T>
    void operator = (FieldVector <T> & y);
    template <typename T>
    void operator = (std::vector <T> & y);
  };

  template <typename T1, typename T2>
  void operator += (FieldVector <T1> & x, T2 & y);

  template <typename T1, typename T2>
  void operator -= (FieldVector <T1> & x, T2 & y);

  template <typename T1, typename T2>
  void operator *= (FieldVector <T1> & x, T2 & y);

  template <typename T1, typename T2>
  void operator /= (FieldVector <T1> & x, T2 & y);

  template <typename T1, typename T2>
  void operator += (FieldVector <T1> & x, FieldVector <T2> & y);

  template <typename T1, typename T2>
  void operator -= (FieldVector <T1> & x, FieldVector <T2> & y);

  template <typename T1, typename T2>
  T1 operator * (FieldVector <T1> & x, FieldVector <T2> & y);

  template <typename T1, typename T2>
  FieldVector <T1> cross (FieldVector <T1> & x, FieldVector <T2> & y);

  template <typename T>
  void operator << (std::ostringstream & in, FieldVector <T> & y);

  template <typename ElementType>
  FieldVector <ElementType>::FieldVector ()
    {
      a[0] = 0.0;
      a[1] = 0.0;
      a[2] = 0.0;
    }

  template <typename ElementType>
  template <typename T>
  FieldVector <ElementType>::FieldVector (T x)
    {
      a[0] = x;
      a[1] = x;
      a[2] = x;
    }

  template <typename ElementType>
  ElementType & FieldVector <ElementType>::operator [] (unsigned int n)
    {
      return this->a[n];
    }

  template <typename ElementType>
  ElementType FieldVector <ElementType>::norm ()
    {
      return ( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );
    }

  template <typename ElementType>
  template <typename T1, typename T2>
  void FieldVector <ElementType>::mv (T1 y, FieldVector <T2> & x)
    {
      a[0] = y * x[0];
      a[1] = y * x[1];
      a[2] = y * x[2];
    }

  template <typename ElementType>
  template <typename T1, typename T2>
  void FieldVector <ElementType>::pmv (T1 y, FieldVector <T2> & x)
    {
      a[0] += y * x[0];
      a[1] += y * x[1];
      a[2] += y * x[2];
    }

  template <typename ElementType>
  template <typename T1, typename T2>
  void FieldVector <ElementType>::mmv (T1 y, FieldVector <T2> & x)
    {
      a[0] -= y * x[0];
      a[1] -= y * x[1];
      a[2] -= y * x[2];
    }

  template <typename ElementType>
  template <typename T1, typename T2>
  void FieldVector <ElementType>::dv (T1 y, FieldVector <T2> & x)
    {
      a[0] = x[0] / y;
      a[1] = x[1] / y;
      a[2] = x[2] / y;
    }

  template <typename ElementType>
  template <typename T1, typename T2>
  void FieldVector <ElementType>::pdv (T1 y, FieldVector <T2> & x)
    {
      a[0] += x[0] / y;
      a[1] += x[1] / y;
      a[2] += x[2] / y;
    }

  template <typename ElementType>
  template <typename T1, typename T2>
  void FieldVector <ElementType>::mdv (T1 y, FieldVector <T2> & x)
    {
      a[0] -= x[0] / y;
      a[1] -= x[1] / y;
      a[2] -= x[2] / y;
    }

  template <typename ElementType>
  void FieldVector <ElementType>::operator = (ElementType & y)
    {
      a[0] = y;
      a[1] = y;
      a[2] = y;
    }

  template <typename ElementType>
  template <typename T>
  void FieldVector <ElementType>::operator = (FieldVector <T> & y)
    {
      a[0] = y[0];
      a[1] = y[1];
      a[2] = y[2];
    }

  template <typename ElementType>
  template <typename T>
  void FieldVector <ElementType>::operator = (std::vector <T> & y)
    {
      if (y.size() != 3)
	{
	  std::cout << std::string("A field-vector with other than three parameters is not accepted.") << std::endl;
	  exit(1);
	}
      a[0] = y[0];
      a[1] = y[1];
      a[2] = y[2];
    }

  template <typename T1, typename T2>
  void operator += (FieldVector <T1> & x, T2 & y)
  {
    x[0] += y;
    x[1] += y;
    x[2] += y;
  }

  template <typename T1, typename T2>
  void operator -= (FieldVector <T1> & x, T2 & y)
  {
    x[0] -= y;
    x[1] -= y;
    x[2] -= y;
  }

  template <typename T1, typename T2>
  void operator *= (FieldVector <T1> & x, T2 & y)
  {
    x[0] *= y;
    x[1] *= y;
    x[2] *= y;
  }

  template <typename T1, typename T2>
  void operator /= (FieldVector <T1> & x, T2 & y)
  {
    x[0] /= y;
    x[1] /= y;
    x[2] /= y;
  }

  template <typename T1, typename T2>
  void operator += (FieldVector <T1> & x, FieldVector <T2> & y)
  {
    x[0] += y[0];
    x[1] += y[1];
    x[2] += y[2];
  }

  template <typename T1, typename T2>
  void operator -= (FieldVector <T1> & x, FieldVector <T2> & y)
  {
    x[0] -= y[0];
    x[1] -= y[1];
    x[2] -= y[2];
  }

  template <typename T1, typename T2>
  T1 operator * (FieldVector <T1> & x, FieldVector <T2> & y)
  {
    return ( x[0]*y[0] + x[1]*y[1] + x[2]*y[2] );
  }

  template <typename T1, typename T2>
  FieldVector <T1> cross (FieldVector <T1> & x, FieldVector <T2> & y)
  {
    FieldVector<T1> z (0.0);
    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = x[2] * y[0] - x[0] * y[2];
    z[2] = x[0] * y[1] - x[1] * y[0];
    return (z);
  }

  template <typename T>
  inline void operator << (std::ostringstream & in, FieldVector <T> & y)
  {
    in << "(" << y[0] << "," << y[1] << "," << y[2] << ")" ;
  }
}
#endif
