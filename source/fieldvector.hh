/********************************************************************************************************
 *  fieldvector.hh : Implementation of the field vector class for the analysis.
 ********************************************************************************************************/

#ifndef FIELDVECTOR_HH_
#define FIELDVECTOR_HH_

namespace Darius
{

  /* Define the types of values used in the analysis.							*/
  typedef double			Double;
  typedef std::complex<double>		Complex;
  MPI_Datatype 				MPI_TYPE = MPI_DOUBLE;

  /** FieldVector is a class constructed to be used for vectors of electric and magnetic fields.     	*/
  template <class ElementType>
  class FieldVector {
  public:
    ElementType a [3];

    FieldVector ()
    {
      a[0] = 0.0;
      a[1] = 0.0;
      a[2] = 0.0;
    };

    template <class T>
    FieldVector (T x)
    {
      a[0] = x;
      a[1] = x;
      a[2] = x;
    }

    ElementType& operator [] (unsigned int n)
    {
      return this->a[n];
    };

    ElementType norm()
    {
      return ( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] );
    };

    template <class T1, class T2>
    void mv(T1 y, FieldVector<T2>& x)
    {
      a[0] = y * x[0];
      a[1] = y * x[1];
      a[2] = y * x[2];
    }

    template <class T1, class T2>
    void pmv(T1 y, FieldVector<T2>& x)
    {
      a[0] += y * x[0];
      a[1] += y * x[1];
      a[2] += y * x[2];
    }

    template <class T1, class T2>
    void mmv(T1 y, FieldVector<T2>& x)
    {
      a[0] -= y * x[0];
      a[1] -= y * x[1];
      a[2] -= y * x[2];
    }

    template <class T1, class T2>
    void dv(T1 y, FieldVector<T2>& x)
    {
      a[0] = x[0] / y;
      a[1] = x[1] / y;
      a[2] = x[2] / y;
    }

    template <class T1, class T2>
    void pdv(T1 y, FieldVector<T2>& x)
    {
      a[0] += x[0] / y;
      a[1] += x[1] / y;
      a[2] += x[2] / y;
    }

    template <class T1, class T2>
    void mdv(T1 y, FieldVector<T2>& x)
    {
      a[0] -= x[0] / y;
      a[1] -= x[1] / y;
      a[2] -= x[2] / y;
    }

    void operator= (ElementType& y)
    {
      a[0] = y;
      a[1] = y;
      a[2] = y;
    };

    template <class T>
    void operator= (FieldVector<T>& y)
    {
      a[0] = y[0];
      a[1] = y[1];
      a[2] = y[2];
    }

    template <class T>
    void operator= (std::vector<T>& y)
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

  };

  /* Define the various operators used for the field vector class.					*/

  template <class T1, class T2>
  void operator+= (FieldVector<T1>& x, T2& y)
  {
    x[0] += y;
    x[1] += y;
    x[2] += y;
  }

  template <class T1, class T2>
  void operator-= (FieldVector<T1>& x, T2& y)
  {
    x[0] -= y;
    x[1] -= y;
    x[2] -= y;
  }

  template <class T1, class T2>
  void operator*= (FieldVector<T1>& x, T2& y)
  {
    x[0] *= y;
    x[1] *= y;
    x[2] *= y;
  }

  template <class T1, class T2>
  void operator/= (FieldVector<T1>& x, T2& y)
  {
    x[0] /= y;
    x[1] /= y;
    x[2] /= y;
  }

  template <class T1, class T2>
  void operator+= (FieldVector<T1>& x, FieldVector<T2>& y)
  {
    x[0] += y[0];
    x[1] += y[1];
    x[2] += y[2];
  }

  template <class T1, class T2>
  void operator-= (FieldVector<T1>& x, FieldVector<T2>& y)
  {
    x[0] -= y[0];
    x[1] -= y[1];
    x[2] -= y[2];
  }

  template <class T1, class T2>
  T1 operator* (FieldVector<T1>& x, FieldVector<T2>& y)
  {
    return ( x[0]*y[0] + x[1]*y[1] + x[2]*y[2] );
  }

  template <class T1, class T2>
  FieldVector<T1> cross (FieldVector<T1>& x, FieldVector<T2>& y)
  {
    FieldVector<T1> z (0.0);
    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = x[2] * y[0] - x[0] * y[2];
    z[2] = x[0] * y[1] - x[1] * y[0];
    return (z);
  }

  template <class T>
  inline void operator<< (std::ostringstream &in, FieldVector<T>& y)
  {
    in << "(" << y[0] << "," << y[1] << "," << y[2] << ")" ;
  }
}

#endif /* FIELDVECTOR_HH_ */
