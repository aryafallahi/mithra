// database.h
//
#ifndef database_h
#define database_h

#include <fstream>
#include <string>
#include <vector>

#include "fieldvector.h"
#include "stdinclude.h"

namespace Darius
{
  struct BunchInitialize
  {
    std::string bunchType_;
    std::string distribution_;
    unsigned int numberOfParticles_;
    Double cloudCharge_;
    Double initialGamma_;
    Double initialBeta_;
    FieldVector <Double> initialDirection_;
    std::vector <FieldVector<Double> > position_;
    FieldVector <unsigned int> numbers_;
    FieldVector <Double> latticeConstants_;
    FieldVector <Double> sigmaPosition_;
    FieldVector <Double> sigmaGammaBeta_;
    Double tranTrun_;
    Double longTrun_;
    std::string fileName_;
    Double lambda_;
    Double bF_;
    Double bFP_;
    bool shotNoise_;
    FieldVector <Double> betaVector_;
    BunchInitialize ();
  };
}
namespace Darius
{
  template <typename T>
  class AdvanceField
  {
  public:
    void advanceMagneticPotentialNSFD (T * v0, T * v1, T * v2, T * v3, T * v31, T * v32, T * v4, T * v41, T * v42, T * v5, T * v51, T * v52, T * v6, T * v61, T * v62, T * v7, T * v8, T * v9);
    void advanceScalarPotentialNSFD (T * v0, T * v1, T * v2, T * v3, T * v31, T * v32, T * v4, T * v41, T * v42, T * v5, T * v51, T * v52, T * v6, T * v61, T * v62, T * v7, T * v8, T * v9);
    void advanceMagneticPotentialFD (T * v0, T * v1, T * v2, T * v3, T * v31, T * v32, T * v4, T * v41, T * v42, T * v5, T * v51, T * v52, T * v6, T * v61, T * v62, T * v7, T * v8, T * v9);
    void advanceScalarPotentialFD (T * v0, T * v1, T * v2, T * v3, T * v31, T * v32, T * v4, T * v41, T * v42, T * v5, T * v51, T * v52, T * v6, T * v61, T * v62, T * v7, T * v8, T * v9);
    void advanceBoundaryF (T * v0, T * v1, T * v2, T * v3, T * v4, T * v5, T * v6, T * v7, T * v8, T * v9, T * v10, T * v11, T * v12, T * v13);
    void advanceBoundaryS (T * v0, T * v1, T * v2, T * v3, T * v4, T * v5, T * v6, T * v7, T * v8, T * v9, T * v10, T * v11, T * v12, T * v13);
    void advanceEdgeF (T * v0, T * v1, T * v2, T * v3, T * v4, T * v5, T * v6, T * v7, T * v8, T * v9, T * v10, T * v11, T * v12, T * v13, T * v14, T * v15, T * v16, T * v17, T * v18, T * v19);
    void advanceEdgeS (T * v0, T * v1, T * v2, T * v3, T * v4, T * v5, T * v6, T * v7, T * v8, T * v9, T * v10, T * v11, T * v12, T * v13, T * v14, T * v15, T * v16, T * v17, T * v18, T * v19);
    void advanceCornerF (T * v0, T * v1, T * v2, T * v3, T * v4, T * v5, T * v6, T * v7, T * v8, T * v9, T * v10, T * v11, T * v12, T * v13, T * v14, T * v15, T * v16, T * v17, T * v18, T * v19, T * v20, T * v21, T * v22, T * v23);
    void advanceCornerS (T * v0, T * v1, T * v2, T * v3, T * v4, T * v5, T * v6, T * v7, T * v8, T * v9, T * v10, T * v11, T * v12, T * v13, T * v14, T * v15, T * v16, T * v17, T * v18, T * v19, T * v20, T * v21, T * v22, T * v23);
    T * ufa_;
    T * ufB_;
    T alpha_;
    T beta_;
  };
}
namespace Darius
{
  template <typename T>
  struct UpdateField
  {
    T jw;
    T jwt;
    T jb;
    Double w;
    Double b;
    Double dt;
    Double dx;
    Double dy;
    Double dz;
    Double dx2;
    Double dy2;
    Double dz2;
    T (a) [6];
    T (bB) [5];
    T (cB) [5];
    T (dB) [5];
    T (eE) [5];
    T (fE) [5];
    T (gE) [5];
    T (hC) [17];
    AdvanceField <Double> af;
    Double * v0;
    Double * v1;
    Double * v2;
    Double * v3;
    Double * v31;
    Double * v32;
    Double * v4;
    Double * v41;
    Double * v42;
    Double * v5;
    Double * v51;
    Double * v52;
    Double * v6;
    Double * v61;
    Double * v62;
    Double * v7;
    Double * v8;
    Double * v9;
    Double * v10;
    Double * v11;
    Double * v12;
    Double * v13;
    Double * v14;
    Double * v15;
    Double * v16;
    Double * v17;
    Double * v18;
    Double * v19;
    Double * anp1;
    Double * an;
    Double * anm1;
    Double * fnp1;
    Double * fn;
    Double * fnm1;
    Double * jn;
    Double * rn;
    Double * en;
    Double * bn;
    unsigned int N0m1;
    unsigned int N1m1;
    unsigned int npm1;
  };
}
namespace Darius
{
  template <typename T>
  struct SampleField
  {
    std::ofstream * file;
    FieldVector <Double> position;
    int i;
    int j;
    int k;
    int m;
    FieldVector <T> et;
    FieldVector <T> at;
    FieldVector <T> bt;
    FieldVector <T> jt;
    T f;
    T q;
    Double dxr;
    Double dyr;
    Double dzr;
    T c2;
    T jw;
    T p;
    double c1;
    unsigned int N;
  };
}
namespace Darius
{
  template <typename T>
  struct VisualizeField
  {
    std::ofstream * file;
    std::string fileName;
    std::string path;
    std::string name;
    unsigned int i;
    unsigned int j;
    unsigned int k;
    unsigned int m;
    std::vector <std::vector<T> > v;
    T jw;
  };
}
namespace Darius
{
  template <typename T>
  struct ProfileField
  {
    unsigned int i;
    unsigned int j;
    unsigned int k;
    unsigned int l;
    unsigned int m;
    std::ofstream * file;
    std::string fileName;
    FieldVector <T> e;
    FieldVector <T> a;
    FieldVector <T> b;
    Double dt;
    T jw;
    int kmin;
    int kmax;
  };
}
namespace Darius
{
  struct UpdateBunchParallel
  {
    int i;
    int j;
    int k;
    int n;
    long int m;
    bool b1;
    bool b2;
    std::vector <Double> qSB;
    std::vector <Double> qSF;
    std::vector <Double> qRB;
    std::vector <Double> qRF;
    FieldVector <Double> et;
    FieldVector <Double> bt;
    Double dxr;
    Double dyr;
    Double dzr;
    Double d2;
    double d1;
    FieldVector <Double> gbm;
    FieldVector <Double> gbp;
    FieldVector <Double> gbpl;
    Double sz;
    Double cy;
    Double lz;
    Double ly;
    Double tsignal;
    Double d;
    Double l;
    Double zRp;
    Double wrp;
    Double zRs;
    Double wrs;
    Double x;
    Double y;
    Double z;
    Double p;
    Double t;
    FieldVector <Double> rv;
    FieldVector <Double> yv;
    FieldVector <Double> rl;
    FieldVector <Double> ex;
    FieldVector <Double> ez;
    FieldVector <Double> eT;
    FieldVector <Double> by;
    FieldVector <Double> bz;
    FieldVector <Double> bT;
    Double tl;
    Double t0;
    Double r0;
    UpdateBunchParallel ();
  };
}
namespace Darius
{
  struct UpdateBunch
  {
    Double dx;
    Double dy;
    Double dz;
    Double dt;
    Double dtb;
    Double ku;
    Double b0;
    Double r1;
    Double r2;
    Double ct;
    Double st;
    int nL;
    int i;
    Charge Q;
  };
}
namespace Darius
{
  struct SampleBunch
  {
    std::ofstream * file;
    Double g;
    Double q;
    Double qT;
    FieldVector <Double> r;
    FieldVector <Double> gb;
    FieldVector <Double> r2;
    FieldVector <Double> gb2;
    FieldVector <Double> rT;
    FieldVector <Double> gbT;
    FieldVector <Double> r2T;
    FieldVector <Double> gb2T;
  };
}
namespace Darius
{
  struct VisualizeBunch
  {
    std::ofstream * file;
    std::string fileName;
    std::string path;
    std::string name;
    unsigned int i;
    unsigned int N;
  };
}
namespace Darius
{
  struct ProfileBunch
  {
    std::ofstream * file;
    std::string fileName;
  };
}
namespace Darius
{
  struct UpdateCurrent
  {
    int ip;
    int jp;
    int kp;
    int im;
    int jm;
    int km;
    long int m;
    Double dx;
    Double dy;
    Double dz;
    Double dxp;
    Double dyp;
    Double dzp;
    Double dxm;
    Double dym;
    Double dzm;
    Double x1;
    Double x2;
    Double y1;
    Double y2;
    Double z1;
    Double z2;
    Double q;
    FieldVector <Double> rp;
    FieldVector <Double> rm;
    FieldVector <Double> r;
    FieldVector <Double> jcp;
    FieldVector <Double> jcm;
    Double rc;
    Double dv;
    double c;
    std::vector <FieldVector<Double> > jt;
    std::vector <Double> rt;
  };
}
namespace Darius
{
  struct SampleRadiationPower
  {
    std::vector <std::ofstream*> file;
    std::vector <Double> w;
    int k;
    int l;
    Double dxr;
    Double dyr;
    Double dzr;
    double c;
    std::vector <std::vector<std::vector<Double> > > pLi;
    std::vector <Double> pL;
    std::vector <Double> pG;
    Double dt;
    Double dx;
    Double dy;
    Double dz;
    Double pc;
    unsigned int N;
    unsigned int Nl;
    unsigned int Nf;
    unsigned int Nz;
    std::vector <std::vector<std::vector<Double> > > fdt;
    unsigned int m;
    std::vector <std::vector<Complex> > ep;
    std::vector <std::vector<Complex> > em;
  };
}
namespace Darius
{
  struct SampleRadiationEnergy
  {
    std::vector <std::ofstream*> file;
    std::vector <Double> w;
    int k;
    int l;
    Double dxr;
    Double dyr;
    Double dzr;
    double c;
    std::vector <std::vector<std::vector<Double> > > pLi;
    std::vector <Double> pL;
    std::vector <Double> pG;
    Double dt;
    Double dx;
    Double dy;
    Double dz;
    Double pc;
    unsigned int N;
    unsigned int Nl;
    unsigned int Nf;
    std::vector <std::vector<std::vector<Double> > > fdt;
  };
}
namespace Darius
{
  template <typename T>
  void AdvanceField <T>::advanceMagneticPotentialNSFD (T * v0, T * v1, T * v2, T * v3, T * v31, T * v32, T * v4, T * v41, T * v42, T * v5, T * v51, T * v52, T * v6, T * v61, T * v62, T * v7, T * v8, T * v9)
    {
      *v0 =
	  *ufa_     * *v2 - *v1 + alpha_ * (
	  *(ufa_+1) * ( *v3 + *v4 + beta_ * ( *v31 + *v32 + *v41 + *v42 ) )   +
	  *(ufa_+2) * ( *v5 + *v6 + beta_ * ( *v51 + *v52 + *v61 + *v62 ) ) ) +
	  *(ufa_+3) * ( *v7 + *v8 ) +
	  *(ufa_+4) * *v9;
      *(v0+1) =
	  *ufa_     * *(v2+1) - *(v1+1) + alpha_ * (
	  *(ufa_+1) * ( *(v3+1) + *(v4+1) + beta_ * ( *(v31+1) + *(v32+1) + *(v41+1) + *(v42+1) ) )   +
	  *(ufa_+2) * ( *(v5+1) + *(v6+1) + beta_ * ( *(v51+1) + *(v52+1) + *(v61+1) + *(v62+1) ) ) ) +
	  *(ufa_+3) * ( *(v7+1) + *(v8+1) ) +
	  *(ufa_+4) * *(v9+1);
      *(v0+2) =
	  *ufa_     * *(v2+2) - *(v1+2) + alpha_ * (
	  *(ufa_+1) * ( *(v3+2) + *(v4+2) + beta_ * ( *(v31+2) + *(v32+2) + *(v41+2) + *(v42+2) ) )   +
	  *(ufa_+2) * ( *(v5+2) + *(v6+2) + beta_ * ( *(v51+2) + *(v52+2) + *(v61+2) + *(v62+2) ) ) ) +
	  *(ufa_+3) * ( *(v7+2) + *(v8+2) ) +
	  *(ufa_+4) * *(v9+2);
    }
}
namespace Darius
{
  template <typename T>
  void AdvanceField <T>::advanceScalarPotentialNSFD (T * v0, T * v1, T * v2, T * v3, T * v31, T * v32, T * v4, T * v41, T * v42, T * v5, T * v51, T * v52, T * v6, T * v61, T * v62, T * v7, T * v8, T * v9)
    {
      *v0 =
	  *ufa_     * *v2 - *v1 + alpha_ * (
	  *(ufa_+1) * ( *v3 + *v4 + beta_ * ( *v31 + *v32 + *v41 + *v42 ) )   +
	  *(ufa_+2) * ( *v5 + *v6 + beta_ * ( *v51 + *v52 + *v61 + *v62 ) ) ) +
	  *(ufa_+3) * ( *v7 + *v8 ) +
	  *(ufa_+5) * *v9;
    }
}
namespace Darius
{
  template <typename T>
  void AdvanceField <T>::advanceMagneticPotentialFD (T * v0, T * v1, T * v2, T * v3, T * v31, T * v32, T * v4, T * v41, T * v42, T * v5, T * v51, T * v52, T * v6, T * v61, T * v62, T * v7, T * v8, T * v9)
    {
      *v0 =
	  *ufa_     * *v2 - *v1 +
	  *(ufa_+1) * ( *v3 + *v4 ) +
	  *(ufa_+2) * ( *v5 + *v6 ) +
	  *(ufa_+3) * ( *v7 + *v8 ) +
	  *(ufa_+4) * *v9;
      *(v0+1) =
	  *ufa_     * *(v2+1) - *(v1+1) +
	  *(ufa_+1) * ( *(v3+1) + *(v4+1) ) +
	  *(ufa_+2) * ( *(v5+1) + *(v6+1) ) +
	  *(ufa_+3) * ( *(v7+1) + *(v8+1) ) +
	  *(ufa_+4) * *(v9+1);
      *(v0+2) =
	  *ufa_     * *(v2+2) - *(v1+2) +
	  *(ufa_+1) * ( *(v3+2) + *(v4+2) ) +
	  *(ufa_+2) * ( *(v5+2) + *(v6+2) ) +
	  *(ufa_+3) * ( *(v7+2) + *(v8+2) ) +
	  *(ufa_+4) * *(v9+2);
    }
}
namespace Darius
{
  template <typename T>
  void AdvanceField <T>::advanceScalarPotentialFD (T * v0, T * v1, T * v2, T * v3, T * v31, T * v32, T * v4, T * v41, T * v42, T * v5, T * v51, T * v52, T * v6, T * v61, T * v62, T * v7, T * v8, T * v9)
    {
      *v0 =
	  *ufa_     * *v2 - *v1 +
	      *(ufa_+1) * ( *v3 + *v4 ) +
	      *(ufa_+2) * ( *v5 + *v6 ) +
	      *(ufa_+3) * ( *v7 + *v8 ) +
	      *(ufa_+5) * *v9;
    }
}
namespace Darius
{
  template <typename T>
  void AdvanceField <T>::advanceBoundaryF (T * v0, T * v1, T * v2, T * v3, T * v4, T * v5, T * v6, T * v7, T * v8, T * v9, T * v10, T * v11, T * v12, T * v13)
    {
      *v0 =
	  *ufB_     * ( *v1 + *v5 ) +
	  *(ufB_+1) * *v3 +
	  *(ufB_+2) * ( *v2 + *v4 ) +
	  *(ufB_+3) * ( *v6 + *v7 + *v10 + *v11 ) +
	  *(ufB_+4) * ( *v8 + *v9 + *v12 + *v13 );
      *(v0+1) =
	  *ufB_     * ( *(v1+1) + *(v5+1) ) +
	  *(ufB_+1) * *(v3+1) +
	  *(ufB_+2) * ( *(v2+1) + *(v4+1) ) +
	  *(ufB_+3) * ( *(v6+1) + *(v7+1) + *(v10+1) + *(v11+1) ) +
	  *(ufB_+4) * ( *(v8+1) + *(v9+1) + *(v12+1) + *(v13+1) );
      *(v0+2) =
	  *ufB_     * ( *(v1+2) + *(v5+2) ) +
	  *(ufB_+1) * *(v3+2) +
	  *(ufB_+2) * ( *(v2+2) + *(v4+2) ) +
	  *(ufB_+3) * ( *(v6+2) + *(v7+2) + *(v10+2) + *(v11+2) ) +
	  *(ufB_+4) * ( *(v8+2) + *(v9+2) + *(v12+2) + *(v13+2) );
    }
}
namespace Darius
{
  template <typename T>
  void AdvanceField <T>::advanceBoundaryS (T * v0, T * v1, T * v2, T * v3, T * v4, T * v5, T * v6, T * v7, T * v8, T * v9, T * v10, T * v11, T * v12, T * v13)
    {
      *v0 =
	  *ufB_     * ( *v1 + *v5 ) +
	  *(ufB_+1) * *v3 +
	  *(ufB_+2) * ( *v2 + *v4 ) +
	  *(ufB_+3) * ( *v6 + *v7 + *v10 + *v11 ) +
	  *(ufB_+4) * ( *v8 + *v9 + *v12 + *v13 );
    }
}
namespace Darius
{
  template <typename T>
  void AdvanceField <T>::advanceEdgeF (T * v0, T * v1, T * v2, T * v3, T * v4, T * v5, T * v6, T * v7, T * v8, T * v9, T * v10, T * v11, T * v12, T * v13, T * v14, T * v15, T * v16, T * v17, T * v18, T * v19)
    {
      *v0     = *ufB_ *     ( *v3 + *v8 ) +
		*(ufB_+1) * ( *v5 + *v6 ) +
		*(ufB_+2) * ( *v2 + *v9 ) +
		*(ufB_+3) * ( *v1 + *v4 + *v7 + *v10 ) +
		*(ufB_+4) * ( *v12 + *v13 + *v14 + *v15 + *v16 + *v17 + *v18 + *v19 ) -
		*v11;

      *(v0+1) = *ufB_ *     ( *(v3+1) + *(v8+1) ) +
		*(ufB_+1) * ( *(v5+1) + *(v6+1) ) +
		*(ufB_+2) * ( *(v2+1) + *(v9+1) ) +
		*(ufB_+3) * ( *(v1+1) + *(v4+1) + *(v7+1) + *(v10+1) ) +
		*(ufB_+4) * ( *(v12+1)+ *(v13+1)+ *(v14+1)+ *(v15+1)+ *(v16+1)+ *(v17+1)+ *(v18+1)+ *(v19+1)) -
		*(v11+1);

      *(v0+2) = *ufB_ *     ( *(v3+2) + *(v8+2) ) +
      		*(ufB_+1) * ( *(v5+2) + *(v6+2) ) +
      		*(ufB_+2) * ( *(v2+2) + *(v9+2) ) +
      		*(ufB_+3) * ( *(v1+2) + *(v4+2) + *(v7+2) + *(v10+2) ) +
      		*(ufB_+4) * ( *(v12+2)+ *(v13+2)+ *(v14+2)+ *(v15+2)+ *(v16+2)+ *(v17+2)+ *(v18+2)+ *(v19+2)) -
      		*(v11+2);

    }
}
namespace Darius
{
  template <typename T>
  void AdvanceField <T>::advanceEdgeS (T * v0, T * v1, T * v2, T * v3, T * v4, T * v5, T * v6, T * v7, T * v8, T * v9, T * v10, T * v11, T * v12, T * v13, T * v14, T * v15, T * v16, T * v17, T * v18, T * v19)
    {
      *v0     = *ufB_ *     ( *v3 + *v8 ) +
		*(ufB_+1) * ( *v5 + *v6 ) +
		*(ufB_+2) * ( *v2 + *v9 ) +
		*(ufB_+3) * ( *v1 + *v4 + *v7 + *v10 ) +
		*(ufB_+4) * ( *v12 + *v13 + *v14 + *v15 + *v16 + *v17 + *v18 + *v19 ) -
		*v11;
    }
}
namespace Darius
{
  template <typename T>
  void AdvanceField <T>::advanceCornerF (T * v0, T * v1, T * v2, T * v3, T * v4, T * v5, T * v6, T * v7, T * v8, T * v9, T * v10, T * v11, T * v12, T * v13, T * v14, T * v15, T * v16, T * v17, T * v18, T * v19, T * v20, T * v21, T * v22, T * v23)
    {
      *v0     = - ( *v1  * *(ufB_+16) + *v2  * *(ufB_+8) +
		    *v3  * *(ufB_+1) + *v4  * *(ufB_+16) + *v5  * *(ufB_+9) +
		    *v6  * *(ufB_+2) + *v7  * *(ufB_+16) + *v8  * *(ufB_+10) +
		    *v9  * *(ufB_+3) + *v10 * *(ufB_+16) + *v11 * *(ufB_+11) +
		    *v12 * *(ufB_+4) + *v13 * *(ufB_+16) + *v14 * *(ufB_+12) +
		    *v15 * *(ufB_+5) + *v16 * *(ufB_+16) + *v17 * *(ufB_+13) +
		    *v18 * *(ufB_+6) + *v19 * *(ufB_+16) + *v20 * *(ufB_+14) +
		    *v21 * *(ufB_+7) + *v22 * *(ufB_+16) + *v23 * *(ufB_+15)) / *ufB_;

      *(v0+1) = - ( *(v1+1)  * *(ufB_+16) + *(v2+1)  * *(ufB_+8) +
      		    *(v3+1)  * *(ufB_+1) + *(v4+1)  * *(ufB_+16) + *(v5+1)  * *(ufB_+9) +
      		    *(v6+1)  * *(ufB_+2) + *(v7+1)  * *(ufB_+16) + *(v8+1)  * *(ufB_+10) +
      		    *(v9+1)  * *(ufB_+3) + *(v10+1) * *(ufB_+16) + *(v11+1) * *(ufB_+11) +
      		    *(v12+1) * *(ufB_+4) + *(v13+1) * *(ufB_+16) + *(v14+1) * *(ufB_+12) +
      		    *(v15+1) * *(ufB_+5) + *(v16+1) * *(ufB_+16) + *(v17+1) * *(ufB_+13) +
      		    *(v18+1) * *(ufB_+6) + *(v19+1) * *(ufB_+16) + *(v20+1) * *(ufB_+14) +
      		    *(v21+1) * *(ufB_+7) + *(v22+1) * *(ufB_+16) + *(v23+1) * *(ufB_+15)) / *ufB_;

      *(v0+2) = - ( *(v1+2)  * *(ufB_+16) + *(v2+2)  * *(ufB_+8) +
		    *(v3+2)  * *(ufB_+1) + *(v4+2)  * *(ufB_+16) + *(v5+2)  * *(ufB_+9) +
		    *(v6+2)  * *(ufB_+2) + *(v7+2)  * *(ufB_+16) + *(v8+2)  * *(ufB_+10) +
		    *(v9+2)  * *(ufB_+3) + *(v10+2) * *(ufB_+16) + *(v11+2) * *(ufB_+11) +
		    *(v12+2) * *(ufB_+4) + *(v13+2) * *(ufB_+16) + *(v14+2) * *(ufB_+12) +
		    *(v15+2) * *(ufB_+5) + *(v16+2) * *(ufB_+16) + *(v17+2) * *(ufB_+13) +
		    *(v18+2) * *(ufB_+6) + *(v19+2) * *(ufB_+16) + *(v20+2) * *(ufB_+14) +
		    *(v21+2) * *(ufB_+7) + *(v22+2) * *(ufB_+16) + *(v23+2) * *(ufB_+15)) / *ufB_;

    }
}
namespace Darius
{
  template <typename T>
  void AdvanceField <T>::advanceCornerS (T * v0, T * v1, T * v2, T * v3, T * v4, T * v5, T * v6, T * v7, T * v8, T * v9, T * v10, T * v11, T * v12, T * v13, T * v14, T * v15, T * v16, T * v17, T * v18, T * v19, T * v20, T * v21, T * v22, T * v23)
    {
      *v0     = - ( *v1  * *(ufB_+16) + *v2  * *(ufB_+8) +
	  *v3  * *(ufB_+1) + *v4  * *(ufB_+16) + *v5  * *(ufB_+9) +
	  *v6  * *(ufB_+2) + *v7  * *(ufB_+16) + *v8  * *(ufB_+10) +
	  *v9  * *(ufB_+3) + *v10 * *(ufB_+16) + *v11 * *(ufB_+11) +
	  *v12 * *(ufB_+4) + *v13 * *(ufB_+16) + *v14 * *(ufB_+12) +
	  *v15 * *(ufB_+5) + *v16 * *(ufB_+16) + *v17 * *(ufB_+13) +
	  *v18 * *(ufB_+6) + *v19 * *(ufB_+16) + *v20 * *(ufB_+14) +
	  *v21 * *(ufB_+7) + *v22 * *(ufB_+16) + *v23 * *(ufB_+15)) / *ufB_;
    }
}
#endif
