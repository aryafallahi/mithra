/********************************************************************************************************
 *  database.hh : Implementation of the look up tables and data bases in darius
 ********************************************************************************************************/

#ifndef	DATABASE_P_HH_
#define	DATABASE_P_HH_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define	tab_N	8192

namespace Darius
{
  /* The structure of data for initializing a bunch.							*/
  struct BunchInitialize
  {
    /* Type of the bunch which is one of the manual, ellipsoid, cylinder, cube, and 3D-crystal. If it is
     * manual the charge at points of the position vector will be produced.    				*/
    std::string     			bunchType_;

    /* Type of the distributions (transverse or longitudinal) in the bunch.				*/
    std::string     			distribution_;

    /* Total number of macroparticles in the bunch.                                                     */
    unsigned int       			numberOfParticles_;

    /* Total charge of the bunch in pC.                                                                 */
    Double  				cloudCharge_;

    /* Initial energy of the bunch in MeV.                                                              */
    Double            			initialGamma_;

    /* Initial normalized speed of the bunch.                                                           */
    Double            			initialBeta_;

    /* Initial movement direction of the bunch, which is a unit vector.                                 */
    FieldVector<Double>			initialDirection_;

    /* Position of the center of the bunch in the unit of length scale.                           	*/
    std::vector<FieldVector<Double> >	position_;

    /* Number of macroparticles in each direction for 3Dcrystal type.                                   */
    FieldVector<unsigned int>		numbers_;

    /* Lattice constant in x, y, and z directions for 3D crystal type.                                  */
    FieldVector<Double>			latticeConstants_;

    /* Spread in position for each of the directions in the unit of length scale. For the 3D crystal
     * type, it will be the spread in position for each micro-bunch of the crystal.			*/
    FieldVector<Double>			sigmaPosition_;

    /* Spread in energy in each direction.                                                              */
    FieldVector<Double>			sigmaGammaBeta_;

    /* Store the truncation transverse distance for the electron generation.				*/
    Double				tranTrun_;

    /* Store the truncation longitudinal distance for the electron generation.				*/
    Double				longTrun_;

    /* Name of the file for reading the electrons distribution from.					*/
    std::string				fileName_;

    /* The radiation wavelength corresponding to the bunch length outside the undulator			*/
    Double				lambda_;

    /* Bunching factor for the initialization of the bunch.						*/
    Double				bF_;

    /* Phase of the bunching factor for the initialization of the bunch.				*/
    Double				bFP_;

    /* Boolean flag determining the activation of shot-noise.						*/
    bool				shotNoise_;

    /* Initial beta vector of the bunch, which is obtained as the product of beta and direction.	*/
    FieldVector<Double>			betaVector_;

    /* Initialize the parameters for the bunch initialization to some first values.                     */
    BunchInitialize()
    {
      bunchType_			= "";
      distribution_			= "";
      numberOfParticles_ 		= 0;
      cloudCharge_			= 0.0;
      initialGamma_			= 0.0;
      initialBeta_			= 0.0;
      initialDirection_			= 0.0;
      position_.clear();
      numbers_				= 0;
      latticeConstants_			= 0.0;
      sigmaPosition_			= 0.0;
      sigmaGammaBeta_			= 0.0;
      tranTrun_				= 0.0;
      longTrun_				= 0.0;
      fileName_				= "";
      bF_				= 0.0;
      bFP_				= 0.0;
      shotNoise_			= false;
      lambda_				= 0.0;
    }
  };

  /* Structure of data required for updating the field.							*/
  template <typename T> class AdvanceField
  {
  public:

    void advanceMagneticPotentialNSFD (T* v0,
				       T* v1,
				       T* v2,
				       T* v3, T* v31, T* v32,
				       T* v4, T* v41, T* v42,
				       T* v5, T* v51, T* v52,
				       T* v6, T* v61, T* v62,
				       T* v7,
				       T* v8,
				       T* v9)
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
    };

    void advanceScalarPotentialNSFD (T* v0,
				     T* v1,
				     T* v2,
				     T* v3, T* v31, T* v32,
				     T* v4, T* v41, T* v42,
				     T* v5, T* v51, T* v52,
				     T* v6, T* v61, T* v62,
				     T* v7,
				     T* v8,
				     T* v9)
    {
      *v0 =
	  *ufa_     * *v2 - *v1 + alpha_ * (
	  *(ufa_+1) * ( *v3 + *v4 + beta_ * ( *v31 + *v32 + *v41 + *v42 ) )   +
	  *(ufa_+2) * ( *v5 + *v6 + beta_ * ( *v51 + *v52 + *v61 + *v62 ) ) ) +
	  *(ufa_+3) * ( *v7 + *v8 ) +
	  *(ufa_+5) * *v9;
    }

    void advanceMagneticPotentialFD (T* v0,
				     T* v1,
				     T* v2,
				     T* v3, T* v31, T* v32,
				     T* v4, T* v41, T* v42,
				     T* v5, T* v51, T* v52,
				     T* v6, T* v61, T* v62,
				     T* v7,
				     T* v8,
				     T* v9)
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
    };

    void advanceScalarPotentialFD (T* v0,
				   T* v1,
				   T* v2,
				   T* v3, T* v31, T* v32,
				   T* v4, T* v41, T* v42,
				   T* v5, T* v51, T* v52,
				   T* v6, T* v61, T* v62,
				   T* v7,
				   T* v8,
				   T* v9)
    {
      *v0 =
	  *ufa_     * *v2 - *v1 +
	      *(ufa_+1) * ( *v3 + *v4 ) +
	      *(ufa_+2) * ( *v5 + *v6 ) +
	      *(ufa_+3) * ( *v7 + *v8 ) +
	      *(ufa_+5) * *v9;
    }

    void advanceBoundaryF(T* v0,
                          T* v1, T* v2, T* v3, T* v4, T* v5, T* v6, T* v7,
                          T* v8, T* v9, T* v10,T* v11,T* v12,T* v13)
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
    };

    void advanceBoundaryS(T* v0,
                          T* v1, T* v2, T* v3, T* v4, T* v5, T* v6, T* v7,
                          T* v8, T* v9, T* v10,T* v11,T* v12,T* v13)
    {
      *v0 =
	  *ufB_     * ( *v1 + *v5 ) +
	  *(ufB_+1) * *v3 +
	  *(ufB_+2) * ( *v2 + *v4 ) +
	  *(ufB_+3) * ( *v6 + *v7 + *v10 + *v11 ) +
	  *(ufB_+4) * ( *v8 + *v9 + *v12 + *v13 );
    };

    void advanceEdgeF (T* v0,  T* v1,  T* v2,
                       T* v3,  T* v4,  T* v5,
                       T* v6,  T* v7,  T* v8,
                       T* v9,  T* v10, T* v11,
                       T* v12, T* v13, T* v14, T* v15,
                       T* v16, T* v17, T* v18, T* v19)
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

    };

    void advanceEdgeS (T* v0,  T* v1,  T* v2,
                       T* v3,  T* v4,  T* v5,
		       T* v6,  T* v7,  T* v8,
                       T* v9,  T* v10, T* v11,
                       T* v12, T* v13, T* v14, T* v15,
                       T* v16, T* v17, T* v18, T* v19)
    {
      *v0     = *ufB_ *     ( *v3 + *v8 ) +
		*(ufB_+1) * ( *v5 + *v6 ) +
		*(ufB_+2) * ( *v2 + *v9 ) +
		*(ufB_+3) * ( *v1 + *v4 + *v7 + *v10 ) +
		*(ufB_+4) * ( *v12 + *v13 + *v14 + *v15 + *v16 + *v17 + *v18 + *v19 ) -
		*v11;
    }

    void advanceCornerF (
                       T* v0,  T* v1,  T* v2,
                       T* v3,  T* v4,  T* v5,
                       T* v6,  T* v7,  T* v8,
                       T* v9,  T* v10, T* v11,
                       T* v12, T* v13, T* v14,
                       T* v15, T* v16, T* v17,
                       T* v18, T* v19, T* v20,
                       T* v21, T* v22, T* v23)
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

    };

    void advanceCornerS (
                         T* v0,  T* v1,  T* v2,
                         T* v3,  T* v4,  T* v5,
                         T* v6,  T* v7,  T* v8,
                         T* v9,  T* v10, T* v11,
                         T* v12, T* v13, T* v14,
                         T* v15, T* v16, T* v17,
                         T* v18, T* v19, T* v20,
                         T* v21, T* v22, T* v23)
    {
      *v0     = - ( *v1  * *(ufB_+16) + *v2  * *(ufB_+8) +
	  *v3  * *(ufB_+1) + *v4  * *(ufB_+16) + *v5  * *(ufB_+9) +
	  *v6  * *(ufB_+2) + *v7  * *(ufB_+16) + *v8  * *(ufB_+10) +
	  *v9  * *(ufB_+3) + *v10 * *(ufB_+16) + *v11 * *(ufB_+11) +
	  *v12 * *(ufB_+4) + *v13 * *(ufB_+16) + *v14 * *(ufB_+12) +
	  *v15 * *(ufB_+5) + *v16 * *(ufB_+16) + *v17 * *(ufB_+13) +
	  *v18 * *(ufB_+6) + *v19 * *(ufB_+16) + *v20 * *(ufB_+14) +
	  *v21 * *(ufB_+7) + *v22 * *(ufB_+16) + *v23 * *(ufB_+15)) / *ufB_;
    };

    T*	ufa_;
    T*  ufB_;
    T	alpha_, beta_;
  };

  /* Structure of data required for updating the field.							*/
  template <typename T> struct UpdateField
  {
    T					jw, jwt, jb;
    Double				w,  b;
    Double				dt,  dx,  dy,  dz;
    Double                              dx2, dy2, dz2;
    T					a 	[6];

    T					bB  	[5];
    T					cB  	[5];
    T					dB  	[5];

    T					eE	[5];
    T					fE	[5];
    T					gE	[5];

    T					hC	[17];

    AdvanceField<Double>		af;

    Double 				*v0, *v1, *v2;
    Double 				*v3, *v31,*v32;
    Double 				*v4, *v41,*v42;
    Double 				*v5, *v51,*v52;
    Double 				*v6, *v61,*v62;
    Double 				*v7, *v8, *v9;
    Double 				*v10,*v11,*v12,*v13;
    Double 				*v14,*v15,*v16,*v17,*v18,*v19;

    Double				*anp1, *an, *anm1;
    Double				*fnp1, *fn, *fnm1;
    Double				*jn,   *rn;
    Double				*en,   *bn;

    unsigned int			N0m1, N1m1, npm1;
  };

  /* Structure of data required for sampling the field.							*/
  template <typename T> struct SampleField
  {
    std::ofstream*			file;
    FieldVector<Double>			position;
    int					i, j, k, m;
    FieldVector<T>			et, at, bt, jt;
    T					f, q;
    Double				dxr, dyr, dzr;
    T					c2, jw, p;
    double				c1;
    unsigned int			N;
  };

  /* Structure of data required for visualizing the field.						*/
  template <typename T> struct VisualizeField
  {
    std::ofstream* 			file;
    std::string				fileName, path, name;
    unsigned int 			i, j, k, m;
    std::vector<std::vector<T> >	v;
    T					jw;
  };

  /* Structure of data required for saving the field profile.						*/
  template <typename T> struct ProfileField
  {
    unsigned int 			i, j, k, l, m;
    std::ofstream* 			file;
    std::string				fileName;
    FieldVector<T>			e, a, b;
    Double				dt;
    T					jw;
    int                                 kmin, kmax;
  };

  /* The parameters needed for parallel operation of update bunch.					*/
  struct UpdateBunchParallel
  {
    int					i, j, k, n;
    long int                            m;
    bool                                b1, b2;
    std::vector<Double>                 qSB, qSF, qRB, qRF;

    FieldVector<Double>			et, bt;
    Double				dxr, dyr, dzr;
    Double				d2;
    double				d1;
    FieldVector<Double>			gbm, gbp, gbpl;
    Double				sz, cy;

    Double                              lz, ly;
    Double				tsignal;
    Double				d, l, zRp, wrp, zRs, wrs, x, y, z, p, t;
    FieldVector<Double>			rv, yv, rl;
    FieldVector<Double>			ex, ez, eT;
    FieldVector<Double>			by, bz, bT;
    Double				tl, t0;
    Double				r0;

    UpdateBunchParallel()
    {
      qSB.clear(); qSF.clear(); qRB.clear(); qRF.clear();
    };
  };

  /* Structure of data required for saving updating the bunch.						*/
  struct UpdateBunch
  {
    Double				dx, dy, dz, dt, dtb;
    Double				ku, b0;
    Double				r1, r2;
    Double				ct, st;
    int                                 nL, i;
    Charge                              Q;
  };

  /* Structure of data required for sampling the bunch.							*/
  struct SampleBunch
  {
    std::ofstream*			file;
    Double				g, q, qT;
    FieldVector<Double>			r,  gb,  r2,  gb2;
    FieldVector<Double>			rT, gbT, r2T, gb2T;
  };

  /* Structure of data required for visualizing the bunch.						*/
  struct VisualizeBunch
  {
    std::ofstream* 			file;
    std::string				fileName, path, name;
    unsigned int			i, N;
  };

  /* structure of the data required for profiling the bunch.						*/
  struct ProfileBunch
  {
    std::ofstream* 			file;
    std::string				fileName;
  };

  /* structure of the data required for updating the current.						*/
  struct UpdateCurrent
  {
    int                                                 ip, jp, kp, im, jm, km;
    long int                                            m;
    Double                                              dx, dy, dz;
    Double                                              dxp, dyp, dzp, dxm, dym, dzm;
    Double                                              x1, x2, y1, y2, z1, z2, q;
    FieldVector<Double>                                 rp, rm, r;
    FieldVector<Double>                                 jcp, jcm;
    Double						rc;
    Double                                 		dv;
    double						c;
    std::vector<FieldVector<Double> >                   jt;
    std::vector<Double>                                 rt;
  };

  /* Structure of data required for sampling the field.							*/
  struct SampleRadiationPower
  {
    std::vector<std::ofstream*>				file;
    std::vector<Double>					w;
    int							k, l;
    Double						dxr, dyr, dzr;
    double						c;
    std::vector<std::vector<std::vector<Double> > >	pLi;
    std::vector<Double>					pL, pG;
    Double						dt, dx, dy, dz;
    Double						pc;
    unsigned int					N, Nl, Nf, Nz;
    std::vector<std::vector<std::vector<Double> > >	fdt;
    unsigned int					m;
    std::vector<std::vector<Complex> >			ep, em;
  };

  /* Structure of data required for sampling the field.							*/
  struct SampleRadiationEnergy
  {
    std::vector<std::ofstream*>                         file;
    std::vector<Double>                                 w;
    int                                                 k, l;
    Double                                              dxr, dyr, dzr;
    double                                              c;
    std::vector<std::vector<std::vector<Double> > >     pLi;
    std::vector<Double>                                 pL, pG;
    Double                                              dt, dx, dy, dz;
    Double                                              pc;
    unsigned int                                        N, Nl, Nf;
    std::vector<std::vector<std::vector<Double> > >     fdt;
  };

}

#endif
