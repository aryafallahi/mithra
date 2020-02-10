/********************************************************************************************************
 *  database.h : Implementation of the header files related to the database
 ********************************************************************************************************/

#ifndef DATABASE_H_
#define DATABASE_H_

#include <fstream>
#include <list>
#include <string>
#include <vector>

#include "fieldvector.h"
#include "stdinclude.h"

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

    /* List of initial particles. Read from file or from OPAL.						*/
    std::list<Charge>   inputVector_;

    /* Bunching factor for the initialization of the bunch.						*/
    Double				bF_;

    /* Phase of the bunching factor for the initialization of the bunch.				*/
    Double				bFP_;

    /* Boolean flag determining the activation of shot-noise.						*/
    bool				shotNoise_;

    /* Initial beta vector of the bunch, which is obtained as the product of beta and direction.	*/
    FieldVector<Double>			betaVector_;

    /* Initialize the parameters for the bunch initialization to some first values.                     */
    BunchInitialize ();
  };

  /* Structure of data required for updating the field.							*/
  class AdvanceField
  {
  public:

    typedef Double* P;

    void advanceMagneticPotentialNSFD 	(P v0 , P v1 , P v2 ,
					 P v3 , P v31, P v32,
					 P v4 , P v41, P v42,
					 P v5 , P v51, P v52,
					 P v6 , P v61, P v62,
					 P v7 , P v8 , P v9 );

    void advanceScalarPotentialNSFD 	(P v0 , P v1 , P v2 ,
					 P v3 , P v31, P v32,
					 P v4 , P v41, P v42,
					 P v5 , P v51, P v52,
					 P v6 , P v61, P v62,
					 P v7 , P v8 , P v9 );

    void advanceMagneticPotentialFD 	(P v0 , P v1 , P v2 ,
					 P v3 , P v31, P v32,
					 P v4 , P v41, P v42,
					 P v5 , P v51, P v52,
					 P v6 , P v61, P v62,
					 P v7 , P v8 , P v9 );

    void advanceScalarPotentialFD 	(P v0 , P v1 , P v2 ,
					 P v3 , P v31, P v32,
					 P v4 , P v41, P v42,
					 P v5 , P v51, P v52,
					 P v6 , P v61, P v62,
					 P v7 , P v8 , P v9 );

    void advanceBoundaryF 		(P v0 , P v1 , P v2 ,
					 P v3 , P v4 , P v5 ,
					 P v6 , P v7 , P v8 ,
					 P v9 , P v10, P v11,
					 P v12, P v13);

    void advanceBoundaryS 		(P v0 , P v1 , P v2 ,
					 P v3 , P v4 , P v5 ,
					 P v6 , P v7 , P v8 ,
					 P v9 , P v10, P v11,
					 P v12, P v13);

    void advanceEdgeF 			(P v0 , P v1 , P v2 ,
					 P v3 , P v4 , P v5 ,
					 P v6 , P v7 , P v8 ,
					 P v9 , P v10, P v11,
					 P v12, P v13, P v14,
					 P v15, P v16, P v17,
					 P v18, P v19);

    void advanceEdgeS 			(P v0 , P v1 , P v2 ,
					 P v3 , P v4 , P v5 ,
					 P v6 , P v7 , P v8 ,
					 P v9 , P v10, P v11,
					 P v12, P v13, P v14,
					 P v15, P v16, P v17,
					 P v18, P v19);

    void advanceCornerF 		(P v0 , P v1 , P v2 ,
					 P v3 , P v4 , P v5 ,
					 P v6 , P v7 , P v8 ,
					 P v9 , P v10, P v11,
					 P v12, P v13, P v14,
					 P v15, P v16, P v17,
					 P v18, P v19, P v20,
					 P v21, P v22, P v23);

    void advanceCornerS 		(P v0 , P v1 , P v2 ,
					 P v3 , P v4 , P v5 ,
					 P v6 , P v7 , P v8 ,
					 P v9 , P v10, P v11,
					 P v12, P v13, P v14,
					 P v15, P v16, P v17,
					 P v18, P v19, P v20,
					 P v21, P v22, P v23);
    P 					ufa_;
    P 					ufB_;
    Double 				alpha_;
    Double 				beta_;
  };

  /* Structure of data required for updating the field.							*/
  struct UpdateField
  {
    Double				jw, jwt, jb;
    Double				w,  b;
    Double				dt,  dx,  dy,  dz;
    Double				dx2, dy2, dz2;
    Double				a 	[6];

    Double				bB  	[5];
    Double				cB  	[5];
    Double				dB  	[5];

    Double				eE	[5];
    Double				fE	[5];
    Double				gE	[5];

    Double				hC	[17];

    AdvanceField			af;

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
  struct SampleField
  {
    std::ofstream*			file;
    FieldVector<Double>			position;
    int					i, j, k, m;
    FieldVector<Double>			et, at, bt, jt;
    Double				f, q;
    Double				dxr, dyr, dzr;
    Double				c2, jw, p;
    double				c1;
    unsigned int			N;
  };

  /* Structure of data required for visualizing the field.						*/
  struct VisualizeField
  {
    std::ofstream* 			file;
    std::string				fileName, path, name;
    unsigned int 			i, j, k, m;
    std::vector<std::vector<Double> >	v;
    Double				jw;
  };

  /* Structure of data required for saving the field profile.						*/
  struct ProfileField
  {
    unsigned int 			i, j, k, l, m;
    std::ofstream* 			file;
    std::string				fileName;
    FieldVector<Double>			e, a, b;
    Double				dt;
    Double				jw;
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

    UpdateBunchParallel();
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
    Double				longTrun, longTrunT;
    
    SampleBunch();
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
    int					ip, jp, kp, im, jm, km;
    long int				m;
    Double				dx, dy, dz;
    Double				dxp, dyp, dzp, dxm, dym, dzm;
    Double				x1, x2, y1, y2, z1, z2, q;
    FieldVector<Double>			rp, rm, r;
    FieldVector<Double>			jcp, jcm;
    Double				rc;
    Double				dv;
    Double				c;
    std::vector<FieldVector<Double> >	jt;
    std::vector<Double>			rt;
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
    std::vector<std::ofstream*>				file;
    std::vector<Double>					w;
    int							k, l;
    Double						dxr, dyr, dzr;
    Double						c;
    std::vector<std::vector<std::vector<Double> > >     pLi;
    std::vector<Double>					pL, pG;
    Double						dt, dx, dy, dz;
    Double						pc;
    unsigned int					N, Nl, Nf;
    std::vector<std::vector<std::vector<Double> > >     fdt;
  };

  /* Structure of data required for bunch profile from screens.							*/
  struct SampleScreenProfile
  {
    std::vector<std::ofstream*> 	files;
    std::vector<std::string>		fileNames;
  };
}

#endif
