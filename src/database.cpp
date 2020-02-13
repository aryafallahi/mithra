/********************************************************************************************************
 *  database.cpp : Implementation of functions related to the database in mithra
 ********************************************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "database.h"
#include "stdinclude.h"

namespace Darius
{
  /* Initialize the parameters for the bunch initialization to some first values.                     	*/
  BunchInitialize::BunchInitialize ()
  {
    bunchType_			= "";
    distribution_		= "";
    numberOfParticles_ 		= 0;
    cloudCharge_		= 0.0;
    initialGamma_		= 0.0;
    initialBeta_		= 0.0;
    initialDirection_		= 0.0;
    position_.clear();
    numbers_			= 0;
    latticeConstants_		= 0.0;
    sigmaPosition_		= 0.0;
    sigmaGammaBeta_		= 0.0;
    tranTrun_			= 0.0;
    longTrun_			= 0.0;
    fileName_			= "";
    bF_				= 0.0;
    bFP_			= 0.0;
    shotNoise_			= false;
    lambda_			= 0.0;
  }

  UpdateBunchParallel::UpdateBunchParallel ()
  {
    qSB.clear(); qSF.clear(); qRB.clear(); qRF.clear();
  }

  /* Advance the magnetic potential using the Non-standard Finite-Difference algorithm.			*/
  void AdvanceField::advanceMagneticPotentialNSFD (P v0 , P v1 , P v2 ,
						   P v3 , P v31, P v32,
						   P v4 , P v41, P v42,
						   P v5 , P v51, P v52,
						   P v6 , P v61, P v62,
						   P v7 , P v8 , P v9 )
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

  /* Advance the scalar potential using the Non-standard Finite-Difference algorithm.			*/
  void AdvanceField::advanceScalarPotentialNSFD (P v0 , P v1 , P v2 ,
						 P v3 , P v31, P v32,
						 P v4 , P v41, P v42,
						 P v5 , P v51, P v52,
						 P v6 , P v61, P v62,
						 P v7 , P v8 , P v9 )
  {
    *v0 =
	*ufa_     * *v2 - *v1 + alpha_ * (
	*(ufa_+1) * ( *v3 + *v4 + beta_ * ( *v31 + *v32 + *v41 + *v42 ) )   +
	*(ufa_+2) * ( *v5 + *v6 + beta_ * ( *v51 + *v52 + *v61 + *v62 ) ) ) +
	*(ufa_+3) * ( *v7 + *v8 ) +
	*(ufa_+5) * *v9;
  };

  /* Advance the magnetic potential using the Standard Finite-Difference algorithm.			*/
  void AdvanceField::advanceMagneticPotentialFD (P v0 , P v1 , P v2 ,
						 P v3 , P v31, P v32,
						 P v4 , P v41, P v42,
						 P v5 , P v51, P v52,
						 P v6 , P v61, P v62,
						 P v7 , P v8 , P v9 )
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

  /* Advance the scalar potential using the Standard Finite-Difference algorithm.			*/
  void AdvanceField::advanceScalarPotentialFD	(P v0 , P v1 , P v2 ,
						 P v3 , P v31, P v32,
						 P v4 , P v41, P v42,
						 P v5 , P v51, P v52,
						 P v6 , P v61, P v62,
						 P v7 , P v8 , P v9 )
  {
    *v0 =
	*ufa_     * *v2 - *v1 +
	*(ufa_+1) * ( *v3 + *v4 ) +
	*(ufa_+2) * ( *v5 + *v6 ) +
	*(ufa_+3) * ( *v7 + *v8 ) +
	*(ufa_+5) * *v9;
  }

  /* Advance the magnetic potential at the boundaries of the computational domain.			*/
  void AdvanceField::advanceBoundaryF 		(P v0 , P v1 , P v2 ,
						 P v3 , P v4 , P v5 ,
						 P v6 , P v7 , P v8 ,
						 P v9 , P v10, P v11,
						 P v12, P v13)
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

  /* Advance the scalar potential at the boundaries of the computational domain.			*/
  void AdvanceField::advanceBoundaryS 		(P v0 , P v1 , P v2 ,
						 P v3 , P v4 , P v5 ,
						 P v6 , P v7 , P v8 ,
						 P v9 , P v10, P v11,
						 P v12, P v13)
  {
    *v0 =
	*ufB_     * ( *v1 + *v5 ) +
	*(ufB_+1) * *v3 +
	*(ufB_+2) * ( *v2 + *v4 ) +
	*(ufB_+3) * ( *v6 + *v7 + *v10 + *v11 ) +
	*(ufB_+4) * ( *v8 + *v9 + *v12 + *v13 );
  }

  /* Advance the magnetic potential at the edges of the computational domain.				*/
  void AdvanceField::advanceEdgeF 		(P v0 , P v1 , P v2 ,
						 P v3 , P v4 , P v5 ,
						 P v6 , P v7 , P v8 ,
						 P v9 , P v10, P v11,
						 P v12, P v13, P v14,
						 P v15, P v16, P v17,
						 P v18, P v19)
  {
    *v0     =
	*ufB_ *     ( *v3 + *v8 ) +
	*(ufB_+1) * ( *v5 + *v6 ) +
	*(ufB_+2) * ( *v2 + *v9 ) +
	*(ufB_+3) * ( *v1 + *v4 + *v7 + *v10 ) +
	*(ufB_+4) * ( *v12 + *v13 + *v14 + *v15 + *v16 + *v17 + *v18 + *v19 ) -
	*v11;

    *(v0+1) =
	*ufB_ *     ( *(v3+1) + *(v8+1) ) +
	*(ufB_+1) * ( *(v5+1) + *(v6+1) ) +
	*(ufB_+2) * ( *(v2+1) + *(v9+1) ) +
	*(ufB_+3) * ( *(v1+1) + *(v4+1) + *(v7+1) + *(v10+1) ) +
	*(ufB_+4) * ( *(v12+1)+ *(v13+1)+ *(v14+1)+ *(v15+1)+ *(v16+1)+ *(v17+1)+ *(v18+1)+ *(v19+1)) -
	*(v11+1);

    *(v0+2) =
	*ufB_ *     ( *(v3+2) + *(v8+2) ) +
	*(ufB_+1) * ( *(v5+2) + *(v6+2) ) +
	*(ufB_+2) * ( *(v2+2) + *(v9+2) ) +
	*(ufB_+3) * ( *(v1+2) + *(v4+2) + *(v7+2) + *(v10+2) ) +
	*(ufB_+4) * ( *(v12+2)+ *(v13+2)+ *(v14+2)+ *(v15+2)+ *(v16+2)+ *(v17+2)+ *(v18+2)+ *(v19+2)) -
	*(v11+2);

  }

  /* Advance the scalar potential at the edges of the computational domain.				*/
  void AdvanceField::advanceEdgeS 		(P v0 , P v1 , P v2 ,
						 P v3 , P v4 , P v5 ,
						 P v6 , P v7 , P v8 ,
						 P v9 , P v10, P v11,
						 P v12, P v13, P v14,
						 P v15, P v16, P v17,
						 P v18, P v19)
  {
    *v0     =
	*ufB_ *     ( *v3 + *v8 ) +
	*(ufB_+1) * ( *v5 + *v6 ) +
	*(ufB_+2) * ( *v2 + *v9 ) +
	*(ufB_+3) * ( *v1 + *v4 + *v7 + *v10 ) +
	*(ufB_+4) * ( *v12 + *v13 + *v14 + *v15 + *v16 + *v17 + *v18 + *v19 ) -
	*v11;
  }

  /* Advance the magnetic potential at the corners of the computational domain.				*/
  void AdvanceField::advanceCornerF 		(P v0 , P v1 , P v2 ,
						 P v3 , P v4 , P v5 ,
						 P v6 , P v7 , P v8 ,
						 P v9 , P v10, P v11,
						 P v12, P v13, P v14,
						 P v15, P v16, P v17,
						 P v18, P v19, P v20,
						 P v21, P v22, P v23)
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

  /* Advance the scalar potential at the corners of the computational domain.				*/
  void AdvanceField::advanceCornerS 		(P v0 , P v1 , P v2 ,
						 P v3 , P v4 , P v5 ,
						 P v6 , P v7 , P v8 ,
						 P v9 , P v10, P v11,
						 P v12, P v13, P v14,
						 P v15, P v16, P v17,
						 P v18, P v19, P v20,
						 P v21, P v22, P v23)
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
