// database.cpp
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "database.h"
#include "stdinclude.h"

namespace Darius
{
  BunchInitialize::BunchInitialize ()
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

  UpdateBunchParallel::UpdateBunchParallel ()
    {
      qSB.clear(); qSF.clear(); qRB.clear(); qRF.clear();
    }
}
