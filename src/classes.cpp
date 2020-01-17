/********************************************************************************************************
 *  classes.hh : Implementation of the functions related to classes used in mithra
 ********************************************************************************************************/

#include <algorithm>
#include <fstream>

#include "classes.h"

namespace Darius
{

  /*** Mesh class ***************************************************************************************/

  /* Show the stored values for the mesh.								*/
  void Mesh::show ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string(" Length scale = ") + stringify(lengthScale_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Total mesh length vector = ") + stringify(meshLength_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Meshing resolution vector = ") + stringify(meshResolution_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Time scale = ") + stringify(timeScale_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Total simulation time = ") + stringify(totalTime_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Mesh truncation order = ") + stringify(truncationOrder_));
    if 	( spaceCharge_ == true )
      printmessage(std::string(__FILE__), __LINE__, std::string(" Space-charge = true "));
    else if 	( spaceCharge_ == false )
      printmessage(std::string(__FILE__), __LINE__, std::string(" Space-charge = false "));
    if 	( solver_ == NSFD )
      printmessage(std::string(__FILE__), __LINE__, std::string(" Solver = Non-standard finite-difference "));
    else if 	( solver_ == FD )
      printmessage(std::string(__FILE__), __LINE__, std::string(" Solver = finite-difference "));

  }

  /* Initialize the parameters in the mesh initializer.							*/
  void Mesh::initialize ()
  {
    spaceCharge_ 		= false;
    solver_			= NSFD;
  }

  /*** Bunch class **************************************************************************************/

  /* The constructor clears and initializes the internal data structure.				*/
  Bunch::Bunch ()
  {
    /* Initialize the parameters for the bunch to some first values.                      		*/
    timeStep_			= 0.0;

    /** Initialize the parameters for the bunch to some first values.                           	*/
    sampling_			= false;
    directory_			= "./";
    basename_             	= "";
    rhythm_			= 0.0;

    /** Initialize the parameters for the vtk visualization to some first values.               	*/
    bunchVTK_             	= false;
    bunchVTKDirectory_    	= "./";
    bunchVTKBasename_    	= "";
    bunchVTKRhythm_      	= 0.0;

    /** Initialize the parameters for saving the bunch profile to some first values.           		*/
    bunchProfile_            	= false;
    bunchProfileDirectory_    	= "./";
    bunchProfileBasename_     	= "";
    bunchProfileTime_.clear();
    bunchProfileRhythm_		= 0.0;
  }

  /* Initialize a bunch with a manual type. This bunch produces one charge equal to the cloudCharge_. 	*/
  void Bunch::initializeManual (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia)
  {
    /* Declare the required parameters for the initialization of charge vectors.                      	*/
    Charge        charge;

    /* Determine the properties of each charge point and add them to the charge vector.               	*/
    charge.q  	= bunchInit.cloudCharge_;
    charge.rnp    	= bunchInit.position_[ia];
    charge.gbnp.mv( bunchInit.initialGamma_, bunchInit.betaVector_ );

    /* Insert this charge to the charge list if and only if it resides in the processor's portion.    	*/
    if ( ( charge.rnp[2] < zp[1] || rank == size - 1 ) && ( charge.rnp[2] >= zp[0] || rank == 0 ) )
      chargeVector.push_back(charge);
  }

  /* Initialize a bunch with an ellipsoid type. This bunch produces a number of charges equal to the
   * numberOfParticles_ with the total charge equal to the cloudCharge_ which are distributed in an
   * ellipsoid with dimensions given by sigmaPosition_ and center given by the position vector. The
   * particles have uniform energy distribution centered at initialEnergy_ with variances determined by
   * sigmaGammaBeta_.                                                         				*/
  void Bunch::initializeEllipsoid (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia)
  {
    /* Save the initially given number of particles.							*/
    unsigned int	Np = bunchInit.numberOfParticles_, i, Np0 = chargeVector.size();

    /* Declare the required parameters for the initialization of charge vectors.                      	*/
    Charge            charge; charge.q  = bunchInit.cloudCharge_ / Np;
    FieldVector<Double> gb (0.0); gb.mv( bunchInit.initialGamma_, bunchInit.betaVector_ );
    FieldVector<Double> r  (0.0);
    FieldVector<Double> t  (0.0);
    Double            t0, g;
    Double		zmin = 1e100;
    Double		Ne, bF, bFi;
    unsigned int	bmi;

    /* Check the bunching factor.                                                                     	*/
    if ( bunchInit.bF_ > 2.0 || bunchInit.bF_ < 0.0 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The bunching factor can not be larger than one or a negative value !!!") );
	exit(1);
      }

    /* Declare the function for injecting the shot noise.						*/
    auto insertCharge = [&] (Charge q) {

      for ( unsigned int ii = 0; ii < 4; ii++ )
	{
	  /* The random modulation is introduced depending on the shot-noise being activated.		*/
	  if ( bunchInit.shotNoise_ )
	    {
	      /* Obtain the number of beamlet.								*/
	      bmi = int( ( charge.rnp[2] - zmin ) / ( bunchInit.lambda_ / 2.0 ) );

	      /* Obtain the phase and amplitude of the modulation.					*/
	      bFi = bF * sqrt( - 2.0 * log( halton( 8 , bmi ) ) );

	      q.rnp[2]  = charge.rnp[2] - ( bunchInit.lambda_ / 2.0 ) / 4 * ii;

	      q.rnp[2] -= bunchInit.lambda_ / ( 2.0 * PI ) * bFi * sin( 2.0 * PI / ( bunchInit.lambda_ / 2.0 ) * q.rnp[2] + 2.0 * PI * halton( 9 , bmi ) );
	    }
	  else if ( bunchInit.lambda_ != 0.0)
	    {
	      q.rnp[2]  = charge.rnp[2] - ( bunchInit.lambda_ / 2.0 ) / 4 * ii;

	      q.rnp[2] -= bunchInit.lambda_ / ( 2.0 * PI ) * bunchInit.bF_ * sin( 2.0 * PI / ( bunchInit.lambda_ / 2.0 ) * q.rnp[2] + bunchInit.bFP_ * PI / 180.0 );
	    }

	  /* Before add the bunch to the global charge vector, a correction on the position of the bunch
	   * should be made. This correction assures that the bunch properties are valid at the entrance
	   * of the undulator.										*/
	  g		= sqrt( 1.0 + charge.gbnp.norm() );
	  q.rnp[0] -= ( charge.gbnp[0] / g - bunchInit.betaVector_[0] ) * ( zu_ - charge.rnp[2] ) / ( bunchInit.betaVector_[2] + beta_ );
	  q.rnp[1] -= ( charge.gbnp[1] / g - bunchInit.betaVector_[1] ) * ( zu_ - charge.rnp[2] ) / ( bunchInit.betaVector_[2] + beta_ );
	  q.rnp[2] -= ( charge.gbnp[2] / g - bunchInit.betaVector_[2] ) * ( zu_ - charge.rnp[2] ) / ( bunchInit.betaVector_[2] + beta_ );

	  /* Insert this charge to the charge list if and only if it resides in the processor's
	   * portion.                                                                               	*/
	  if ( ( q.rnp[2] < zp[1] || rank == size - 1 ) && ( q.rnp[2] >= zp[0] || rank == 0 ) )
	    chargeVector.push_back(q);
	}
    };

    /* If the shot noise is on, we need the minimum value of the bunch z coordinate to be able to
     * calculate the FEL bucket number.									*/
    if ( bunchInit.shotNoise_ )
      {
	for (i = 0; i < Np / 4; i++)
	  {
	    if ( bunchInit.distribution_ == "uniform" )
	      zmin = std::min(   ( 2.0 * halton(2, i + Np0) - 1.0 ) * bunchInit.sigmaPosition_[2] , zmin );
	    else if ( bunchInit.distribution_ == "gaussian" )
	      zmin = std::min(   bunchInit.sigmaPosition_[2] * sqrt( - 2.0 * log( halton(2, i + Np0) ) ) * sin( 2.0 * PI * halton(3, i + Np0) ) , zmin );
	    else
	      {
		printmessage(std::string(__FILE__), __LINE__, std::string("The longitudinal type is not correctly given to the code !!!") );
		exit(1);
	      }
	  }

	if ( bunchInit.distribution_ == "uniform" )
	  for ( ; i < unsigned( Np / 4 * ( 1.0 + bunchInit.lambda_ * sqrt( 2.0 * PI ) / ( 2.0 * bunchInit.sigmaPosition_[2] ) ) ); i++)
	    {
	      t0  = bunchInit.lambda_ * sqrt( - 2.0 * log( halton( 2, i + Np0 ) ) ) * sin( 2.0 * PI * halton( 3, i + Np0 ) );
	      t0 += ( t0 < 0.0 ) ? ( - bunchInit.sigmaPosition_[2] ) : ( bunchInit.sigmaPosition_[2] );

	      zmin = std::min(   t0 , zmin );
	    }

	zmin = zmin + bunchInit.position_[ia][2];

	/* Obtain the average number of electrons per FEL beamlet.					*/
	Ne = bunchInit.cloudCharge_ * ( bunchInit.lambda_ / 2.0 ) / ( 2.0 * bunchInit.sigmaPosition_[2] );

	/* Set the bunching factor level for the shot noise depending on the given values.		*/
	bF = ( bunchInit.bF_ == 0.0 ) ? 1.0 / sqrt(Ne) : bunchInit.bF_;

	printmessage(std::string(__FILE__), __LINE__, std::string("The standard deviation of the bunching factor for the shot noise implementation is set to ") + stringify(bF) );
      }

    /* Determine the properties of each charge point and add them to the charge vector.               	*/
    for (i = 0; i < Np / 4; i++)
      {
	/* Determine the transverse coordinate.							*/
	r[0] = bunchInit.sigmaPosition_[0] * sqrt( - 2.0 * log( halton(0, i + Np0) ) ) * cos( 2.0 * PI * halton(1, i + Np0) );
	r[1] = bunchInit.sigmaPosition_[1] * sqrt( - 2.0 * log( halton(0, i + Np0) ) ) * sin( 2.0 * PI * halton(1, i + Np0) );

	/* Determine the longitudinal coordinate.							*/
	if ( bunchInit.distribution_ == "uniform" )
	  r[2] = ( 2.0 * halton(2, i + Np0) - 1.0 ) * bunchInit.sigmaPosition_[2];
	else if ( bunchInit.distribution_ == "gaussian" )
	  r[2] = bunchInit.sigmaPosition_[2] * sqrt( - 2.0 * log( halton(2, i + Np0) ) ) * sin( 2.0 * PI * halton(3, i + Np0) );
	else
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The longitudinal type is not correctly given to the code !!!") );
	    exit(1);
	  }

	/* Determine the transverse momentum.								*/
	t[0] = bunchInit.sigmaGammaBeta_[0] * sqrt( - 2.0 * log( halton(4, i + Np0) ) ) * cos( 2.0 * PI * halton(5, i + Np0) );
	t[1] = bunchInit.sigmaGammaBeta_[1] * sqrt( - 2.0 * log( halton(4, i + Np0) ) ) * sin( 2.0 * PI * halton(5, i + Np0) );
	t[2] = bunchInit.sigmaGammaBeta_[2] * sqrt( - 2.0 * log( halton(6, i + Np0) ) ) * cos( 2.0 * PI * halton(7, i + Np0) );

	if ( sqrt( pow(r[0],2) + pow(r[1],2) ) < bunchInit.tranTrun_ && fabs(r[2]) < bunchInit.longTrun_)
	  {
	    /* Shift the generated charge to the center position and momentum space.			*/
	    charge.rnp    = bunchInit.position_[ia];
	    charge.rnp   += r;

	    charge.gbnp   = gb;
	    charge.gbnp  += t;

	    /* Insert this charge and the mirrored ones into the charge vector.				*/
	    insertCharge(charge);
	  }
      }

    /* If the longitudinal type of the bunch is uniform a tapered part needs to be added to remove the
     * CSE from the tail of the bunch.									*/
    if ( bunchInit.distribution_ == "uniform" )
      for ( ; i < unsigned( Np / 4 * ( 1.0 + bunchInit.lambda_ * sqrt( 2.0 * PI ) / ( 2.0 * bunchInit.sigmaPosition_[2] ) ) ); i++)
	{
	  r[0] = bunchInit.sigmaPosition_[0] * sqrt( - 2.0 * log( halton(0, i + Np0) ) ) * cos( 2.0 * PI * halton(1, i + Np0) );
	  r[1] = bunchInit.sigmaPosition_[1] * sqrt( - 2.0 * log( halton(0, i + Np0) ) ) * sin( 2.0 * PI * halton(1, i + Np0) );

	  /* Determine the longitudinal coordinate.							*/
	  r[2] = bunchInit.lambda_ * sqrt( - 2.0 * log( halton(2, i + Np0) ) ) * sin( 2.0 * PI * halton(3, i + Np0) );
	  r[2] += ( r[2] < 0.0 ) ? ( - bunchInit.sigmaPosition_[2] ) : ( bunchInit.sigmaPosition_[2] );

	  /* Determine the transverse momentum.								*/
	  t[0] = bunchInit.sigmaGammaBeta_[0] * sqrt( - 2.0 * log( halton(4, i + Np0) ) ) * cos( 2.0 * PI * halton(5, i + Np0) );
	  t[1] = bunchInit.sigmaGammaBeta_[1] * sqrt( - 2.0 * log( halton(4, i + Np0) ) ) * sin( 2.0 * PI * halton(5, i + Np0) );
	  t[2] = bunchInit.sigmaGammaBeta_[2] * sqrt( - 2.0 * log( halton(6, i + Np0) ) ) * cos( 2.0 * PI * halton(7, i + Np0) );

	  if ( sqrt( pow(r[0],2) + pow(r[1],2) ) < bunchInit.tranTrun_ && fabs(r[2]) < bunchInit.longTrun_)
	    {
	      /* Shift the generated charge to the center position and momentum space.			*/
	      charge.rnp   = bunchInit.position_[ia];
	      charge.rnp  += r;

	      charge.gbnp  = gb;
	      charge.gbnp += t;

	      /* Insert this charge and the mirrored ones into the charge vector.			*/
	      insertCharge(charge);
	    }
	}

    /* Reset the value for the number of particle variable according to the installed number of
     * macro-particles and perform the corresponding changes.                                         	*/
    bunchInit.numberOfParticles_ = chargeVector.size();
  }

  /* Initialize a bunch with a 3D-crystal type. This bunch produces a number of charges equal to the
   * numberOfParticles_ with the total charge equal to the cloudCharge_ which are arranged in a 3D crystal.
   * The number of particles in each direction is given by numbers_. Therefore, numberOfParticles_ should
   * be a multiple of the product of all these three numbers. The ratio gives the number of particles in
   * each crystal point. Position of each particle is determined by the lattice constants and the crystal
   * is centered at the position_ vector. At each point, the charges have a small Gaussian distribution
   * around the crsytal.                                           					*/
  void Bunch::initialize3DCrystal (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia)
  {
    /* Check if the numberOfParticles_ is a multiple of the product of values in numbers_.            	*/
    if ( bunchInit.numberOfParticles_ % (bunchInit.numbers_[0] * bunchInit.numbers_[1] * bunchInit.numbers_[2]) != 0 )
      {
	printmessage(std::string(__FILE__), __LINE__,
		     std::string("The number of the particles and their lattice numbers do not match !!!") );
	exit(1);
      }

    /* Declare the required parameters for the initialization of charge vectors.                      	*/
    Charge            	charge;
    FieldVector<Double> 	gb (0.0);
    gb.mv( bunchInit.initialGamma_, bunchInit.betaVector_ );
    unsigned int np = bunchInit.numberOfParticles_ / (bunchInit.numbers_[0] * bunchInit.numbers_[1] * bunchInit.numbers_[2]);

    /* Clear the charge vector for adding the charges.							*/
    chargeVector.clear();

    /* Determine the properties of each charge point and add them to the charge vector.               	*/
    for (unsigned int i = 0; i < bunchInit.numbers_[0]; i++)
      {
	for (unsigned int j = 0; j < bunchInit.numbers_[1]; j++)
	  {
	    for (unsigned int k = 0; k < bunchInit.numbers_[2]; k++)
	      {
		for (unsigned int l = 0; l < np; l++)
		  {
		    charge.q  = bunchInit.cloudCharge_ / bunchInit.numberOfParticles_;

		    charge.rnp[0]  = bunchInit.position_[ia][0] + ( i + 1.0 - 0.5 * bunchInit.numbers_[0] ) * bunchInit.latticeConstants_[0];
		    charge.rnp[1]  = bunchInit.position_[ia][1] + ( j + 1.0 - 0.5 * bunchInit.numbers_[1] ) * bunchInit.latticeConstants_[1];
		    charge.rnp[2]  = bunchInit.position_[ia][2] + ( k + 1.0 - 0.5 * bunchInit.numbers_[2] ) * bunchInit.latticeConstants_[2];
		    charge.rnp[0] += 0.5 * bunchInit.sigmaPosition_[0] * sqrt( - 2.0 * log( halton(0,i) ) ) * sin( 2.0 * PI * halton(1,i) );
		    charge.rnp[1] += 0.5 * bunchInit.sigmaPosition_[1] * sqrt( - 2.0 * log( halton(2,i) ) ) * sin( 2.0 * PI * halton(3,i) );
		    charge.rnp[2] += 0.5 * bunchInit.sigmaPosition_[2] * sqrt( - 2.0 * log( halton(4,i) ) ) * sin( 2.0 * PI * halton(5,i) );

		    charge.gbnp     = gb;
		    charge.gbnp[0] += bunchInit.sigmaGammaBeta_[0] * sqrt( - 2.0 * log( halton(6,i) ) ) * sin( 2.0 * PI * halton(7,i) );
		    charge.gbnp[1] += bunchInit.sigmaGammaBeta_[1] * sqrt( - 2.0 * log( halton(8,i) ) ) * sin( 2.0 * PI * halton(9,i) );
		    charge.gbnp[2] += bunchInit.sigmaGammaBeta_[2] * sqrt( - 2.0 * log( halton(10,i)) ) * sin( 2.0 * PI * halton(11,i));

		    /* Insert this charge to the charge list if and only if it resides in the processor
		     * portion.    									*/
		    if ( ( charge.rnp[2] < zp[1] || rank == size - 1 ) && ( charge.rnp[2] >= zp[0] || rank == 0 ) )
		      chargeVector.push_back(charge);
		  }
	      }
	  }
      }
  }

  /* Initialize a bunch with a file type. This bunch produces a number of charges read from a given file.
   * The number of initialized charge is equal to the vertical length of the table in the text file. The
   * file format should contain the charge value, 3 position coordinates and 3 momentum coordinates of
   * of the charge distribution.							                */
  void Bunch::initializeFile (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia)
  {
    /* Impose bF unless bF = 0.									*/
    if ( bunchInit.bF_ != 0 ){
      bunchInit.numberOfParticles_ *= 4;
      for (std::list<Charge>::iterator it = bunchInit.inputVector_.begin(); it != bunchInit.inputVector_.end(); ++it)
	{
	  it->q /= 4;
	  for ( unsigned int ii = 1; ii < 4; ii++ )
	    {
	      Charge charge = *it;
	      /* Do particle mirroring to suppress initial bunching factor                                */
	      charge.rnp[2]  -=  ( bunchInit.lambda_ / 2.0 ) / 4 * ii;
	      /* Impose the given bF                                                                      */
	      charge.rnp[2] -= bunchInit.lambda_ / ( 2.0 * PI ) * bunchInit.bF_ * sin( 2.0 * PI / ( bunchInit.lambda_ / 2.0 ) * charge.rnp[2]);
	      bunchInit.inputVector_.insert(it,charge);
	    }
	}
      printmessage(std::string(__FILE__), __LINE__, std::string("Added bunching factor, and increased number of particles by a factor of 4.") );
    }

    /* Fill in the charge vector.										 */
    distributeParticles(bunchInit.inputVector_, zp, rank, size);
    chargeVector.clear();
    chargeVector = bunchInit.inputVector_;
    bunchInit.inputVector_.clear();

    /* Calculate the total amount of installed particles.						*/
    unsigned int NqL = chargeVector.size(), NqG = 0;
    MPI_Reduce(&NqL,&NqG,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    
    /* Check the size of the charge vector with the number of particles.                            	*/
    if ( bunchInit.numberOfParticles_ != NqG && rank == 0 )
      {
        printmessage(std::string(__FILE__), __LINE__, std::string("The number of the particles and the input vector size do not match !!!") );
        exit(1);
      }
    printmessage(std::string(__FILE__), __LINE__, std::string("The number of electrons per macro particle is " + stringify(chargeVector.begin()->q) ) );
  }

  /* Redistribute particles among processors.								*/
  void Bunch::distributeParticles (std::list<Charge>& chargeVector, Double (zp) [2], int rank, int size)
  {
    /* Note that this function only redistributes q, rnp, gbnp, but NOT rnm and gbnm.			*/

    std::vector<Double> sendCV;
    std::list<Charge>::iterator it = chargeVector.begin();
    while(it != chargeVector.end())
      {
	if ( ( it->rnp[2] < zp[1] || rank == size - 1 ) && ( it->rnp[2] >= zp[0] || rank == 0 ) )
	  it++;
	else
	  {
	    sendCV.push_back(it->q);
	    sendCV.push_back(it->rnp[0]);
	    sendCV.push_back(it->rnp[1]);
	    sendCV.push_back(it->rnp[2]);
	    sendCV.push_back(it->gbnp[0]);
	    sendCV.push_back(it->gbnp[1]);
	    sendCV.push_back(it->gbnp[2]);
	    it = chargeVector.erase(it);
	  }
      }

    /* Send the charges that were erased from the charge vector.					*/
    int sizeSend = sendCV.size();
    int *counts = new int [size], *disps = new int [size];
    MPI_Allgather( &sizeSend, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD );
    for (int i = 0; i < size; i++)
      disps[i] = (i > 0) ? (disps[i-1] + counts[i-1]) : 0;
    std::vector<Double> recvCV (disps[size-1] + counts[size-1], 0);
    MPI_Allgatherv(&sendCV[0], sizeSend, MPI_DOUBLE,
		   &recvCV[0], counts, disps, MPI_DOUBLE, MPI_COMM_WORLD);

    /* And now replace all the charges in the charge vector of the corresponding processor.  */
    unsigned i = 0;
    Charge charge;
    while ( i < recvCV.size() )
      {
	if (( recvCV[i+3] < zp[1] || rank == size - 1 ) && ( recvCV[i+3] >= zp[0] || rank == 0 ))
	  {
	    charge.q = recvCV[i++];
	    charge.rnp[0] = recvCV[i++];
	    charge.rnp[1] = recvCV[i++];
	    charge.rnp[2] = recvCV[i++];
	    charge.gbnp[0] = recvCV[i++];
	    charge.gbnp[1] = recvCV[i++];
	    charge.gbnp[2] = recvCV[i++];
	    chargeVector.push_back(charge);	      
	  }
	else
	  i += 7;
      }
  }

  /* Compute bunch statistics.							                	*/
  SampleBunch Bunch::computeBunchSample (std::list<Charge> chargeVector, int size)
  {
    SampleBunch sb;
    sb.longTrun = -1e3;

    for (std::list<Charge>::iterator iter = chargeVector.begin(); iter != chargeVector.end(); iter++)
      {
	sb.q	 += iter->q;
	sb.r .pmv(   iter->q, iter->rnp );
	sb.gb.pmv(   iter->q, iter->gbnp);	

	for (int l = 0; l < 3; l++)
	  {
	    sb.r2 [l] += iter->rnp[l]  * iter->rnp[l]  * iter->q;
	    sb.gb2[l] += iter->gbnp[l] * iter->gbnp[l] * iter->q;
	  }
	if (iter->rnp[2] > sb.longTrun)
	  sb.longTrun = iter->rnp[2];
      }

    /* Add the contribution from each processor.							*/
    MPI_Allreduce(&sb.q    , &sb.qT    , 1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&sb.r[0] , &sb.rT[0] , 3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&sb.r2[0], &sb.r2T[0], 3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&sb.gb[0], &sb.gbT[0], 3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&sb.gb2[0],&sb.gb2T[0],3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    double *lTs = new double [size];
    MPI_Allgather(&sb.longTrun, 1, MPI_DOUBLE, lTs, 1, MPI_DOUBLE, MPI_COMM_WORLD );

    /* Divide the obtained values by the number of particles to calculate the mean value.		*/
    sb.rT   /= sb.qT;
    sb.r2T  /= sb.qT;
    sb.gbT  /= sb.qT;
    sb.gb2T /= sb.qT;
    sb.longTrunT = *std::max_element(lTs, lTs + size);

    return sb;
  }

  /* Show the stored values for the bunch.                                                          	*/
  void Bunch::show ()
  {
    for (unsigned int i = 0; i < bunchInit_.size(); i++)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" Type of the bunch initialization = ") + bunchInit_[i].bunchType_ );
	if (bunchInit_[i].bunchType_ == "ellipsoid")
	  printmessage(std::string(__FILE__), __LINE__, std::string(" Type of the current profile = ") + bunchInit_[i].distribution_ );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Number of macro-particles = ") + stringify(bunchInit_[i].numberOfParticles_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Total charge of the cloud = ") + stringify(bunchInit_[i].cloudCharge_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Initial mean gamma of the bunch = ") + stringify(bunchInit_[i].initialGamma_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Initial mean speed of the bunch = ") + stringify(bunchInit_[i].initialBeta_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Initial direction of the bunch speed = ") + stringify(bunchInit_[i].initialDirection_) );
	for ( unsigned int ia = 0; ia < bunchInit_[i].position_.size(); ia++)
	  printmessage(std::string(__FILE__), __LINE__, std::string(" Initial position of the bunch = ") + stringify(bunchInit_[i].position_[ia]));
	printmessage(std::string(__FILE__), __LINE__, std::string(" Position spread of the bunch = ") + stringify(bunchInit_[i].sigmaPosition_));
	printmessage(std::string(__FILE__), __LINE__, std::string(" Momentum spread of the bunch = ") + stringify(bunchInit_[i].sigmaGammaBeta_));
	printmessage(std::string(__FILE__), __LINE__, std::string(" Longitudinal truncation of the bunch = ") + stringify(bunchInit_[i].longTrun_));
	if ( bunchInit_[i].bunchType_ == "3D-crystal" )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Number of crystal points = ") + stringify(bunchInit_[i].numbers_) );
	    printmessage(std::string(__FILE__), __LINE__, std::string(" Lattice constant of the crystal = ") + stringify(bunchInit_[i].latticeConstants_) );
	  }
	if ( bunchInit_[i].bunchType_ == "file" )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string(" File name of the bunch = ") + bunchInit_[i].fileName_ );
	  }
	printmessage(std::string(__FILE__), __LINE__, std::string(" Initial bunching factor in the bunch = ") + stringify(bunchInit_[i].bF_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Time step for the bunch calculation = ") + stringify(timeStep_) );
      }
    if ( sampling_ )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" Sampling of the bunch data is enabled. ") );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Directory for the bunch sampling = ") + stringify(directory_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Base name for the bunch sampling = ") + stringify(basename_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Rhythm the bunch sampling = ") + stringify(rhythm_) );
      }
    if ( bunchVTK_ )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" vtk visualization of the bunch is enabled. ") );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Directory for the bunch visualization = ") + stringify(bunchVTKDirectory_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Base name for the bunch visualization = ") + stringify(bunchVTKBasename_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Rhythm the bunch visualization = ") + stringify(bunchVTKRhythm_) );
      }
    if ( bunchProfile_ )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" Save the bunch profile is enabled. ") );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Directory for the bunch profile = ") + stringify(bunchProfileDirectory_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Base name for the bunch profile = ") + stringify(bunchProfileBasename_) );
	for (unsigned int i = 0; i < bunchProfileTime_.size(); i++)
	  printmessage(std::string(__FILE__), __LINE__, std::string(" Bunch profile will be saved at = ") + stringify(bunchProfileTime_[i]) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Rhythm the bunch profiling = ") + stringify(bunchProfileRhythm_) );
      }
  }

  /*** Signal class *************************************************************************************/

  /* Initialize the values of the parameters.								*/
  Signal::Signal ()
  {
    t0_ 		= 0.0;
    s_  		= 0.0;
    f0_ 		= 1.0;
    cep_		= 0.0;
    signalType_	= GAUSSIAN;
  }

  /* Initializer with signal type, time offset, variance, frequency and carrier-envelope-phase.       	*/
  void Signal::initialize (std::string type, Double l0, Double s, Double l, Double cep)
  {
    /* Initialize the signal type.                                                                    	*/
    if      ( type.compare("neumann") == 0 )             signalType_ = NEUMANN;
    else if ( type.compare("gaussian") == 0 )            signalType_ = GAUSSIAN;
    else if ( type.compare("secant-hyperbolic") == 0 )   signalType_ = SECANT;
    else if ( type.compare("flat-top") == 0 )   	   signalType_ = FLATTOP;
    else { std::cout << type << " is an unknown signal type for the given set of parameters." << std::endl; exit(1); }

    /* Initialize the time delay, variance and the frequency of the carrier.                          	*/
    t0_ = l0;
    s_  = s;
    f0_ = 1 / l;

    /* Initialize the carrier envelope phase.                                                         	*/
    cep_ = cep * PI / 180;

    /* Check if variance is unequal zero.                                                             	*/
    if (s_ == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" Variance of signal is set to zero. "));
	printmessage(std::string(__FILE__), __LINE__, std::string(" This is not allowed because we divide through the variance. "));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
  }

  Double Signal::self (Double & t, Double & phase)
  {
    /* If the signal is out of the 20*s range around the center, consider it to be zero.		*/
    if ( fabs(t - t0_) > 10.0 * s_ )	return ( 0.0 );
    else
      {
	if  ( signalType_ == NEUMANN )
	  return ( - cos( 2 * PI * f0_ * (t-t0_) + cep_ + phase) * 2.7724 * (t-t0_) / (s_*s_) * exp( -1.3862 * (t-t0_) * (t-t0_) / (s_*s_) ) );

	else if  ( signalType_ == GAUSSIAN )
	  return ( cos( 2*PI*f0_ * (t-t0_) + cep_ + phase) * exp( -1.3862 * (t-t0_) * (t-t0_) / (s_*s_) ) );

	else if  ( signalType_ == SECANT )
	  return ( cos( 2*PI*f0_ * (t-t0_) + cep_ + phase) / cosh( (t-t0_) / s_ ) );

	else if  ( signalType_ == FLATTOP )
	  {
	    if      ( t - t0_ <= - s_ / 2.0 )
	      return ( cos( 2*PI*f0_ * (t-t0_) + cep_ + phase ) * exp( -1.3862 * pow((t-t0_+s_/2.0)*f0_/2.0,2) ) );
	    else if ( t - t0_ <= s_ / 2.0   )
	      return ( cos( 2*PI*f0_ * (t-t0_) + cep_ + phase ) );
	    else
	      return ( cos( 2*PI*f0_ * (t-t0_) + cep_ + phase ) * exp( -1.3862 * pow((t-t0_-s_/2.0)*f0_/2.0,2) ) );
	  }
      }
    return (0.0);
  }

  /* Show the stored values for this signal.                                                          	*/
  void Signal::show ()
  {
    if      (signalType_ == NEUMANN)
      printmessage(std::string(__FILE__), __LINE__, std::string(" Signal type = neumann pulse"));
    else if (signalType_ == GAUSSIAN)
      printmessage(std::string(__FILE__), __LINE__, std::string(" Signal type = gaussian pulse"));
    else if (signalType_ == SECANT)
      printmessage(std::string(__FILE__), __LINE__, std::string(" Signal type = secant hyperbolic pulse"));
    else if (signalType_ == FLATTOP)
      printmessage(std::string(__FILE__), __LINE__, std::string(" Signal type = flat-top pulse"));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Signal offset         = ") + stringify(t0_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Signal variance       = ") + stringify(s_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Signal frequency  = ") + stringify(f0_));
    printmessage(std::string(__FILE__), __LINE__, std::string(" Signal carrier envelope phase = ") + stringify(cep_));
  }

  /*** Seed class ***************************************************************************************/

  Seed::Seed ()
  {
    seedType_ 			= PLANEWAVE;
    c0_				= 0.0;
    amplitude_			= 0.0;
    radius_.resize(2,0.0);

    sampling_			= false;
    samplingType_		= ATPOINT;
    samplingField_.clear();
    samplingDirectory_		= "";
    samplingBasename_		= "";
    samplingRhythm_		= 0.0;
    samplingPosition_.clear();
    samplingRes_		= 0.0;

    vtk_.clear();

    profile_			= false;
    profileField_.clear();
    profileDirectory_		= "";
    profileBasename_		= "";
    profileTime_.clear();
    profileRhythm_		= 0.0;
  }

  void Seed::initialize (std::string        	type,
			 std::vector<Double>    position,
			 std::vector<Double>    direction,
			 std::vector<Double>    polarization,
			 Double                 amplitude,
			 std::vector<Double>	radius,
			 Signal                 signal)
  {
    /* Set the seed type according to the returned string for seedType.                   		*/
    if      ( type.compare("plane-wave"         ) == 0 ) seedType_ = PLANEWAVE;
    else if ( type.compare("confined-plane-wave") == 0 ) seedType_ = PLANEWAVECONFINED;
    else if ( type.compare("gaussian-beam"      ) == 0 ) seedType_ = GAUSSIANBEAM;
    else    { std::cout << type << " is an unknown type." << std::endl; exit(1); }

    /* Set the vectors position, direction and polarization for the seed class.                 	*/
    position_ 	= position;
    polarization_ 	= polarization;
    direction_ 	= direction;

    /* check if length of diection vector is zero and normalize the vector.				*/
    if ( direction_.norm() == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The direction vector has length zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
    else
      {
	Double vl = sqrt( direction_.norm() );
	direction_ /= vl;
      }

    /* check if length of polarization vector is zero and normalize the vector.				*/
    if ( polarization_.norm() == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The polarization vector has length zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
    else
      {
	Double vl = sqrt( polarization_.norm() );
	polarization_ /= vl;
      }

    /* check if polarization and direction are normal to each other.                       		*/
    if ( fabs( polarization_ * direction_ ) > 1.0e-50 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("Polarization is not normal to the direction."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }

    /* Initialize the amplitude of the seed.                                                    	*/
    amplitude_ 		= amplitude;

    /* Initialize the Rayleigh radius of the Gaussian beam.                                           	*/
    radius_ 		= radius;

    /* check if length of polarization vector is zero and normalize the vector.                       	*/
    if ( seedType_ == GAUSSIANBEAM && radius_[0] * radius_[1] == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("One of the radii of the Gaussian beam is set to zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }

    /* Initialize the signal of the seed.                                                       	*/
    signal_ 		= signal;
  }

  /* Return the potentials at any desired location and time.                    			*/
  void Seed::fields (FieldVector <Double> & aufpunkt, Double & time, FieldVector <Double> & a)
  {
    /* Transfer the coordinate from the bunch rest frame to the lab frame.				*/
    rl[0] = aufpunkt[0]; rl[1] = aufpunkt[1];
    rl[2] = gamma_ * ( aufpunkt[2] + beta_ * c0_ * ( time + dt_ ) );
    tl    = gamma_ * ( time + dt_  + beta_ / c0_ * aufpunkt[2]    );

    /* Calculate the distance to the reference position along the propagation direction.   		*/
    rv = rl; rv -= position_;
    z  = rv * direction_ ;

    /* Compute propagation delay and subtract it from the time.                         		*/
    tl -= z / c0_;

    /* Reset the carrier envelope phase of the pulse.							*/
    p = 0.0;

    /* Now manipulate the electric field vector depending on the specific seed given.          		*/
    if ( seedType_ == PLANEWAVE )
      {
	/* Retrieve signal value at corrected time.                                                   	*/
	tsignal = signal_.self(tl, p);

	/* Calculate the field only if the signal value is larger than a limit.				*/
	if ( fabs(tsignal) < 1.0e-100 )
	  a = 0.0;
	else
	  a.mv( amplitude_ * tsignal , polarization_ );
      }
    else if ( seedType_ == PLANEWAVECONFINED )
      {
	/* Retrieve signal value at corrected time.                                                   	*/
	tsignal = signal_.self(tl, p);

	/* Calculate the transverse distance to the center line.   					*/
	x  = rv * polarization_;
	yv = cross(direction_, polarization_);
	y  = rv * yv;

	/* Calculate the field only if the signal value is larger than a limit.				*/
	if ( fabs(tsignal) < 1.0e-100  || fabs(x) > radius_[0] || fabs(y) > radius_[1] )
	  a = 0.0;
	else
	  a.mv( amplitude_ * tsignal , polarization_ );
      }
    else if ( seedType_ == GAUSSIANBEAM )
      {
	/* Retrieve signal value at corrected time.                                         		*/
	tsignal = signal_.self(tl, p);

	if ( fabs(tsignal) < 1.0e-100 ) a = 0.0;
	else
	  {
	    /* Provide vector to store transverse, longitudinal and total  electric field.      	*/
	    a = polarization_;

	    /* Calculate the transverse distance to the center line.   					*/
	    x  = rv * polarization_;
	    yv = cross(direction_, polarization_);
	    y  = rv * yv;

	    /* Calculate the wavelength corresponding to the given central frequency.              	*/
	    l = c0_ / signal_.f0_;

	    /* Calculate the Rayleigh length and the relative radius of the beam.                	*/
	    zRp = PI * radius_[0] * radius_[0] / l;
	    wrp = sqrt(1.0 + z * z / ( zRp * zRp ));
	    zRs = PI * radius_[1] * radius_[1] / l;
	    wrs = sqrt(1.0 + z * z / ( zRs * zRs ));

	    /* Compute the transverse vector between the point and the reference point.          	*/
	    p         = 0.5 * ( atan(z/zRp) + atan(z/zRs) - PI ) - PI*z/l * ( pow(x/(zRp*wrp),2) + pow(y/(zRs*wrs),2) );
	    tsignal   = signal_.self(tl, p);
	    t         = exp( - pow(x/(radius_[0]*wrp),2) - pow(y/(radius_[1]*wrs),2) ) / sqrt(wrs*wrp) * amplitude_;
	    a.mv( t * tsignal, a);
	  }
      }

    /* Now transfer the computed magnetic vector potential into the bunch rest frame.			*/
    a[2] *= gamma_;
  }

  /* Initialize the data-base for field visualization.							*/
  Seed::vtk::vtk ()
  {
    sample_			= false;
    field_.clear();
    directory_			= "";
    basename_			= "";
    rhythm_			= 0.0;
    type_			= ALLDOMAIN;
    plane_			= ZNORMAL;
    position_			= 0.0;
  }

  /* Set the sampling type of the seed.                                                               	*/
  SamplingType Seed::samplingType (std::string samplingType)
  {
    if      ( samplingType.compare("at-point")   == 0 )	return( ATPOINT  );
    else if ( samplingType.compare("over-line")  == 0 )	return( OVERLINE );
    else { std::cout << samplingType << " is an unknown sampling type." << std::endl; exit(1); }
  }

  /* Set the sampling type of the seed.                                                               	*/
  SamplingType Seed::vtkType (std::string vtkType)
  {
    if      ( vtkType.compare("in-plane") == 0 )	return( INPLANE   );
    else if ( vtkType.compare("all-domain") == 0 )	return( ALLDOMAIN );
    else { std::cout << vtkType << " is an unknown vtk type." << std::endl; exit(1); }
  }

  /* Set the plane type for vtk in plane visualization.                                           	*/
  PlaneType Seed::planeType (std::string planeType)
  {
    if      ( planeType.compare("yz")     == 0 )	return( XNORMAL );
    else if ( planeType.compare("xz")     == 0 )	return( YNORMAL );
    else if ( planeType.compare("xy")     == 0 )	return( ZNORMAL );
    else { std::cout << planeType << " is an unknown vtk plane type." << std::endl; exit(1); }
  }

  /* Set the field sampling type of the seed.                                                         */
  FieldType Seed::fieldType (std::string fieldType)
  {
    if      ( fieldType.compare("Ex") == 0 )  return(Ex);
    else if ( fieldType.compare("Ey") == 0 )  return(Ey);
    else if ( fieldType.compare("Ez") == 0 )  return(Ez);
    else if ( fieldType.compare("Bx") == 0 )  return(Bx);
    else if ( fieldType.compare("By") == 0 )  return(By);
    else if ( fieldType.compare("Bz") == 0 )  return(Bz);
    else if ( fieldType.compare("Ax") == 0 )  return(Ax);
    else if ( fieldType.compare("Ay") == 0 )  return(Ay);
    else if ( fieldType.compare("Az") == 0 )  return(Az);
    else if ( fieldType.compare("Jx") == 0 )  return(Jx);
    else if ( fieldType.compare("Jy") == 0 )  return(Jy);
    else if ( fieldType.compare("Jz") == 0 )  return(Jz);
    else if ( fieldType.compare("Q") == 0 )   return(Q);
    else if ( fieldType.compare("F") == 0 )   return(F);
    else { std::cout << fieldType << " is an unknown sampling field." << std::endl; exit(1); }
  }

  /* Show the stored values for this signal.                                                          */
  void Seed::show ()
  {
    if        (seedType_ == PLANEWAVE)      		printmessage(std::string(__FILE__), __LINE__, std::string("Seed type = plane-wave"));
    else if   (seedType_ == PLANEWAVECONFINED)   	printmessage(std::string(__FILE__), __LINE__, std::string("Seed type = confined-plane-wave"));
    else if   (seedType_ == GAUSSIANBEAM)   		printmessage(std::string(__FILE__), __LINE__, std::string("Seed type = gaussian-beam"));
    printmessage(std::string(__FILE__), __LINE__, std::string("Seed position [")
    + stringify(position_[0]) + std::string("; ")
    + stringify(position_[1]) + std::string("; ")
    + stringify(position_[2]) + std::string("]"));
    printmessage(std::string(__FILE__), __LINE__, std::string("Seed direction [")
    + stringify(direction_[0]) + std::string("; ")
    + stringify(direction_[1]) + std::string("; ")
    + stringify(direction_[2]) + std::string("]"));
    printmessage(std::string(__FILE__), __LINE__, std::string("Seed amplitude = ") + stringify(amplitude_));
    if      ( seedType_ == PLANEWAVE || seedType_ == GAUSSIANBEAM )
      printmessage(std::string(__FILE__), __LINE__, std::string("Seed polarization [")
    + stringify(polarization_[0]) + std::string("; ")
    + stringify(polarization_[1]) + std::string("; ")
    + stringify(polarization_[2]) + std::string("]"));
    if      ( seedType_ == GAUSSIANBEAM )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed parallel Rayleigh-radius = ") + stringify(radius_[0]));
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed perpendicular Rayleigh-radius = ") + stringify(radius_[1]));
      }
    if (sampling_)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed sampling at a position is enabled.") );
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed sampling directory = ") +  samplingDirectory_);
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed sampling base-name = ") +  samplingBasename_);
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed sampling rhythm = ") +  stringify(samplingRhythm_) );
	for (unsigned int i=0; i < samplingPosition_.size(); i++)
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("Seed sampling position = ") +  stringify(samplingPosition_[i]) );
	  }
      }
    for (unsigned int i = 0; i < vtk_.size(); i++)
      {
	if (vtk_[i].sample_)
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("Seed visualization in all domain is enabled.") );
	    printmessage(std::string(__FILE__), __LINE__, std::string("Seed visualization directory = ") +  vtk_[i].directory_);
	    printmessage(std::string(__FILE__), __LINE__, std::string("Seed visualization base-name = ") +  vtk_[i].basename_);
	    printmessage(std::string(__FILE__), __LINE__, std::string("Seed visualization rhythm = ") +  stringify(vtk_[i].rhythm_) );
	  }
      }
    if (profile_)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed profile saving at a time is enabled.") );
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed profile directory = ") +  profileDirectory_);
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed profile base-name = ") +  profileBasename_);
	printmessage(std::string(__FILE__), __LINE__, std::string("Seed profile rhythm = ") +  stringify(profileRhythm_) );
	for (unsigned int i=0; i < profileTime_.size(); i++)
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("Seed profile time = ") +  stringify(profileTime_[i]) );
	  }
      }
    signal_.show();
  }

  /*** Undulator class **********************************************************************************/

  Undulator::Undulator ()
  {
    k_		= 0.0;
    lu_		= 0.0;
    rb_		= 0.0;
    length_		= 0.0;
    beta_		= 0.0;
    gamma_		= 1.0;
    dt_		= 0.0;
    theta_		= 0.0;
    type_		= STATIC;
    seedType_ 	= PLANEWAVE;
    c0_		= 0.0;
    amplitude_	= 0.0;
    radius_.resize(2,0.0);
  }

  /* Set the type of the undulator.                                                                   	*/
  UndulatorType Undulator::undulatorType (std::string undulatorType)
  {
    if      ( undulatorType.compare("static")   == 0 )   return( STATIC  );
    else if ( undulatorType.compare("optical")  == 0 )   return( OPTICAL );
    else { std::cout << undulatorType << " is an unknown sampling type." << std::endl; exit(1); }
  }

  /* Initialize the data of the undulator according to the input parameters.				*/
  void Undulator::initialize (std::string        	type,
			      std::vector<Double>    	position,
			      std::vector<Double>    	direction,
			      std::vector<Double>    	polarization,
			      Double                 	amplitude,
			      std::vector<Double>	radius,
			      Double			wavelength,
			      Signal                   	signal)
  {
    /* Set the seed type according to the returned string for seedType.                   		*/
    if      ( type.compare("plane-wave"         ) == 0 ) seedType_ = PLANEWAVE;
    else if ( type.compare("plane-wave-confined") == 0 ) seedType_ = PLANEWAVECONFINED;
    else if ( type.compare("gaussian-beam"      ) == 0 ) seedType_ = GAUSSIANBEAM;
    else    { std::cout << type << " is an unknown type." << std::endl; exit(1); }

    /* Set the vectors position, direction and polarization for the seed class.                 	*/
    position_ 		= position;
    polarization_ 	= polarization;
    direction_ 		= direction;

    /* check if length of diection vector is zero and normalize the vector.                           	*/
    if ( direction_.norm() == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The direction vector has length zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
    else
      {
	Double vl = sqrt( direction_.norm() );
	direction_ /= vl;
      }

    /* check if length of polarization vector is zero and normalize the vector.                       	*/
    if ( polarization_.norm() == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The polarization vector has length zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
    else
      {
	Double vl = sqrt( polarization_.norm() );
	polarization_ /= vl;
      }

    /* check if polarization and direction are normal to each other.                       		*/
    if ( fabs( polarization_ * direction_ ) > 1.0e-50 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("Polarization is not normal to the direction."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }

    /* Initialize the amplitude of the seed.                                                    	*/
    amplitude_ 		= amplitude;

    /* Initialize the Rayleigh radius of the Gaussian beam.                                           	*/
    radius_ 		= radius;

    /* Initialize the undulator period according to the given wavelength for the signal.		*/
    lu_			= wavelength;

    /* check if the given radius of the gaussian beam makes sense.                       		*/
    if ( seedType_ == GAUSSIANBEAM && radius_[0] * radius_[1] == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("One of the radii of the Gaussian beam is set to zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }

    /* Initialize the signal of the seed.                                                       	*/
    signal_ 		= signal;
  }

  /* Show the stored values for the undulator.                                                        	*/
  void Undulator::show ()
  {
    if ( type_ == STATIC )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" Undulator type is static.") );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Undulator parameter = ") + stringify(k_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Undulator period = ") + stringify(lu_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Undulator begin = ") + stringify(rb_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Undulator length = ") + stringify(length_) );
	printmessage(std::string(__FILE__), __LINE__, std::string(" Magnetic field angle with respect to x = ") + stringify(theta_ * 180 / PI));
      }
    else if ( type_ == OPTICAL )
      {
	if 		(seedType_ == PLANEWAVE)      printmessage(std::string(__FILE__), __LINE__, std::string("Beam type = plane-wave"));
	else if 	(seedType_ == GAUSSIANBEAM)   printmessage(std::string(__FILE__), __LINE__, std::string("Beam type = gaussian-beam"));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam position [")
	+ stringify(position_[0]) + std::string("; ")
	+ stringify(position_[1]) + std::string("; ")
	+ stringify(position_[2]) + std::string("]"));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam direction [")
	+ stringify(direction_[0]) + std::string("; ")
	+ stringify(direction_[1]) + std::string("; ")
	+ stringify(direction_[2]) + std::string("]"));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam amplitude = ")
	+ stringify(amplitude_));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam polarization [")
	+ stringify(polarization_[0]) + std::string("; ")
	+ stringify(polarization_[1]) + std::string("; ")
	+ stringify(polarization_[2]) + std::string("]"));
	if      ( seedType_ == GAUSSIANBEAM )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("Beam parallel Rayleigh-radius = ") + stringify(radius_[0]));
	    printmessage(std::string(__FILE__), __LINE__, std::string("Beam perpendicular Rayleigh-radius = ") + stringify(radius_[1]));
	  }
	signal_.show();
      }
  }

  /*** ExtField class ***********************************************************************************/

  ExtField::ExtField ()
  {}

  /* Initialize the data of the undulator according to the input parameters.				*/
  void ExtField::initialize (std::string                type,
			     std::vector<Double>        position,
			     std::vector<Double>        direction,
			     std::vector<Double>        polarization,
			     Double                     amplitude,
			     std::vector<Double>        radius,
			     Double                     wavelength,
			     Signal                     signal)
  {
    /* Set the seed type according to the returned string for seedType.                               	*/
    if      ( type.compare("plane-wave"         ) == 0 ) seedType_ = PLANEWAVE;
    else if ( type.compare("plane-wave-confined") == 0 ) seedType_ = PLANEWAVECONFINED;
    else if ( type.compare("gaussian-beam"      ) == 0 ) seedType_ = GAUSSIANBEAM;
    else    { std::cout << type << " is an unknown type." << std::endl; exit(1); }

    /* Set the vectors position, direction and polarization for the seed class.                       	*/
    position_         = position;
    polarization_     = polarization;
    direction_        = direction;

    /* check if length of diection vector is zero and normalize the vector.                           	*/
    if ( direction_.norm() == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The direction vector has length zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
    else
      {
	Double vl = sqrt( direction_.norm() );
	direction_ /= vl;
      }

    /* check if length of polarization vector is zero and normalize the vector.                       	*/
    if ( polarization_.norm() == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The polarization vector has length zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }
    else
      {
	Double vl = sqrt( polarization_.norm() );
	polarization_ /= vl;
      }

    /* check if polarization and direction are normal to each other.                       		*/
    if ( fabs( polarization_ * direction_ ) > 1.0e-50 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("Polarization is not normal to the direction."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }

    /* Initialize the amplitude of the seed.                                                          	*/
    amplitude_        = amplitude;

    /* Initialize the Rayleigh radius of the Gaussian beam.                                           	*/
    radius_           = radius;

    /* check if length of polarization vector is zero and normalize the vector.                       	*/
    if ( seedType_ == GAUSSIANBEAM && radius_[0] * radius_[1] == 0.0)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("One of the radii of the Gaussian beam is set to zero."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	exit(1);
      }

    /* Initialize the signal of the seed.                                                             	*/
    signal_           = signal;
  }

  /* Show the stored values for the external field.							*/
  void ExtField::show ()
  {
    if ( type_ == EMWAVE )
      {
	if            (seedType_ == PLANEWAVE)      printmessage(std::string(__FILE__), __LINE__, std::string("Beam type = plane-wave"));
	else if       (seedType_ == GAUSSIANBEAM)   printmessage(std::string(__FILE__), __LINE__, std::string("Beam type = gaussian-beam"));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam position [")
	+ stringify(position_[0]) + std::string("; ")
	+ stringify(position_[1]) + std::string("; ")
	+ stringify(position_[2]) + std::string("]"));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam direction [")
	+ stringify(direction_[0]) + std::string("; ")
	+ stringify(direction_[1]) + std::string("; ")
	+ stringify(direction_[2]) + std::string("]"));
	printmessage(std::string(__FILE__), __LINE__, std::string("Beam amplitude = ") + stringify(amplitude_));
	if      ( seedType_ == PLANEWAVE || seedType_ == GAUSSIANBEAM )
	  printmessage(std::string(__FILE__), __LINE__, std::string("Beam polarization [")
	+ stringify(polarization_[0]) + std::string("; ")
	+ stringify(polarization_[1]) + std::string("; ")
	+ stringify(polarization_[2]) + std::string("]"));
	if      ( seedType_ == GAUSSIANBEAM )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("Beam parallel Rayleigh-radius = ") + stringify(radius_[0]));
	    printmessage(std::string(__FILE__), __LINE__, std::string("Beam perpendicular Rayleigh-radius = ") + stringify(radius_[1]));
	  }
	signal_.show();
      }
  }

  /*** FreeElectronLaser class **************************************************************************/

  /* Set the sampling type of the radiation power.                                                	*/
  void FreeElectronLaser::RadiationPower::samplingType (std::string samplingType)
  {
    if      ( samplingType.compare("at-point")   == 0 )   samplingType_ = ATPOINT;
    else if ( samplingType.compare("over-line")  == 0 )   samplingType_ = OVERLINE;
    else { std::cout << samplingType << " is an unknown sampling type." << std::endl; exit(1); }
  }

  /* Initialize the values for initializing the radiation power.                     			*/
  FreeElectronLaser::RadiationPower::RadiationPower ()
  {
    z_.clear();
    sampling_		= false;
    directory_		= "";
    basename_		= "";
    lineBegin_		= 0.0;
    lineEnd_		= 0.0;
    res_		= 0.0;
    samplingType_	= ATPOINT;
    lambda_.clear();
    lambdaMin_		= 0.0;
    lambdaMax_		= 0.0;
    lambdaRes_		= 0.0;
  }

  /* Initialize the values for initializing the radiation power.                     			*/
  FreeElectronLaser::vtk::vtk ()
  {
    z_			= 0.0;
    sampling_		= false;
    directory_		= "";
    basename_		= "";
    rhythm_		= 0.0;
  }

  /* Set the sampling type of the radiation power.                                                	*/
  void FreeElectronLaser::RadiationEnergy::samplingType (std::string samplingType)
  {
    if      ( samplingType.compare("at-point")   == 0 )   samplingType_ = ATPOINT;
    else if ( samplingType.compare("over-line")  == 0 )   samplingType_ = OVERLINE;
    else { std::cout << samplingType << " is an unknown sampling type." << std::endl; exit(1); }
  }

  /* Initialize the values for initializing the radiation energy.                     			*/
  FreeElectronLaser::RadiationEnergy::RadiationEnergy ()
  {
    z_.clear();
    sampling_		= false;
    directory_		= "";
    basename_		= "";
    lineBegin_		= 0.0;
    lineEnd_		= 0.0;
    res_		= 0.0;
    samplingType_	= ATPOINT;
    lambda_.clear();
    lambdaMin_		= 0.0;
    lambdaMax_		= 0.0;
    lambdaRes_		= 0.0;
  }
}
