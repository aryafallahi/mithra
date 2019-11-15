/********************************************************************************************************
 *  classes.hh : Implementation of the classes used in the darius code
 ********************************************************************************************************/

#ifndef CLASSES_HH_
#define CLASSES_HH_

namespace Darius
{

  /* Structure containing all the parsed parameters for mesh.                                           */
  struct Mesh
  {

    /* The parsed data related to the mesh (ls).                                                       	*/
    Double		        lengthScale_;

    /* Coordinate of the center for the computational domain (x0, y0, z0).				*/
    FieldVector<Double>		meshCenter_;

    /* Length of the mesh box in the three dimensions (lx, ly, lz).					*/
    FieldVector<Double>		meshLength_;

    /* Mesh resolution in the three dimensions (dx, dy, dz).						*/
    FieldVector<Double>		meshResolution_;

    /* The parsed data related to the time marching scheme.                                             */
    Double    		        timeScale_;
    Double		        timeStep_;
    Double		        totalTime_;

    /* Truncation order for the finite difference mesh. It can be either one or two.			*/
    unsigned int	        truncationOrder_;

    /* Boolean flag returning the status of space charge assumption.					*/
    bool			spaceCharge_;

    /* Solver type that will be used to update the field values.					*/
    SolverType			solver_;

    /* Show the stored values for the mesh.                                                      	*/
    void show()
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

    /* Initialize the parameters in the mesh initializer.						*/
    void initialize()
    {
      spaceCharge_ 	= false;
      solver_		= NSFD;
    }
  };

  /* Bunch class includes the data for the bunch properties and the functions for initializing the
   * bunches and also evaluating the bunch properties.                                              	*/
  class Bunch
  {

  public:

    typedef std::list<Charge> ChargeVector;

    /** The constructor clears the internal data structure.                                             */
    Bunch()
    {
      /* Initialize the parameters for the bunch to some first values.                      		*/
      timeStart_			= 0.0;
      timeStep_				= 0.0;

      /** Initialize the parameters for the bunch to some first values.                           	*/
      sampling_				= false;
      directory_			= "./";
      basename_             		= "";
      rhythm_				= 0.0;

      /** Initialize the parameters for the vtk visualization to some first values.               	*/
      bunchVTK_             		= false;
      bunchVTKDirectory_    		= "./";
      bunchVTKBasename_    		= "";
      bunchVTKRhythm_      		= 0.0;

      /** Initialize the parameters for saving the bunch profile to some first values.           	*/
      bunchProfile_            		= false;
      bunchProfileDirectory_    	= "./";
      bunchProfileBasename_     	= "";
      bunchProfileTime_.clear();
      bunchProfileRhythm_		= 0.0;
    };

    /****************************************************************************************************/

    /* Initialize a bunch with a manual type. This bunch produces one charge equal to the cloudCharge_. */

    void initializeManual (BunchInitialize bunchInit, ChargeVector& chargeVector, Double zp [2], int rank, int size, int ia)
    {
      /* Declare the required parameters for the initialization of charge vectors.                      */
      Charge        charge;
      Double		gb;

      /* Determine the properties of each charge point and add them to the charge vector.               */
      charge.q  	= bunchInit.cloudCharge_;
      charge.rnp    	= bunchInit.position_[ia];
      gb            	= bunchInit.initialBeta_ / sqrt( 1.0 - bunchInit.initialBeta_ * bunchInit.initialBeta_ );
      (charge.gbnp).mv(gb,bunchInit.initialDirection_);

      /* Insert this charge to the charge list if and only if it resides in the processor's portion.    */
      if ( ( charge.rnp[2] < zp[1] || rank == size - 1 ) && ( charge.rnp[2] >= zp[0] || rank == 0 ) )
	chargeVector.push_back(charge);
    }

    /****************************************************************************************************/

    /* Initialize a bunch with an ellipsoid type. This bunch produces a number of charges equal to the
     * numberOfParticles_ with the total charge equal to the cloudCharge_ which are distributed in an
     * ellipsoid with dimensions given by sigmaPosition_ and center given by the position vector. The
     * particles have uniform energy distribution centered at initialEnergy_ with variances determined by
     * sigmaGammaBeta_.                                                         			*/

    void initializeEllipsoid (BunchInitialize bunchInit, ChargeVector& chargeVector, Double zp [2], int rank, int size, int ia)
    {
      /* Save the initially given number of particles.							*/
      unsigned int	Np = bunchInit.numberOfParticles_, i, Np0 = chargeVector.size();

      /* Declare the required parameters for the initialization of charge vectors.                      */
      Charge            charge; charge.q  = bunchInit.cloudCharge_ / Np;
      Double gb 	= bunchInit.initialBeta_ / sqrt( 1.0 - bunchInit.initialBeta_ * bunchInit.initialBeta_);
      FieldVector<Double> r (0.0);
      FieldVector<Double> t (0.0);
      Double            t0, t1, t2 = -1.0;
      Double		zmin = 1e100;
      Double		Ne, bF, bFi;
      unsigned int	bmi;

      /* Declare the function for injecting the shot noise.						*/
      auto insertCharge = [&] (Charge q) {

	for ( unsigned int ii = 0; ii < 4; ii++ )
	  {
	    /* The random modulation is introduced depending on the shot-noise being activated.		*/
	    if ( bunchInit.shotNoise_ )
	      {
		/* Obtain the number of beamlet.							*/
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

	    /* Insert this charge to the charge list if and only if it resides in the processor's
	     * portion.                                                                               	*/
	    if ( ( q.rnp[2] < zp[1] || rank == size - 1 ) && ( q.rnp[2] >= zp[0] || rank == 0 ) )
	      chargeVector.push_back(q);
	  }
      };

      /* Check the bunching factor.                                                                     */
      if ( bunchInit.bF_ > 2.0 || bunchInit.bF_ < 0.0 )
	{
	  printmessage(std::string(__FILE__), __LINE__, std::string("The bunching factor can not be larger than one or a negative value !!!") );
	  exit(1);
	}

      /* If the shot noise is on, we need the minimum value of the bunch z coordinate to be able to
       * calculate the FEL bucket number.								*/
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
	    for ( ; i < int( Np / 4 * ( 1.0 + bunchInit.lambda_ * sqrt( 2.0 * PI ) / ( 2.0 * bunchInit.sigmaPosition_[2] ) ) ); i++)
	      {
		t0  = bunchInit.lambda_ * sqrt( - 2.0 * log( halton( 2, i + Np0 ) ) ) * sin( 2.0 * PI * halton( 3, i + Np0 ) );
		t0 += ( t0 < 0.0 ) ? ( - bunchInit.sigmaPosition_[2] ) : ( bunchInit.sigmaPosition_[2] );

		zmin = std::min(   t0 , zmin );
	      }

	  zmin = zmin + bunchInit.position_[ia][2];

	  /* Obtain the average number of electrons per FEL beamlet.					*/
	  Ne = bunchInit.cloudCharge_ * ( bunchInit.lambda_ / 2.0 ) / ( 2.0 * bunchInit.sigmaPosition_[2] );

	  /* Set the bunching factor level for the shot noise depending on the given values.		*/
	  bF = ( bunchInit.bF_ == 0.0 ) ? 1.0 / sqrt(Ne) : bunchInit.bF_; //bF /= sqrt(2.0);

	  printmessage(std::string(__FILE__), __LINE__, std::string("The standard deviation of the bunching factor for the shot noise implementation is set to ") + stringify(bF) );
	}

      /* Determine the properties of each charge point and add them to the charge vector.               */
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

	  /* Determine the longitudinal momentum.							*/
	  if ( bunchInit.distribution_ == "uniform" )
	    t[2] = ( 2.0 * halton(6, i + Np0) - 1.0 ) * bunchInit.sigmaGammaBeta_[2];
	  else if ( bunchInit.distribution_ == "gaussian" )
	    t[2] = bunchInit.sigmaGammaBeta_[2] * sqrt( - 2.0 * log( halton(6, i + Np0) ) ) * cos( 2.0 * PI * halton(7, i + Np0) );

	  if ( fabs(r[0]) < bunchInit.tranTrun_ && fabs(r[1]) < bunchInit.tranTrun_ && fabs(r[2]) < bunchInit.longTrun_)
	    {
	      /* Shift the generated charge to the center position and momentum space.			*/
	      charge.rnp    = bunchInit.position_[ia];
	      charge.rnp   += r;

	      (charge.gbnp).mv(gb, bunchInit.initialDirection_);
	      charge.gbnp  += t;

	      /* Insert this charge and the mirrored ones into the charge vector.			*/
	      insertCharge(charge);
	    }
	}

      /* If the longitudinal type of the bunch is uniform a tapered part needs to be added to remove the
       * CSE from the tail of the bunch.								*/
      if ( bunchInit.distribution_ == "uniform" )
	for ( ; i < int( Np / 4 * ( 1.0 + bunchInit.lambda_ * sqrt( 2.0 * PI ) / ( 2.0 * bunchInit.sigmaPosition_[2] ) ) ); i++)
	  {
	    r[0] = bunchInit.sigmaPosition_[0] * sqrt( - 2.0 * log( halton(0, i + Np0) ) ) * cos( 2.0 * PI * halton(1, i + Np0) );
	    r[1] = bunchInit.sigmaPosition_[1] * sqrt( - 2.0 * log( halton(0, i + Np0) ) ) * sin( 2.0 * PI * halton(1, i + Np0) );

	    /* Determine the longitudinal coordinate.							*/
	    r[2] = bunchInit.lambda_ * sqrt( - 2.0 * log( halton(2, i + Np0) ) ) * sin( 2.0 * PI * halton(3, i + Np0) );
	    r[2] += ( r[2] < 0.0 ) ? ( - bunchInit.sigmaPosition_[2] ) : ( bunchInit.sigmaPosition_[2] );

	    /* Determine the transverse momentum.							*/
	    t[0] = bunchInit.sigmaGammaBeta_[0] * sqrt( - 2.0 * log( halton(4, i + Np0) ) ) * cos( 2.0 * PI * halton(5, i + Np0) );
	    t[1] = bunchInit.sigmaGammaBeta_[1] * sqrt( - 2.0 * log( halton(4, i + Np0) ) ) * sin( 2.0 * PI * halton(5, i + Np0) );

	    /* Determine the longitudinal momentum.							*/
	    if ( bunchInit.distribution_ == "uniform" )
	      t[2] = ( 2.0 * halton(6, i + Np0) - 1.0 ) * bunchInit.sigmaGammaBeta_[2];
	    else if ( bunchInit.distribution_ == "gaussian" )
	      t[2] = bunchInit.sigmaGammaBeta_[2] * sqrt( - 2.0 * log( halton(6, i + Np0) ) ) * cos( 2.0 * PI * halton(7, i + Np0) );

	    if ( fabs(r[0]) < bunchInit.tranTrun_ && fabs(r[1]) < bunchInit.tranTrun_ && fabs(r[2]) < bunchInit.longTrun_)
	      {
		/* Shift the generated charge to the center position and momentum space.		*/
		charge.rnp   = bunchInit.position_[ia];
		charge.rnp  += r;

		(charge.gbnp).mv(gb, bunchInit.initialDirection_);
		charge.gbnp += t;

		/* Insert this charge and the mirrored ones into the charge vector.			*/
		insertCharge(charge);
	      }
	  }

      /* Reset the value for the number of particle variable according to the installed number of
       * macro-particles and perform the corresponding changes.                                         */
      bunchInit.numberOfParticles_ = chargeVector.size();
    }

    /****************************************************************************************************/

    /* Initialize a bunch with a 3D-crystal type. This bunch produces a number of charges equal
     * to the numberOfParticles_ with the total charge equal to the cloudCharge_ which are arranged in a
     * 3D crystal. The number of particles in each direction is given by numbers_. Therefore,
     * numberOfParticles_ should be a multiple of the product of all these three numbers. The ratio gives
     * the number of particles in each crystal point. Position of each particle is determined by the
     * lattice constants and the crystal is centered at the position_ vector. At each point, the charges
     * have a small Gaussian distribution around the crsytal.                                           */

    void initialize3DCrystal (BunchInitialize bunchInit, ChargeVector& chargeVector, Double zp [2], int rank, int size, int ia)
    {
      /* Check if the numberOfParticles_ is a multiple of the product of values in numbers_.            */
      if ( bunchInit.numberOfParticles_ % (bunchInit.numbers_[0] * bunchInit.numbers_[1] * bunchInit.numbers_[2]) != 0 )
	{
	  printmessage(std::string(__FILE__), __LINE__,
		       std::string("The number of the particles and their lattice numbers do not match !!!") );
	  exit(1);
	}

      /* Declare the required parameters for the initialization of charge vectors.                      */
      Charge            charge;
      Double gb       = bunchInit.initialBeta_ / sqrt( 1.0 - bunchInit.initialBeta_ * bunchInit.initialBeta_);
      unsigned int np = bunchInit.numberOfParticles_ / (bunchInit.numbers_[0] * bunchInit.numbers_[1] * bunchInit.numbers_[2]);

      /* Clear the charge vector for adding the charges.						*/
      chargeVector.clear();

      /* Determine the properties of each charge point and add them to the charge vector.               */
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

		      (charge.gbnp).mv(gb, bunchInit.initialDirection_);
		      charge.gbnp[0] += bunchInit.sigmaGammaBeta_[0] * sqrt( - 2.0 * log( halton(6,i) ) ) * sin( 2.0 * PI * halton(7,i) );
		      charge.gbnp[1] += bunchInit.sigmaGammaBeta_[1] * sqrt( - 2.0 * log( halton(8,i) ) ) * sin( 2.0 * PI * halton(9,i) );
		      charge.gbnp[2] += bunchInit.sigmaGammaBeta_[2] * sqrt( - 2.0 * log( halton(10,i)) ) * sin( 2.0 * PI * halton(11,i));

		      /* Insert this charge to the charge list if and only if it resides in the processor's portion.    */
		      if ( ( charge.rnp[2] < zp[1] || rank == size - 1 ) && ( charge.rnp[2] >= zp[0] || rank == 0 ) )
			chargeVector.push_back(charge);
		    }
		}
	    }
	}
    }

    /****************************************************************************************************/

    /* Initialize a bunch with a file type. This bunch produces a number of charges read from a given file.
     * The number of initialized charge is equal to the vertical length of the table in the text file. The
     * file format should contain the charge value, 3 position coordinates and 3 momentum coordinates of
     * of the charge distribution.							                */

    void initializeFile (BunchInitialize bunchInit, ChargeVector& chargeVector, Double zp [2], int rank, int size, int ia)
    {

      /* Declare the required parameters for the initialization of charge vectors.                	*/
      Charge                    charge;
      FieldVector<Double>	gb (0.0);
      FieldVector<Double>	r  (0.0);
      Double 			gbm = bunchInit.initialBeta_ / sqrt( 1.0 - bunchInit.initialBeta_ * bunchInit.initialBeta_);

      /* Clear the charge vector for adding the charges.						*/
      chargeVector.clear();

      /* Read the file and fill the position and electric field vectors according to the saved values.	*/
      std::ifstream myfile ( bunchInit.fileName_.c_str() );

      while (myfile.good())
	{
	  charge.q  = bunchInit.cloudCharge_ / bunchInit.numberOfParticles_;

	  myfile >> r[0];
	  myfile >> r[1];
	  myfile >> r[2];

	  myfile >> gb[0];
	  myfile >> gb[1];
	  myfile >> gb[2];

	  charge.rnp   = bunchInit.position_[ia];
	  charge.rnp  += r;

	  (charge.gbnp).mv(gbm, bunchInit.initialDirection_);
	  charge.gbnp += gb;

	  /* Insert this charge to the charge list if and only if it resides in the processor's portion.*/
	  if ( ( charge.rnp[2] < zp[1] || rank == size - 1 ) && ( charge.rnp[2] >= zp[0] || rank == 0 ) )
	    chargeVector.push_back(charge);

	  if ( fabs(charge.rnp[0]) > bunchInit.tranTrun_ || fabs(charge.rnp[1]) > bunchInit.tranTrun_ )
	    printmessage(std::string(__FILE__), __LINE__, std::string("Warning: The particle coordinate is out of the truncation bunch. "
		"The results may be inaccurate !!!") );
	}

      /* Calculate the total amount of installed particles.						*/
      unsigned int NqL = chargeVector.size(), NqG = 0;
      MPI_Reduce(&NqL,&NqG,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);

      /* Check the size of the charge vector with the number of particles.                            	*/
      if ( bunchInit.numberOfParticles_ != NqG && rank == 0 )
	{
	  printmessage(std::string(__FILE__), __LINE__, std::string("The number of the particles and the file size do not match !!!") );
	  exit(1);
	}
    }

    /****************************************************************************************************/

    /* The data containing the bunch initialization parameters.						*/
    std::vector<BunchInitialize>	bunchInit_;

    /* The directory parsed for the whole project.                                                      */
    std::string				directory_;

    /* Base name for writing the outputs of the electron acceleration analysis.                         */
    std::string				basename_;

    /* Boolean variable that determines if the bunch sampling should be done or not.			*/
    bool				sampling_;

    /* Store the time step for updating the electron motion.						*/
    Double				timeStep_;

    /* Rhythm of writing the bunch macroscopic values in the output file.                              	*/
    Double         			rhythm_;

    /* Time points at which the bunch is produced.                                                      */
    Double                    		timeStart_;

    /* Boolean parameter that determines if the vtk visualization should be done.                      	*/
    bool				bunchVTK_;

    /* The directory in which the bunch vtk files should be saved.                                   	*/
    std::string				bunchVTKDirectory_;

    /* Name of the files in which the vtk visualization should be saved.                               	*/
    std::string				bunchVTKBasename_;

    /* Rhythm of producing the vtk files. It should be double value bigger than the time step.         	*/
    Double				bunchVTKRhythm_;

    /* Boolean parameter that determines if the bunch profile should be saved.                         	*/
    bool				bunchProfile_;

    /* The directory in which the bunch profile should be saved.                                   	*/
    std::string				bunchProfileDirectory_;

    /* Name of the files in which the bunch profile should be saved.                                   	*/
    std::string				bunchProfileBasename_;

    /* Vector of time points at which the bunch profile should be saved.                               	*/
    std::vector<Double>      		bunchProfileTime_;

    /* Rhythm of saving the bunch profile. It should be a double value bigger than the time step.	*/
    Double				bunchProfileRhythm_;

    /* Show the stored values for the bunch.                                                          	*/
    void show()
    {
      for (unsigned int i = 0; i < bunchInit_.size(); i++)
	{
	  printmessage(std::string(__FILE__), __LINE__, std::string(" Type of the bunch initialization = ") + bunchInit_[i].bunchType_ );
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
	  printmessage(std::string(__FILE__), __LINE__, std::string(" Initialization time of the bunch = ") + stringify(timeStart_) );
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

  };

  /* Define the main signal class.                                                                      */
  class Signal
  {
  public:
    Signal()
  {
      t0_ 		= 0.0;
      s_  		= 0.0;
      f0_ 		= 0.0;
      cep_		= 0.0;
      signalType_	= GAUSSIAN;
  };

    /* Initializer with signal type, time offset, variance, frequency and carrier-envelope-phase.       */
    void initialize (std::string type, Double l0, Double s, Double l, Double cep)
    {
      /* Initialize the signal type.                                                                    */
      if      ( type.compare("neumann") == 0 )             signalType_ = NEUMANN;
      else if ( type.compare("gaussian") == 0 )            signalType_ = GAUSSIAN;
      else if ( type.compare("secant-hyperbolic") == 0 )   signalType_ = SECANT;
      else if ( type.compare("flat-top") == 0 )   	   signalType_ = FLATTOP;
      else { std::cout << type << " is an unknown signal type for the given set of parameters." << std::endl; exit(1); }

      /* Initialize the time delay, variance and the frequency of the carrier.                          */
      t0_ = l0;
      s_  = s;
      f0_ = 1 / l;

      /* Initialize the carrier envelope phase.                                                         */
      cep_ = cep * PI / 180;

      /* Check if variance is unequal zero.                                                             */
      if (s_ == 0.0)
	{
	  printmessage(std::string(__FILE__), __LINE__, std::string(" Variance of signal is set to zero. "));
	  printmessage(std::string(__FILE__), __LINE__, std::string(" This is not allowed because we divide through the variance. "));
	  printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	  exit(1);
	}
    };

  public:

    /* Store type of the signal.                                                                        */
    SignalType			signalType_;

    /* Time offset of the signal.                                                                       */
    Double     			t0_;

    /* Variance of the signal. Variance is defined according to the point where the intensity of the
     * signal, i.e. signal squared is half of the maximum.                                              */
    Double     			s_;

    /* Frequency of the modulation.                                                                     */
    Double     			f0_;

    /* Carrier envelope phase of the modulation.                                                        */
    Double     			cep_;

    /* Provide the signal at time t.                                                                    */
    Double self (Double& t, Double& phase)
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

    /* Show the stored values for this signal.                                                          */
    void show()
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
  };

  /* Define the main seed class.                                                                  	*/
  class Seed
  {

  public:

    Seed()
  {
      seedType_ 		= PLANEWAVE;
      c0_			= 0.0;
      amplitude_		= 0.0;
      radius_.resize(2,0.0);

      sampling_			= false;
      samplingType_		= ATPOINT;
      samplingField_.clear();
      samplingDirectory_	= "";
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
  };

    void initialize (std::string        	type,
		     std::vector<Double>    	position,
		     std::vector<Double>    	direction,
		     std::vector<Double>    	polarization,
		     Double                 	amplitude,
		     std::vector<Double>	radius,
		     Signal                   	signal)
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

      /* check if length of diection vector is zero and normalize the vector.                           */
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

      /* check if length of polarization vector is zero and normalize the vector.                       */
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
      amplitude_ 	= amplitude;

      /* Initialize the Rayleigh radius of the Gaussian beam.                                           */
      radius_ 		= radius;

      /* check if length of polarization vector is zero and normalize the vector.                       */
      if ( seedType_ == GAUSSIANBEAM && radius_[0] * radius_[1] == 0.0)
	{
	  printmessage(std::string(__FILE__), __LINE__, std::string("One of the radii of the Gaussian beam is set to zero."));
	  printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	  exit(1);
	}

      /* Initialize the signal of the seed.                                                       	*/
      signal_ 		= signal;
    };

  public:

    /* Store type of the seed.                                                                    	*/
    SeedType           			seedType_;

    /* Speed of light value in terms of the given length-scale and time-scale.				*/
    Double				c0_;

    /* Store reference position of the seed. For waves, it is a reference position and for
     * hertzian dipole, it is the position of the dipole.                                               */
    FieldVector<Double>      		position_;

    /* Store the direction of the seed. For waves, it is the propagation direction and for a
     * hertzian dipole, it is the direction of the dipole.                                              */
    FieldVector<Double>       		direction_;

    /* Store the polarization of the wave in the seed.                                            	*/
    FieldVector<Double>       		polarization_;

    /* Store the amplitude of the seed.                                                           	*/
    Double                  		amplitude_;

    /* Store the Rayleigh radius of the Gaussian beam in the parallel and perpendicular directions.	*/
    std::vector<Double>      		radius_;

    /* Store the signal class for this seed.                                                      	*/
    Signal                   		signal_;

    /* Parameters for Lorentz transformation.								*/
    Double				beta_;
    Double				gamma_;
    Double				dt_;

  private:

    /* Store all the required variables in the computations.                                            */
    Double                    		gamma;
    Double				tsignal;
    Double				d, l, zRp, wrp, zRs, wrs, x, y, z, p, t;
    FieldVector<Double>			rv, yv, ax, az;
    FieldVector<Double>			rl;
    Double				tl;

  public:

    /* Return the potentials at any desired location and time.                    			*/
    void fields (FieldVector<Double>& aufpunkt, Double& time, FieldVector<Double>& a)
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

      /* Now manipulate the electric field vector depending on the specific seed given.          	*/
      if ( seedType_ == PLANEWAVE )
	{
	  /* Retrieve signal value at corrected time.                                                   */
	  tsignal = signal_.self(tl, p);

	  /* Calculate the field only if the signal value is larger than a limit.			*/
	  if ( fabs(tsignal) < 1.0e-100 )
	    a = 0.0;
	  else
	    a.mv( amplitude_ * tsignal , polarization_ );
	}
      else if ( seedType_ == PLANEWAVECONFINED )
	{
	  /* Retrieve signal value at corrected time.                                                   */
	  tsignal = signal_.self(tl, p);

	  /* Calculate the transverse distance to the center line.   					*/
	  x  = rv * polarization_;
	  yv = cross(direction_, polarization_);
	  y  = rv * yv;

	  /* Calculate the field only if the signal value is larger than a limit.			*/
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

	      /* Calculate the transverse distance to the center line.   				*/
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

    /* Store the data required for sampling the radiated field in a point.                              */
    bool                                sampling_;
    SamplingType                        samplingType_;
    std::vector<FieldType>              samplingField_;
    std::string                         samplingDirectory_;
    std::string                         samplingBasename_;
    Double                        	samplingRhythm_;
    std::vector<FieldVector<Double> >   samplingPosition_;
    FieldVector<Double>                 samplingLineBegin_;
    FieldVector<Double>                 samplingLineEnd_;
    FieldVector<Double>                 samplingSurfaceBegin_;
    FieldVector<Double>                 samplingSurfaceEnd_;
    unsigned int                        samplingRes_;

    /* Store data required for visualizing the radiated field in all-domain.                            */
    struct vtk
    {
      bool                              sample_;
      std::vector<FieldType>            field_;
      std::string                       directory_;
      SamplingType			type_;
      std::string                       basename_;
      Double                      	rhythm_;
      PlaneType				plane_;
      FieldVector<Double>		position_;

      vtk()
      {
	sample_				= false;
	field_.clear();
	directory_			= "";
	basename_			= "";
	rhythm_				= 0.0;
	type_				= ALLDOMAIN;
	plane_				= ZNORMAL;
	position_			= 0.0;
      }
    };
    std::vector<vtk> 			vtk_;


    /* Store data required for writing the profile of the field in all-domain.                          */
    bool                                profile_;
    std::vector<FieldType>              profileField_;
    std::string                         profileDirectory_;
    std::string                         profileBasename_;
    std::vector<Double>                 profileTime_;
    Double                        	profileRhythm_;

    /* Set the sampling type of the seed.                                                               */
    SamplingType samplingType(std::string samplingType)
    {
      if      ( samplingType.compare("at-point")   == 0 )	return( ATPOINT  );
      else if ( samplingType.compare("over-line")  == 0 )	return( OVERLINE );
      else { std::cout << samplingType << " is an unknown sampling type." << std::endl; exit(1); }
    }

    /* Set the sampling type of the seed.                                                               */
    SamplingType vtkType(std::string vtkType)
    {
      if      ( vtkType.compare("in-plane") == 0 )	return( INPLANE   );
      else if ( vtkType.compare("all-domain") == 0 )	return( ALLDOMAIN );
      else { std::cout << vtkType << " is an unknown vtk type." << std::endl; exit(1); }
    }

    /* Set the plane type for vtk in plane visualization.                                           	*/
    PlaneType planeType(std::string planeType)
    {
      if      ( planeType.compare("yz")     == 0 )	return( XNORMAL );
      else if ( planeType.compare("xz")     == 0 )	return( YNORMAL );
      else if ( planeType.compare("xy")     == 0 )	return( ZNORMAL );
      else { std::cout << planeType << " is an unknown vtk plane type." << std::endl; exit(1); }
    }

    /* Set the field sampling type of the seed.                                                         */
    FieldType fieldType(std::string fieldType)
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
    void show()
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
  };

  /* Define the structure containing the main parameters for the undulator.				*/
  class Undulator
  {

  public:

    Undulator()
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
  };

    /* The magnetic field of the undulator.                                                       	*/
    Double				k_;

    /* The period of the undulator.                                                       		*/
    Double				lu_;

    /* The start position of the undulator.								*/
    Double				rb_;

    /* The length of the undulator.									*/
    unsigned int			length_;

    /* The normalized velocity and the equivalent gamma of the undulator movement.			*/
    Double				beta_;
    Double				gamma_;
    Double				dt_;

    /* Angle of the undulator polarization with respect to x axis.					*/
    Double				theta_;

    /* Type of the undulator, it can be an optical or a static undulator.				*/
    UndulatorType			type_;

    /* Store type of the seed.                                                                    	*/
    SeedType           			seedType_;

    /* Speed of light value in terms of the given length-scale and time-scale.				*/
    Double				c0_;

    /* Store reference position of the seed. For waves, it is a reference position and for
     * hertzian dipole, it is the position of the dipole.                                               */
    FieldVector<Double>      		position_;

    /* Store the direction of the seed. For waves, it is the propagation direction and for a
     * hertzian dipole, it is the direction of the dipole.                                              */
    FieldVector<Double>       		direction_;

    /* Store the polarization of the wave in the seed.                                            	*/
    FieldVector<Double>       		polarization_;

    /* Store the amplitude of the seed.                                                           	*/
    Double                  		amplitude_;

    /* Store the Rayleigh radius of the Gaussian beam in the parallel and perpendicular directions.	*/
    std::vector<Double>      		radius_;

    /* Store the signal class for this seed.                                                      	*/
    Signal                   		signal_;

    /* Set the type of the undulator.                                                                   */
    UndulatorType undulatorType(std::string undulatorType)
    {
      if      ( undulatorType.compare("static")   == 0 )   return( STATIC  );
      else if ( undulatorType.compare("optical")  == 0 )   return( OPTICAL );
      else { std::cout << undulatorType << " is an unknown sampling type." << std::endl; exit(1); }
    }

    void initialize (std::string        	type,
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
      position_ 	= position;
      polarization_ 	= polarization;
      direction_ 	= direction;

      /* check if length of diection vector is zero and normalize the vector.                           */
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

      /* check if length of polarization vector is zero and normalize the vector.                       */
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
      amplitude_ 	= amplitude;

      /* Initialize the Rayleigh radius of the Gaussian beam.                                           */
      radius_ 		= radius;

      /* Initialize the undulator period according to the given wavelength for the signal.		*/
      lu_		= wavelength;

      /* check if the given radius of the gaussian beam makes sense.                       		*/
      if ( seedType_ == GAUSSIANBEAM && radius_[0] * radius_[1] == 0.0)
	{
	  printmessage(std::string(__FILE__), __LINE__, std::string("One of the radii of the Gaussian beam is set to zero."));
	  printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	  exit(1);
	}

      /* Initialize the signal of the seed.                                                       	*/
      signal_ 		= signal;
    };

    /* Show the stored values for the undulator.                                                        */
    void show()
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
  };

  /* Define the structure containing the main parameters for the undulator.                           */
  class ExtField
  {

  public:

    ExtField()
  {
  };

    /* Type of the undulator, it can be an optical or a static undulator.                               */
    ExtFieldType                        type_;

    /* Store type of the seed.                                                                          */
    SeedType                            seedType_;

    /* Speed of light value in terms of the given length-scale and time-scale.                          */
    Double                              c0_;

    /* Store reference position of the seed. For waves, it is a reference position and for
     * hertzian dipole, it is the position of the dipole.                                               */
    FieldVector<Double>                 position_;

    /* Store the direction of the seed. For waves, it is the propagation direction and for a
     * hertzian dipole, it is the direction of the dipole.                                              */
    FieldVector<Double>                 direction_;

    /* Store the polarization of the wave in the seed.                                                  */
    FieldVector<Double>                 polarization_;

    /* Store the amplitude of the seed.                                                                 */
    Double                              amplitude_;

    /* Store the Rayleigh radius of the Gaussian beam in the parallel and perpendicular directions.     */
    std::vector<Double>                 radius_;

    /* Store the signal class for this seed.                                                            */
    Signal                              signal_;

    void initialize (std::string                type,
		     std::vector<Double>        position,
		     std::vector<Double>        direction,
		     std::vector<Double>        polarization,
		     Double                     amplitude,
		     std::vector<Double>        radius,
		     Double                     wavelength,
		     Signal                     signal)
    {
      /* Set the seed type according to the returned string for seedType.                               */
      if      ( type.compare("plane-wave"         ) == 0 ) seedType_ = PLANEWAVE;
      else if ( type.compare("plane-wave-confined") == 0 ) seedType_ = PLANEWAVECONFINED;
      else if ( type.compare("gaussian-beam"      ) == 0 ) seedType_ = GAUSSIANBEAM;
      else    { std::cout << type << " is an unknown type." << std::endl; exit(1); }

      /* Set the vectors position, direction and polarization for the seed class.                       */
      position_         = position;
      polarization_     = polarization;
      direction_        = direction;

      /* check if length of diection vector is zero and normalize the vector.                           */
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

      /* check if length of polarization vector is zero and normalize the vector.                       */
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

      /* Initialize the amplitude of the seed.                                                          */
      amplitude_        = amplitude;

      /* Initialize the Rayleigh radius of the Gaussian beam.                                           */
      radius_           = radius;

      /* check if length of polarization vector is zero and normalize the vector.                       */
      if ( seedType_ == GAUSSIANBEAM && radius_[0] * radius_[1] == 0.0)
	{
	  printmessage(std::string(__FILE__), __LINE__, std::string("One of the radii of the Gaussian beam is set to zero."));
	  printmessage(std::string(__FILE__), __LINE__, std::string("Exit!"));
	  exit(1);
	}

      /* Initialize the signal of the seed.                                                             */
      signal_           = signal;
    };

    /* Show the stored values for the undulator.                                                        */
    void show()
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
  };

  /* Structure containing all the parsed parameters for the free electron laser.			*/
  struct FreeElectronLaser
  {

    /* Structure contianing the parsed parameters for the radiation power.				*/
    struct RadiationPower
    {
      /* The distance of the planes from the bunch to sample the radiation power.			*/
      std::vector<Double>			z_;

      /* Store the data required for sampling the radiation power.					*/
      bool					sampling_;

      /* Store the directory in which the radiation data is saved.					*/
      std::string				directory_;

      /* Store the base-name of the file in which the radiation data is saved.				*/
      std::string				basename_;

      /* The begin and end of the line as well as the resolution on which the radiation data is saved.	*/
      Double                 			lineBegin_;
      Double                 			lineEnd_;
      unsigned int                              res_;

      /* Store the sampling type of the radiation power.						*/
      SamplingType                        	samplingType_;

      /* The wavelength of the harmonic whose power should be plotted.					*/
      std::vector<Double>			lambda_;

      /* The wavelength sweep data for the power computation.						*/
      Double					lambdaMin_;
      Double					lambdaMax_;
      unsigned int				lambdaRes_;


      /* Set the sampling type of the radiation power.                                                	*/
      void samplingType(std::string samplingType)
      {
	if      ( samplingType.compare("at-point")   == 0 )   samplingType_ = ATPOINT;
	else if ( samplingType.compare("over-line")  == 0 )   samplingType_ = OVERLINE;
	else { std::cout << samplingType << " is an unknown sampling type." << std::endl; exit(1); }
      }

      /* Initialize the values for initializing the radiation power.                     		*/
      RadiationPower()
      {
	z_.clear();
	sampling_	= false;
	directory_	= "";
	basename_	= "";
	lineBegin_	= 0.0;
	lineEnd_	= 0.0;
	res_		= 0.0;
	samplingType_	= ATPOINT;
	lambda_.clear();
	lambdaMin_	= 0.0;
	lambdaMax_	= 0.0;
	lambdaRes_	= 0.0;
      }
    } 						radiationPower_;

    /* Structure containing the parsed parameters for the power-visualization.				*/
    struct vtk
    {
      /* The distance of the plane from the bunch to sample the radiation power.			*/
      Double					z_;

      /* Store the data required for sampling the radiation power.					*/
      bool					sampling_;

      /* Store the directory in which the radiation data is saved.					*/
      std::string				directory_;

      /* Store the base-name of the file in which the radiation data is saved.				*/
      std::string				basename_;

      /* Store the rhythm for calculating and saving the radiation power.				*/
      Double					rhythm_;

      /* The wavelength of the harmonic whose power should be plotted.					*/
      Double					lambda_;

      /* Initialize the values for initializing the radiation power.                     		*/
      vtk()
      {
	z_		= 0.0;
	sampling_	= false;
	directory_	= "";
	basename_	= "";
	rhythm_		= 0.0;
      }
    } 						vtk_;

    /* Structure contianing the parsed parameters for the radiation energy.				*/
    struct RadiationEnergy
    {
      /* The distance of the planes from the bunch to sample the radiation power.			*/
      std::vector<Double>			z_;

      /* Store the data required for sampling the radiation power.					*/
      bool					sampling_;

      /* Store the directory in which the radiation data is saved.					*/
      std::string				directory_;

      /* Store the base-name of the file in which the radiation data is saved.				*/
      std::string				basename_;

      /* Store the rhythm for calculating and saving the radiation power.				*/
      Double					rhythm_;

      /* The begin and end of the line as well as the resolution on which the radiation data is saved.	*/
      Double                 			lineBegin_;
      Double                 			lineEnd_;
      Double                              	res_;

      /* Store the sampling type of the radiation power.						*/
      SamplingType                        	samplingType_;

      /* The wavelength of the harmonic whose power should be plotted.					*/
      std::vector<Double>			lambda_;

      /* The wavelength sweep data for the power computation.						*/
      Double					lambdaMin_;
      Double					lambdaMax_;
      Double					lambdaRes_;

      /* Set the sampling type of the radiation power.                                                	*/
      void samplingType(std::string samplingType)
      {
	if      ( samplingType.compare("at-point")   == 0 )   samplingType_ = ATPOINT;
	else if ( samplingType.compare("over-line")  == 0 )   samplingType_ = OVERLINE;
	else { std::cout << samplingType << " is an unknown sampling type." << std::endl; exit(1); }
      }

      /* Initialize the values for initializing the radiation energy.                     		*/
      RadiationEnergy()
      {
	z_.clear();
	sampling_	= false;
	directory_	= "";
	basename_	= "";
	lineBegin_	= 0.0;
	lineEnd_	= 0.0;
	res_		= 0.0;
	samplingType_	= ATPOINT;
	lambda_.clear();
	lambdaMin_	= 0.0;
	lambdaMax_	= 0.0;
	lambdaRes_	= 0.0;
      }
    } 						radiationEnergy_;


  };
}       /* End of namespace Cyrus3d.                                                                    */
#endif
