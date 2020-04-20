/********************************************************************************************************
 *  solver.cpp : Implementation of the functions for the solver class in the darius code
 ********************************************************************************************************/

#include <algorithm>

#include "solver.h"

namespace MITHRA
{
  Solver::Solver (Mesh& 				mesh,
		  Bunch& 				bunch,
		  Seed& 				seed,
		  std::vector<Undulator>&		undulator,
		  std::vector<ExtField>& 		extField,
		  std::vector<FreeElectronLaser>& 	FEL)
  : mesh_ 	( mesh ),
    bunch_ 	( bunch ),
    seed_ 	( seed ),
    undulator_ 	( undulator ),
    extField_ 	( extField ),
    FEL_ 	( FEL )
  {
    /* Clear the vectors in the FdTd class.								*/
    anp1_ = new std::vector<FieldVector<Double> > ();
    an_   = new std::vector<FieldVector<Double> > ();
    anm1_ = new std::vector<FieldVector<Double> > ();

    fnp1_ = new std::vector<Double> ();
    fn_   = new std::vector<Double> ();
    fnm1_ = new std::vector<Double> ();

    r_.clear();

    /* Reset the number of nodes to zero.								*/
    N0_ = N1_ = N2_ = 0;

    /* Reset the time and the time number of the FdTd class.						*/
    timep1_	   =  mesh_.timeStep_;
    time_  	   =  0.0;
    timem1_	   = -mesh_.timeStep_;
    nTime_ 	   =  0;
    nTimeBunch_    =  0;

    /* Initialize the value of MPI variables.								*/
    MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
    MPI_Comm_size(MPI_COMM_WORLD,&size_);

    /* Initialize the speed of light value according to the given length scale and time scale.		*/
    c0_ = C0 / mesh_.lengthScale_ * mesh_.timeScale_;
    m0_ = MU_ZERO / mesh_.lengthScale_;
    e0_ = 1.0 / ( c0_ * c0_ * m0_ );
  }

  /******************************************************************************************************
   * Using the parsed data, set the required parameters for simulation.
   ******************************************************************************************************/

  void Solver::setSimulationParameters ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string("::: Boosting the given mesh parameters into the electron rest frame ") );

    /****************************************************************************************************/

    /* Initialize the corresponding length and time scales in other parts of the solver.		*/
    seed_.c0_ 			 = c0_;
    seed_.signal_.t0_ 		/= c0_;
    seed_.signal_.f0_ 		*= c0_;
    seed_.signal_.s_          	/= c0_;
    for (std::vector<Undulator>::iterator iter = undulator_.begin(); iter != undulator_.end(); iter++)
      {
	iter->c0_ 	       	 = c0_;
	iter->signal_.t0_     	/= c0_;
	iter->signal_.f0_ 	*= c0_;
	iter->signal_.s_ 	/= c0_;
      }
    for (std::vector<ExtField>::iterator iter = extField_.begin(); iter != extField_.end(); iter++)
      {
	iter->c0_ 	         = c0_;
	iter->signal_.t0_     	/= c0_;
	iter->signal_.f0_ 	*= c0_;
	iter->signal_.s_ 	/= c0_;
      }

    /****************************************************************************************************/

    /* According to the set parameters for length and time scales correct the amplitudes of the seed
     * and undulator. Note that the seed amplitude is the amplitude of the vector potential, whereas the
     * undulator amplitude is the amplitude of the fields.						*/
    seed_.amplitude_ 	        = seed_.a0_ * EM * c0_ / EC;
    for (std::vector<Undulator>::iterator iter = undulator_.begin(); iter != undulator_.end(); iter++)
      iter->amplitude_ 		= iter->a0_ * EM * c0_ * 2 * PI * iter->signal_.f0_ / EC;
    for (std::vector<ExtField>::iterator iter = extField_.begin(); iter != extField_.end(); iter++)
      iter->amplitude_        	= iter->a0_ * EM * c0_ * 2 * PI * iter->signal_.f0_ / EC;

    /****************************************************************************************************/

    /* Calculate the average gamma of the input bunches.						*/
    Double gamma = 0.0;
    for (unsigned int i = 0; i < bunch_.bunchInit_.size(); i++)
      gamma += bunch_.bunchInit_[i].initialGamma_ / bunch_.bunchInit_.size();

    /* Now, depending on the undulator type determine the maximum and minimum gamma of the bunch
     * travelling through the undulator.								*/
    Double gmin = 1.0e100, gmax = -1.0e100, g = 0.0;
    for (std::vector<Undulator>::iterator iter = undulator_.begin(); iter != undulator_.end(); iter++)
      {
	if ( iter->type_ == STATIC )
	  {
	    g  	 = gamma / sqrt( 1.0 + iter->k_ * iter->k_ / 2.0 );
	    gmin = ( gmin < g ) ? gmin : g;
	    gmax = ( gmax > g ) ? gmax : g;
	  }
	else if ( iter->signal_.signalType_ == FLATTOP )
	  {
	    g  = gamma / sqrt( 1.0 + iter->a0_ * iter->a0_ / 2.0 );
	    gmin = ( gmin < g ) ? gmin : g;
	    gmax = ( gmax > g ) ? gmax : g;
	  }
	else
	  {
	    /* For optical undulator, when a pulse other than flattop is the pulse format, the gamma of
	     * the electrons change throughout the interaction.						*/
	    g  = gamma / sqrt( 1.0 + iter->a0_ * iter->a0_ / 2.0 );
	    gmin = ( gmin < g ) ? gmin : g;
	    g  = gamma;
	    gmax = ( gmax > g ) ? gmax : g;
	  }

      }

    /* Set the gamma of the moving frame equal to the average of the maximum and minimum gamma.		*/
    gamma_ = ( gmin + gmax ) / 2.0;

    /****************************************************************************************************/

    /* Boost and initialize the undulator related parameters.						*/
    beta_ = sqrt( 1.0 - 1.0 / ( gamma_ * gamma_ ) );
    for (std::vector<Undulator>::iterator iter = undulator_.begin(); iter != undulator_.end(); iter++)
      if ( iter->type_ == OPTICAL ) iter->lu_ /= ( 1 + beta_ );

    /* Sort the undulators according to their beginning point.						*/
    std::sort(undulator_.begin(), undulator_.end(), undulatorCompare);

    /* Now shift all the undulator modules so that the first module starts at zero.			*/
    for (std::vector<Undulator>::reverse_iterator iter = undulator_.rbegin(); iter != undulator_.rend(); iter++)
      iter->rb_ -= undulator_[0].rb_;

    /****************************************************************************************************/

    /* Set the modulation wavelength for each bunch using the given gamma of the bunch. 		*/
    for (unsigned int i = 0; i < bunch_.bunchInit_.size(); i++)
      {
	/* First determine the beta vector of the bunch.						*/
	bunch_.bunchInit_[i].initialBeta_	= sqrt( 1.0 - 1.0 / pow( bunch_.bunchInit_[i].initialGamma_ , 2 ) );
	bunch_.bunchInit_[i].betaVector_.mv( bunch_.bunchInit_[i].initialBeta_, bunch_.bunchInit_[i].initialDirection_);

	/* Calculate the modulation wavelength of the bunch for the initial bunching factor or the shot
	 * noise implementation.									*/
	Double lu				= ( undulator_.size() > 0 ) ? undulator_[0].lu_ : 0.0;
	bunch_.bunchInit_[i].lambda_		= lu / ( 2.0 * gamma_ * gamma_ ) * bunch_.bunchInit_[i].betaVector_[2] / beta_;

	printmessage(std::string(__FILE__), __LINE__, std::string("Modulation wavelength of the bunch outside the undulator is set to " + stringify( bunch_.bunchInit_[i].lambda_ ) ) );
      }

    printmessage(std::string(__FILE__), __LINE__, std::string("The given mesh parameters are boosted into the electron rest frame :::") );

    /****************************************************************************************************/

    /* The Lorentz boost parameters should also be transfered to the seed class in order to correctly
     * compute the fields within the computational domain.						*/
    seed_.beta_    	= beta_;
    seed_.gamma_   	= gamma_;
  }

  /******************************************************************************************************
   * Boost the mesh into the electron rest frame.
   ******************************************************************************************************/

  void Solver::lorentzBoostMesh ()
  {
    /****************************************************************************************************/

    /* Boost the mesh data into the electron rest frame.						*/
    mesh_.meshLength_[2] 	*= gamma_;
    mesh_.meshResolution_[2] 	*= gamma_;
    mesh_.meshCenter_[2] 	*= gamma_;
    mesh_.totalTime_		/= gamma_;

    /****************************************************************************************************/

    /* Set the mesh parameters according to the given simulation.					*/
    if ( mesh_.solver_ == NSFD )
      {
	/* Adjust the transverse mesh resolution to match the stability criterion.			*/
	Double t = 1.0 / sqrt( pow( mesh_.meshResolution_[2] / mesh_.meshResolution_[0], 2.0 ) + pow( mesh_.meshResolution_[2] / mesh_.meshResolution_[1], 2.0 ) );
	if ( t < 1.02 )
	  {
	    mesh_.meshResolution_[0] *= 1.02 / t;
	    mesh_.meshResolution_[1] *= 1.02 / t;

	    printmessage(std::string(__FILE__), __LINE__, std::string("Transverse discretization along x is set to " + stringify(mesh_.meshResolution_[0]) ) );
	    printmessage(std::string(__FILE__), __LINE__, std::string("Transverse discretization along y is set to " + stringify(mesh_.meshResolution_[1]) ) );
	  }

	/* Based on the dispersion condition, the field time step can be obtained.			*/
	mesh_.timeStep_		 = mesh_.meshResolution_[2] / c0_;
	printmessage(std::string(__FILE__), __LINE__, std::string("Time step for the field update is set to " + stringify(mesh_.timeStep_ * gamma_) ) );
      }
    else if ( mesh_.solver_ == FD )
      {
	/* Based on the dispersion condition, the field time step can be obtained.			*/
	mesh_.timeStep_		 = 1.0 / ( c0_ * sqrt( 1.0 / pow(mesh_.meshResolution_[0], 2.0) + 1.0 / pow(mesh_.meshResolution_[1], 2.0) + 1.0 / pow(mesh_.meshResolution_[2], 2.0) ) );
	printmessage(std::string(__FILE__), __LINE__, std::string("Time step for the field update is set to " + stringify(mesh_.timeStep_ * gamma_) ) );
      }
  }

  /******************************************************************************************************
   * Boost particles into the electron rest frame.
   ******************************************************************************************************/

  void Solver::lorentzBoostBunch ()
  {
    /****************************************************************************************************/

    /* Set the bunch update time step if it is given, otherwise set it according to the MITHRA rules.	*/
    bunch_.timeStep_			/= gamma_;

    /* Adjust the given bunch time step according to the given field time step.				*/
    bunch_.timeStep_ 			 = mesh_.timeStep_ / ceil(mesh_.timeStep_ / bunch_.timeStep_);
    nUpdateBunch_    			 = mesh_.timeStep_ / bunch_.timeStep_;
    printmessage(std::string(__FILE__), __LINE__, std::string("Time step for the bunch update is set to " + stringify(bunch_.timeStep_ * gamma_) ) );

    /****************************************************************************************************/

    /* Boost the bunch sampling parameters into the electron rest frame.				*/
    bunch_.rhythm_			/= gamma_;
    bunch_.bunchVTKRhythm_		/= gamma_;
    for (unsigned int i = 0; i < bunch_.bunchProfileTime_.size(); i++)
      bunch_.bunchProfileTime_[i] 	/= gamma_;
    bunch_.bunchProfileRhythm_		/= gamma_;

    /****************************************************************************************************/

    /* Boost the coordinates of the macro-particles into the bunch rest frame. During the boosting, the
     * maximum value of the z-coordinate is important, since it needs to be used in the shift that will
     * be introduced to the bunch.									*/

    Double zmaxL = -1.0e100, zmaxG;
    for (auto iterQ = chargeVectorn_.begin(); iterQ != chargeVectorn_.end(); iterQ++ )
      {
	Double g  	= std::sqrt(1.0 + iterQ->gbnp.norm2());
	Double bz 	= iterQ->gbnp[2] / g;
	iterQ->rnp[2]  *= gamma_;
	iterQ->gbnp[2] 	= gamma_ * g * ( bz - beta_ );

	zmaxL 		= std::max( zmaxL , iterQ->rnp[2] );
      }
    MPI_Allreduce(&zmaxL, &zmaxG, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    /****************************************************************************************************/

    /* Here, we define the shift in time such that the bunch end is at the begin of fringing field
     * section. For optical undulator such a separation should be considered in the input parameters
     * where offset is given.										*/

    /* This shift in time makes sure that the maximum z in the bunch at time zero is at undulator[0].dist_
     * away from the undulator begin ( the default distance is 2 undulator periods for static undulators
     * and 10 undulator periods for optical ones, due to different fring field formats used in the two
     * cases.)												*/
    if ( undulator_.size() > 0 )
      {
	Double nl = ( undulator_[0].type_ == STATIC ) ? 2.0 : 10.0;
	if (undulator_[0].dist_ == 0.0) undulator_[0].dist_ = nl * undulator_[0].lu_;
	else if (undulator_[0].dist_ < nl * undulator_[0].lu_)
	  printmessage(std::string(__FILE__), __LINE__, std::string("Warning: the undulator is set very close to the bunch, the results may be inaccurate.") );
	dt_ 		= - 1.0 / ( beta_ * undulator_[0].c0_ ) * ( zmaxG + undulator_[0].dist_ / gamma_ );
	printmessage(std::string(__FILE__), __LINE__, std::string("Initial distance from bunch head to undulator is ") + stringify(undulator_[0].dist_) );
      }

    /* The same shift in time should also be done for the seed field.					*/
    seed_.dt_		= dt_;

    /* With the above definition in the time begin the entrance of the undulator, i.e. z = 0 in the lab
     * frame, corresponds to the z = zmax + undulator_[0].dist_ / gamma_ in the bunch rest frame at
     * the initialization instant, i.e. timeBunch = 0.0.
     * In MITHRA, we assume that the particle move on a straight line before reaching the start point.
     * This forces the bunch properties to be as given at the start point of the undulator. The start
     * point is here the undulator entrance minus the undulator distance.				*/
    bunch_.zu_ 		= zmaxG;
    bunch_.beta_ 	= beta_;

    /****************************************************************************************************/

    for (auto iterQ = chargeVectorn_.begin(); iterQ != chargeVectorn_.end(); iterQ++ )
      {
	Double g	= std::sqrt(1.0 + iterQ->gbnp.norm2());
	iterQ->rnp[0]  += iterQ->gbnp[0] / g * ( iterQ->rnp[2] - bunch_.zu_ ) * beta_;
	iterQ->rnp[1]  += iterQ->gbnp[1] / g * ( iterQ->rnp[2] - bunch_.zu_ ) * beta_;
	iterQ->rnp[2]  += iterQ->gbnp[2] / g * ( iterQ->rnp[2] - bunch_.zu_ ) * beta_;
      }

    /* Distribute particles in their respective processor, depending on their longituinal coordinate.	*/
    distributeParticles(chargeVectorn_);

    /* Print the total number of macro-particles for the user.						*/
    unsigned int NqL = chargeVectorn_.size(), NqG = 0;
    MPI_Reduce(&NqL,&NqG,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    printmessage(std::string(__FILE__), __LINE__, std::string("The total number of macro-particles is equal to ") + stringify(NqG) + std::string(" .") );

    /* Initialize the total number of charges.								*/
    Nc_ = chargeVectorn_.size();
  }

  /******************************************************************************************************
   * Distribute particles in their respective processor, depending on their longituinal coordinate.
   ******************************************************************************************************/

  void Solver::distributeParticles (std::list<Charge>& chargeVector)
  {
    /* Note that this function only redistributes q, rnp, gbnp, but NOT rnm and gbnm.			*/
    std::vector<Double> sendCV;
    std::list<Charge>::iterator it = chargeVector.begin();
    while(it != chargeVector.end())
      {
	if ( ( it->rnp[2] < zp_[1] || rank_ == size_ - 1 ) && ( it->rnp[2] >= zp_[0] || rank_ == 0 ) )
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
    int *counts = new int [size_], *disps = new int [size_];
    MPI_Allgather( &sizeSend, 1, MPI_INT, counts, 1, MPI_INT, MPI_COMM_WORLD );
    for (int i = 0; i < size_; i++)
      disps[i] = (i > 0) ? (disps[i-1] + counts[i-1]) : 0;
    std::vector<Double> recvCV (disps[size_-1] + counts[size_-1], 0);
    MPI_Allgatherv(&sendCV[0], sizeSend, MPI_DOUBLE, &recvCV[0], counts, disps, MPI_DOUBLE, MPI_COMM_WORLD);

    /* And now replace all the charges in the charge vector of the corresponding processor.  		*/
    unsigned i = 0;
    Charge charge;
    while (i < recvCV.size() )
      {
	if ( (recvCV[i+3] < zp_[1] || rank_ == size_ - 1) && (recvCV[i+3] >= zp_[0] || rank_ == 0) )
	  {
	    charge.q 		= recvCV[i++];
	    charge.rnp[0] 	= recvCV[i++];
	    charge.rnp[1] 	= recvCV[i++];
	    charge.rnp[2] 	= recvCV[i++];
	    charge.gbnp[0] 	= recvCV[i++];
	    charge.gbnp[1] 	= recvCV[i++];
	    charge.gbnp[2] 	= recvCV[i++];
	    chargeVector.push_back(charge);
	  }
	else
	  i += 7;
      }
  }

  /******************************************************************************************************
   * Initialize the matrix for the field values and the coordinates.
   ******************************************************************************************************/

  void Solver::initialize ()
  {
    /* Using the parsed data, set the required parameters for simulation.				*/
    setSimulationParameters();

    /* As the first step, initialize the given bunch according to to the parameters given in the
     * datainput file.											*/
    initializeBunch();
    timeBunch_ = time_;

    /* Transfer the whole quantities to the electron rest frame.					*/
    lorentzBoostMesh();

    /* Initialize the spatial and temporal mesh of the problem.						*/
    initializeMesh();

    /* Transfer the bunch to the lab frame and adapt the bunch for the FEL simulation, i.e. add the shot
     * noise, distribute the mirrored macro-particles, and add the bunch tail.				*/
    lorentzBoostBunch();

    /* Initialize the update data for the field.							*/
    initializeField();

    /* If sampling is enabled initialize the required data for sampling the field and saving it.	*/
    if (seed_.sampling_)			initializeSeedSampling();

    /* Initialize the required data for visualizing and saving the field.				*/
    initializeSeedVTK();

    /* If profiling is enabled initialize the required data for profiling the field and saving it.	*/
    if (seed_.profile_)				initializeSeedProfile();

    /* Initialize the data needed for updating the bunches.						*/
    initializeBunchUpdate();

    /* If sampling or visualizing the radiation power is enabled initialize the required data for
     * calculating the radiation power and saving it.							*/
    initializePowerSample(); initializePowerVisualize();

    /* If sampling the radiation energy is enabled initialize the required data for calculating the
     * radiation energy and saving it.									*/
    initializeEnergySample();

    /* Initialize required data for saving particles hitting screens.					*/
    initializeScreenProfile();
  }

  /******************************************************************************************************
   * Initialize the temporal and spatial mesh of the problem.
   ******************************************************************************************************/

  void Solver::initializeMesh ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the temporal and spatial mesh of the problem.") );

    /* Declare the required variables in the calculations to avoid redundant data declaration.		*/
    unsigned int 	m = 0;

    /* First the mesh-length needs to be adjusted to contain a multiple number of the mesh-resolution in
     * each direction.											*/
    N0_ = (int) ( mesh_.meshLength_[0] / mesh_.meshResolution_[0] ) + 2;
    N1_ = (int) ( mesh_.meshLength_[1] / mesh_.meshResolution_[1] ) + 2;
    N2_ = (int) ( mesh_.meshLength_[2] / mesh_.meshResolution_[2] ) + 2;
    N1N0_ = N1_*N0_;
    mesh_.meshLength_[0] = ( N0_ - 1 ) * mesh_.meshResolution_[0];
    mesh_.meshLength_[1] = ( N1_ - 1 ) * mesh_.meshResolution_[1];
    mesh_.meshLength_[2] = ( N2_ - 1 ) * mesh_.meshResolution_[2];

    /* Evaluate the number of nodes in each processor.							*/
    if ( size_ > 1 )
      {
	if ( rank_ == 0 )
	  {
	    np_ = N2_ / size_ + 1;
	    k0_ = 0;
	  }
	else if ( rank_ == size_ - 1 )
	  {
	    np_ = N2_ - ( size_ - 1 ) * ( N2_ / size_ ) + 3;
	    k0_ = ( size_ - 1 ) * ( N2_ / size_ ) - 1;
	  }
	else
	  {
	    np_ = N2_ / size_ + 2;
	    k0_ = rank_ * ( N2_ / size_ ) - 1;
	  }
      }
    else
      {
	np_ = N2_;
	k0_ = 0;
      }

    /* Now according to the number of nodes, initialize each of the field and coordinate matrices.	*/
    FieldVector<Double> ZERO_VECTOR (0.0);
    anp1_ = new std::vector<FieldVector<Double> > (N1N0_*np_, ZERO_VECTOR);
    an_   = new std::vector<FieldVector<Double> > (N1N0_*np_, ZERO_VECTOR);
    anm1_ = new std::vector<FieldVector<Double> > (N1N0_*np_, ZERO_VECTOR);
    jn_  .resize(N1N0_*np_, ZERO_VECTOR);
    r_   .resize(N1N0_*np_, ZERO_VECTOR);
    en_  .resize(N1N0_*np_, ZERO_VECTOR);
    bn_  .resize(N1N0_*np_, ZERO_VECTOR);
    pic_ .resize(N1N0_*np_, false);

    if ( mesh_.spaceCharge_ )
      {
	fnp1_ = new std::vector<Double> (N1N0_*np_, 0.0 );
	fn_   = new std::vector<Double> (N1N0_*np_, 0.0 );
	fnm1_ = new std::vector<Double> (N1N0_*np_, 0.0 );
	rn_   .resize(N1N0_*np_, 0.0);
      }

    /* Set the borders of the computational mesh.                                                     	*/
    xmin_ = mesh_.meshCenter_[0] - mesh_.meshLength_[0] / 2.0;
    xmax_ = mesh_.meshCenter_[0] + mesh_.meshLength_[0] / 2.0;
    ymin_ = mesh_.meshCenter_[1] - mesh_.meshLength_[1] / 2.0;
    ymax_ = mesh_.meshCenter_[1] + mesh_.meshLength_[1] / 2.0;
    zmin_ = mesh_.meshCenter_[2] - mesh_.meshLength_[2] / 2.0;
    zmax_ = mesh_.meshCenter_[2] + mesh_.meshLength_[2] / 2.0;

    /* Now set up the coordinate of nodes according to the number of nodes and mesh-length.		*/
    for (int i = 0; i < N0_; i++)
      for (int j = 0; j < N1_; j++)
	for (int k = 0; k < np_; k++)
	  {
	    m = N1N0_*k+N1_*i+j;
	    r_[m][0] = xmin_ + i         * mesh_.meshResolution_[0];
	    r_[m][1] = ymin_ + j         * mesh_.meshResolution_[1];
	    r_[m][2] = zmin_ + (k + k0_) * mesh_.meshResolution_[2];

	    if 		( k == 0 ) 		zp_[0] = r_[m][2];
	    else if 	( k == np_ - 2 ) 	zp_[1] = r_[m][2];
	  }

    /* Initialize the time points for the fields in the domain.						*/
    timep1_ += mesh_.timeStep_;
    timem1_ -= mesh_.timeStep_;

    printmessage(std::string(__FILE__), __LINE__, std::string("The temporal and spatial mesh of the problem is initialized. :::") );
  }

  /******************************************************************************************************
   * Initialize the data for updating the field in the fdtd algorithm.
   ******************************************************************************************************/

  void Solver::initializeField ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string(" ::: Initializing the field update data") );

    /* Declare the required variables in the calculations to avoid redundant data declaration.		*/
    unsigned int      m = 0;

    /* Initialize the data for updating the currents.                                                 	*/
    uc_.dx = mesh_.meshResolution_[0];
    uc_.dy = mesh_.meshResolution_[1];
    uc_.dz = mesh_.meshResolution_[2];
    uc_.dv = - m0_ * EC / mesh_.timeStep_ /   ( uc_.dx * uc_.dy * uc_.dz );
    uc_.rc = - EC / e0_ / 2.0 /               ( uc_.dx * uc_.dy * uc_.dz );
    FieldVector<Double> ZERO_VECTOR (0.0);
    uc_.jt.resize(N1N0_,ZERO_VECTOR);
    if ( mesh_.spaceCharge_ ) uc_.rt.resize(N1N0_,0.0);

    /* Now calculate all the coefficients needed for updating the fields with time.			*/

    uf_.dt	 = mesh_.timeStep_;

    uf_.dx	 = mesh_.meshResolution_[0];
    uf_.dy	 = mesh_.meshResolution_[1];
    uf_.dz	 = mesh_.meshResolution_[2];

    uf_.dx2    = 2.0 * uf_.dx;
    uf_.dy2    = 2.0 * uf_.dy;
    uf_.dz2    = 2.0 * uf_.dz;

    uf_.N0m1   = N0_ - 1;
    uf_.N1m1   = N1_ - 1;
    uf_.npm1   = np_ - 1;

    /* Coefficients of the Nonstandard finite difference method.					*/
    Double beta   = ( 1.0 + 0.02 / ( pow(uf_.dz/uf_.dx,2.0) + pow(uf_.dz/uf_.dy,2.0) ) ) / 4.0;
    Double alpha  = 1.0 - 2.0 * beta;
    uf_.af.alpha_ = alpha;
    uf_.af.beta_  = beta / alpha;

    uf_.a[0]	 = 2.0 * ( 1.0 - alpha * pow(c0_*uf_.dt/uf_.dx,2) - alpha * pow(c0_*uf_.dt/uf_.dy,2) - pow(c0_*uf_.dt/uf_.dz,2) );
    uf_.a[1]	 = pow(c0_*uf_.dt/uf_.dx,2.0);
    uf_.a[2]	 = pow(c0_*uf_.dt/uf_.dy,2.0);
    uf_.a[3]	 = pow(c0_*uf_.dt/uf_.dz,2.0) - 2.0 * ( beta * pow(c0_*uf_.dt/uf_.dx,2) + beta * pow(c0_*uf_.dt/uf_.dy,2) );
    uf_.a[4]	 = pow(c0_*uf_.dt,2.0) * uc_.dv;
    uf_.a[5]   = pow(c0_*uf_.dt,2.0) * uc_.rc;

    uf_.af.ufa_ = &uf_.a[0];

    /* x = 0 and x = h boundary coefficients.								*/

    Double alpha1 = 0.0;
    Double alpha2 = 0.0;
    Double p      = ( 1.0 + cos(alpha1) * cos(alpha2) ) / ( cos(alpha1) + cos(alpha2) );
    Double q      = - 1.0 / ( cos(alpha1) + cos(alpha2) );

    Double d   = 1.0 / ( 2.0 * uf_.dt * uf_.dx ) + p / ( 2.0 * c0_ * uf_.dt * uf_.dt );

    uf_.bB[0]	 = (   1.0 / ( 2.0 * uf_.dt * uf_.dx ) - p / ( 2.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.bB[1]	 = ( - 1.0 / ( 2.0 * uf_.dt * uf_.dx ) - p / ( 2.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.bB[2]  = (   p   / ( c0_ * uf_.dt * uf_.dt ) + q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( uf_.dy * uf_.dy ) + c0_ / ( uf_.dz * uf_.dz ) ) ) / d;
    uf_.bB[3]  = - q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( 2.0 * uf_.dy * uf_.dy ) ) / d ;
    uf_.bB[4]  = - q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( 2.0 * uf_.dz * uf_.dz ) ) / d ;

    /* y = 0 and y = h boundary coefficients.								*/

    d  	 = 1.0 / ( 2.0 * uf_.dt * uf_.dy ) + p / ( 2.0 * c0_ * uf_.dt * uf_.dt );

    uf_.cB[0]	 = (   1.0 / ( 2.0 * uf_.dt * uf_.dy ) - p / ( 2.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.cB[1]	 = ( - 1.0 / ( 2.0 * uf_.dt * uf_.dy ) - p / ( 2.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.cB[2]  = (   p / ( c0_ * uf_.dt * uf_.dt ) + q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( uf_.dx * uf_.dx ) + c0_ / ( uf_.dz * uf_.dz ) ) ) / d;
    uf_.cB[3]  = - q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( 2.0 * uf_.dx * uf_.dx ) ) / d ;
    uf_.cB[4]  = - q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( 2.0 * uf_.dz * uf_.dz ) ) / d ;

    /* z = 0 and z = h boundary coefficients.								*/

    alpha1 = 0.0;
    alpha2 = 0.0;
    p      = ( 1.0 + cos(alpha1) * cos(alpha2) ) / ( cos(alpha1) + cos(alpha2) );
    q      = - 1.0 / ( cos(alpha1) + cos(alpha2) );

    d  	 = 1.0 / ( 2.0 * uf_.dt * uf_.dz ) + p / ( 2.0 * c0_ * uf_.dt * uf_.dt );

    uf_.dB[0]	 = (   1.0 / ( 2.0 * uf_.dt * uf_.dz ) - p / ( 2.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.dB[1]	 = ( - 1.0 / ( 2.0 * uf_.dt * uf_.dz ) - p / ( 2.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.dB[2]  = (   p / ( c0_ * uf_.dt * uf_.dt ) + q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( uf_.dx * uf_.dx ) + c0_ / ( uf_.dy * uf_.dy ) ) ) / d;
    uf_.dB[3]  = - q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( 2.0 * uf_.dx * uf_.dx ) ) / d ;
    uf_.dB[4]  = - q * ( mesh_.truncationOrder_ - 1.0 ) * ( c0_ / ( 2.0 * uf_.dy * uf_.dy ) ) / d ;

    /* Edge coefficients for edges along z.								*/
    d  	= ( 1.0 / uf_.dy + 1.0 / uf_.dx ) / ( 4.0 * uf_.dt ) + 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt );
    uf_.eE[0] = ( - ( 1.0 / uf_.dy - 1.0 / uf_.dx ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.eE[1] = (   ( 1.0 / uf_.dy - 1.0 / uf_.dx ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.eE[2] = (   ( 1.0 / uf_.dy + 1.0 / uf_.dx ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.eE[3] = ( 3.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt ) - c0_ / ( 4.0 * uf_.dz * uf_.dz ) ) / d;
    uf_.eE[4] = c0_ / ( 8.0 * uf_.dz * uf_.dz ) / d;

    /* Edge coefficients for edges along x.								*/
    d  	= ( 1.0 / uf_.dz + 1.0 / uf_.dy ) / ( 4.0 * uf_.dt ) + 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt );
    uf_.fE[0] = ( - ( 1.0 / uf_.dz - 1.0 / uf_.dy ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.fE[1] = (   ( 1.0 / uf_.dz - 1.0 / uf_.dy ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.fE[2] = (   ( 1.0 / uf_.dz + 1.0 / uf_.dy ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.fE[3] = ( 3.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt ) - c0_ / ( 4.0 * uf_.dx * uf_.dx ) ) / d;
    uf_.fE[4] = c0_ / ( 8.0 * uf_.dx * uf_.dx ) / d;

    /* Edge coefficients for edges along y.								*/
    d  	= ( 1.0 / uf_.dx + 1.0 / uf_.dz ) / ( 4.0 * uf_.dt ) + 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt );
    uf_.gE[0] = ( - ( 1.0 / uf_.dx - 1.0 / uf_.dz ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.gE[1] = (   ( 1.0 / uf_.dx - 1.0 / uf_.dz ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.gE[2] = (   ( 1.0 / uf_.dx + 1.0 / uf_.dz ) / ( 4.0 * uf_.dt ) - 3.0 / ( 8.0 * c0_ * uf_.dt * uf_.dt ) ) / d;
    uf_.gE[3] = ( 3.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt ) - c0_ / ( 4.0 * uf_.dy * uf_.dy ) ) / d;
    uf_.gE[4] = c0_ / ( 8.0 * uf_.dy * uf_.dy ) / d;

    /* Corner coefficients.										*/
    uf_.hC[0]  =   ( - 1.0 / uf_.dx - 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[1]  =   (   1.0 / uf_.dx - 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[2]  =   ( - 1.0 / uf_.dx + 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[3]  =   ( - 1.0 / uf_.dx - 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[4]  =   (   1.0 / uf_.dx + 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[5]  =   (   1.0 / uf_.dx - 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[6]  =   ( - 1.0 / uf_.dx + 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[7]  =   (   1.0 / uf_.dx + 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[8]  = - ( - 1.0 / uf_.dx - 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[9]  = - (   1.0 / uf_.dx - 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[10] = - ( - 1.0 / uf_.dx + 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[11] = - ( - 1.0 / uf_.dx - 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[12] = - (   1.0 / uf_.dx + 1.0 / uf_.dy - 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[13] = - (   1.0 / uf_.dx - 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[14] = - ( - 1.0 / uf_.dx + 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[15] = - (   1.0 / uf_.dx + 1.0 / uf_.dy + 1.0 / uf_.dz ) / ( 8.0 * uf_.dt ) - 1.0 / ( 4.0 * c0_ * uf_.dt * uf_.dt );
    uf_.hC[16] = 1.0 / ( 2.0 * c0_ * uf_.dt * uf_.dt );

    /* The seed field should be initialized in the computational mesh as well. This should be done for
     * the mesh points within the TF domain.								*/
    if (seed_.amplitude_ > 1.0e-50)
      {
	for (int i = 2; i < N0_ - 2; i++)
	  for (int j = 2; j < N1_ - 2; j++)
	    for (int k = 2 * abs(signof(rank_) - 1 ) ; k < np_ - 2 * abs(signof(rank_ - size_ + 1) + 1 ); k++)
	      {
		m = N1N0_ * k + N1_ * i + j;
		seed_.fields(r_[m], time_,   (*an_)[m]   );
		seed_.fields(r_[m], timem1_, (*anm1_)[m] );
	      }
      }

    printmessage(std::string(__FILE__), __LINE__, std::string(" The field update data are initialized. :::") );
  }

  /******************************************************************************************************
   * Initialize the data required for sampling seed or the total field in the computational domain.
   ******************************************************************************************************/

  void Solver::initializeSeedSampling ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string(" ::: Initializing the field sampling data") );

    /* Return an error if the seed sampling rhythm is still zero.					*/
    if ( seed_.samplingRhythm_ == 0 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The sampling rhythm of the field is zero although sampling is activated !!!") );
	exit(1);
      }

    /* Lorentz boost the seed sampling rhythm to the electron rest frame.				*/
    seed_.samplingRhythm_		/= gamma_;

    /* Perform the lorentz boost for the sampling data.							*/
    for (unsigned i = 0; i < seed_.samplingPosition_.size(); i++)
      seed_.samplingPosition_[i][2] 	*= gamma_;
    seed_.samplingLineBegin_[2]	*= gamma_;
    seed_.samplingLineEnd_  [2]	*= gamma_;

    /* If sampling type is plot over line initialize the positions according to the line begin and line
     * end.                                                                              		*/
    if ( seed_.samplingType_ == OVERLINE )
      {
	FieldVector<Double> l = seed_.samplingLineEnd_;
	l       -= seed_.samplingLineBegin_;
	Double n = seed_.samplingRes_ ;
	l       /= n;

	FieldVector<Double> position, r;

	for (unsigned i = 0; i < n; i++)
	  {
	    position[0] = seed_.samplingLineBegin_[0] + i * l[0];
	    position[1] = seed_.samplingLineBegin_[1] + i * l[1];
	    position[2] = seed_.samplingLineBegin_[2] + i * l[2];
	    seed_.samplingPosition_.push_back(position);
	  }

      }

    /* Remove the point from the sampling locations if the point does not correspond to the corresponding
     * processor.											*/
    std::vector<FieldVector<Double> >	samplingPosition; samplingPosition.clear();
    for (unsigned int n = 0; n < seed_.samplingPosition_.size(); ++n)
      {
	/* Get the position of the sampling point.							*/
	sf_.position 	= seed_.samplingPosition_[n];

	/* Check if the sampling point resides in the computational mesh.				*/
	if ( sf_.position[0] < xmax_ - ub_.dx && sf_.position[0] > xmin_ + ub_.dx &&
	    sf_.position[1] < ymax_ - ub_.dy && sf_.position[1] > ymin_ + ub_.dy &&
	    sf_.position[2] < zmax_ - ub_.dz && sf_.position[2] > zmin_ + ub_.dz )
	  {
	    if ( sf_.position[2] < zp_[1] && sf_.position[2] >= zp_[0] )
	      samplingPosition.push_back(seed_.samplingPosition_[n]);
	  }
	else
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The sampling point does not reside in the grid. No data is saved.") );
	  }
      }
    seed_.samplingPosition_ = samplingPosition;
    sf_.N  = seed_.samplingPosition_.size();

    /* Create the file stream and directories for saving the seed sampling data.			*/
    if ( sf_.N > 0 )
      {
	std::string baseFilename = "";
	if (!(isabsolute(seed_.samplingBasename_))) baseFilename = seed_.samplingDirectory_;
	baseFilename += seed_.samplingBasename_ + "-" + stringify(rank_) + TXT_FILE_SUFFIX;

	createDirectory(baseFilename, 0);
	sf_.file = new std::ofstream(baseFilename.c_str(),std::ios::trunc);
      }

    /* Set the constants to change the units in the output file.					*/
    sf_.Ce = mesh_.lengthScale_ / pow( mesh_.timeScale_ , 2 );
    sf_.Cb = 1.0 / ( mesh_.lengthScale_ * mesh_.timeScale_  );
    sf_.Ca = 1.0 / mesh_.timeScale_;
    sf_.Cj = 1.0 / ( pow( mesh_.lengthScale_ , 2 ) * mesh_.timeScale_ );
    sf_.Cf = pow( mesh_.lengthScale_ / mesh_.timeScale_ , 2 );

    printmessage(std::string(__FILE__), __LINE__, std::string(" The field sampling data are initialized. :::") );
  }

  /******************************************************************************************************
   * Initialize the required data for visualizing and saving the field.
   ******************************************************************************************************/

  void Solver::initializeSeedVTK ()
  {
    /* Resize the visualization field database according to the number of defined visualizations.	*/
    vf_.resize( seed_.vtk_.size() );

    /* Perform the lorentz boost for the visualization data.						*/
    for (unsigned int i = 0; i < seed_.vtk_.size(); i++)
      {
	/* Continue the loop if the sampling option for this vtk file is not enabled.			*/
	if ( !seed_.vtk_[i].sample_ ) continue;

	/* Return an error if the seed visualization rhythm is still zero.				*/
	if ( seed_.vtk_[i].rhythm_ == 0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The visualization rhythm of the field is zero although visualization is activated !!!") );
	    exit(1);
	  }

	/* Lorentz boost the seed visualization rhythm to the electron rest frame.			*/
	seed_.vtk_[i].rhythm_			/= gamma_;

	/* Lorentz boost the plane-position to the electron rest frame.				*/
	seed_.vtk_[i].position_[2]		*= gamma_;

	if (!(isabsolute(seed_.vtk_[i].basename_)))
	  seed_.vtk_[i].basename_ = seed_.vtk_[i].directory_ + seed_.vtk_[i].basename_;
	splitFilename(seed_.vtk_[i].basename_, vf_[i].path, vf_[i].name);

	/* If the directory of the baseFilename does not exist create this directory.			*/
	createDirectory(seed_.vtk_[i].basename_, rank_);

	/* Initialize the vector to save the data and write to the vtk file.				*/
	std::vector<Double> ZERO_VECTOR ( (seed_.vtk_[i].field_).size(), 0.0);
	if 	  ( seed_.vtk_[i].type_ == ALLDOMAIN )
	  vf_[i].v.resize( N1N0_*np_, ZERO_VECTOR);
	else if ( seed_.vtk_[i].type_ == INPLANE )
	  {
	    if      ( seed_.vtk_[i].plane_ == XNORMAL )
	      {
		/* Check if the given position for the plane resides in the computational domain.	*/
		if ( seed_.vtk_[i].position_[0] > xmax_ - ub_.dx || seed_.vtk_[i].position_[0] < xmin_ + ub_.dx )
		  {
		    printmessage(std::string(__FILE__), __LINE__, std::string("The plane does not reside in the grid. No data is saved.") );
		    seed_.vtk_.erase( seed_.vtk_.begin() + i );
		    vf_.erase( vf_.begin() + i );
		  }
		else
		  vf_[i].v.resize( N1_*np_, ZERO_VECTOR);
	      }
	    else if ( seed_.vtk_[i].plane_ == YNORMAL )
	      {
		/* Check if the given position for the plane resides in the computational domain.	*/
		if ( seed_.vtk_[i].position_[1] > ymax_ - ub_.dy || seed_.vtk_[i].position_[1] < ymin_ + ub_.dy )
		  {
		    printmessage(std::string(__FILE__), __LINE__, std::string("The plane does not reside in the grid. No data is saved.") );
		    seed_.vtk_.erase( seed_.vtk_.begin() + i );
		    vf_.erase( vf_.begin() + i );
		  }
		else
		  vf_[i].v.resize( N0_*np_, ZERO_VECTOR);
	      }
	    else if ( seed_.vtk_[i].plane_ == ZNORMAL )
	      {
		/* Check if the given position for the plane resides in the computational domain.	*/
		if ( seed_.vtk_[i].position_[2] > zmax_ - ub_.dz || seed_.vtk_[i].position_[2] < zmin_ + ub_.dz )
		  {
		    printmessage(std::string(__FILE__), __LINE__, std::string("The plane does not reside in the grid. No data is saved.") );
		    seed_.vtk_.erase( seed_.vtk_.begin() + i );
		    vf_.erase( vf_.begin() + i );
		  }
		else
		  {
		    if ( seed_.vtk_[i].position_[2] < zp_[1] && seed_.vtk_[i].position_[2] >= zp_[0] )
		      vf_[i].v.resize( N1N0_,   ZERO_VECTOR);
		  }
	      }
	  }
      }
  }

  /******************************************************************************************************
   * Initialize the required data for profiling and saving the field.
   ******************************************************************************************************/

  void Solver::initializeSeedProfile ()
  {
    /* Return an error if the bunch profiling rhythm is still zero and no time is set.			*/
    if ( seed_.profileRhythm_ == 0 && seed_.profileTime_.size() == 0 )
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("The profiling rhythm of the field is zero and no time is set although profiling of the field is activated !!!") );
	exit(1);
      }

    /* Lorentz boost the seed profiling rhythm to the electron rest frame.				*/
    seed_.profileRhythm_	/= gamma_;

    /* Perform the lorentz boost for the profiling data.						*/
    for (unsigned i = 0; i < seed_.profileTime_.size(); i++)
      seed_.profileTime_[i] 	/= gamma_;

    if (!(isabsolute(seed_.profileBasename_))) seed_.profileBasename_ = seed_.profileDirectory_ + seed_.profileBasename_;

    /* If the directory of the baseFilename does not exist create this directory.			*/
    createDirectory(seed_.profileBasename_, rank_);

    pf_.dt = 2.0 * mesh_.timeStep_;
  }

  /******************************************************************************************************
   * Initialize the required data for updating the bunch.
   ******************************************************************************************************/

  void Solver::initializeBunchUpdate ()
  {
    /* Calculate the absolute values for the time step and the mesh resolution.				*/
    ub_.dt	= mesh_.timeStep_;
    ub_.dtb	= c0_ * bunch_.timeStep_;
    ub_.dx	= mesh_.meshResolution_[0];
    ub_.dy	= mesh_.meshResolution_[1];
    ub_.dz	= mesh_.meshResolution_[2];
    ub_.r1    = - EC / ( EM * c0_ ) * bunch_.timeStep_ / 2.0;
    ub_.r2    = - EC / EM * bunch_.timeStep_ / 2.0;

    /* If bunch sampling is enabled, initialize the required data for sampling and saving the bunch.	*/
    if (bunch_.sampling_)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the bunch sampling data.") );

	std::string baseFilename = "";
	if (!(isabsolute(bunch_.basename_))) baseFilename = bunch_.directory_;
	baseFilename += bunch_.basename_ + TXT_FILE_SUFFIX;

	/* If the directory of the baseFilename does not exist create this directory.			*/
	createDirectory(baseFilename, rank_);

	sb_.file = new std::ofstream(baseFilename.c_str(),std::ios::trunc);

	/* Return an error if the bunch sampling rhythm is still zero.					*/
	if ( bunch_.rhythm_ == 0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The sampling rhythm of the bunch is zero although sampling is activated !!!") );
	    exit(1);
	  }

	printmessage(std::string(__FILE__), __LINE__, std::string(" The sampling data are initialized. :::") );
      }

    /* If bunch visualization is enabled, initialize the required data for visualizing the bunch.	*/
    if (bunch_.bunchVTK_)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the bunch visualization data.") );

	if (!(isabsolute(bunch_.bunchVTKBasename_))) bunch_.bunchVTKBasename_ = bunch_.bunchVTKDirectory_ + bunch_.bunchVTKBasename_;

	/* If the directory of the baseFilename does not exist create this directory.			*/
	createDirectory(bunch_.bunchVTKBasename_, rank_);

	/* Return an error if the bunch visualization rhythm is still zero.				*/
	if ( bunch_.bunchVTKRhythm_ == 0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The visualization rhythm of the bunch is zero although visualization is activated !!!") );
	    exit(1);
	  }

	printmessage(std::string(__FILE__), __LINE__, std::string(" The bunch visualization data are initialized. :::") );
      }

    /* If writing the bunch profile is enabled, initialize the required data for profiling the bunch.	*/
    if (bunch_.bunchProfile_)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the bunch profiling data.") );

	if (!(isabsolute(bunch_.bunchProfileBasename_))) bunch_.bunchProfileBasename_ = bunch_.bunchProfileDirectory_ + bunch_.bunchProfileBasename_;

	/* If the directory of the baseFilename does not exist create this directory.			*/
	createDirectory(bunch_.bunchProfileBasename_, rank_);

	/* Return an error if the bunch profiling rhythm is still zero and no time is set.		*/
	if ( bunch_.bunchProfileRhythm_ == 0 && bunch_.bunchProfileTime_.size() == 0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("The profiling rhythm of the bunch is zero and no time is set although profiling of the bunch is activated !!!") );
	    exit(1);
	  }

	printmessage(std::string(__FILE__), __LINE__, std::string(" The bunch profiling data are initialized. :::") );
      }
  }

  /******************************************************************************************************
   * Initialize the charge vector containing the bunch given by the user.
   ******************************************************************************************************/

  void Solver::initializeBunch ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string("[[[ Initializing the bunch and prepare the charge vector ") );

    /* Clear the global charge vector and define a temporary one.					*/
    chargeVectorn_.clear();
    std::list<Charge> qv;

    for (unsigned int i = 0; i < bunch_.bunchInit_.size(); i++)
      {
	/* Clear the temprary charge vector.								*/
	qv.clear();

	/* Initialize the bunch in the code.								*/
	if 	  ( bunch_.bunchInit_[i].bunchType_ == "manual" )
	  {
	    for ( unsigned int ia = 0; ia < bunch_.bunchInit_[i].position_.size(); ia++)
	      bunch_.initializeManual(	bunch_.bunchInit_[i], qv, zp_, rank_, size_, ia);
	  }
	else if ( bunch_.bunchInit_[i].bunchType_ == "ellipsoid" )
	  {
	    for ( unsigned int ia = 0; ia < bunch_.bunchInit_[i].position_.size(); ia++)
	      bunch_.initializeEllipsoid( bunch_.bunchInit_[i], qv, rank_, size_, ia);
	  }
	else if ( bunch_.bunchInit_[i].bunchType_ == "3D-crystal" )
	  {
	    for ( unsigned int ia = 0; ia < bunch_.bunchInit_[i].position_.size(); ia++)
	      bunch_.initialize3DCrystal(	bunch_.bunchInit_[i], qv, zp_, rank_, size_, ia);
	  }
	else if ( bunch_.bunchInit_[i].bunchType_ == "file" )
	  {
	    for ( unsigned int ia = 0; ia < bunch_.bunchInit_[i].position_.size(); ia++)
	      bunch_.initializeFile(		bunch_.bunchInit_[i], qv, zp_, rank_, size_, ia);
	  }

	/* Add the bunch distribution to the global charge vector.					*/
	chargeVectorn_.splice(chargeVectorn_.end(),qv);
      }

    printmessage(std::string(__FILE__), __LINE__, std::string("The bunch is initialized and the charge vector is prepared. ]]]") );
  }

  /******************************************************************************************************
   * Update the fields for one time-step
   ******************************************************************************************************/

  void Solver::bunchUpdate ()
  {
    /* First define a parameter for the processor number.						*/
    UpdateBunchParallel	        ubp;
    MPI_Status                  status;
    int                         msgtag1 = 1, msgtag2 = 2;

    /* Set the count for particles outside the transverse domain of the bunch to zero.			*/
    ubp.nt = 0;

    /* Loop over the charge points in the bunch, extract the real field values of the seed at their point,
     * superpose with the undulator field and eventually accelerate the particles within the field.	*/
    auto iter = chargeVectorn_.begin();
    while ( iter != chargeVectorn_.end() )
      {
	/* If the particle does not belong to this processor continue the loop over particles         	*/
	if ( !( ( iter->rnp[2] < zp_[1] || rank_ == size_ - 1 ) && ( iter->rnp[2] >= zp_[0] || rank_ == 0 ) ) )
	  {
	    ++iter;
	    continue;
	  }

	/* By default the charge stays in the computational domain after update.			*/
	ubp.dq = false;

	/* Get the boolean flag determining if the particle resides in the computational domain.      	*/
	ubp.b1x = ( iter->rnp[0] < xmax_ - ub_.dx && iter->rnp[0] > xmin_ + ub_.dx );
	ubp.b1y = ( iter->rnp[1] < ymax_ - ub_.dy && iter->rnp[1] > ymin_ + ub_.dy );
	ubp.b1z = ( iter->rnp[2] < zp_[1]         && iter->rnp[2] >= zp_[0] );

	/* Initialize the undulator fields.								*/
	ubp.bt = 0.0;
	ubp.et = 0.0;

	/* Compute the fields if the particle has passed the entrance zone of the undulator.		*/
	if ( iter->e == 1.0 )
	  {
	    /* Calculate the undulator field at the particle position.				        */
	    undulatorField(ubp, iter->rnp);

	    /* Calculate the external field at the particle position and add to the undulator field.	*/
	    externalField(ubp, iter->rnp);

	    if ( ubp.b1x && ubp.b1y && ubp.b1z )
	      {
		ubp.dxr = modf( ( iter->rnp[0] - xmin_ ) / ub_.dx , &ubp.d1 ); ubp.i = (int) ubp.d1;
		ubp.dyr = modf( ( iter->rnp[1] - ymin_ ) / ub_.dy , &ubp.d1 ); ubp.j = (int) ubp.d1;
		ubp.dzr = modf( ( iter->rnp[2] - zmin_ ) / ub_.dz , &ubp.d1 ); ubp.k = (int) ubp.d1;

		/* Get the index of the cell.								*/
		ubp.m   = ( ubp.k - k0_) * N1N0_ + ubp.i * N1_ + ubp.j;

		if (!pic_[ubp.m            ])     fieldEvaluate(ubp.m            );
		if (!pic_[ubp.m+N1_        ])     fieldEvaluate(ubp.m+N1_        );
		if (!pic_[ubp.m+1          ])     fieldEvaluate(ubp.m+1          );
		if (!pic_[ubp.m+N1_+1      ])     fieldEvaluate(ubp.m+N1_+1      );
		if (!pic_[ubp.m+N1N0_      ])     fieldEvaluate(ubp.m+N1N0_      );
		if (!pic_[ubp.m+N1N0_+N1_  ])     fieldEvaluate(ubp.m+N1N0_+N1_  );
		if (!pic_[ubp.m+N1N0_+1    ])     fieldEvaluate(ubp.m+N1N0_+1    );
		if (!pic_[ubp.m+N1N0_+N1_+1])     fieldEvaluate(ubp.m+N1N0_+N1_+1);

		/* Calculate and interpolate the electric field to find the value at the bunch point.	*/
		ubp.et.pmv( ( 1.0 - ubp.dxr ) * ( 1.0 - ubp.dyr ) * ( 1.0 - ubp.dzr) , en_[ubp.m		]);
		ubp.et.pmv(         ubp.dxr   * ( 1.0 - ubp.dyr ) * ( 1.0 - ubp.dzr) , en_[ubp.m+N1_	  	]);
		ubp.et.pmv( ( 1.0 - ubp.dxr ) *         ubp.dyr   * ( 1.0 - ubp.dzr) , en_[ubp.m+1	  	]);
		ubp.et.pmv(         ubp.dxr   *         ubp.dyr   * ( 1.0 - ubp.dzr) , en_[ubp.m+N1_+1	  	]);
		ubp.et.pmv( ( 1.0 - ubp.dxr ) * ( 1.0 - ubp.dyr ) *         ubp.dzr  , en_[ubp.m+N1N0_	  	]);
		ubp.et.pmv(         ubp.dxr   * ( 1.0 - ubp.dyr ) *         ubp.dzr  , en_[ubp.m+N1N0_+N1_  	]);
		ubp.et.pmv( ( 1.0 - ubp.dxr ) *         ubp.dyr   *         ubp.dzr  , en_[ubp.m+N1N0_+1	]);
		ubp.et.pmv(         ubp.dxr   *         ubp.dyr   *         ubp.dzr  , en_[ubp.m+N1N0_+N1_+1	]);

		/* Calculate and interpolate the magnetic field to find the value at the bunch point.	*/
		ubp.bt.pmv( ( 1.0 - ubp.dxr ) * ( 1.0 - ubp.dyr ) * ( 1.0 - ubp.dzr) , bn_[ubp.m		]);
		ubp.bt.pmv(         ubp.dxr   * ( 1.0 - ubp.dyr ) * ( 1.0 - ubp.dzr) , bn_[ubp.m+N1_		]);
		ubp.bt.pmv( ( 1.0 - ubp.dxr ) *         ubp.dyr   * ( 1.0 - ubp.dzr) , bn_[ubp.m+1		]);
		ubp.bt.pmv(         ubp.dxr   *         ubp.dyr   * ( 1.0 - ubp.dzr) , bn_[ubp.m+N1_+1		]);
		ubp.bt.pmv( ( 1.0 - ubp.dxr ) * ( 1.0 - ubp.dyr ) *         ubp.dzr  , bn_[ubp.m+N1N0_		]);
		ubp.bt.pmv(         ubp.dxr   * ( 1.0 - ubp.dyr ) *         ubp.dzr  , bn_[ubp.m+N1N0_+N1_	]);
		ubp.bt.pmv( ( 1.0 - ubp.dxr ) *         ubp.dyr   *         ubp.dzr  , bn_[ubp.m+N1N0_+1	]);
		ubp.bt.pmv(         ubp.dxr   *         ubp.dyr   *         ubp.dzr  , bn_[ubp.m+N1N0_+N1_+1	]);
	      }
	    else if ( !(ubp.b1x) && !(ubp.b1y) && ubp.b1z )
	      {
		/* Add one to the number of particles that do not reside in the transverse size of the
		 * computational domain.								*/
		ubp.nt++;
	      }
	  }
	else if ( undulator_.size() > 0 )
	  {
	    /* Update the emission vector flag based on the particle position in lab frame.		*/
	    ubp.lz 	= gamma_ * ( iter->rnp[2] + beta_ * c0_ * ( timeBunch_ + dt_ ) );
	    iter->e 	= ( ubp.lz > - undulator_[0].dist_ ) ? 1.0 : 0.0;
	  }
	else
	  iter->e	= 1.0;

	/* Update the velocity of the particle according to the calculated electric and magnetic field.	*/

	/* First, calculate the value of (gamma*beta)-.							*/
	ubp.gbm = iter->gbnp;
	ubp.gbm.pmv( ub_.r1 , ubp.et );

	/* Second, calculate the value of (gamma*beta)'.				                */
	ubp.gbp   = cross( ubp.gbm , ubp.bt );
	ubp.d1    = sqrt( 1.0 + ubp.gbm.norm2() );
	ubp.gbp.mv( ub_.r2 / ubp.d1, ubp.gbp);
	ubp.gbp  += ubp.gbm;

	/* Third, calculate the (gamma*beta)+.								*/
	ubp.gbpl  = cross( ubp.gbp , ubp.bt );
	ubp.gbpl.mv( 2.0 / ( ubp.d1 / ub_.r2 + ub_.r2 / ubp.d1 * ubp.bt.norm2() ), ubp.gbpl);
	ubp.gbpl += ubp.gbm;

	/* Fourth, update the (gamma*beta) vector.						        */
	iter->gbnp = ubp.gbpl;
	iter->gbnp.pmv( ub_.r1 , ubp.et );

	/* Determine the final position of the particle.				                */
	iter->rnp.pmv( ub_.dtb / sqrt (1.0 + iter->gbnp.norm2()) , iter->gbnp );

	/* If the particle enters the adjacent computational domain, save it to communication buffer. 	*/
	if ( iter->rnp[2] < zp_[0] && rank_ != 0 )
	  {
	    ubp.qSB.push_back( iter->q       );
	    ubp.qSB.push_back( iter->rnp [0] );
	    ubp.qSB.push_back( iter->rnp [1] );
	    ubp.qSB.push_back( iter->rnp [2] );
	    ubp.qSB.push_back( iter->gbnp[0] );
	    ubp.qSB.push_back( iter->gbnp[1] );
	    ubp.qSB.push_back( iter->gbnp[2] );
	    ubp.qSB.push_back( iter->rnm [0] );
	    ubp.qSB.push_back( iter->rnm [1] );
	    ubp.qSB.push_back( iter->rnm [2] );
	    ubp.qSB.push_back( iter->gbnm[0] );
	    ubp.qSB.push_back( iter->gbnm[1] );
	    ubp.qSB.push_back( iter->gbnm[2] );
	    ubp.qSB.push_back( iter->e 	     );
	    ubp.dq = true;
	  }
	else if ( iter->rnp[2] >= zp_[1] && rank_ != size_ - 1 )
	  {
	    ubp.qSF.push_back( iter->q       );
	    ubp.qSF.push_back( iter->rnp [0] );
	    ubp.qSF.push_back( iter->rnp [1] );
	    ubp.qSF.push_back( iter->rnp [2] );
	    ubp.qSF.push_back( iter->gbnp[0] );
	    ubp.qSF.push_back( iter->gbnp[1] );
	    ubp.qSF.push_back( iter->gbnp[2] );
	    ubp.qSF.push_back( iter->rnm [0] );
	    ubp.qSF.push_back( iter->rnm [1] );
	    ubp.qSF.push_back( iter->rnm [2] );
	    ubp.qSF.push_back( iter->gbnm[0] );
	    ubp.qSF.push_back( iter->gbnm[1] );
	    ubp.qSF.push_back( iter->gbnm[2] );
	    ubp.qSF.push_back( iter->e 	     );
	    ubp.dq = true;
	  }

	/* Delete the charge if it leaves the computational domain.					*/
	if (ubp.dq)
	  chargeVectorn_.erase(iter);
	else
	  ++iter;
      }

    /* Now communicate the charges which propagate throughout the borders to other processors.		*/
    if (rank_ != 0)
      MPI_Send(&ubp.qSB[0],ubp.qSB.size(),MPI_DOUBLE,rank_-1,msgtag1,MPI_COMM_WORLD);

    if (rank_ != size_ - 1)
      {
	MPI_Probe(rank_+1,msgtag1,MPI_COMM_WORLD,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&ub_.nL);
	ubp.qRF.resize(ub_.nL);
	MPI_Recv(&ubp.qRF[0],ub_.nL,MPI_DOUBLE,rank_+1,msgtag1,MPI_COMM_WORLD,&status);
      }

    if (rank_ != size_ - 1)
      MPI_Send(&ubp.qSF[0],ubp.qSF.size(),MPI_DOUBLE,rank_+1,msgtag2,MPI_COMM_WORLD);

    if (rank_ != 0)
      {
	MPI_Probe(rank_-1,msgtag2,MPI_COMM_WORLD,&status);
	MPI_Get_count(&status,MPI_DOUBLE,&ub_.nL);
	ubp.qRB.resize(ub_.nL);
	MPI_Recv(&ubp.qRB[0],ub_.nL,MPI_DOUBLE,rank_-1,msgtag2,MPI_COMM_WORLD,&status);
      }

    /* Now insert the newly incoming particles in this processor to the list of particles.            */
    unsigned i = 0;
    while ( i < ubp.qRF.size() )
      {
	ub_.Q.q       = ubp.qRF[i++];
	ub_.Q.rnp [0] = ubp.qRF[i++];
	ub_.Q.rnp [1] = ubp.qRF[i++];
	ub_.Q.rnp [2] = ubp.qRF[i++];
	ub_.Q.gbnp[0] = ubp.qRF[i++];
	ub_.Q.gbnp[1] = ubp.qRF[i++];
	ub_.Q.gbnp[2] = ubp.qRF[i++];
	ub_.Q.rnm [0] = ubp.qRF[i++];
	ub_.Q.rnm [1] = ubp.qRF[i++];
	ub_.Q.rnm [2] = ubp.qRF[i++];
	ub_.Q.gbnm[0] = ubp.qRF[i++];
	ub_.Q.gbnm[1] = ubp.qRF[i++];
	ub_.Q.gbnm[2] = ubp.qRF[i++];
	ub_.Q.e       = ubp.qRF[i++];

	chargeVectorn_.push_back(ub_.Q);
      }
    i = 0;
    while ( i < ubp.qRB.size() )
      {
	ub_.Q.q       = ubp.qRB[i++];
	ub_.Q.rnp [0] = ubp.qRB[i++];
	ub_.Q.rnp [1] = ubp.qRB[i++];
	ub_.Q.rnp [2] = ubp.qRB[i++];
	ub_.Q.gbnp[0] = ubp.qRB[i++];
	ub_.Q.gbnp[1] = ubp.qRB[i++];
	ub_.Q.gbnp[2] = ubp.qRB[i++];
	ub_.Q.rnm [0] = ubp.qRB[i++];
	ub_.Q.rnm [1] = ubp.qRB[i++];
	ub_.Q.rnm [2] = ubp.qRB[i++];
	ub_.Q.gbnm[0] = ubp.qRB[i++];
	ub_.Q.gbnm[1] = ubp.qRB[i++];
	ub_.Q.gbnm[2] = ubp.qRB[i++];
	ub_.Q.e       = ubp.qRB[i++];

	chargeVectorn_.push_back(ub_.Q);
      }

    /* Add the number of charge points that do not reside in the transverse size of the domain from
     * different processors.										*/
    MPI_Reduce(&ubp.nt, &ub_.nt,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if ( rank_ == 0 && ub_.nt > 0 )
      printmessage(std::string(__FILE__), __LINE__, std::string("Warning: " + stringify(ub_.nt) +
								" particles have transverse dimensions that are larger than the computational domain size.") );
  }

  /******************************************************************************************************
   * Sample the bunch data and save it to the given file.
   ******************************************************************************************************/

  void Solver::bunchSample ()
  {
    /* First, we need to evaluate the charge cloud properties and for that the corresponding data should
     * be initialized.                                                                     		*/
    sb_.q   = 0.0;
    sb_.r   = 0.0;
    sb_.r2  = 0.0;
    sb_.gb  = 0.0;
    sb_.gb2 = 0.0;

    /* Now, we perform an addition of all the charges. Since the point charges are equal, we do not need
     * to do weighted additions.                                            				*/
    for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
      {
	if ( ( iter->rnp[2] >= zp_[0] || rank_ == 0 ) && ( iter->rnp[2] < zp_[1] || rank_ == size_ - 1 ) )
	  {
	    sb_.q	 += iter->q;
	    sb_.r .pmv(   iter->q, iter->rnp );
	    sb_.gb.pmv(   iter->q, iter->gbnp);

	    for (int l = 0; l < 3; l++)
	      {
		sb_.r2 [l] += iter->rnp[l]  * iter->rnp[l]  * iter->q;
		sb_.gb2[l] += iter->gbnp[l] * iter->gbnp[l] * iter->q;
	      }
	  }
      }

    /* Add the contribution from each processor to the first element.					*/
    MPI_Reduce(&sb_.q    , &sb_.qT    , 1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sb_.r[0] , &sb_.rT[0] , 3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sb_.r2[0], &sb_.r2T[0], 3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sb_.gb[0], &sb_.gbT[0], 3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sb_.gb2[0],&sb_.gb2T[0],3,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

    if ( rank_ == 0 )
      {
	/* Divide the obtained values by the number of particles to calculate the mean value.		*/
	sb_.rT   /= sb_.qT;
	sb_.r2T  /= sb_.qT;
	sb_.gbT  /= sb_.qT;
	sb_.gb2T /= sb_.qT;

	(*sb_.file).setf(std::ios::scientific);
	(*sb_.file).precision(4);

	/** Write time into the first column.                                                 		*/
	*sb_.file << timeBunch_ << "\t";

	/** Now write the calculated values for the charge distribution in this row.          		*/
	*sb_.file << sb_.rT[0]   << "\t" << sb_.rT[1]   << "\t" << sb_.rT[2]   << "\t";
	*sb_.file << sb_.gbT[0]  << "\t" << sb_.gbT[1]  << "\t" << sb_.gbT[2]  << "\t";
	*sb_.file << sqrt( sb_.r2T[0]  - sb_.rT[0]  * sb_.rT[0]  ) << "\t" ;
	*sb_.file << sqrt( sb_.r2T[1]  - sb_.rT[1]  * sb_.rT[1]  ) << "\t" ;
	*sb_.file << sqrt( sb_.r2T[2]  - sb_.rT[2]  * sb_.rT[2]  ) << "\t" ;
	*sb_.file << sqrt( sb_.gb2T[0] - sb_.gbT[0] * sb_.gbT[0] ) << "\t" ;
	*sb_.file << sqrt( sb_.gb2T[1] - sb_.gbT[1] * sb_.gbT[1] ) << "\t" ;
	*sb_.file << sqrt( sb_.gb2T[2] - sb_.gbT[2] * sb_.gbT[2] ) << std::endl ;
      }
  }

  /******************************************************************************************************
   * Visualize the bunch as vtk files and save them to the file with given name.
   ******************************************************************************************************/

  void Solver::bunchVisualize ()
  {
    /* First define a parameters.                         						*/
    Double                            gamma, beta;

    /* The old files if existing should be deleted.                                               	*/
    vb_.fileName = bunch_.bunchVTKBasename_ + "-p" + stringify(rank_) + "-" + stringify(nTimeBunch_) + VTU_FILE_SUFFIX;
    vb_.file = new std::ofstream(vb_.fileName.c_str(),std::ios::trunc);

    vb_.file->setf(std::ios::scientific);
    vb_.file->precision(4);

    /* Store the number of particles in the simulation.							*/
    vb_.N = 0;
    for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
      if ( ( iter->rnp[2] >= zp_[0] || rank_ == 0 ) && ( iter->rnp[2] < zp_[1] || rank_ == size_ - 1 ) )
	vb_.N++;

    /* Write the initial data for the vtk file.                                                     	*/
    *vb_.file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
	<< std::endl;
    *vb_.file << "<UnstructuredGrid>"                                                     	<< std::endl;
    *vb_.file << "<Piece NumberOfPoints=\"" << vb_.N + 1 << "\" NumberOfCells=\"" << 1 << "\">"
	<< std::endl;

    /* Insert the coordinates of the grid for the charge points.                                      	*/
    *vb_.file << "<Points>"                                                             	<< std::endl;
    *vb_.file << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" 	<< std::endl;

    for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
      {
	if ( ( iter->rnp[2] >= zp_[0] || rank_ == 0 ) && ( iter->rnp[2] < zp_[1] || rank_ == size_ - 1 ) )
	  *vb_.file << iter->rnp[0] << " " << iter->rnp[1] << " " << iter->rnp[2] 		<< std::endl;
      }

    *vb_.file << r_[0][0] << " " << r_[0][1] << " " << r_[0][2]         			<< std::endl;
    *vb_.file << "</DataArray>"                                                          	<< std::endl;
    *vb_.file << "</Points>"                                                              	<< std::endl;

    /* Insert each cell vertices number into the vtk file.                                            	*/
    *vb_.file << "<Cells>"                                                               	<< std::endl;
    *vb_.file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"       	<< std::endl;
    for (unsigned i = 0; i < vb_.N + 1 ; ++i) *vb_.file << i << " ";
    *vb_.file                                                                         		<< std::endl;
    *vb_.file << "</DataArray>"                                                               	<< std::endl;
    *vb_.file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"         	<< std::endl;
    *vb_.file << vb_.N + 1                                                    			<< std::endl;
    *vb_.file << "</DataArray>"                                                            	<< std::endl;
    *vb_.file << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"              	<< std::endl;
    *vb_.file << 2                                                                       	<< std::endl;
    *vb_.file << "</DataArray>"                                                           	<< std::endl;
    *vb_.file << "</Cells>"                                                                	<< std::endl;

    *vb_.file << "<PointData Vectors = \"charge\">"                                    		<< std::endl;
    *vb_.file << "<DataArray type=\"Float64\" Name=\"charge\" NumberOfComponents=\"3\" format=\"ascii\">"
	<< std::endl;
    for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
      {
	if ( ( iter->rnp[2] >= zp_[0] || rank_ == 0 ) && ( iter->rnp[2] < zp_[1] || rank_ == size_ - 1 ) )
	  {
	    gamma = sqrt( 1.0 + iter->gbnp.norm2() );
	    beta  = iter->gbnp[2] / gamma;
	    *vb_.file << iter->q << " " <<  gamma * gamma_ * ( 1.0 + beta_ * beta )
            						    << " " << gamma * gamma_ * ( 1.0 + beta_ * beta ) * 0.512   << std::endl;
	  }
      }
    *vb_.file << 0.0 << " " << 0.0 << " " << 0.0						<< std::endl;
    *vb_.file << 0.0 << " " << 0.0 << " " << 0.0						<< std::endl;
    *vb_.file << "</DataArray>"                                                             	<< std::endl;
    *vb_.file << "</PointData>"                                                             	<< std::endl;
    *vb_.file << "</Piece>"                                                                 	<< std::endl;
    *vb_.file << "</UnstructuredGrid>"                                                      	<< std::endl;
    *vb_.file << "</VTKFile>"                                                              	<< std::endl;

    /* Close the file.                                                                          	*/
    (*vb_.file).close();

    /* Connect the vtk files in the root processor.							*/
    if ( rank_ == 0 )
      {
	/* Write the parallel vtk file for combining the files.						*/
	vb_.fileName = bunch_.bunchVTKBasename_ + "-" + stringify(nTimeBunch_) + PTU_FILE_SUFFIX;
	vb_.file = new std::ofstream(vb_.fileName.c_str(),std::ios::trunc);

	/* Obtain the path and name of the files.							*/
	splitFilename(bunch_.bunchVTKBasename_, vb_.path, vb_.name);

	*vb_.file << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"<< std::endl;
	*vb_.file << "<PUnstructuredGrid> GhostLevel = \"0\""                                        	<< std::endl;

	/* Insert the coordinates of the grid for the charge cloud.                                   	*/
	*vb_.file << "<PPoints>"                                                                    	<< std::endl;
	*vb_.file << "<PDataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\" />"
	    << std::endl;
	*vb_.file << "</PPoints>"                                                                  	<< std::endl;

	*vb_.file << "<PPointData>"                                                                	<< std::endl;
	*vb_.file << "<PDataArray type=\"Float64\" Name=\"charge\" NumberOfComponents=\"3\" format=\"ascii\" />"
	    << std::endl;
	*vb_.file << "</PPointData>"                                                                	<< std::endl;
	for (int i = 0; i < size_; ++i )
	  {
	    vb_.fileName = vb_.name + "-p" + stringify(i) + "-" + stringify(nTimeBunch_) + VTU_FILE_SUFFIX;
	    *vb_.file << "<Piece  Source=\"" << vb_.fileName << "\"/>"                         		<< std::endl;
	  }
	*vb_.file << "</PUnstructuredGrid>"                                                         	<< std::endl;
	*vb_.file << "</VTKFile>"                                                                  	<< std::endl;

	(*vb_.file).close();
      }
  }

  /******************************************************************************************************
   * Write the total profile of the field into the given file name.
   ******************************************************************************************************/

  void Solver::bunchProfile ()
  {
    /* The old files if existing should be deleted.                                           		*/
    pb_.fileName = bunch_.bunchProfileBasename_ + "-p" + stringify(rank_) + "-" + stringify(nTime_) + TXT_FILE_SUFFIX;
    pb_.file = new std::ofstream(pb_.fileName.c_str(),std::ios::trunc);

    (*pb_.file).setf(std::ios::scientific);
    (*pb_.file).precision(15);
    (*pb_.file).width(40);

    for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
      {
	if ( ( iter->rnp[2] >= zp_[0] || rank_ == 0 ) && ( iter->rnp[2] < zp_[1] || rank_ == size_ - 1 ) )
	  {
	    /* Loop over the particles and print the data of each particle into the file.		*/
	    *pb_.file << iter->q  	<< "\t";
	    *pb_.file << iter->rnp[0]  	<< "\t";
	    *pb_.file << iter->rnp[1]  	<< "\t";
	    *pb_.file << iter->rnp[2]  	<< "\t";
	    *pb_.file << iter->gbnp[0] 	<< "\t";
	    *pb_.file << iter->gbnp[1] 	<< "\t";
	    *pb_.file << iter->gbnp[2] 	<< std::endl;
	  }
      }

    /* Close the file.											*/
    (*pb_.file).close();
  }

  /******************************************************************************************************
   * Calculate the magnetic field of undulator and add it to the magnetic field of the seed.
   ******************************************************************************************************/

  void Solver::undulatorField (UpdateBunchParallel & ubp, FieldVector <Double> & r)
  {
    /* For each external field given by the user, add the external fields to the undulator field.     	*/
    for (std::vector<Undulator>::iterator iter = undulator_.begin(); iter != undulator_.end(); iter++)
      {

	/* Calculate the undulator magnetic field.							*/
	ub_.b0	= (iter->lu_ != 0.0 ) ? EM * c0_ * 2 * PI / iter->lu_ * iter->k_ / EC : 0.0;

	/* Calculate the undulator wave number.								*/
	ub_.ku 	= (iter->lu_ != 0.0 ) ? 2 * PI / iter->lu_ : 0.0;

	/* Calculate the sine and cosine functions of the undulator angle.				*/
	ub_.ct	= cos( iter->theta_ );
	ub_.st	= sin( iter->theta_ );

	if ( iter->type_ == STATIC )
	  {
	    /* First find the position with respect to the undulator begin point. The equation below
	     * assumes that the bunch at time zero resides in a distance gamma*rb_ from the first
	     * undulator.										*/
	    ubp.lz = gamma_ * ( r[2] + beta_ * c0_ * ( timeBunch_ + dt_ ) ) - iter->rb_;
	    ubp.ly = r[0] * ub_.ct + r[1] * ub_.st;

	    /* Now, calculate the undulator field according to the obtained position.               	*/
	    if ( ubp.lz >= 0.0 && ubp.lz <= iter->length_ * iter->lu_ )
	      {
		ubp.d1     = ub_.b0 * cosh(ub_.ku * ubp.ly) * sin(ub_.ku * ubp.lz) * gamma_;
		ubp.bt[0] += ubp.d1 * ub_.ct;
		ubp.bt[1] += ubp.d1 * ub_.st;
		ubp.bt[2] += ub_.b0 * sinh(ub_.ku * ubp.ly) * cos(ub_.ku * ubp.lz);

		ubp.d1    *= c0_ * beta_;
		ubp.et[1] +=   ubp.d1 * ub_.ct;
		ubp.et[0] += - ubp.d1 * ub_.st;
		ubp.et[2] += 0.0;
	      }
	    else if ( ubp.lz < 0.0 )
	      {
		ubp.sz = exp( - pow( ub_.ku * ubp.lz , 2 ) / 2.0 );

		if ( iter != undulator_.begin() )
		  {
		    ubp.i  = iter - undulator_.begin() - 1;
		    ubp.r0 = undulator_[ubp.i].rb_ + undulator_[ubp.i].length_ * undulator_[ubp.i].lu_ - iter->rb_;
		    if ( ubp.lz < ubp.r0 || ubp.r0 == 0.0 ) ubp.sz = 0.0;
		    else
		      ubp.sz *= 0.35875 + 0.48829 * cos( PI * ubp.lz / ubp.r0 ) + 0.14128 * cos( 2.0 * PI * ubp.lz / ubp.r0 ) + 0.01168 * cos( 3.0 * PI * ubp.lz / ubp.r0 );
		  }

		ubp.d1     = ub_.b0 * cosh(ub_.ku * ubp.ly ) * ubp.sz * ub_.ku * ubp.lz * gamma_;
		ubp.bt[0] += ubp.d1 * ub_.ct;
		ubp.bt[1] += ubp.d1 * ub_.st;
		ubp.bt[2] += ub_.b0 * sinh(ub_.ku * ubp.ly) * ubp.sz;

		ubp.d1    *= c0_ * beta_;
		ubp.et[1] +=   ubp.d1 * ub_.ct;
		ubp.et[0] += - ubp.d1 * ub_.st;
		ubp.et[2] += 0.0;
	      }
	    else if ( ubp.lz > iter->length_ * iter->lu_ )
	      {
		ubp.t0 = ubp.lz - iter->length_ * iter->lu_;

		ubp.sz = exp( - pow( ub_.ku *  ubp.t0 , 2 ) / 2.0 );

		if ( iter != undulator_.end() )
		  {
		    ubp.i  = iter - undulator_.begin() + 1;
		    ubp.r0 = undulator_[ubp.i].rb_ - iter->rb_ - iter->length_ * iter->lu_;
		    if ( ubp.t0 > ubp.r0 || ubp.r0 == 0.0 ) ubp.sz = 0.0;
		    else
		      ubp.sz *= 0.35875 + 0.48829 * cos( PI * ubp.t0 / ubp.r0 ) + 0.14128 * cos( 2.0 * PI * ubp.t0 / ubp.r0 ) + 0.01168 * cos( 3.0 * PI * ubp.t0 / ubp.r0 );
		  }

		ubp.d1     = ub_.b0 * cosh(ub_.ku * ubp.ly ) * ubp.sz * ub_.ku * ubp.t0 * gamma_;
		ubp.bt[0] += ubp.d1 * ub_.ct;
		ubp.bt[1] += ubp.d1 * ub_.st;
		ubp.bt[2] += ub_.b0 * sinh(ub_.ku * ubp.ly ) * ubp.sz;

		ubp.d1    *= c0_ * beta_;
		ubp.et[1] +=   ubp.d1 * ub_.ct;
		ubp.et[0] += - ubp.d1 * ub_.st;
		ubp.et[2] += 0.0;
	      }
	  }
	else if ( iter->type_ == OPTICAL )
	  {
	    /* Transfer the coordinate from the bunch rest frame to the lab frame.			*/
	    ubp.rl[0] = r[0]; ubp.rl[1] = r[1];
	    ubp.rl[2] = gamma_ * ( r[2] + beta_ * c0_ * ( timeBunch_ + dt_ ) );
	    ubp.t0    = gamma_ * ( timeBunch_ + dt_  + beta_ / c0_ * r[2] );

	    /* Calculate the distance to the reference position along the propagation direction.   	*/
	    ubp.rv = ubp.rl; ubp.rv -= iter->position_;
	    ubp.z  = ubp.rv * iter->direction_ ;

	    /* Compute propagation delay and subtract it from the time.                         	*/
	    ubp.tl = ubp.t0 - ubp.z / c0_;

	    /* Reset the carrier envelope phase of the pulse.						*/
	    ubp.p = 0.0;

	    /* Now manipulate the electric field vector depending on the specific seed given.		*/
	    if ( iter->seedType_ == PLANEWAVE )
	      {
		/* Retrieve signal value at corrected time.                                   		*/
		ubp.tsignal = iter->signal_.self(ubp.tl, ubp.p);

		/* Calculate the field only if the signal value is larger than a limit.			*/
		if ( fabs(ubp.tsignal) < 1.0e-100 ) { ubp.eT = 0.0; ubp.bT = 0.0; }
		else
		  {
		    /* Provide vector to store the electric field of the undulator.              	*/
		    ubp.eT.mv( iter->amplitude_ * ubp.tsignal, iter->polarization_ );

		    /* Provide vector to store the magnetic field of the undulator.         		*/
		    ubp.bT = cross( iter->direction_, iter->polarization_ );
		    ubp.bT.mv( iter->amplitude_ * ubp.tsignal / c0_, ubp.bT );
		  }
	      }
	    else if ( iter->seedType_ == PLANEWAVECONFINED )
	      {
		/* Retrieve signal value at corrected time.                                     	*/
		ubp.tsignal = iter->signal_.self(ubp.tl, ubp.p);

		/* Calculate the transverse distance to the center line.   				*/
		ubp.x  = ubp.rv * iter->polarization_;
		ubp.yv = cross( iter->direction_, iter->polarization_ );
		ubp.y  = ubp.rv * ubp.yv;

		/* Calculate the field only if the signal value is larger than a limit.			*/
		if ( fabs(ubp.tsignal) < 1.0e-100 || sqrt( pow(ubp.x/iter->radius_[0],2) + pow(ubp.y/iter->radius_[1],2) ) > 1.0 )
		  {
		    ubp.eT = 0.0;
		    ubp.bT = 0.0;
		  }
		else
		  {
		    /* Provide vector to store the electric field of the undulator.                 	*/
		    ubp.eT.mv( iter->amplitude_ * ubp.tsignal, iter->polarization_ );

		    /* Provide vector to store the magnetic field of the undulator.              	*/
		    ubp.bT = cross( iter->direction_, iter->polarization_ );
		    ubp.bT.mv( iter->amplitude_ * ubp.tsignal / c0_, ubp.bT );
		  }
	      }
	    else if ( iter->seedType_ == GAUSSIANBEAM )
	      {
		/* Retrieve signal value at corrected time.                                     	*/
		ubp.tsignal = iter->signal_.self(ubp.tl, ubp.p);

		if ( fabs(ubp.tsignal) < 1.0e-100 ) { ubp.eT = 0.0; ubp.bT = 0.0; }
		else
		  {
		    /* Provide vector to store transverse, longitudinal and total  electric field. 	*/
		    ubp.ex = iter->polarization_;
		    ubp.ez = iter->direction_;

		    /* Calculate the transverse distance to the center line.   				*/
		    ubp.x  = ubp.rv * iter->polarization_;
		    ubp.yv = cross( iter->direction_, iter->polarization_ );
		    ubp.y  = ubp.rv * ubp.yv;

		    /* Calculate the wavelength corresponding to the given central frequency.		*/
		    ubp.l = c0_ / iter->signal_.f0_;

		    /* Calculate the Rayleigh length and the relative radius of the beam.     		*/
		    ubp.zRp = PI * pow( iter->radius_[0] , 2 ) / ubp.l;
		    ubp.zRs = PI * pow( iter->radius_[1] , 2 ) / ubp.l;
		    ubp.wrp = sqrt( 1.0 + pow( ubp.z / ubp.zRp , 2 ) );
		    ubp.wrs = sqrt( 1.0 + pow( ubp.z / ubp.zRs , 2 ) );

		    /* Compute the transverse vector between the point and the reference point.   	*/
		    ubp.p   = 0.5 * ( atan( ubp.z / ubp.zRp ) + atan( ubp.z / ubp.zRs ) ) - PI * ubp.z / ubp.l *
			( pow( ubp.x / ( ubp.zRp * ubp.wrp) , 2 ) + pow( ubp.y / ( ubp.zRs * ubp.wrs ) , 2 ) );
		    ubp.t   = exp( - pow( ubp.x/(iter->radius_[0]*ubp.wrp), 2) - pow( ubp.y/(iter->radius_[1]*ubp.wrs), 2) ) / sqrt(ubp.wrs*ubp.wrp);

		    ubp.ex.mv( ubp.t * iter->amplitude_, 						iter->polarization_ );
		    ubp.ez.mv( ubp.t * iter->amplitude_ * ( - ubp.x / ( ubp.wrp * ubp.zRp ) ), 	iter->direction_    );
		    ubp.bz.mv( ubp.t * iter->amplitude_ * ( - ubp.y / ( ubp.wrs * ubp.zRs ) ) / c0_, iter->direction_    );
		    ubp.by = cross( iter->direction_, ubp.ex); ubp.by /= c0_;

		    /* Retrieve signal value at corrected time.                                       	*/
		    ubp.p 	-= PI/2.0;
		    ubp.tsignal	 = iter->signal_.self(ubp.tl, ubp.p);
		    ubp.ex 	*= ubp.tsignal; ubp.by *= ubp.tsignal;

		    ubp.p 	+= PI/2.0 + atan(ubp.z/ubp.zRp);
		    ubp.tsignal	 = iter->signal_.self(ubp.tl, ubp.p);
		    ubp.ez	*= ubp.tsignal;

		    ubp.p 	+= atan(ubp.z/ubp.zRs) - atan(ubp.z/ubp.zRp);
		    ubp.tsignal	 = iter->signal_.self(ubp.tl, ubp.p);
		    ubp.bz	*= ubp.tsignal;

		    /* Calculate the total electric and magnetic field.					*/
		    ubp.eT = ubp.ex; ubp.eT += ubp.ez;
		    ubp.bT = ubp.by; ubp.bT += ubp.bz;
		  }
	      }

	    /* Now transfer the computed magnetic vector potential into the bunch rest frame.		*/
	    ubp.bt[0] += gamma_ * ( ubp.bT[0] + beta_ / c0_ * ubp.eT[1] );
	    ubp.bt[1] += gamma_ * ( ubp.bT[1] - beta_ / c0_ * ubp.eT[0] );
	    ubp.bt[2] += ubp.bT[2];

	    ubp.et[0] += gamma_ * ( ubp.eT[0] - beta_ * c0_ * ubp.bT[1] );
	    ubp.et[1] += gamma_ * ( ubp.eT[1] + beta_ * c0_ * ubp.bT[0] );
	    ubp.et[2] += ubp.eT[2];
	  }
      }
  }

  /******************************************************************************************************
   * Calculate the field of external field and add it to the field of the seed.
   ******************************************************************************************************/

  void Solver::externalField (UpdateBunchParallel & ubp, FieldVector <Double> & r)
  {

    /* Transfer the coordinate from the bunch rest frame to the lab frame.                            	*/
    ubp.rl[0] = r[0];
    ubp.rl[1] = r[1];
    ubp.rl[2] = gamma_ * ( r[2] + beta_ * c0_ * ( timeBunch_ + dt_ ) );
    ubp.t0    = gamma_ * ( timeBunch_ + dt_  + beta_ / c0_ * r[2] );

    /* For each external field given by the user, add the external fields to the undulator field.     	*/
    for (std::vector<ExtField>::iterator iter = extField_.begin(); iter != extField_.end(); iter++)
      {
	/* Calculate the distance to the reference position along the propagation direction.          	*/
	ubp.rv  = ubp.rl; ubp.rv -= iter->position_;
	ubp.z   = ubp.rv * iter->direction_ ;

	/* Compute propagation delay and subtract it from the time.                                   	*/
	ubp.tl  = ubp.t0 - ubp.z / c0_;

	/* Reset the carrier envelope phase of the pulse.                                             	*/
	ubp.p   = 0.0;

	/* Now manipulate the electric field vector depending on the specific seed given.             	*/
	if ( iter->seedType_ == PLANEWAVE )
	  {
	    /* Retrieve signal value at corrected time.                                               	*/
	    ubp.tsignal = iter->signal_.self(ubp.tl, ubp.p);

	    /* Calculate the field only if the signal value is larger than a limit.                   	*/
	    if ( fabs(ubp.tsignal) < 1.0e-100 ) { ubp.eT = 0.0; ubp.bT = 0.0; }
	    else
	      {
		/* Provide vector to store the electric field of the external field.                  	*/
		ubp.eT.mv( iter->amplitude_ * ubp.tsignal, iter->polarization_ );

		/* Provide vector to store the magnetic field of the external field.                  	*/
		ubp.bT = cross( iter->direction_, iter->polarization_ );
		ubp.bT.mv( iter->amplitude_ * ubp.tsignal / c0_, ubp.bT );
	      }
	  }
	else if ( iter->seedType_ == PLANEWAVECONFINED )
	  {
	    /* Retrieve signal value at corrected time.                                               	*/
	    ubp.tsignal = iter->signal_.self(ubp.tl, ubp.p);

	    /* Calculate the transverse distance to the center line.   					*/
	    ubp.x  = ubp.rv * iter->polarization_;
	    ubp.yv = cross( iter->direction_, iter->polarization_ );
	    ubp.y  = ubp.rv * ubp.yv;

	    /* Calculate the field only if the signal value is larger than a limit.			*/
	    if ( fabs(ubp.tsignal) < 1.0e-100 || sqrt( pow(ubp.x/iter->radius_[0],2) + pow(ubp.y/iter->radius_[1],2) ) > 1.0 )
	      {
		ubp.eT = 0.0;
		ubp.bT = 0.0;
	      }
	    else
	      {
		/* Provide vector to store the electric field of the undulator.                       	*/
		ubp.eT.mv( iter->amplitude_ * ubp.tsignal, iter->polarization_ );

		/* Provide vector to store the magnetic field of the undulator.                       	*/
		ubp.bT = cross( iter->direction_, iter->polarization_ );
		ubp.bT.mv( iter->amplitude_ * ubp.tsignal / c0_, ubp.bT );
	      }
	  }
	else if ( iter->seedType_ == GAUSSIANBEAM )
	  {
	    /* Retrieve signal value at corrected time.                                               	*/
	    ubp.tsignal = iter->signal_.self(ubp.tl, ubp.p);

	    if ( fabs(ubp.tsignal) < 1.0e-100 ) { ubp.eT = 0.0; ubp.bT = 0.0; }
	    else
	      {
		/* Provide vector to store transverse, longitudinal and total  electric field.        	*/
		ubp.ex = iter->polarization_;
		ubp.ez = iter->direction_;

		/* Calculate the transverse distance to the center line.                              	*/
		ubp.x  = ubp.rv * iter->polarization_;
		ubp.yv = cross( iter->direction_, iter->polarization_ );
		ubp.y  = ubp.rv * ubp.yv;

		/* Calculate the wavelength corresponding to the given central frequency.             	*/
		ubp.l = c0_ / iter->signal_.f0_;

		/* Calculate the Rayleigh length and the relative radius of the beam.                 	*/
		ubp.zRp = PI * pow( iter->radius_[0] , 2 ) / ubp.l;
		ubp.zRs = PI * pow( iter->radius_[1] , 2 ) / ubp.l;
		ubp.wrp = sqrt( 1.0 + pow( ubp.z / ubp.zRp , 2 ) );
		ubp.wrs = sqrt( 1.0 + pow( ubp.z / ubp.zRs , 2 ) );

		/* Compute the transverse vector between the point and the reference point.           	*/
		ubp.p   = 0.5 * ( atan( ubp.z / ubp.zRp ) + atan( ubp.z / ubp.zRs ) ) - PI * ubp.z / ubp.l *
		    ( pow( ubp.x / ( ubp.zRp * ubp.wrp) , 2 ) + pow( ubp.y / ( ubp.zRs * ubp.wrs ) , 2 ) );
		ubp.t   = exp( - pow( ubp.x/(iter->radius_[0]*ubp.wrp), 2) - pow( ubp.y/(iter->radius_[1]*ubp.wrs), 2) ) / sqrt(ubp.wrs*ubp.wrp);

		ubp.ex.mv( ubp.t * iter->amplitude_,                                             iter->polarization_ );
		ubp.ez.mv( ubp.t * iter->amplitude_ * ( - ubp.x / ( ubp.wrp * ubp.zRp ) ),       iter->direction_    );
		ubp.bz.mv( ubp.t * iter->amplitude_ * ( - ubp.y / ( ubp.wrs * ubp.zRs ) ) / c0_, iter->direction_    );
		ubp.by  = cross( iter->direction_, ubp.ex);
		ubp.by /= c0_;

		/* Retrieve signal value at corrected time.                                           	*/
		ubp.p         -= PI/2.0;
		ubp.tsignal    = iter->signal_.self(ubp.tl, ubp.p);
		ubp.ex        *= ubp.tsignal; ubp.by *= ubp.tsignal;

		ubp.p         += PI/2.0 + atan(ubp.z/ubp.zRp);
		ubp.tsignal    = iter->signal_.self(ubp.tl, ubp.p);
		ubp.ez        *= ubp.tsignal;

		ubp.p         += atan(ubp.z/ubp.zRs) - atan(ubp.z/ubp.zRp);
		ubp.tsignal    = iter->signal_.self(ubp.tl, ubp.p);
		ubp.bz        *= ubp.tsignal;

		/* Calculate the total electric and magnetic field.                                   	*/
		ubp.eT = ubp.ex; ubp.eT += ubp.ez;
		ubp.bT = ubp.by; ubp.bT += ubp.bz;
	      }
	  }

	/* Now transfer the computed magnetic vector potential into the bunch rest frame.             	*/
	ubp.bt[0] += gamma_ * ( ubp.bT[0] + beta_ / c0_ * ubp.eT[1] );
	ubp.bt[1] += gamma_ * ( ubp.bT[1] - beta_ / c0_ * ubp.eT[0] );
	ubp.bt[2] += ubp.bT[2];

	ubp.et[0] += gamma_ * ( ubp.eT[0] - beta_ * c0_ * ubp.bT[1] );
	ubp.et[1] += gamma_ * ( ubp.eT[1] + beta_ * c0_ * ubp.bT[0] );
	ubp.et[2] += ubp.eT[2];
      }
  }

  /******************************************************************************************************
   * Initialize the data required for sampling and saving the radiation energy at the given position.
   ******************************************************************************************************/

  void Solver::initializeEnergySample ()
  {
    printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the FEL radiation energy data.") );
    re_.clear(); re_.resize(FEL_.size());

    /* Loop over the different FEL output parameters and calculate the radiation energy if the energy
     * calculation is activated.                                                                      	*/
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {
	/* Initialize if and only if the sampling of the power is enabled.                            	*/
	if (!FEL_[jf].radiationEnergy_.sampling_) continue;

	/* Perform the lorentz boost for the sampling data.                                           	*/
	for (unsigned int i = 0; i < FEL_[jf].radiationEnergy_.z_.size(); i++)
	  FEL_[jf].radiationEnergy_.z_[i] = FEL_[jf].radiationEnergy_.z_[i] * gamma_;
	FEL_[jf].radiationEnergy_.lineBegin_ = FEL_[jf].radiationEnergy_.lineBegin_ * gamma_;
	FEL_[jf].radiationEnergy_.lineEnd_   = FEL_[jf].radiationEnergy_.lineEnd_   * gamma_;

	/* If sampling type is plot over line initialize the positions according to the line begin and line
	 * end.                                                                                       	*/
	if ( FEL_[jf].radiationEnergy_.samplingType_ == OVERLINE )
	  {
	    Double l = 0.0;
	    FieldVector<Double> position;
	    while ( fabs(l) < fabs(FEL_[jf].radiationEnergy_.lineEnd_ - FEL_[jf].radiationEnergy_.lineBegin_) )
	      {
		FEL_[jf].radiationEnergy_.z_.push_back( FEL_[jf].radiationEnergy_.lineBegin_ + l);
		l += FEL_[jf].radiationEnergy_.res_ * gamma_;
	      }
	  }

	re_[jf].N  = FEL_[jf].radiationEnergy_.z_.size();

	/* Add the normalized wavelength sweep to the vector of wavelengths.                          	*/
	for (Double rw = FEL_[jf].radiationEnergy_.lambdaMin_; rw < FEL_[jf].radiationEnergy_.lambdaMax_; rw += FEL_[jf].radiationEnergy_.lambdaRes_)
	  FEL_[jf].radiationEnergy_.lambda_.push_back(rw);
	re_[jf].Nl = FEL_[jf].radiationEnergy_.lambda_.size();

	/* Based on the number of points to calculate and the size of the threads, allocate memory for
	 * saving the powers.                                                                         	*/
	re_[jf].pL.resize(re_[jf].Nl * re_[jf].N, 0.0);
	re_[jf].pG.resize(re_[jf].Nl * re_[jf].N, 0.0);

	re_[jf].Nf = 0;

	/* Initialize the file streams to save the data.                                       		*/
	re_[jf].file.resize(re_[jf].Nl);
	re_[jf].w.resize(re_[jf].Nl);
	for (unsigned int i = 0; i < re_[jf].Nl; i++)
	  {
	    std::string baseFilename = "";
	    if (!(isabsolute(FEL_[jf].radiationEnergy_.basename_))) baseFilename = FEL_[jf].radiationEnergy_.directory_;
	    baseFilename += FEL_[jf].radiationEnergy_.basename_ + "-" + stringify(i) + TXT_FILE_SUFFIX;

	    /* If the directory of the baseFilename does not exist create this directory.		*/
	    createDirectory(baseFilename, rank_);

	    re_[jf].file[i] = new std::ofstream(baseFilename.c_str(),std::ios::trunc);
	    ( *(re_[jf].file[i]) ).setf(std::ios::scientific);
	    ( *(re_[jf].file[i]) ).precision(15);
	    ( *(re_[jf].file[i]) ).width(40);

	    /* Determine the number of time points needed to calculate the amplitude of each
	     * radiation harmonic.                                                                    	*/
	    Double dt = undulator_[0].lu_ / FEL_[jf].radiationEnergy_.lambda_[i] / ( gamma_ * c0_ );
	    re_[jf].Nf = ( unsigned( dt / mesh_.timeStep_ ) > re_[jf].Nf ) ? unsigned( dt/mesh_.timeStep_ ) : re_[jf].Nf;

	    /* Calculate the angular frequency for each wavelength.                                   	*/
	    re_[jf].w[i] = 2 * PI / dt;
	  }

	/* Based on the obtained Nf, resize the vectors for saving the time domain data.          	*/
	re_[jf].fdt.resize(re_[jf].Nf, std::vector<std::vector<Double> > (re_[jf].N * N1_ * N0_, std::vector<Double> (4,0.0) ) );

	re_[jf].dt = mesh_.timeStep_;
	re_[jf].dx = mesh_.meshResolution_[0];
	re_[jf].dy = mesh_.meshResolution_[1];
	re_[jf].dz = mesh_.meshResolution_[2];

	re_[jf].pc  = 2.0 * re_[jf].dx * re_[jf].dy / ( m0_ * re_[jf].Nf * re_[jf].Nf ) * pow(mesh_.lengthScale_,2) / pow(mesh_.timeScale_,3);
      }

    printmessage(std::string(__FILE__), __LINE__, std::string(" The FEL radiation energy data is initialized. :::") );
  }

  /******************************************************************************************************
   * Sample the radiation energy at the given position and save it to the file.
   ******************************************************************************************************/

  void Solver::energySample ()
  {
    /* Declare the temporary parameters needed for calculating the radiated power.                    	*/
    unsigned int              mi, ni;
    FieldVector<Double>       et, bt;
    Complex                   ew1, bw1, ew2, bw2, ex;

    /* Loop over the different FEL output parameters and calculate the radiation energy if the energy
     * calculation is activated.                                                                      	*/
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {
	/* Initialize if and only if the sampling of the power is enabled.                            	*/
	if (!FEL_[jf].radiationEnergy_.sampling_) continue;

	/* First reset all the previously calculated powers.                                         	*/
	for (unsigned k = 0; k < re_[jf].N; ++k)
	  for (unsigned l = 0; l < re_[jf].Nl; ++l)
	    re_[jf].pL[k * re_[jf].Nl + l] = 0.0;

	/* Loop over the sampling positions, transverse dicretizations, and frequency to calculate the
	 * radiated power at the specific point and frequency.                                        	*/
	for (unsigned k = 0; k < re_[jf].N; ++k)
	  {
	    /* Obtain the z index of the cell containing the point.                                   	*/
	    re_[jf].dzr = modf( ( FEL_[jf].radiationEnergy_.z_[k] - ( mesh_.meshCenter_[2] - mesh_.meshLength_[2] / 2.0 ) ) / mesh_.meshResolution_[2] , &re_[jf].c);
	    re_[jf].k   = (int) re_[jf].c;

	    /* Do not continue if this index is not supported by the processor.                       	*/
	    if ( re_[jf].k < k0_ || re_[jf].k > k0_ + np_ - 1) continue;

	    /* Loop over the transverse indices.                                                      	*/
	    for (int i = 2; i < N0_ - 2; i += 1)
	      for (int j = 2; j < N1_ - 2; j += 1)
		{
		  /* Get the index in the computation grid as well as the field storage grid. 		*/
		  mi = ( re_[jf].k - k0_) * N1_* N0_ + i*N1_ + j;
		  ni = k * N1_* N0_ + i*N1_ + j;

		  /* Calculate and interpolate the electric field to find the value at the bunch point.	*/
		  et[0] = ( 1.0 - re_[jf].dzr ) * en_[mi][0] + re_[jf].dzr * en_[mi+N1N0_][0];
		  et[1] = ( 1.0 - re_[jf].dzr ) * en_[mi][1] + re_[jf].dzr * en_[mi+N1N0_][1];

		  /* Calculate and interpolate the magnetic field to find its value at the bunch point.	*/
		  bt[0] = ( 1.0 - re_[jf].dzr ) * bn_[mi][0] + re_[jf].dzr * bn_[mi+N1N0_][0];
		  bt[1] = ( 1.0 - re_[jf].dzr ) * bn_[mi][1] + re_[jf].dzr * bn_[mi+N1N0_][1];

		  /* Transform the fields to the lab frame.                                   */
		  re_[jf].fdt[nTime_ % re_[jf].Nf][ni][0] = gamma_ * ( et[0] + c0_ * beta_ * bt[1] );
		  re_[jf].fdt[nTime_ % re_[jf].Nf][ni][1] = gamma_ * ( et[1] - c0_ * beta_ * bt[0] );

		  re_[jf].fdt[nTime_ % re_[jf].Nf][ni][2] = gamma_ * ( bt[0] - beta_ / c0_ * et[1] );
		  re_[jf].fdt[nTime_ % re_[jf].Nf][ni][3] = gamma_ * ( bt[1] + beta_ / c0_ * et[0] );

		  /* Add the contribution of this field to the radiation power.                       	*/
		  for (unsigned l = 0; l < re_[jf].Nl; l++)
		    {
		      ew1 = Complex (0.0, 0.0);
		      bw1 = Complex (0.0, 0.0);
		      for (unsigned m = 0; m < re_[jf].Nf; m++)
			{
			  ex  = exp( I * ( re_[jf].w[l] * m * mesh_.timeStep_ ) );
			  ew1 += re_[jf].fdt[m][ni][0] * ex;
			  bw1 += re_[jf].fdt[m][ni][3] / ex;
			}

		      ew2 = Complex (0.0, 0.0);
		      bw2 = Complex (0.0, 0.0);
		      for (unsigned m = 0; m < re_[jf].Nf; m++)
			{
			  ex  = exp( I * ( re_[jf].w[l] * m * mesh_.timeStep_ ) );
			  ew2 += re_[jf].fdt[m][ni][1] * ex;
			  bw2 += re_[jf].fdt[m][ni][2] / ex;
			}

		      re_[jf].pL[k * re_[jf].Nl + l] += re_[jf].pc * ( std::real( ew1 * bw1 ) - std::real( ew2 * bw2 ) );
		    }
		}
	  }

	/* Add the data from each processor together at the root processor.                           	*/
	MPI_Allreduce(&re_[jf].pL[0],&re_[jf].pG[0],re_[jf].N*re_[jf].Nl,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	/* If the rank of the processor is equal to zero, i.e. root processor save the fields into the
	 * given file.                                                                                	*/
	for (unsigned l = 0; l < re_[jf].Nl; l++)
	  {
	    if ( rank_ == int( l % size_ ) )
	      {
		for (unsigned k = 0; k < re_[jf].N; ++k)
		  *(re_[jf].file[l]) << gamma_ * ( FEL_[jf].radiationEnergy_.z_[k] + beta_ * c0_ * timeBunch_ )
		  - pow( gamma_ * beta_ , 2 ) * undulator_[0].rb_ << "\t" << re_[jf].pG[k * re_[jf].Nl + l] << "\t";
		*(re_[jf].file[l]) << std::endl;
	      }
	  }
      }
  }

  /******************************************************************************************************
   * Initialize the data required for storing particles hitting a screen at the given position.
   ******************************************************************************************************/

  void Solver::initializeScreenProfile ()
  {
    scrp_.clear(); scrp_.resize(FEL_.size());
    printmessage(std::string(__FILE__), __LINE__, std::string("::: Initializing the screen to measure bunch profile.") );
    for ( unsigned int jf = 0; jf < FEL_.size(); jf++)
      {
	/* Initialize if and only if screens have been enabled.		                            	*/
	if (!FEL_[jf].screenProfile_.sampling_) continue;

	if (!(isabsolute(FEL_[jf].screenProfile_.basename_))) FEL_[jf].screenProfile_.basename_ = FEL_[jf].screenProfile_.directory_ + FEL_[jf].screenProfile_.basename_;

	/* According to the given rhythm add to the position vector.					*/
	Double	z = 0.0;
	while ( z < ( undulator_.end()->rb_ + undulator_.end()->length_ * undulator_.end()->lu_ ) )
	  {
	    FEL_[jf].screenProfile_.pos_.push_back(z);
	    z += FEL_[jf].screenProfile_.rhythm_;
	  }

	/*  resize vectors storing filenames                   */
	scrp_[jf].fileNames.resize((FEL_[jf].screenProfile_.pos_).size());
	scrp_[jf].files.resize((FEL_[jf].screenProfile_.pos_).size());

	/* If the directory of the baseFilename does not exist create this directory.			*/
	createDirectory(FEL_[jf].screenProfile_.basename_, rank_);

	/* Order screens                          							*/
	std::sort(FEL_[jf].screenProfile_.pos_.begin(), FEL_[jf].screenProfile_.pos_.end());

	MPI_Barrier(MPI_COMM_WORLD);

	/* Open files for each screen and write its position in first line.                    	       	*/
	for ( unsigned int i = 0; i < (FEL_[jf].screenProfile_.pos_).size(); i++ ){
	    printmessage(std::string(__FILE__), __LINE__, std::string("Screen ") + stringify(i) + std::string(" is at distance ") + stringify(FEL_[jf].screenProfile_.pos_[i]) + std::string(" from the undulator beginning.") );
	    scrp_[jf].fileNames[i] = FEL_[jf].screenProfile_.basename_ + "-p" + stringify(rank_) + "-screen" + stringify(i) + TXT_FILE_SUFFIX;
	    scrp_[jf].files[i] = new std::ofstream(scrp_[jf].fileNames[i].c_str(),std::ios::trunc);
	    (*scrp_[jf].files[i]).setf(std::ios::scientific);
	    (*scrp_[jf].files[i]).precision(15);
	    (*scrp_[jf].files[i]).width(40);
	    //*scrp_[jf].files[i] << "Screen distance from the beginning of the undulator = " << stringify(FEL_[jf].screenProfile_.pos_[i]) << std::endl;
	}

	/* Return an error if the screen sampling is activated but no screen is given.	       		*/
	if ( FEL_[jf].screenProfile_.pos_.size() == 0 )
	  {
	    printmessage(std::string(__FILE__), __LINE__, std::string("No position is set for the screen although the screen sampling is activated !!!") );
	    exit(1);
	  }

	printmessage(std::string(__FILE__), __LINE__, std::string(" The screen profiling data is initialized. :::") );
      }
  }

  /******************************************************************************************************
   * Store the bunch profile from particles hitting a screen at the given position and save it to the file.
   ******************************************************************************************************/

  void Solver::screenProfile ()
  {
    for ( unsigned jf = 0; jf < FEL_.size(); jf++)
      {
	/* Continued past this function if screen profile sampling is not enabled.			*/
	if (!FEL_[jf].screenProfile_.sampling_) continue;

	/* Iterate through all the screens.								*/
	for (unsigned i = 0; i < (FEL_[jf].screenProfile_.pos_).size(); i++ )
	  {
	    Double lzScreen = FEL_[jf].screenProfile_.pos_[i];
	    for (auto iter = chargeVectorn_.begin(); iter != chargeVectorn_.end(); iter++)
	      {
		/* Only look at particles which belong to the domain of this processor.			*/
		if ( ( iter->rnp[2] >= zp_[0] || rank_ == 0 ) && ( iter->rnp[2] < zp_[1] || rank_ == size_ - 1 ) )
		  {
		    Double lzm = gamma_ * ( iter->rnm[2] + beta_ * c0_ * ( timeBunch_- mesh_.timeStep_ + dt_ ) );
		    if ( lzm >= lzScreen )	continue;

		    Double lzp = gamma_ * ( iter->rnp[2] + beta_ * c0_ * ( timeBunch_ + dt_) );
		    if ( lzp <  lzScreen )	continue;

		    /* Interpolate quantities and write them in file.					*/
		    // *scrp_[jf].files[i] << iter->q  	        << "\t";
		    *scrp_[jf].files[i] << interp( lzm, lzp, iter->rnm[0], iter->rnp[0], lzScreen ) << "\t";
		    *scrp_[jf].files[i] << interp( lzm, lzp, iter->rnm[1], iter->rnp[1], lzScreen ) << "\t";
		    Double tm  = gamma_ * ( timeBunch_ + dt_ - mesh_.timeStep_	+ beta_ / c0_ * iter->rnm[2] );
		    Double tp  = gamma_ * ( timeBunch_ + dt_	  		+ beta_ / c0_ * iter->rnp[2] );
		    *scrp_[jf].files[i] << interp( lzm, lzp, tm, tp, lzScreen ) << "\t";

		    /* Use different positions since momenta are found at (time - 1/2*bunchTimeStep).	*/
		    Double gm = sqrt( 1 + pow(iter->gbnm[0], 2) + pow(iter->gbnm[1], 2) + pow(iter->gbnm[2], 2) );
		    lzm -= gamma_ * c0_ *  .5 * bunch_.timeStep_ * ( iter->gbnm[2] / gm + beta_ );
		    Double gp = sqrt( 1 + pow(iter->gbnp[0], 2) + pow(iter->gbnp[1], 2) + pow(iter->gbnp[2], 2) );
		    lzp -= gamma_ * c0_ *  .5 * bunch_.timeStep_ * ( iter->gbnp[2] / gp + beta_ );
		    Double gbzm = gamma_ * ( iter->gbnm[2] + beta_ * gm );
		    Double gbzp = gamma_ * ( iter->gbnp[2] + beta_ * gp );

		    /* Write the momentum data into the file.						*/
		    *scrp_[jf].files[i] << interp( lzm, lzp, iter->gbnm[0], iter->gbnp[0], lzScreen ) << "\t";
		    *scrp_[jf].files[i] << interp( lzm, lzp, iter->gbnm[1], iter->gbnp[1], lzScreen ) << "\t";
		    *scrp_[jf].files[i] << interp( lzm, lzp, gbzm, gbzp, lzScreen ) << std::endl;
		  }
	      }
	  }
      }
  }

  /******************************************************************************************************
   * Finalize the field calculations.
   ******************************************************************************************************/

  void Solver::finalize ()
  {
    /* If sampling was active, close the file since the simulation is finished now.			*/
    if (seed_.sampling_ && sf_.N > 0) 	(*sf_.file).close();

    /* If bunch sampling was active, close the file since the simulation is finished now.		*/
    if (bunch_.sampling_) 	(*sb_.file).close();
  }

  /******************************************************************************************************
   * Define the boolean function for comparing undulator begins.
   ******************************************************************************************************/

  bool Solver::undulatorCompare (Undulator i, Undulator j)
  { return ( i.rb_ < j.rb_ ); }

  /****************************************************************************************************
   * Define the function for linear interpolation.
   ****************************************************************************************************/

  Double Solver::interp( Double x0, Double x1, Double y0, Double y1, Double x )
  {
    return y0 + ( x - x0 ) / ( x1 - x0 ) * ( y1 - y0 );
  }

}
