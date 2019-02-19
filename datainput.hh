/********************************************************************************************************
 *  datainput.hh : Implementation of the parameter parser for the code
 ******************************************************************************************************** */

#ifndef DATAINPUT_HH_
#define DATAINPUT_HH_

namespace Darius
{
  /* The class of functions used for reading the text file of parameters and parsing them to the darius
   * solver.                                                                                            */
  class ParseDarius
  {
  private:

    /* The parameters and data files needed for parsing the values.                                     */
    std::list<std::string>&		jobFile_;
    Mesh& 				mesh_;
    Bunch&				bunch_;
    Seed&				seed_;
    std::vector<Undulator>&		undulator_;
    std::vector<ExtField>&              extField_;
    std::vector<FreeElectronLaser>&	FEL_;

  public:

    ParseDarius (std::list<std::string>& jobFile, Mesh& mesh, Bunch& bunch, Seed& seed,
		 std::vector<Undulator>& undulator, std::vector<ExtField>& extField,
		 std::vector<FreeElectronLaser>& FEL)
  :   jobFile_ ( jobFile ), mesh_ ( mesh ), bunch_ ( bunch ), seed_ ( seed ), undulator_ ( undulator ),
      extField_ ( extField ), FEL_ ( FEL )
  {
  };

    /* Read the parameters from the file and set all the parsed parameters for FEL simulation.         	*/
    void setJobParameters()
    {
      std::list<std::string>::iterator iter = jobFile_.begin();
      do
	{
	  /* Parse the solver parameters.                                                           	*/
	  if      (*iter == "MESH")		this->readMesh		( iter );

	  /* Parse the materials parameters.                                                         	*/
	  else if (*iter == "BUNCH")            this->readBunch		( iter );

	  /* Parse the surfaces parameters.                                                         	*/
	  else if (*iter == "FIELD")        	this->readField		( iter );

	  /* Parse the undulator parameters.                                                        	*/
	  else if (*iter == "UNDULATOR")        this->readUndulator	( iter );

	  /* Parse the external field parameters.                                                       */
	  else if (*iter == "EXTERNAL-FIELD")   this->readExtField      ( iter );

	  /* Parse the FEL output parameters.                                                        	*/
	  else if (*iter == "FEL-OUTPUT")       this->readFEL		( iter );

	  else { std::cout << (*iter) << " is not a defined group." << std::endl; exit(1); }

	  ++iter;
	}
      while (iter != jobFile_.end());
    };

    /* Read the parameters parsed for the mesh in the solver.                                    	*/
    void readMesh ( std::list<std::string>::iterator& iter )
    {
      ++iter;
      if (*iter != "{") { std::cout << "The solver directory is empty" << std::endl; exit(1); }
      else ++iter;
      do
	{
	  if (parameterName(*iter) == "length-scale")
	    {
	      if          (stringValue(*iter) == "METER")       mesh_.lengthScale_    	= 1.0;
	      else if     (stringValue(*iter) == "DECIMETER")   mesh_.lengthScale_    	= 1.0e-1;
	      else if     (stringValue(*iter) == "CENTIMETER")  mesh_.lengthScale_    	= 1.0e-2;
	      else if     (stringValue(*iter) == "MILLIMETER")  mesh_.lengthScale_    	= 1.0e-3;
	      else if     (stringValue(*iter) == "MICROMETER")  mesh_.lengthScale_    	= 1.0e-6;
	      else if     (stringValue(*iter) == "NANOMETER")   mesh_.lengthScale_    	= 1.0e-9;
	      else if     (stringValue(*iter) == "ANGSTROM")    mesh_.lengthScale_	= 1.0e-10;
	      else                                              mesh_.lengthScale_	= doubleValue(*iter);
	    }
	  else if (parameterName(*iter) == "mesh-lengths")
	    {
	      std::vector<Double> meshLength = vectorDoubleValue(*iter);
	      mesh_.meshLength_ = meshLength;
	    }
	  else if (parameterName(*iter) == "mesh-resolution")
	    {
	      std::vector<Double> meshResolution = vectorDoubleValue(*iter);
	      mesh_.meshResolution_ = meshResolution;
	    }
	  else if (parameterName(*iter) == "mesh-center")
	    {
	      std::vector<Double> meshCenter = vectorDoubleValue(*iter);
	      mesh_.meshCenter_ = meshCenter;
	    }
	  else if (parameterName(*iter) == "total-time")        	mesh_.totalTime_ 	= doubleValue(*iter);
	  else if (parameterName(*iter) == "bunch-time-step")   	bunch_.timeStep_ 	= doubleValue(*iter);
	  else if (parameterName(*iter) == "bunch-time-start")        	bunch_.timeStart_	= doubleValue(*iter);
	  else if (parameterName(*iter) == "mesh-truncation-order")
	    {
	      mesh_.truncationOrder_ 	= intValue(*iter);
	      if ( mesh_.truncationOrder_ != 1 && mesh_.truncationOrder_ != 2 )
		{
		  printmessage(std::string(__FILE__), __LINE__, std::string("Mesh truncation order can not be different from one or two.") );
		  exit(1);
		}
	    }
	  else if (parameterName(*iter) == "time-scale")
	    {
	      if          (stringValue(*iter) == "SECOND")      	mesh_.timeScale_      	= 1.0;
	      else if     (stringValue(*iter) == "MILLISECOND") 	mesh_.timeScale_      	= 1.0e-3;
	      else if     (stringValue(*iter) == "MICROSECOND") 	mesh_.timeScale_      	= 1.0e-6;
	      else if     (stringValue(*iter) == "NANOSECOND")  	mesh_.timeScale_      	= 1.0e-9;
	      else if     (stringValue(*iter) == "PICOSECOND")  	mesh_.timeScale_      	= 1.0e-12;
	      else if     (stringValue(*iter) == "FEMTOSECOND") 	mesh_.timeScale_      	= 1.0e-15;
	      else if     (stringValue(*iter) == "ATTOSECOND")  	mesh_.timeScale_      	= 1.0e-18;
	      else                                              	mesh_.timeScale_      	= doubleValue(*iter);
	    }
	  else if (parameterName(*iter) == "space-charge") 		mesh_.spaceCharge_ 	= boolValue(*iter);
	  else { std::cout << parameterName(*iter) << " is not defined in solver group." << std::endl; exit(1); }
	  ++iter;
	}
      while (*iter != "}");
    }

    /* Read the parameters parsed for the bunch in the darius solver.                     		*/
    void readBunch ( std::list<std::string>::iterator& iter )
    {
      ++iter;
      if (*iter != "{") { std::cout << "The bunch directory is empty" << std::endl; exit(1); }
      else ++iter;
      do
	{
	  /* The data related to the bunch initialization are to be parsed.                            */
	  if (*iter == "bunch-initialization")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The bunch-initialization directory is empty" << std::endl; exit(1); }
	      else ++iter;

	      BunchInitialize	bunchInit;

	      do
		{
		  if      (parameterName(*iter) == "type")              	bunchInit.bunchType_		= stringValue(*iter);
		  else if (parameterName(*iter) == "distribution") 		bunchInit.distribution_		= stringValue(*iter);
		  else if (parameterName(*iter) == "charge")            	bunchInit.cloudCharge_		= doubleValue(*iter);
		  else if (parameterName(*iter) == "number-of-particles")       bunchInit.numberOfParticles_       = intValue(*iter);
		  else if (parameterName(*iter) == "gamma")
		    {
		      bunchInit.initialGamma_ 	= doubleValue(*iter);
		      bunchInit.initialBeta_ 	= sqrt( 1.0 - 1.0 / ( bunchInit.initialGamma_ * bunchInit.initialGamma_ ) );
		    }
		  else if (parameterName(*iter) == "direction")
		    {
		      std::vector<Double> direction 	= vectorDoubleValue(*iter);
		      bunchInit.initialDirection_	= direction;
		    }
		  else if (parameterName(*iter) == "position")
		    {
		      std::vector<Double> position 	= vectorDoubleValue(*iter);
		      FieldVector<Double> p (0.0); p	= position;
		      bunchInit.position_.push_back(p);
		    }
		  else if (parameterName(*iter) == "numbers")
		    {
		      std::vector<unsigned int> numbers = vectorIntValue(*iter);
		      bunchInit.numbers_		= numbers;
		    }
		  else if (parameterName(*iter) == "lattice-constants")
		    {
		      std::vector<Double> lattice	= vectorDoubleValue(*iter);
		      bunchInit.latticeConstants_	= lattice;
		    }
		  else if (parameterName(*iter) == "sigma-position")
		    {
		      std::vector<Double> sigma		= vectorDoubleValue(*iter);
		      bunchInit.sigmaPosition_		= sigma;
		    }
		  else if (parameterName(*iter) == "sigma-momentum")
		    {
		      std::vector<Double> sigma		= vectorDoubleValue(*iter);
		      bunchInit.sigmaGammaBeta_		= sigma;
		    }
		  else if (parameterName(*iter) == "transverse-truncation")     	bunchInit.tranTrun_		= doubleValue(*iter);
		  else if (parameterName(*iter) == "longitudinal-truncation")   	bunchInit.longTrun_		= doubleValue(*iter);
		  else if (parameterName(*iter) == "file")              		bunchInit.fileName_		= stringValue(*iter);
		  else if (parameterName(*iter) == "bunching-factor")			bunchInit.bF_			= doubleValue(*iter);
		  else if (parameterName(*iter) == "bunching-factor-phase")		bunchInit.bFP_			= doubleValue(*iter);
		  else if (parameterName(*iter) == "shot-noise")			bunchInit.shotNoise_		= boolValue(*iter);
		  else { std::cout << parameterName(*iter) << " is not defined in bunch-initialization group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");

	      /* Add the bunch to the set of bunches in the database.				*/
	      (bunch_.bunchInit_).push_back(bunchInit);
	    }

	  /* The data related to the bunch sampling are to be parsed.                          	*/
	  else if (*iter == "bunch-sampling")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The bunch-sampling directory is empty" << std::endl; exit(1); }
	      else ++iter;
	      do
		{
		  if 	  (parameterName(*iter) == "sample") 		bunch_.sampling_	= boolValue(*iter);
		  else if (parameterName(*iter) == "directory")         bunch_.directory_	= stringValue(*iter);
		  else if (parameterName(*iter) == "base-name")         bunch_.basename_	= stringValue(*iter);
		  else if (parameterName(*iter) == "rhythm")            bunch_.rhythm_   	= doubleValue(*iter);
		  else { std::cout << parameterName(*iter) << " is not defined in the bunch-sampling group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");
	    }

	  /* The data related to the electron acceleration visualization are to be parsed.              */
	  else if (*iter == "bunch-visualization")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The bunch-visualization directory is empty" << std::endl; exit(1); }
	      else ++iter;
	      do
		{
		  if      (parameterName(*iter) == "sample")    		bunch_.bunchVTK_                = boolValue(*iter);
		  else if (parameterName(*iter) == "directory")         	bunch_.bunchVTKDirectory_	= stringValue(*iter);
		  else if (parameterName(*iter) == "base-name") 		bunch_.bunchVTKBasename_        = stringValue(*iter);
		  else if (parameterName(*iter) == "rhythm")    		bunch_.bunchVTKRhythm_          = doubleValue(*iter);
		  else { std::cout << parameterName(*iter) << " is not defined in the bunch-visualization group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");
	    }
	  /* The data related to the field emission bunch report are to be parsed.                      */
	  else if (*iter == "bunch-profile")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The bunch-profile directory is empty" << std::endl; exit(1); }
	      else ++iter;
	      do
		{
		  if      (parameterName(*iter) == "sample")    bunch_.bunchProfile_            = boolValue(*iter);
		  else if (parameterName(*iter) == "directory") bunch_.bunchProfileDirectory_	= stringValue(*iter);
		  else if (parameterName(*iter) == "base-name") bunch_.bunchProfileBasename_    = stringValue(*iter);
		  else if (parameterName(*iter) == "time")      (bunch_.bunchProfileTime_).push_back(doubleValue(*iter));
		  else if (parameterName(*iter) == "rhythm")    bunch_.bunchProfileRhythm_      = doubleValue(*iter);
		  else { std::cout << parameterName(*iter) << " is not defined in the bunch-profile group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");
	    }
	  else { std::cout << parameterName(*iter) << " is not defined in the bunch-profile group." << std::endl; exit(1); }
	  ++iter;
	}
      while (*iter != "}");
    }

    /* Read the parameters parsed for the seed in the darius solver.                               	*/
    void readField ( std::list<std::string>::iterator& iter )
    {
      ++iter;
      if (*iter != "{") { std::cout << "The seed directory is empty" << std::endl; exit(1); }
      else ++iter;
      do
	{
	  /* The data related to the seed initialization are to be parsed.                            	*/
	  if (*iter == "field-initialization")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The seed-initialization directory is empty" << std::endl; exit(1); }
	      else ++iter;

	      std::string		type, signalType;
	      std::vector<Double>	position (3,0.0), direction (3,0.0), polarization (3,0.0);
	      Double                	amplitude, offset, variance, wavelength, cep;
	      Signal                	signal;
	      std::vector<Double>   	radius (2,0.0);

	      do
		{
		  if      (parameterName(*iter) == "type")                              type            = stringValue(*iter);
		  else if (parameterName(*iter) == "position")                          position        = vectorDoubleValue(*iter);
		  else if (parameterName(*iter) == "direction")                      	direction    	= vectorDoubleValue(*iter);
		  else if (parameterName(*iter) == "polarization")                      polarization    = vectorDoubleValue(*iter);
		  else if (parameterName(*iter) == "strength-parameter")                amplitude       = doubleValue(*iter);
		  else if (parameterName(*iter) == "radius-parallel")          		radius[0]       = doubleValue(*iter);
		  else if (parameterName(*iter) == "radius-perpendicular")     		radius[1]       = doubleValue(*iter);
		  else if (parameterName(*iter) == "signal-type")                       signalType      = stringValue(*iter);
		  else if (parameterName(*iter) == "offset")                            offset          = doubleValue(*iter);
		  else if (parameterName(*iter) == "variance")                          variance        = doubleValue(*iter);
		  else if (parameterName(*iter) == "wavelength")                        wavelength      = doubleValue(*iter);
		  else if (parameterName(*iter) == "CEP")                               cep             = doubleValue(*iter);
		  else { std::cout << parameterName(*iter) << " is not defined in seed-initialization group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");

	      if          (signalType.compare("neumann") == 0 || signalType.compare("gaussian") == 0 || signalType.compare("secant-hyperbolic") == 0 ||
		  signalType.compare("flat-top") == 0 )
		signal.initialize(signalType, offset, variance, wavelength, cep);
	      else { std::cout << signalType << " is an unknown signal type." << std::endl; exit(1); }

	      seed_.initialize(type, position, direction, polarization, amplitude, radius, signal);
	    }

	  /* The data related to the seed sampling are to be parsed.                            	*/
	  else if (*iter == "field-sampling")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The seed-sampling directory is empty" << std::endl; exit(1); }
	      else ++iter;

	      do
		{
		  if      (parameterName(*iter) == "sample")         	seed_.sampling_			= boolValue(*iter);
		  else if (parameterName(*iter) == "type")              seed_.samplingType_ 		= seed_.samplingType(stringValue(*iter));
		  else if (parameterName(*iter) == "field")             seed_.samplingField_.push_back(seed_.fieldType(stringValue(*iter)));
		  else if (parameterName(*iter) == "directory")         seed_.samplingDirectory_       	= stringValue(*iter);
		  else if (parameterName(*iter) == "base-name")         seed_.samplingBasename_		= stringValue(*iter);
		  else if (parameterName(*iter) == "rhythm")            seed_.samplingRhythm_          	= doubleValue(*iter);
		  else if (parameterName(*iter) == "position")
		    {
		      std::vector<Double> position = vectorDoubleValue(*iter);
		      FieldVector<Double> positionV; positionV = position;
		      seed_.samplingPosition_.push_back( positionV );
		    }
		  else if (parameterName(*iter) == "line-begin")
		    {
		      std::vector<Double> position = vectorDoubleValue(*iter);
		      seed_.samplingLineBegin_ =  position;
		    }
		  else if (parameterName(*iter) == "line-end")
		    {
		      std::vector<Double> position = vectorDoubleValue(*iter);
		      seed_.samplingLineEnd_ =  position;
		    }
		  else if (parameterName(*iter) == "resolution")   seed_.samplingRes_ = intValue(*iter);
		  else { std::cout << parameterName(*iter) << " is not defined in seed-sampling group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");

	    }

	  /* The data related to the seed visualization are to be parsed.                            	*/
	  else if (*iter == "field-visualization")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The seed-visualization directory is empty" << std::endl; exit(1); }
	      else ++iter;

	      /* Increase the size of the vtk vector in the ssed class by one;				*/
	      unsigned int i = seed_.vtk_.size(); seed_.vtk_.resize(i+1);

	      do
		{
		  if      (parameterName(*iter) == "sample")         	seed_.vtk_[i].sample_		= boolValue(*iter);
		  else if (parameterName(*iter) == "directory")         seed_.vtk_[i].directory_	= stringValue(*iter);
		  else if (parameterName(*iter) == "type")         	seed_.vtk_[i].type_		= seed_.vtkType(stringValue(*iter));
		  else if (parameterName(*iter) == "plane")         	seed_.vtk_[i].plane_ 		= seed_.planeType(stringValue(*iter));
		  else if (parameterName(*iter) == "base-name")         seed_.vtk_[i].basename_		= stringValue(*iter);
		  else if (parameterName(*iter) == "field")             seed_.vtk_[i].field_.push_back(seed_.fieldType(stringValue(*iter)));
		  else if (parameterName(*iter) == "rhythm")            seed_.vtk_[i].rhythm_   	= doubleValue(*iter);
		  else if (parameterName(*iter) == "position")
		    {
		      std::vector<Double> position = vectorDoubleValue(*iter);
		      seed_.vtk_[i].position_ = position;
		    }
		  else { std::cout << parameterName(*iter) << " is not defined in seed-visualization group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");
	    }

	  /* The data related to the seed visualization are to be parsed.                            	*/
	  else if (*iter == "field-profile")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The seed-profile directory is empty" << std::endl; exit(1); }
	      else ++iter;

	      do
		{
		  if      (parameterName(*iter) == "sample")         	seed_.profile_		= boolValue(*iter);
		  else if (parameterName(*iter) == "directory")         seed_.profileDirectory_	= stringValue(*iter);
		  else if (parameterName(*iter) == "base-name")         seed_.profileBasename_	= stringValue(*iter);
		  else if (parameterName(*iter) == "time")            	seed_.profileTime_.push_back(doubleValue(*iter));
		  else if (parameterName(*iter) == "field")             seed_.profileField_.push_back(seed_.fieldType(stringValue(*iter)));
		  else if (parameterName(*iter) == "rhythm")            seed_.profileRhythm_    = doubleValue(*iter);
		  else { std::cout << parameterName(*iter) << " is not defined in seed-profile group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");
	    }

	  ++iter;
	}
      while (*iter != "}");
    }

    /* Read the parameters parsed for the mesh in the solver.                                    	*/
    void readUndulator ( std::list<std::string>::iterator& iter )
    {
      ++iter;
      if (*iter != "{") { std::cout << "The undulator directory is empty" << std::endl; exit(1); }
      else ++iter;
      do
	{
	  /* The data related to the static undulator are to be parsed.                            	*/
	  if (*iter == "static-undulator")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The static undulator directory is empty" << std::endl; exit(1); }
	      else ++iter;

	      Undulator undulator; undulator.type_ = STATIC;

	      do
		{
		  if      (parameterName(*iter) == "undulator-parameter")               undulator.k_            = doubleValue(*iter);
		  else if (parameterName(*iter) == "period")                            undulator.lu_           = doubleValue(*iter);
		  else if (parameterName(*iter) == "polarization-angle")                undulator.theta_        = PI / 180.0 * doubleValue(*iter);
		  else if (parameterName(*iter) == "length")                            undulator.length_       = intValue(*iter);
		  else if (parameterName(*iter) == "offset")                            undulator.rb_       	= doubleValue(*iter);
		  else { std::cout << parameterName(*iter) << " is not defined in the static-undulator group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");

	      /* Add the given undulator to the undulator vector.					*/
	      undulator_.push_back(undulator);
	    }

	  /* The data related to the static undulator are to be parsed.                            	*/
	  if (*iter == "static-undulator-array")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The static undulator directory is empty" << std::endl; exit(1); }
	      else ++iter;

	      Double k, lu, theta, g, t;
	      unsigned int l, N;

	      do
		{
		  if      (parameterName(*iter) == "undulator-parameter")            k		= doubleValue(*iter);
		  else if (parameterName(*iter) == "period")                         lu		= doubleValue(*iter);
		  else if (parameterName(*iter) == "polarization-angle")             theta	= PI / 180.0 * doubleValue(*iter);
		  else if (parameterName(*iter) == "length")                         l   	= intValue(*iter);
		  else if (parameterName(*iter) == "gap")                            g		= doubleValue(*iter);
		  else if (parameterName(*iter) == "number")                         N		= intValue(*iter);
		  else if (parameterName(*iter) == "tapering-parameter")             t		= doubleValue(*iter);
		  else { std::cout << parameterName(*iter) << " is not defined in the static-undulator-array group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");

	      /* Now add each undulator module to the array of undulators.				*/
	      for (unsigned int i = 0; i < N; i++)
		{

		  /* First, calculate the values of the undulator module.				*/
		  Undulator undulator;
		  undulator.type_ 	= STATIC;
		  undulator.k_  	= k + i * t;
		  undulator.lu_ 	= lu;
		  undulator.theta_ 	= theta;
		  undulator.length_	= l;
		  undulator.rb_		= i * ( l * lu + g );

		  /* Now, add the undulator module ot the array of undulators.				*/
		  undulator_.push_back(undulator);
		}
	    }

	  /* The data related to the optical undulator are to be parsed.                            	*/
	  else if (*iter == "optical-undulator")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The optical undulator directory is empty" << std::endl; exit(1); }
	      else ++iter;

	      std::string		type, signalType;
	      std::vector<Double>	position (3,0.0), direction (3,0.0), polarization (3,0.0);
	      Double                	amplitude, offset, variance, wavelength, cep;
	      Signal                	signal;
	      std::vector<Double>   	radius (2,0.0);
	      Undulator 		undulator; undulator.type_ = OPTICAL;

	      do
		{
		  if 	  (parameterName(*iter) == "beam-type")                         type                    = stringValue(*iter);
		  else if (parameterName(*iter) == "position")                          position                = vectorDoubleValue(*iter);
		  else if (parameterName(*iter) == "direction")                         direction               = vectorDoubleValue(*iter);
		  else if (parameterName(*iter) == "polarization")                      polarization            = vectorDoubleValue(*iter);
		  else if (parameterName(*iter) == "strength-parameter")                amplitude               = doubleValue(*iter);
		  else if (parameterName(*iter) == "radius-parallel")          		radius[0]               = doubleValue(*iter);
		  else if (parameterName(*iter) == "radius-perpendicular")     		radius[1]               = doubleValue(*iter);
		  else if (parameterName(*iter) == "signal-type")                       signalType              = stringValue(*iter);
		  else if (parameterName(*iter) == "offset")                            offset                  = doubleValue(*iter);
		  else if (parameterName(*iter) == "variance")                          variance                = doubleValue(*iter);
		  else if (parameterName(*iter) == "wavelength")                        wavelength              = doubleValue(*iter);
		  else if (parameterName(*iter) == "CEP")                               cep                     = doubleValue(*iter);
		  else { std::cout << parameterName(*iter) << " is not defined in the optical-undulator group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");

	      if (signalType.compare("neumann") == 0 || signalType.compare("gaussian") == 0 || signalType.compare("secant-hyperbolic") == 0 ||
		  signalType.compare("flat-top") == 0 )
		signal.initialize(signalType, offset, variance, wavelength, cep);
	      else { std::cout << signalType << " is an unknown signal type." << std::endl; exit(1); }

	      undulator.initialize(type, position, direction, polarization, amplitude, radius, wavelength, signal);

	      undulator_.push_back(undulator);
	    }

	  ++iter;
	}
      while (*iter != "}");
    }

    /* Read the parameters parsed for the seed in the darius solver.                                    */
    void readExtField ( std::list<std::string>::iterator& iter )
    {
      ++iter;
      if (*iter != "{") { std::cout << "The EXTERNAL-FIELD directory is empty" << std::endl; exit(1); }
      else ++iter;

      ExtField          extField;

      do
	{
	  /* The data related to the seed initialization are to be parsed.                              */
	  if (*iter == "electromagnetic-wave")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The electromagnetic-field directory is empty" << std::endl; exit(1); }
	      else ++iter;

	      extField.type_            = EMWAVE;

	      std::string               type, signalType;
	      std::vector<Double>       position (3,0.0), direction (3,0.0), polarization (3,0.0);
	      Double                    amplitude, offset, variance, wavelength, cep;
	      Signal                    signal;
	      std::vector<Double>       radius (2,0.0);

	      do
		{
		  if      (parameterName(*iter) == "beam-type")                         type                    = stringValue(*iter);
		  else if (parameterName(*iter) == "position")                          position                = vectorDoubleValue(*iter);
		  else if (parameterName(*iter) == "direction")                         direction               = vectorDoubleValue(*iter);
		  else if (parameterName(*iter) == "polarization")                      polarization            = vectorDoubleValue(*iter);
		  else if (parameterName(*iter) == "strength-parameter")                amplitude               = doubleValue(*iter);
		  else if (parameterName(*iter) == "radius-parallel")          radius[0]               = doubleValue(*iter);
		  else if (parameterName(*iter) == "radius-perpendicular")     radius[1]               = doubleValue(*iter);
		  else if (parameterName(*iter) == "signal-type")                       signalType              = stringValue(*iter);
		  else if (parameterName(*iter) == "offset")                            offset                  = doubleValue(*iter);
		  else if (parameterName(*iter) == "variance")                          variance                = doubleValue(*iter);
		  else if (parameterName(*iter) == "wavelength")                        wavelength              = doubleValue(*iter);
		  else if (parameterName(*iter) == "CEP")                               cep                     = doubleValue(*iter);
		  else { std::cout << parameterName(*iter) << " is not defined in the electromagnetic external field group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");

	      if (signalType.compare("neumann") == 0 || signalType.compare("gaussian") == 0 || signalType.compare("secant-hyperbolic") == 0 ||
		  signalType.compare("flat-top") == 0 )
		signal.initialize(signalType, offset, variance, wavelength, cep);
	      else { std::cout << signalType << " is an unknown signal type." << std::endl; exit(1); }

	      extField.initialize(type, position, direction, polarization, amplitude, radius, wavelength, signal);

	      extField_.push_back(extField);
	    }

	  ++iter;
	}
      while (*iter != "}");
    }

    /* Read the parameters parsed for the FEL output in the darius solver.                              */
    void readFEL ( std::list<std::string>::iterator& iter )
    {
      ++iter;
      if (*iter != "{") { std::cout << "The FEL-OUTPUT directory is empty" << std::endl; exit(1); }
      else ++iter;
      do
	{
	  /* The data related to the seed initialization are to be parsed.                              */
	  if (*iter == "radiation-power")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The radiation-power directory is empty" << std::endl; exit(1); }
	      else ++iter;

	      FreeElectronLaser		FEL;

	      do
		{
		  if      (parameterName(*iter) == "sample")         		        FEL.radiationPower_.sampling_	= boolValue(*iter);
		  else if (parameterName(*iter) == "plane-position")	        	FEL.radiationPower_.z_.push_back( doubleValue(*iter) );
		  else if (parameterName(*iter) == "directory")         	        FEL.radiationPower_.directory_ 	= stringValue(*iter);
		  else if (parameterName(*iter) == "base-name")         	        FEL.radiationPower_.basename_	= stringValue(*iter);
		  else if (parameterName(*iter) == "line-begin")		        FEL.radiationPower_.lineBegin_ 	= doubleValue(*iter);
		  else if (parameterName(*iter) == "line-end")			        FEL.radiationPower_.lineEnd_ 	= doubleValue(*iter);
		  else if (parameterName(*iter) == "resolution")		        FEL.radiationPower_.res_ 	= doubleValue(*iter);
		  else if (parameterName(*iter) == "normalized-frequency")	        FEL.radiationPower_.lambda_.push_back(doubleValue(*iter));
		  else if (parameterName(*iter) == "minimum-normalized-frequency")	FEL.radiationPower_.lambdaMin_	= doubleValue(*iter);
		  else if (parameterName(*iter) == "maximum-normalized-frequency")	FEL.radiationPower_.lambdaMax_	= doubleValue(*iter);
		  else if (parameterName(*iter) == "normalized-frequency-resolution")	FEL.radiationPower_.lambdaRes_	= doubleValue(*iter);
		  else if (parameterName(*iter) == "type")              	        FEL.radiationPower_.samplingType(stringValue(*iter));
		  else { std::cout << parameterName(*iter) << " is not defined in radiation-power group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");

	      FEL_.push_back(FEL);
	    }

	  /* The data related to the power visualization are to be parsed.                            	*/
	  if (*iter == "power-visualization")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The power-visualization directory is empty" << std::endl; exit(1); }
	      else ++iter;

	      FreeElectronLaser		FEL;

	      do
		{
		  if      (parameterName(*iter) == "sample")         			FEL.vtk_.sampling_	= boolValue(*iter);
		  else if (parameterName(*iter) == "directory")         		FEL.vtk_.directory_	= stringValue(*iter);
		  else if (parameterName(*iter) == "base-name")         		FEL.vtk_.basename_	= stringValue(*iter);
		  else if (parameterName(*iter) == "plane-position")	        	FEL.vtk_.z_		= doubleValue(*iter);
		  else if (parameterName(*iter) == "rhythm")            		FEL.vtk_.rhythm_   	= doubleValue(*iter);
		  else if (parameterName(*iter) == "normalized-frequency")	        FEL.vtk_.lambda_        = doubleValue(*iter);
		  else { std::cout << parameterName(*iter) << " is not defined in power-visualization group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");

	      FEL_.push_back(FEL);
	    }

	  /* The data related to the seed initialization are to be parsed.                            	*/
	  if (*iter == "radiation-energy")
	    {
	      ++iter;
	      if (*iter != "{") { std::cout << "The radiation-energy directory is empty" << std::endl; exit(1); }
	      else ++iter;

	      FreeElectronLaser		FEL;

	      do
		{
		  if      (parameterName(*iter) == "sample")                            FEL.radiationPower_.sampling_   = boolValue(*iter);
		  else if (parameterName(*iter) == "distance-from-bunch")               FEL.radiationPower_.z_.push_back( doubleValue(*iter) );
		  else if (parameterName(*iter) == "directory")                         FEL.radiationPower_.directory_  = stringValue(*iter);
		  else if (parameterName(*iter) == "base-name")                         FEL.radiationPower_.basename_   = stringValue(*iter);
		  else if (parameterName(*iter) == "line-begin")                        FEL.radiationPower_.lineBegin_  = doubleValue(*iter);
		  else if (parameterName(*iter) == "line-end")                          FEL.radiationPower_.lineEnd_    = doubleValue(*iter);
		  else if (parameterName(*iter) == "resolution")                        FEL.radiationPower_.res_        = doubleValue(*iter);
		  else if (parameterName(*iter) == "normalized-wavelength")             FEL.radiationPower_.lambda_.push_back(doubleValue(*iter));
		  else if (parameterName(*iter) == "minimum-normalized-wavelength")     FEL.radiationPower_.lambdaMin_  = doubleValue(*iter);
		  else if (parameterName(*iter) == "maximum-normalized-wavelength")     FEL.radiationPower_.lambdaMax_  = doubleValue(*iter);
		  else if (parameterName(*iter) == "normalized-wavelength-resolution")  FEL.radiationPower_.lambdaRes_  = doubleValue(*iter);
		  else if (parameterName(*iter) == "type")                              FEL.radiationPower_.samplingType(stringValue(*iter));
		  else { std::cout << parameterName(*iter) << " is not defined in radiation-power group." << std::endl; exit(1); }
		  ++iter;
		}
	      while (*iter != "}");

	      FEL_.push_back(FEL);
	    }

	  ++iter;
	}
      while (*iter != "}");
    }

  };    /* End of parseDarius class.                                                                	*/

}       /* End of namespace Cyrus3d.                                                          		*/
#endif
