/********************************************************************************************************
 *  classes.hh : Implementation of the classes used in mithra
 ********************************************************************************************************/

#ifndef CLASSES_H_
#define CLASSES_H_

#include <list>

#include "database.h"
#include "fieldvector.h"
#include "stdinclude.h"

namespace MITHRA
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
    void show ();

    /* Initialize the parameters in the mesh initializer.						*/
    void initialize ();
  };

  /* Bunch class includes the data for the bunch properties and the functions for initializing the
   * bunches and also evaluating the bunch properties.                                              	*/
  class Bunch
  {

  public:
    typedef std::list <Charge> ChargeVector;

    /* The constructor clears and initializes the internal data structure.				*/
    Bunch ();

    /* Initialize a bunch with a manual type. This bunch produces one charge equal to the cloudCharge_. */
    void initializeManual (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia);

    /* Initialize a bunch with an ellipsoid type. This bunch produces a number of charges equal to the
     * numberOfParticles_ with the total charge equal to the cloudCharge_ which are distributed in an
     * ellipsoid with dimensions given by sigmaPosition_ and center given by the position vector. The
     * particles have uniform energy distribution centered at initialEnergy_ with variances determined by
     * sigmaGammaBeta_.                                                         			*/
    void initializeEllipsoid (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia);

    /* Initialize a bunch with a 3D-crystal type. This bunch produces a number of charges equal
     * to the numberOfParticles_ with the total charge equal to the cloudCharge_ which are arranged in a
     * 3D crystal. The number of particles in each direction is given by numbers_. Therefore,
     * numberOfParticles_ should be a multiple of the product of all these three numbers. The ratio gives
     * the number of particles in each crystal point. Position of each particle is determined by the
     * lattice constants and the crystal is centered at the position_ vector. At each point, the charges
     * have a small Gaussian distribution around the crsytal.                                           */
    void initialize3DCrystal (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia);

    /* Initialize a bunch with a file type. This bunch produces a number of charges read from a given file.
     * The number of initialized charge is equal to the vertical length of the table in the text file. The
     * file format should contain the charge value, 3 position coordinates and 3 momentum coordinates of
     * of the charge distribution.							                */
    void initializeFile (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia);

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

    /* Position of the undulator begin at the instance of bunch initialization.				*/
    Double				zu_;

    /* Beta of the moving frame in the stationary lab frame.						*/
    Double				beta_;

    /* Show the stored values for the bunch.                                                          	*/
    void show ();
  };

  /* Define the main signal class.                                                                      */
  class Signal
  {
  public:

    /* Initialize the values of the parameters.								*/
    Signal ();

    /* Initializer with signal type, time offset, variance, frequency and carrier-envelope-phase.       */
    void initialize (std::string type, Double l0, Double s, Double l, Double cep);

  public:

    /* Store type of the signal.                                                                        */
    SignalType				signalType_;

    /* Time offset of the signal.                                                                       */
    Double     				t0_;

    /* Variance of the signal. Variance is defined according to the point where the intensity of the
     * signal, i.e. signal squared is half of the maximum.                                              */
    Double     				s_;

    /* Frequency of the modulation.                                                                     */
    Double     				f0_;

    /* Carrier envelope phase of the modulation.                                                        */
    Double     				cep_;

    /* Provide the signal at time t.                                                                    */
    Double self (Double& t, Double& phase);

    /* Show the stored values for this signal.                                                          */
    void show ();
  };

  /* Define the main seed class.                                                                  	*/
  class Seed
  {
  public:
    Seed ();

    void initialize (std::string        	type,
		     std::vector<Double>    	position,
		     std::vector<Double>    	direction,
		     std::vector<Double>    	polarization,
		     Double                 	amplitude,
		     std::vector<Double>	radius,
		     Signal                   	signal);

    /* Store data required for visualizing the radiated field in all-domain.                            */
    struct vtk
    {
      bool sample_;
      std::vector <FieldType> field_;
      std::string directory_;
      SamplingType type_;
      std::string basename_;
      Double rhythm_;
      PlaneType plane_;
      FieldVector <Double> position_;

      /* Initialize the data-base for field visualization.						*/
      vtk ();
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
    void fields (FieldVector <Double> & aufpunkt, Double & time, FieldVector <Double> & a);

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

    std::vector <vtk> vtk_;

    /* Store data required for writing the profile of the field in all-domain.                          */
    bool                                profile_;
    std::vector<FieldType>              profileField_;
    std::string                         profileDirectory_;
    std::string                         profileBasename_;
    std::vector<Double>                 profileTime_;
    Double                        	profileRhythm_;

    /* Set the sampling type of the seed.                                                               */
    SamplingType 	samplingType 	(std::string samplingType);

    /* Set the sampling type of the seed.                                                               */
    SamplingType 	vtkType 	(std::string vtkType);

    /* Set the plane type for vtk in plane visualization.                                           	*/
    PlaneType 		planeType 	(std::string planeType);

    /* Set the field sampling type of the seed.                                                         */
    FieldType 		fieldType 	(std::string fieldType);

    /* Show the stored values for this signal.                                                          */
    void show ();
  };

  /* Define the structure containing the main parameters for the undulator.				*/
  class Undulator
  {
  public:
    Undulator ();

    /* The magnetic field of the undulator.                                                       	*/
    Double				k_;

    /* The period of the undulator.                                                       		*/
    Double				lu_;

    /* The start position of the undulator.								*/
    Double				rb_;

    /* The length of the undulator.									*/
    unsigned int			length_;

    /* The initial distance between the bunch head and the undulator begin.						*/
    Double				dist_;
    
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
    UndulatorType undulatorType (std::string undulatorType);

    /* Initialize the data of the undulator according to the input parameters.				*/
    void initialize (std::string        	type,
		     std::vector<Double>    	position,
		     std::vector<Double>    	direction,
		     std::vector<Double>    	polarization,
		     Double                 	amplitude,
		     std::vector<Double>	radius,
		     Double			wavelength,
		     Signal                   	signal);

    /* Show the stored values for the undulator.                                                        */
    void show ();
  };

  class ExtField
  {
  public:

    ExtField();

    /* Type of the undulator, it can be an optical or a static undulator.                   		*/
    ExtFieldType                        type_;

    /* Store type of the seed.                                                                    	*/
    SeedType                            seedType_;

    /* Speed of light value in terms of the given length-scale and time-scale.                      	*/
    Double                              c0_;

    /* Store reference position of the seed. For waves, it is a reference position and for hertzian
     * dipole, it is the position of the dipole.                                            		*/
    FieldVector<Double>                 position_;

    /* Store the direction of the seed. For waves, it is the propagation direction and for hertzian
     * dipole, it is the direction of the dipole.                                              		*/
    FieldVector<Double>                 direction_;

    /* Store the polarization of the wave in the seed.                                          	*/
    FieldVector<Double>                 polarization_;

    /* Store the amplitude of the seed.                                                            	*/
    Double                              amplitude_;

    /* Store the Rayleigh radius of the Gaussian beam in the parallel and perpendicular directions.  	*/
    std::vector<Double>                 radius_;

    /* Store the signal class for this seed.                                                         	*/
    Signal                              signal_;

    /* Initialize the data of the undulator according to the input parameters.				*/
    void initialize (std::string                type,
		     std::vector<Double>        position,
		     std::vector<Double>        direction,
		     std::vector<Double>        polarization,
		     Double                     amplitude,
		     std::vector<Double>        radius,
		     Double                     wavelength,
		     Signal                     signal);

    /* Show the stored values for the external field.							*/
    void show ();
  };

  struct FreeElectronLaser
  {
    /* Structure contianing the parsed parameters for the radiation power.				*/
    struct RadiationSampling
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
      void samplingType(std::string samplingType);

      /* Initialize the values for initializing the radiation power.                     		*/
      RadiationSampling();
    };

    /* Structure containing the parsed parameters for the power or energy visualization.		*/
    struct RadiationVisualization
    {
      /* The distance of the plane from the bunch to sample the radiation power or energy.		*/
      Double					z_;

      /* Store the data required for sampling the radiation power or energy.				*/
      bool					sampling_;

      /* Store the directory in which the radiation data is saved.					*/
      std::string				directory_;

      /* Store the base-name of the file in which the radiation data is saved.				*/
      std::string				basename_;

      /* Store the rhythm for calculating and saving the radiation power or energy.			*/
      Double					rhythm_;

      /* The wavelength of the harmonic whose power or energy should be plotted.			*/
      Double					lambda_;

      /* Initialize the values for initializing the radiation power or energy.                     	*/
      RadiationVisualization ();
    };

    /* Define the variables containing the required structures for power sampling and visualization.	*/
    RadiationSampling 				radiationPower_;
    RadiationVisualization 			vtkPower_;

    /* Define the variables containing the required structures for energy sampling and visualization.	*/
    RadiationSampling 				radiationEnergy_;
    RadiationVisualization 			vtkEnergy_;
    
    /* Parsed parameters for the screens that record the bunch profile. This produces the bunch profile
     * in the lab frame.										*/
    struct ScreenProfile
    {
      /* Flag that activates storing the data required for sampling the bunch on a screen.		*/
      bool					sampling_;

      /* Store the directory in which the bunch profile in the lab frame is saved.			*/
      std::string				directory_;

      /* Store the base-name of the file in which this bunch profile is saved.				*/
      std::string				basename_;

      /* Store the positions in the undulator where the saving of the bunch profile is done.		*/
      std::vector<Double>      			pos_;

      /* Store the rhythm in position for saving the bunch profile.					*/
      Double					rhythm_;

      ScreenProfile ();
    };

    /* Define the variables containing the required structure for profiling the bunch in the lab frame.	*/
    ScreenProfile 				screenProfile_;
  };
}
#endif
