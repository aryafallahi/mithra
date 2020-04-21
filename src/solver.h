/********************************************************************************************************
 *  solver.h : Implementation of the solver class for the darius code
 ********************************************************************************************************/

#ifndef SOLVER_H_
#define SOLVER_H_

#include <iomanip>
#include <list>
#include <vector>

#include "classes.h"
#include "database.h"
#include "fieldvector.h"

namespace MITHRA
{

  /* Class of functions used for the solution of the fields in time domain using FDTD.			*/
  class Solver
  {

  public:

    Solver( Mesh& 				mesh,
	    Bunch& 				bunch,
	    Seed& 				seed,
	    std::vector<Undulator>&		undulator,
	    std::vector<ExtField>& 		extField,
	    std::vector<FreeElectronLaser>& 	FEL);

    /* Using the parsed data, set the required parameters for simulation.				*/
    void 		setSimulationParameters 	();

    /* Boost the mesh into the electron rest frame.							*/
    void 		lorentzBoostMesh		();

    /* Boost particles into the electron rest frame.							*/
    void 		lorentzBoostBunch 		();

    /* Distribute particles in their respective processor, depending on their longituinal coordinate.	*/
    void 		distributeParticles 		(std::list<Charge>& chargeVector);

    /* Get the average gamma and average direction of a bunch read in from a file.					*/
    void 		computeFileGamma 		(BunchInitialize & bunchInit);

    /* Initialize the matrix for the field values and the coordinates.					*/
    void 		initialize			();

    /* Initialize the temporal and spatial mesh of the problem.						*/
    void 		initializeMesh			();

    /* Initialize the data for updating the field in the fdtd algorithm.				*/
    void 		initializeField			();

    /* Initialize the data required for sampling seed or the total field in the computational domain.	*/
    void 		initializeSeedSampling		();

    /* Initialize the required data for visualizing and saving the field.				*/
    void 		initializeSeedVTK		();

    /* Initialize the required data for profiling and saving the field.					*/
    void 		initializeSeedProfile		();

    /* Initialize the required data for updating the bunch.						*/
    void 		initializeBunchUpdate		();

    /* Initialize the charge vector containing the bunch given by the user.				*/
    void 		initializeBunch			();

    /* Update the fields for one time-step.								*/
    void 		bunchUpdate			();

    /* Sample the bunch data and save it to the given file.						*/
    void 		bunchSample			();

    /* Visualize the bunch as vtk files and save them to the file with given name.			*/
    void 		bunchVisualize			();

    /* Write the total profile of the field into the given file name.					*/
    void 		bunchProfile			();

    /* Calculate the magnetic field of undulator and add it to the magnetic field of the seed.		*/
    void 		undulatorField 			(UpdateBunchParallel& ubp,  FieldVector<Double>& r);

    /* Calculate the field of external field and add it to the field of the seed.			*/
    void 		externalField			(UpdateBunchParallel& ubp,  FieldVector<Double>& r);

    /* Initialize the data required for sampling and saving the radiation power at the given position.	*/
    void 		initializePowerSample		();

    /* Sample the radiation power at the given position and save it to the file.			*/
    void 		powerSample			();

    /* Initialize the data required for visualizing the radiation power at the given position.		*/
    void 		initializePowerVisualize	();

    /* Visualize the radiation power at the given position and save it to the file.			*/
    void 		powerVisualize			();

    /* Initialize the data required for sampling and saving the radiation energy at the given position.	*/
    void 		initializeEnergySample		();

    /* Sample the radiation energy at the given position and save it to the file. 			*/
    void 		energySample			();

    /* Initialize the data required for storing particles hitting a screen at the given position.	*/
    void 		initializeScreenProfile		();

    /* Store the bunch profile from particles hitting a screen at the given position and save it to the
     * file. 												*/
    void 		screenProfile			();

    /* Finalize the field calculations. 								*/
    void 		finalize			();

    /* Define the virtual function for field evaluation.						*/
    virtual void 	fieldEvaluate 			(long int m) = 0;

    /* Define the boolean function for comparing undulator begins.					*/
    static bool 	undulatorCompare 		(Undulator i, Undulator j);

    /* Define the function for linear interpolation.							*/
    Double	 	interp				(Double x0, Double x1, Double y0, Double y1, Double x);


    /****************************************************************************************************
     * List of required parameters in the FdTd code.
     ****************************************************************************************************/

    /* The parsed parameters for the simulation are written in four classes: mesh, bunch, seed, and
     * undulator.											*/
    Mesh& 								mesh_;
    Bunch&								bunch_;
    Seed&								seed_;
    std::vector<Undulator>&						undulator_;
    std::vector<ExtField>&                                              extField_;
    std::vector<FreeElectronLaser>&					FEL_;

    /* The vector potential at the nodes in the computational mesh at three different time points.	*/
    std::vector<FieldVector<Double> >* 					anp1_;
    std::vector<FieldVector<Double> >* 					an_;
    std::vector<FieldVector<Double> >* 					anm1_;

    /* The static potential at the nodes in the computational mesh at three different time points.	*/
    std::vector<Double>* 						fnp1_;
    std::vector<Double>* 						fn_;
    std::vector<Double>* 						fnm1_;

    /* Boolean vector determining the inclusion of particles.                                           */
    std::vector<bool>                                                   pic_;

    /* The vector potential at the nodes in the computational mesh at three different time points.	*/
    std::vector<FieldVector<Double> > 					en_;
    std::vector<FieldVector<Double> > 					bn_;

    /* The current density at the nodes in the computational mesh at three different time points.	*/
    std::vector<FieldVector<Double> > 					jn_;

    /* The charge density at the nodes in the computational mesh at three different time points.	*/
    std::vector<Double> 						rn_;

    /* The coordinate of the nodes in the computational mesh.						*/
    std::vector<FieldVector<Double> > 					r_;

    /* Number of nodes in each direction.								*/
    int									N0_, N1_, N2_, N1N0_;

    /* total number of charges in the simulation.							*/
    unsigned int							Nc_;

    /* Number of nodes along z in the specific processor and the index of the first column.		*/
    int									np_, k0_;

    /* z coordinates of the critical points in the mesh partitioning.					*/
    Double								zp_ [2];

    /* borders of the computational mesh.                                                               */
    Double                                                              xmin_, xmax_;
    Double                                                              ymin_, ymax_;
    Double                                                              zmin_, zmax_;

    /* Time and the time step number for the field calculations.					*/
    Double								timep1_;
    Double								time_;
    Double								timem1_;
    unsigned int 							nTime_;

    /* vector of charge structures containing the position and momentum of the charge points.		*/
    std::list<Charge>							chargeVectorn_;

    /* Time and the time step number for the bunch calculations.					*/
    Double								timeBunch_;
    unsigned int 							nTimeBunch_;

    /* Number of bunch updates within each field update.						*/
    Double 								nUpdateBunch_;

    /* The gamma, beta and dt factor for the moving frame derived from the first undulator parameter.	*/
    Double								gamma_;
    Double								beta_;
    Double								dt_;

    /* Define a structure containing the parameters needed to update the values. These parameters are
     * defined once in the class to avoid declaring them every time a field is updated.			*/
    UpdateField								uf_;

    /* Define a structure containing the parameters needed to sample the fields the values. These
     * parameters are defined once in the class to avoid declaring them every time a field is updated.	*/
    SampleField								sf_;

    /* Define a structure containing the parameters needed to visualize the field profile. These
     * parameters are defined once in the class to avoid declaring them every time a field is updated.	*/
    std::vector<VisualizeField>						vf_;

    /* Define a structure containing the parameters needed to save the field profile. These
     * parameters are defined once in the class to avoid declaring them every time a field is updated.	*/
    ProfileField							pf_;

    /* Define a structure containing the parameters needed to update the bunch distribution. These
     * parameters are defined once in the class to avoid declaring them every time a field is updated.	*/
    UpdateBunch						                ub_;

    /* Define a structure containing the parameters needed to sample the bunch values. These parameters
     * are defined once in the class to avoid declaring them every time a field is updated.		*/
    SampleBunch								sb_;

    /* Define a structure containing the parameters needed to visualize the bunch distribution. These
     * parameters are defined once in the class to avoid declaring them every time a field is updated.	*/
    VisualizeBunch							vb_;

    /* Define a structure containing the parameters needed to save the bunch profile. These parameters
     * are defined once in the class to avoid declaring them every time a field is updated.		*/
    ProfileBunch							pb_;

    /* Define a structure containing the parameters needed to update the current. These parameters
     * are defined once in the class to avoid declaring them every time a field is updated.		*/
    UpdateCurrent						        uc_;

    /* Define a structure containing the parameters needed to sample the FEL radiation power the values.
     * These parameters are defined once in the class to avoid declaring them every time a field is
     * updated.												*/
    std::vector<SampleRadiationPower>				        rp_;

    /* Define a structure containing the parameters needed to sample the FEL radiation power the values.
     * These parameters are defined once in the class to avoid declaring them every time a field is
     * updated.												*/
    std::vector<SampleRadiationEnergy>				        re_;

    /* Define a structure containing the parameters needed to store the particles hitting screens
     * These parameters are defined once in the class to avoid declaring them every time a bunch is
     * updated.												*/
    std::vector<SampleScreenProfile>				        scrp_;

    /* Define the value of MPI variables.								*/
    int									rank_, size_;

    /* Speed of light value in terms of the given length-scale and time-scale.				*/
    Double								c0_, m0_, e0_;
  };

}

#endif
