/********************************************************************************************************
 *  fdtd.hh : Implementation of the real fdtd time marching solution class for the mithra code
 ********************************************************************************************************/

#ifndef FDTD_HH_
#define FDTD_HH_

#include "classes.h"
#include "solver.h"

namespace Mithra
{

  /* Class of functions used for the solution of the fields in time domain using FDTD.			*/
  class FdTd : public Solver
  {
  public:
    FdTd (Mesh& 				mesh,
	  Bunch& 				bunch,
	  Seed& 				seed,
	  std::vector<Undulator>&		undulator,
	  std::vector<ExtField>& 		extField,
	  std::vector<FreeElectronLaser>& 	FEL);

    /* The function which is called for solving the fields in time domain.				*/
    void solve ();

    /* Reset the currents to zero.									*/
    void currentReset ();

    /* Update the currents at cell points for the filed update.						*/
    void currentUpdate ();

    /* Communicate the currents among different processors.						*/
    void currentCommunicate ();

    /* Update the fields for one time-step								*/
    void fieldUpdate ();

    /* Evaluate the field of the m'th pixel from the potentials.					*/
    void fieldEvaluate (long int m);

    /* Sample the field and save it to the given file.							*/
    void fieldSample ();

    /* Visualize the field as vtk files on the whole domain and save them to the file with given name.	*/
    void fieldVisualizeAllDomain 	(unsigned int ivtk);

    /* Visualize the field as vtk files in plane and save them to the file with the given name.		*/
    void fieldVisualizeInPlane 		(unsigned int ivtk);

    /* Visualize the field as vtk files in a plane normal to x axis and save them to the file with the
     * given name.											*/
    void fieldVisualizeInPlaneXNormal 	(unsigned int ivtk);

    /* Visualize the field as vtk files in a plane normal to y axis and save them to the file with the
     * given name.											*/
    void fieldVisualizeInPlaneYNormal 	(unsigned int ivtk);

    /* Visualize the field as vtk files in a plane normal to z axis and save them to the file with the
     * given name.											*/
    void fieldVisualizeInPlaneZNormal 	(unsigned int ivtk);

    /* Write the total profile of the field into the given file name.					*/
    void fieldProfile ();
  };
}
#endif
