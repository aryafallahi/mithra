/********************************************************************************************************
 Darius-1.0 (Completely Numerical Calculation of Free Electron Laser Radiation)
 Version 1.0, copyright 2014, A. Fallahi
 ********************************************************************************************************
 darius.cpp: Main program file
 MPI version should be compiled in the mithra folder using:
 mpiCC -O3 source/darius.cpp -o project/darius
 to commit to git use
 git add --all - git commit -m "message" - git push
 on Mac:
 mpiCC -lstdc++ -O3 source/darius.cpp -o project/darius
 ********************************************************************************************************/

/* standard C++ header files 										*/
#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <map>
#include <sys/time.h>
#include <cmath>
#include <complex>
#include <algorithm>
#include <utility>
#include <sstream>
#include <list>
#include <limits>
#include <mpi.h>
#include <algorithm>
#include <sys/stat.h>

/** Include darius header files.                                                                        */
#include "fieldvector.hh"
#include "stdinclude.hh"
#include "readdata.hh"
#include "database.hh"
#include "classes.hh"
#include "datainput.hh"
#include "readdata.hh"
#include "solver.hh"
#include "fdtd.hh"
#include "fdtdSC.hh"
#include "radiation.hh"

int main (int argc, char* argv[])
{

  /* initialize MPI, finalize is done automatically on exit                                             */
  MPI_Init(&argc,&argv);

  /* Activate namespaces                                                                                */
  using namespace Darius;
  using namespace std;

  /* Retrieve time when we start the simulation.                                                        */
  timeval simulationStart, simulationEnd;
  gettimeofday(&simulationStart, NULL);

  /* Hello message                                                                                      */
  printmessage(std::string(__FILE__), __LINE__, std::string(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::") );
  printmessage(std::string(__FILE__), __LINE__, std::string("MITHRA-2.0: Completely Numerical Calculation of Free Electron Laser Radiation)") );
  printmessage(std::string(__FILE__), __LINE__, std::string("Version 2.0, Copyright 2019, Arya Fallahi") );
  printmessage(std::string(__FILE__), __LINE__, std::string("Written by Arya Fallahi, IT'IS Foundation, Zurich, Switzerland") );
  printmessage(std::string(__FILE__), __LINE__, std::string(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::") );

  /* Parse the command line options                                                                     */
  std::list<std::string> jobFile = read_file(argv[1]);
  cleanJobFile(jobFile);

  /* Create the solver database.                                                                        */
  Mesh                                  mesh;
  mesh.initialize();

  /* Create the bunch database.                                                                         */
  Bunch                                 bunch;

  /* Create the seed database.                                                                          */
  Seed                                  seed;

  /* Create the undulator database.                                                                     */
  std::vector<Undulator>                undulator;
  undulator.clear();

  /* Create the external field database.                                                                */
  std::vector<ExtField>                 extField;
  extField.clear();

  /* Create the free electron laser database.                                                           */
  std::vector<FreeElectronLaser>        FEL;
  FEL.clear();

  /* Open input parameter parser and instantiate the databases.                                         */
  ParseDarius parser (jobFile, mesh, bunch, seed, undulator, extField, FEL);
  parser.setJobParameters();

  /* Show the parameters for the simulation.                                                            */
  mesh.show();
  bunch.show();
  seed.show();
  for (unsigned int i = 0; i < undulator.size(); i++) 	undulator[i].show();
  for (unsigned int i = 0; i < extField.size();  i++) 	extField[i] .show();

  /* Initialize the class for the FDTD computations.                                                    */
  FdTd   fdtd   (mesh, bunch, seed, undulator, extField, FEL);
  FdTdSC fdtdsc (mesh, bunch, seed, undulator, extField, FEL);

  /* Solve for the fields and the bunch distribution over the specified time.                           */
  if ( mesh.spaceCharge_ )
    fdtdsc.solve();
  else
    fdtd.solve();

  /* Calculate the total simulation time.                                                               */
  gettimeofday(&simulationEnd, NULL);
  Double deltaTime = ( simulationEnd.tv_usec - simulationStart.tv_usec ) / 1.0e6;
  deltaTime += ( simulationEnd.tv_sec - simulationStart.tv_sec );
  printmessage(std::string(__FILE__), __LINE__, std::string("::: total simulation time [seconds] = ") + stringify(deltaTime) );

  MPI_Finalize();
}
