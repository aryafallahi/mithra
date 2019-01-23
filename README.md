# mithra
MITHRA is a full-wave numerical solver for free -electron lasers.

********************************************************************************************************
MITHRA-2.0 (Completely Numerical Calculation of Free Electron Laser Radiation)
Version 2.0, copyright 2019, Arya Fallahi
********************************************************************************************************
 
darius.cpp: Main program file

stdinclude.hh: Set of different standard functions and constants used in the code.

readdata.hh: Implementation of the functions reading the lines of the job file.

fieldvector.hh: Implementation of the field vector class for the analysis.

datainput.hh: Implementation of the parameter parser

classes.hh: Implementation of the classes used in the darius code

fdtd.hh: Implementation of the real fdtd time marching solution class for the MITHRA code without space-charge

fdtdSC.hh: Implementation of the real fdtd time marching solution class for the MITHRA code with space-charge

database.hh: Implementation of various data structures used for simulations

solver.hh: Mother class for solving the free-electron laser problem
  
Compiling options used for gcc 4.7.3:

MPI version should be compiled using: " mpiCC -O3 darius.cpp -o ../darius "

********************************************************************************************************
