# mithra
MITHRA is a full-wave numerical solver for free-electron lasers.

Compiling options used for gcc 7.4.0:
- `make` or `make full` to compile binary
- `make install` to compile library and header directory

********************************************************************************************************
MITHRA-2.0 (Completely Numerical Calculation of Free Electron Laser Radiation)
Version 2.0, copyright 2019, Arya Fallahi
********************************************************************************************************
 
mithra.cpp: Main program file

stdinclude: Set of different standard functions and constants used in the code.

readdata: Implementation of the functions reading the lines of the job file.

radiation: Implementation of the functions for the radiation analysis.

fieldvector: Implementation of the field vector class for the analysis.

datainput: Implementation of the parameter parser

classes: Implementation of the classes used in the darius code

fdtd: Implementation of the real fdtd time marching solution class for the MITHRA code without space-charge

fdtdSC: Implementation of the real fdtd time marching solution class for the MITHRA code with space-charge

database: Implementation of various data structures used for simulations

solver: Mother class for solving the free-electron laser problem

********************************************************************************************************
