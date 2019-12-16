/********************************************************************************************************
 *  stdinclude.hh, set of different standard functions and constants used in the code.
 ********************************************************************************************************/


#ifndef STDINCLUDE_HH_
#define STDINCLUDE_HH_

namespace Darius
{
  /* Define the types of the time domain signals.                                                       */
  enum SignalType       	{NEUMANN, GAUSSIAN, SECANT, FLATTOP};

  /* Define the types of the excitations.                                                               */
  enum SeedType   		{PLANEWAVE, PLANEWAVECONFINED, GAUSSIANBEAM};

  /* Define the type of the external field.                                                             */
  enum ExtFieldType             {EMWAVE};

  /* Define the types of the excitations.                                                               */
  enum SamplingType     	{ATPOINT, OVERLINE, INPLANE, ALLDOMAIN};

  /* Define the type of the plane normal to the axis.                                                  	*/
  enum PlaneType     		{XNORMAL, YNORMAL, ZNORMAL};

  /* Define the types of the fields to be sampled.                                                      */
  enum FieldType		{Ex, Ey, Ez, Bx, By, Bz, Ax, Ay, Az, Jx, Jy, Jz, Q, F};

  /* Define the types of the undulator supported by the code.                                           */
  enum UndulatorType		{STATIC, OPTICAL};

  /* Define the type of the solver to be used for the FEL interaction.					*/
  enum SolverType		{FD, NSFD};

  /* Scientific constants.                                                                              */
  const Double PI           	= 3.1415926535;
  const Double EPSILON_ZERO  	= 8.85418782e-12;
  const Double MU_ZERO        	= 4.0 * PI * 1.0e-7;
  const Double C0             	= 1.0 / sqrt(EPSILON_ZERO * MU_ZERO);
  const Double Z0             	= sqrt(MU_ZERO / EPSILON_ZERO);
  const Double EC             	= 1.602e-19;
  const Double EM             	= 9.109e-31;
  const Double HB              	= 1.054e-34;
  const Double KB             	= 1.381e-23;

  /* File suffixes.                                                                                     */
  const std::string VTS_FILE_SUFFIX = ".vts";
  const std::string PTS_FILE_SUFFIX = ".pvts";
  const std::string PTU_FILE_SUFFIX = ".pvtu";
  const std::string TXT_FILE_SUFFIX = ".txt";
  const std::string VTU_FILE_SUFFIX = ".vtu";

  /* Unit imaginary number.										*/
  const Complex I = Complex (0.0,1.0);

  /* Check if a path is an absolute path (beginning with an '/').                                       */
  inline bool isabsolute(std::string filename)
  {
    return (filename.compare(0,1,"/") == 0);
  }

  /* Check if a directory to save the data exists.							*/
  bool pathExist(const std::string &s)
  {
    struct stat buffer;
    return (stat (s.c_str(), &buffer) == 0);
  }

  /* Split the file name to two strings including its path and file name.                               */
  void splitFilename (const std::string& str, std::string& path, std::string& file)
  {
    unsigned found = str.find_last_of("/");
    path = str.substr(0,found+1);
    file = str.substr(found+1);
  }

  /* Check if the directory referred to by the file-name exists. If not create the directory.		*/
  void createDirectory(std::string filename, unsigned int rank)
  {
    std::string path, file;
    splitFilename(filename, path, file);
    if ( !(pathExist(path)) && rank == 0 )
      if ( mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1 )
	{
	  std::cout << "Could not create the directory " << path << ". Probably the given address does not exist." << std::endl;
	  exit(1);
	}
  }

  /* Return the sign of the parameter.									*/
  template <typename T> inline int signof(T x)
  {
    return ( (x > 0) ? 1 : ( (x < 0) ? -1 : 0 ) );
  }

  /* Print a message on the terminal window.                                                            */
  inline void printmessage(std::string filename, unsigned int linenumber, std::string message)
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if (rank == 0)
      {
	/* Set the stream to create for the message.                                                	*/
	std::stringstream printedMessage;
	printedMessage.clear();

	/* Get the time.                                                                               	*/
	time_t rawtime; time(&rawtime);
	std::string timeStr = ctime(&rawtime);

	/* Print the time into the string.                                                             	*/
	printedMessage << timeStr.substr(0,timeStr.size()-1) << " ::: ";

	/* Print the filename and linenumber into the string.                                          	*/
	printedMessage  << filename << ":" << linenumber << " ::: \t \t " << message;

	/* Print on the terminal if there is something written in the stream.                        	*/
	if (printedMessage.str().length() > 0 )     std::cout << printedMessage.str() << std::endl;
      }
  }

  /* Convert any number to a string.                                                                    */
  template<typename numbertype>
  inline std::string stringify(numbertype value)
  {
    std::ostringstream oStream;
    try
    {
	oStream << value;
    }
    catch(std::exception& error)
    {
	printmessage(std::string(__FILE__), __LINE__, std::string("Cannot convert this variable to a string."));
	printmessage(std::string(__FILE__), __LINE__, std::string("Error:") + error.what());
	exit(1);
    }

    return oStream.str();
  }

  /* Define the structure for each charge point.							*/
  struct Charge
  {

    Double			q;		/* Charge of the point in the unit of electron charge.	*/
    FieldVector<Double>		rnp, rnm;	/* Position vector of the charge.			*/
    FieldVector<Double>		gbnp, gbnm;	/* Normalized velocity vector of the charge.		*/

    Charge()
    {
      q = 0.0; rnp = rnm = 0.0; gbnp = gbnm = 0.0;
    };

  };

  /* Function creating a set of halton sequence for the random particle
   * generation.								                        */
  Double halton (unsigned int i, unsigned int j)
  {
    if (i > 20)
      {
	printmessage(std::string(__FILE__), __LINE__, std::string(" dimension can not be larger than 20. ") );
	exit(1);
      }

    unsigned int prime [20] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    int p0, p, k, k0, a;
    Double x = 0.0;

    k0 = j;

    p = prime[i];

    p0 = p;
    k  = k0;
    x  = 0.0;
    while (k > 0)
      {
	a   = k % p;
	x  += a / (double) p0;
	k   = int (k/p);
	p0 *= p;
      }

    return 1.0 - x;
  }

}

#endif
