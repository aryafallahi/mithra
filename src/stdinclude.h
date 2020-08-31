/********************************************************************************************************
 *  stdinclude.hh, set of different standard databases used in the code.
 ********************************************************************************************************/

#ifndef STDINCLUDE_HH_
#define STDINCLUDE_HH_

#include <mpi.h>
#include <sstream>
#include <string>
#include <sys/stat.h>

#include "fieldvector.h"

namespace MITHRA
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
  const Double TPI		= 2.0 * PI;
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

	/* Print the shortened filename and linenumber into the string.                                	*/
	const char * elem = filename.c_str();
	const char * shortfn = elem;
	while ( *elem != '\0' ){
	    if ( *elem == '/' )
	      shortfn = elem + 1;
	    elem = elem + 1;
	}
	printedMessage  << shortfn << ":" << linenumber << " ::: \t \t " << message;

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

    /* Double flag determining if the particle is passing the entrance point of the undulator. This flag
     * can be used for better boosting the bunch to the moving frame. We need to consider it to be double,
     * because this flag needs to be communicated during bunch update.					*/
    Double			e;

    Charge();

  };

  /* Check if a directory to save the data exists.							*/
  bool pathExist (std::string const & s);

  /* Split the file name to two strings including its path and file name.                               */
  void splitFilename (std::string const & str, std::string & path, std::string & file);

  /* Check if the directory referred to by the file-name exists. If not create the directory.		*/
  void createDirectory (std::string filename, unsigned int rank);

  /* Function creating a set of halton sequence for the random particle generation.			*/
  Double halton (unsigned int i, unsigned int j);

  /* Hello message when initialising mithra.														*/
  void helloMessage();

  /* Positive remainder between two double values.							*/
  Double pmod ( const Double& a, const Double& b);

}
#endif
