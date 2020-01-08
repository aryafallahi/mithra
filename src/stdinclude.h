// stdinclude.h
//
#ifndef stdinclude_h
#define stdinclude_h

#include <mpi.h>
#include <sstream>
#include <string>
#include <sys/stat.h>

#include "fieldvector.h"

namespace Darius
{
  Double const PI = 3.1415926535;
  Double const EPSILON_ZERO = 8.85418782e-12;
  Double const MU_ZERO = 4.0 * PI * 1.0e-7;
  Double const C0 = 1.0 / sqrt(EPSILON_ZERO * MU_ZERO);
  Double const Z0 = sqrt(MU_ZERO / EPSILON_ZERO);
  Double const EC = 1.602e-19;
  Double const EM = 9.109e-31;
  Double const HB = 1.054e-34;
  Double const KB = 1.381e-23;
  std::string const VTS_FILE_SUFFIX = ".vts";
  std::string const PTS_FILE_SUFFIX = ".pvts";
  std::string const PTU_FILE_SUFFIX = ".pvtu";
  std::string const TXT_FILE_SUFFIX = ".txt";
  std::string const VTU_FILE_SUFFIX = ".vtu";
  Complex const I = Complex (0.0,1.0);
}

namespace Darius
{
  enum SignalType
  {
    NEUMANN,
    GAUSSIAN,
    SECANT,
    FLATTOP
  };
}
namespace Darius
{
  enum SeedType
  {
    PLANEWAVE,
    PLANEWAVECONFINED,
    GAUSSIANBEAM
  };
}
namespace Darius
{
  enum ExtFieldType
  {
    EMWAVE
  };
}
namespace Darius
{
  enum SamplingType
  {
    ATPOINT,
    OVERLINE,
    INPLANE,
    ALLDOMAIN
  };
}
namespace Darius
{
  enum PlaneType
  {
    XNORMAL,
    YNORMAL,
    ZNORMAL
  };
}
namespace Darius
{
  enum FieldType
  {
    Ex,
    Ey,
    Ez,
    Bx,
    By,
    Bz,
    Ax,
    Ay,
    Az,
    Jx,
    Jy,
    Jz,
    Q,
    F
  };
}
namespace Darius
{
  enum UndulatorType
  {
    STATIC,
    OPTICAL
  };
}
namespace Darius
{
  enum SolverType
  {
    FD,
    NSFD
  };
}
namespace Darius
{
  bool isabsolute (std::string filename);
}
namespace Darius
{
  bool pathExist (std::string const & s);
}
namespace Darius
{
  void splitFilename (std::string const & str, std::string & path, std::string & file);
}
namespace Darius
{
  void createDirectory (std::string filename, unsigned int rank);
}
namespace Darius
{
  template <typename T>
  int signof (T x);
}
namespace Darius
{
  void printmessage (std::string filename, unsigned int linenumber, std::string message);
}
namespace Darius
{
  template <typename numbertype>
  std::string stringify (numbertype value);
}
namespace Darius
{
  struct Charge
  {
    Double q;
    FieldVector <Double> rnp;
    FieldVector <Double> rnm;
    FieldVector <Double> gbnp;
    FieldVector <Double> gbnm;
    Charge ();
  };
}
namespace Darius
{
  Double halton (unsigned int i, unsigned int j);
}
namespace Darius
{
  inline bool isabsolute (std::string filename)
  {
    return (filename.compare(0,1,"/") == 0);
  }
}
namespace Darius
{
  template <typename T>
  inline int signof (T x)
  {
    return ( (x > 0) ? 1 : ( (x < 0) ? -1 : 0 ) );
  }
}
namespace Darius
{
  inline void printmessage (std::string filename, unsigned int linenumber, std::string message)
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
}
namespace Darius
{
  template <typename numbertype>
  inline std::string stringify (numbertype value)
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
}
#endif
