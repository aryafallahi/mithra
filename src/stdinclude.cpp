/********************************************************************************************************
 *  stdinclude.cpp, set of different standard functions used in the code.
 ********************************************************************************************************/

#include "stdinclude.h"

namespace MITHRA
{

  /* Check if a directory to save the data exists.							*/
  bool pathExist (std::string const & s)
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

  Charge::Charge ()
  {
    q = 0.0; rnp = rnm = 0.0; gbnp = gbnm = 0.0;
    e = 0.0;
  }

  /* Function creating a set of halton sequence for the random particle generation.			*/
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
