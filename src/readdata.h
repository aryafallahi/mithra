// readdata.h
//

#ifndef readdata_h
#define readdata_h

#include <list>
#include <string>
#include <vector>

#include "datainput.h"
#include "fieldvector.h"
#include "stdinclude.h"

namespace Darius
{
  std::list <std::string> read_file (char const * filename);

  void cleanJobFile (std::list <std::string> & jobFile);

  std::string parameterName (std::string line);

  std::string stringValue (std::string line);

  Double doubleValue (std::string line);

  int intValue (std::string line);

  bool boolValue (std::string line);

  std::vector <Double> vectorDoubleValue (std::string line);

  std::vector <unsigned int> vectorIntValue (std::string line);

  void mapValue (std::string line, unsigned int & tag, std::string & model);
}
#endif
