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
}
namespace Darius
{
  void cleanJobFile (std::list <std::string> & jobFile);
}
namespace Darius
{
  std::string parameterName (std::string line);
}
namespace Darius
{
  std::string stringValue (std::string line);
}
namespace Darius
{
  Double doubleValue (std::string line);
}
namespace Darius
{
  int intValue (std::string line);
}
namespace Darius
{
  bool boolValue (std::string line);
}
namespace Darius
{
  std::vector <Double> vectorDoubleValue (std::string line);
}
namespace Darius
{
  std::vector <unsigned int> vectorIntValue (std::string line);
}
namespace Darius
{
  void mapValue (std::string line, unsigned int & tag, std::string & model);
}
#endif
