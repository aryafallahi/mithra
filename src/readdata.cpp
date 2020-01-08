// readdata.cpp
//

#include "readdata.h"

namespace Darius
{
  std::list <std::string> read_file (char const * filename)
  {
    std::string                 line;
    std::ifstream               myfile (filename);
    std::list<std::string>      jobFile;

    /* Each line will be a string of the jobFile list.                                                  */
    if (myfile.is_open())
      {
        while (myfile.good())
          {
            std::getline (myfile, line);
            jobFile.push_back(line);
          }
        myfile.close();
      }

    else
      {
        std::cout << "Unable to open file" << std::endl;
        exit(1);
      }

    return jobFile;
  }
}
namespace Darius
{
  void cleanJobFile (std::list <std::string> & jobFile)
  {
    std::list<std::string> cleanedJobFile;
    for (std::list<std::string>::iterator iter = jobFile.begin(); iter != jobFile.end(); ++iter)
      {
        for (std::string::iterator it = (*iter).begin(); it != (*iter).end(); ++it)
          {
            if (*it == '#') { (*iter).erase(it,(*iter).end()); break; }
          }
        for (std::string::iterator it = (*iter).begin(); it != (*iter).end(); ++it)
          {
            if ( (*it == ' ') || (*it == '\t') ) { (*iter).erase(it); --it; }
          }
        if (!(*iter).size() == 0) { cleanedJobFile.push_back(*iter); }
      }
    jobFile.swap(cleanedJobFile);
  }
}
namespace Darius
{
  std::string parameterName (std::string line)
  {
    size_t posEqual = line.find("=");
    std::string name = line.substr(0,posEqual);
    return name;
  }
}
namespace Darius
{
  std::string stringValue (std::string line)
  {
    size_t posEqual = line.find("=");
    std::string value = line.substr(posEqual+1,line.size()-1);
    if ( value.compare(0,1,"\"") == 0 )
      {
        value.erase(value.find("\""),1);
        value.erase(value.find("\""),1);
      }
    return value;
  }
}
namespace Darius
{
  Double doubleValue (std::string line)
  {
    size_t posEqual = line.find("=");
    std::string doubleStr = line.substr(posEqual+1,line.size()-1);
    Double value = std::atof(doubleStr.c_str());
    return value;
  }
}
namespace Darius
{
  int intValue (std::string line)
  {
    size_t posEqual = line.find("=");
    std::string doubleStr = line.substr(posEqual+1,line.size()-1);
    Double value = std::atof(doubleStr.c_str());
    int intvalue = int ( value );
    return intvalue;
  }
}
namespace Darius
{
  bool boolValue (std::string line)
  {
    bool value;
    size_t posEqual = line.find("=");
    std::string boolStr = line.substr(posEqual+1,line.size()-1);
    if ( boolStr.compare("true") == 0 )         value = true;
    else if ( boolStr.compare("false") == 0 )   value = false;
    else{
      printmessage(std::string(__FILE__), __LINE__, std::string("boolValue(std::string line) got unexpected input. Input should be \"true\" or \"false\" ") );
      exit(1);
    }
    return value;
  }
}
namespace Darius
{
  std::vector <Double> vectorDoubleValue (std::string line)
  {
    size_t posEqual = line.find("=");
    line = line.substr(posEqual+1,line.size()-1);
    line.erase(0,1);
    line.replace(line.size()-1,1,",");
    std::vector<Double> doubleVector;
    while ( line.find(",") < line.size() )
      {
        size_t posComma = line.find(",");
        std::string doubleStr = line.substr(0,posComma);
        line.erase(0,posComma+1);
        doubleVector.push_back(std::atof(doubleStr.c_str()));
      }
    return doubleVector;
  }
}
namespace Darius
{
  std::vector <unsigned int> vectorIntValue (std::string line)
  {
    size_t posEqual = line.find("=");
    line = line.substr(posEqual+1,line.size()-1);
    line.erase(0,1);
    line.replace(line.size()-1,1,",");
    std::vector<unsigned int> intVector;
    while ( line.find(",") < line.size() )
      {
        size_t posComma = line.find(",");
        std::string doubleStr = line.substr(0,posComma);
        line.erase(0,posComma+1);
        intVector.push_back(std::atoi(doubleStr.c_str()));
      }
    return intVector;
  }
}
namespace Darius
{
  void mapValue (std::string line, unsigned int & tag, std::string & model)
  {
    size_t posEqual = line.find("=");
    line = line.substr(posEqual+1,line.size()-1);
    line.erase(0,1);
    line.erase(line.size()-1,1);
    size_t posComma = line.find(",");
    std::string tagStr = line.substr(0,posComma);
    tag = std::atof(tagStr.c_str());
    model = line.substr(posComma+2,line.size()-posComma-3);
  }
}
