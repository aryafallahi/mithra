/********************************************************************************************************
 *  readdata.hh : Implementation of the functions reading the lines of the job file.
 ********************************************************************************************************/

#ifndef READDATA_HH_
#define READDATA_HH_

#include <list>
#include <string>
#include <vector>

#include "datainput.h"
#include "fieldvector.h"
#include "stdinclude.h"

namespace Darius
{

  /* Read the data from the input file and store them into a string list.                               */
  std::list <std::string> 	read_file 		(char const * filename);

  /* Clean the stored string vector and make it organized.                                              */
  void 				cleanJobFile 		(std::list <std::string> & jobFile);

  /* Read the ParamaterName at the line.                                                                */
  std::string 			parameterName 		(std::string line);

  /* Read value of a string parameter.                                                                  */
  std::string 			stringValue 		(std::string line);

  /* Read value of a double parameter.                                                                  */
  Double 			doubleValue 		(std::string line);

  /* Read value of an integer parameter.                                                                */
  int 				intValue 		(std::string line);

  /* Read value of a boolean parameter.                                                                 */
  bool 				boolValue 		(std::string line);

  /* Read value of a vector parameter.                                                                  */
  std::vector <Double> 		vectorDoubleValue 	(std::string line);

  /* Read value of a vector parameter.                                                                  */
  std::vector <unsigned int> 	vectorIntValue 		(std::string line);

  /* Read value of a map parameter.                                                                     */
  void 				mapValue 		(std::string line, unsigned int & tag, std::string & model);
}
#endif
