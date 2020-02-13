/********************************************************************************************************
 *  datainput.hh : Implementation of the parameter parser for the code
 *********************************************************************************************************/

#ifndef DATAINPUT_H_
#define DATAINPUT_H_

#include <iostream>
#include <iterator>
#include <list>
#include <string>
#include <vector>

#include "classes.h"
#include "datainput.h"

namespace MITHRA
{

  /* The class of functions used for reading the text file of parameters and parsing them to the darius
   * solver.                                                                                            */
  class ParseDarius
  {

  private:

    /* The parameters and data files needed for parsing the values.                                     */
    std::list<std::string>&		jobFile_;
    Mesh& 				mesh_;
    Bunch&				bunch_;
    Seed&				seed_;
    std::vector<Undulator>&		undulator_;
    std::vector<ExtField>&              extField_;
    std::vector<FreeElectronLaser>&	FEL_;

  public:

    ParseDarius 	(std::list<std::string>& jobFile, Mesh& mesh, Bunch& bunch, Seed& seed,
			 std::vector<Undulator>& undulator, std::vector<ExtField>& extField,
			 std::vector<FreeElectronLaser>& FEL);

    /* Read the parameters from the file and set all the parsed parameters for FEL simulation.         	*/
    void setJobParameters ();

    /* Read the parameters parsed for the mesh in the solver.                                    	*/
    void readMesh 	(std::list <std::string>::iterator & iter);

    /* Read the parameters parsed for the bunch in the darius solver.                     		*/
    void readBunch 	(std::list <std::string>::iterator & iter);

    /* Read the parameters parsed for the seed in the darius solver.                               	*/
    void readField 	(std::list <std::string>::iterator & iter);

    /* Read the parameters parsed for the mesh in the solver.                                    	*/
    void readUndulator 	(std::list <std::string>::iterator & iter);

    /* Read the parameters parsed for the seed in the darius solver.                                    */
    void readExtField 	(std::list <std::string>::iterator & iter);

    /* Read the parameters parsed for the FEL output in the darius solver.                              */
    void readFEL 	(std::list <std::string>::iterator & iter);
  };
}
#endif
