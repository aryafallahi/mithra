// datainput.h
//

#ifndef datainput_h
#define datainput_h

#include <iostream>
#include <iterator>
#include <list>
#include <string>
#include <vector>

#include "classes.h"
#include "datainput.h"

namespace Darius
{
  class ParseDarius
  {
  private:
    std::list <std::string> & jobFile_;
    Mesh & mesh_;
    Bunch & bunch_;
    Seed & seed_;
    std::vector <Undulator> & undulator_;
    std::vector <ExtField> & extField_;
    std::vector <FreeElectronLaser> & FEL_;
  public:
    ParseDarius (std::list <std::string> & jobFile, Mesh & mesh, Bunch & bunch, Seed & seed, std::vector <Undulator> & undulator, std::vector <ExtField> & extField, std::vector <FreeElectronLaser> & FEL);
    void setJobParameters ();
    void readMesh (std::list <std::string>::iterator & iter);
    void readBunch (std::list <std::string>::iterator & iter);
    void readField (std::list <std::string>::iterator & iter);
    void readUndulator (std::list <std::string>::iterator & iter);
    void readExtField (std::list <std::string>::iterator & iter);
    void readFEL (std::list <std::string>::iterator & iter);
  };
}
#endif
