// classes.h
//

#ifndef classes_h
#define classes_h

#include <list>

#include "database.h"
#include "fieldvector.h"
#include "stdinclude.h"


namespace Darius
{
  struct Mesh
  {
    Double lengthScale_;
    FieldVector <Double> meshCenter_;
    FieldVector <Double> meshLength_;
    FieldVector <Double> meshResolution_;
    Double timeScale_;
    Double timeStep_;
    Double totalTime_;
    unsigned int truncationOrder_;
    bool spaceCharge_;
    SolverType solver_;
    void show ();
    void initialize ();
  };
}
namespace Darius
{
  class Bunch
  {
  public:
    typedef std::list <Charge> ChargeVector;
    Bunch ();
    void initializeManual (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia);
    void initializeEllipsoid (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia);
    void initialize3DCrystal (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia);
    void initializeFile (BunchInitialize bunchInit, ChargeVector & chargeVector, Double (zp) [2], int rank, int size, int ia);
    std::vector <BunchInitialize> bunchInit_;
    std::string directory_;
    std::string basename_;
    bool sampling_;
    Double timeStep_;
    Double rhythm_;
    bool bunchVTK_;
    std::string bunchVTKDirectory_;
    std::string bunchVTKBasename_;
    Double bunchVTKRhythm_;
    bool bunchProfile_;
    std::string bunchProfileDirectory_;
    std::string bunchProfileBasename_;
    std::vector <Double> bunchProfileTime_;
    Double bunchProfileRhythm_;
    Double zu_;
    Double beta_;
    void show ();
  };
}
namespace Darius
{
  class Signal
  {
  public:
    Signal ();
    void initialize (std::string type, Double l0, Double s, Double l, Double cep);
  public:
    SignalType signalType_;
    Double t0_;
    Double s_;
    Double f0_;
    Double cep_;
    Double self (Double & t, Double & phase);
    void show ();
  };
}
namespace Darius
{
  class Seed
  {
  public:
    Seed ();
    void initialize (std::string type, std::vector <Double> position, std::vector <Double> direction, std::vector <Double> polarization, Double amplitude, std::vector <Double> radius, Signal signal);
  public:
    SeedType seedType_;
    Double c0_;
    FieldVector <Double> position_;
    FieldVector <Double> direction_;
    FieldVector <Double> polarization_;
    Double amplitude_;
    std::vector <Double> radius_;
    Signal signal_;
    Double beta_;
    Double gamma_;
    Double dt_;
  private:
    Double gamma;
    Double tsignal;
    Double d;
    Double l;
    Double zRp;
    Double wrp;
    Double zRs;
    Double wrs;
    Double x;
    Double y;
    Double z;
    Double p;
    Double t;
    FieldVector <Double> rv;
    FieldVector <Double> yv;
    FieldVector <Double> ax;
    FieldVector <Double> az;
    FieldVector <Double> rl;
    Double tl;
  public:
    void fields (FieldVector <Double> & aufpunkt, Double & time, FieldVector <Double> & a);
    bool sampling_;
    SamplingType samplingType_;
    std::vector <FieldType> samplingField_;
    std::string samplingDirectory_;
    std::string samplingBasename_;
    Double samplingRhythm_;
    std::vector <FieldVector<Double> > samplingPosition_;
    FieldVector <Double> samplingLineBegin_;
    FieldVector <Double> samplingLineEnd_;
    FieldVector <Double> samplingSurfaceBegin_;
    FieldVector <Double> samplingSurfaceEnd_;
    unsigned int samplingRes_;
    struct vtk
    {
      bool sample_;
      std::vector <FieldType> field_;
      std::string directory_;
      SamplingType type_;
      std::string basename_;
      Double rhythm_;
      PlaneType plane_;
      FieldVector <Double> position_;
      vtk ();
    };
    std::vector <vtk> vtk_;
    bool profile_;
    std::vector <FieldType> profileField_;
    std::string profileDirectory_;
    std::string profileBasename_;
    std::vector <Double> profileTime_;
    Double profileRhythm_;
    SamplingType samplingType (std::string samplingType);
    SamplingType vtkType (std::string vtkType);
    PlaneType planeType (std::string planeType);
    FieldType fieldType (std::string fieldType);
    void show ();
  };
}
namespace Darius
{
  class Undulator
  {
  public:
    Undulator ();
    Double k_;
    Double lu_;
    Double rb_;
    unsigned int length_;
    Double beta_;
    Double gamma_;
    Double dt_;
    Double theta_;
    UndulatorType type_;
    SeedType seedType_;
    Double c0_;
    FieldVector <Double> position_;
    FieldVector <Double> direction_;
    FieldVector <Double> polarization_;
    Double amplitude_;
    std::vector <Double> radius_;
    Signal signal_;
    UndulatorType undulatorType (std::string undulatorType);
    void initialize (std::string type, std::vector <Double> position, std::vector <Double> direction, std::vector <Double> polarization, Double amplitude, std::vector <Double> radius, Double wavelength, Signal signal);
    void show ();
  };
}
namespace Darius
{
  class ExtField
  {
  public:
    ExtField ();
    ExtFieldType type_;
    SeedType seedType_;
    Double c0_;
    FieldVector <Double> position_;
    FieldVector <Double> direction_;
    FieldVector <Double> polarization_;
    Double amplitude_;
    std::vector <Double> radius_;
    Signal signal_;
    void initialize (std::string type, std::vector <Double> position, std::vector <Double> direction, std::vector <Double> polarization, Double amplitude, std::vector <Double> radius, Double wavelength, Signal signal);
    void show ();
  };
}
namespace Darius
{
  struct FreeElectronLaser
  {
    struct RadiationPower
    {
      std::vector <Double> z_;
      bool sampling_;
      std::string directory_;
      std::string basename_;
      Double lineBegin_;
      Double lineEnd_;
      unsigned int res_;
      SamplingType samplingType_;
      std::vector <Double> lambda_;
      Double lambdaMin_;
      Double lambdaMax_;
      unsigned int lambdaRes_;
      void samplingType (std::string samplingType);
      RadiationPower ();
    };
    RadiationPower radiationPower_;
    struct vtk
    {
      Double z_;
      bool sampling_;
      std::string directory_;
      std::string basename_;
      Double rhythm_;
      Double lambda_;
      vtk ();
    };
    vtk vtk_;
    struct RadiationEnergy
    {
      std::vector <Double> z_;
      bool sampling_;
      std::string directory_;
      std::string basename_;
      Double rhythm_;
      Double lineBegin_;
      Double lineEnd_;
      Double res_;
      SamplingType samplingType_;
      std::vector <Double> lambda_;
      Double lambdaMin_;
      Double lambdaMax_;
      Double lambdaRes_;
      void samplingType (std::string samplingType);
      RadiationEnergy ();
    };
    RadiationEnergy radiationEnergy_;
  };
}
#endif
