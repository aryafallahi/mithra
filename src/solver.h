// solver.h
//

#ifndef solver_h
#define solver_h

#include <iomanip>
#include <list>
#include <vector>

#include "classes.h"
#include "database.h"
#include "fieldvector.h"

namespace Darius
{
  class Solver
  {
  public:
    Solver (Mesh & mesh, Bunch & bunch, Seed & seed, std::vector <Undulator> & undulator, std::vector <ExtField> & extField, std::vector <FreeElectronLaser> & FEL);
    void lorentzBoost ();
    void initialize ();
    void initializeMesh ();
    void initializeField ();
    void initializeSeedSampling ();
    void initializeSeedVTK ();
    void initializeSeedProfile ();
    void initializeBunchUpdate ();
    void initializeBunch ();
    void bunchUpdate ();
    void bunchSample ();
    void bunchVisualize ();
    void bunchProfile ();
    void undulatorField (UpdateBunchParallel & ubp, FieldVector <Double> & r);
    void externalField (UpdateBunchParallel & ubp, FieldVector <Double> & r);
    void initializePowerSample ();
    void initializePowerVisualize ();
    void powerSample ();
    void powerVisualize ();
    void initializeRadiationEnergySample ();
    void radiationEnergySample ();
    void finalize ();
    virtual void fieldEvaluate (long int m) = 0;
    static bool undulatorCompare (Undulator i, Undulator j);
    Mesh & mesh_;
    Bunch & bunch_;
    Seed & seed_;
    std::vector <Undulator> & undulator_;
    std::vector <ExtField> & extField_;
    std::vector <FreeElectronLaser> & FEL_;
    std::vector <FieldVector<Double> > * anp1_;
    std::vector <FieldVector<Double> > * an_;
    std::vector <FieldVector<Double> > * anm1_;
    std::vector <Double> * fnp1_;
    std::vector <Double> * fn_;
    std::vector <Double> * fnm1_;
    std::vector <bool> pic_;
    std::vector <FieldVector<Double> > en_;
    std::vector <FieldVector<Double> > bn_;
    std::vector <FieldVector<Double> > jn_;
    std::vector <Double> rn_;
    std::vector <FieldVector<Double> > r_;
    int N0_;
    int N1_;
    int N2_;
    int N1N0_;
    unsigned int Nc_;
    std::list <Charge>::iterator iterQB_;
    std::list <Charge>::iterator iterQE_;
    int np_;
    int k0_;
    Double (zp_) [2];
    Double xmin_;
    Double xmax_;
    Double ymin_;
    Double ymax_;
    Double zmin_;
    Double zmax_;
    Double timep1_;
    Double time_;
    Double timem1_;
    unsigned int nTime_;
    std::list <Charge> chargeVectorn_;
    Double timeBunch_;
    unsigned int nTimeBunch_;
    Double nUpdateBunch_;
    Double gamma_;
    Double beta_;
    Double dt_;
    UpdateField <Double> uf_;
    SampleField <Double> sf_;
    std::vector <VisualizeField<Double> > vf_;
    ProfileField <Double> pf_;
    UpdateBunch ub_;
    SampleBunch sb_;
    VisualizeBunch vb_;
    ProfileBunch pb_;
    UpdateCurrent uc_;
    std::vector <SampleRadiationPower> rp_;
    std::vector <SampleRadiationEnergy> re_;
    int rank_;
    int size_;
    Double c0_;
    Double m0_;
    Double e0_;
  };
}
#endif
