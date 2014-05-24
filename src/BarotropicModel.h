#ifndef __BarotropicModel__
#define __BarotropicModel__

#include "barotropic_model_commons.h"

namespace barotropic_model {

/**
 *  This is the base class for several barotropic model variants, e.g.,
 *  different variable stagger configuration and time integrator.
 *
 *  The barotropic equations are
 *
 */
class BarotropicModel {
protected:
    Domain *domain;
    Mesh *mesh;
    IOManager io;
    Field u, v, gd;
    SingleLevelField dut, dvt, dgd;
    Field ut, vt, ght;
    SingleLevelField ghu, ghv;
    bool firstRun;
public:
    BarotropicModel() { firstRun = true; }
    virtual ~BarotropicModel() {}

    virtual void init(int numLon, int numLat) = 0;

    virtual void input(const std::string &fileName) = 0;

    virtual void run(TimeManager &timeManager) = 0;

    virtual void integrate(const TimeLevelIndex &oldTimeIdx, double dt) = 0;

    const Domain& getDomain() const { return *domain; }

    const Mesh& getMesh() const { return *mesh; }

    Field& getZonalWind() { return u; }

    Field& getMeridionalWind() { return v; }

    Field& getGeopotentialDepth() { return gd; }
    
};

}

#endif
