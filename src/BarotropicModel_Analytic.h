#ifndef __BarotropicModel_Analytic__
#define __BarotropicModel_Analytic__

#include "BarotropicModel.h"

namespace barotropic_model {

/**
 *  This barotropic model uses A-grid variable stagger configuration and
 *  implicit midpoint time integration method. The underlying numerical
 *  method is a finite difference, which can conserve the total energy
 *  and total mass exactly.
 */
class BarotropicModel_Analytic : public BarotropicModel {
protected:
    SingleLevelField fu, fv;
    SingleLevelField L1U, L1V, L1H;
    SingleLevelField L2U, L2V, L2H;
    
    double dlon, dlat;
    
    vec cosLat, tanLat;
    vec factorCor;  //>! Coriolis factor: 2*OMEGA*sin(lat)
    vec factorCur;  //>! Curvature factor: tan(lat)/R
    vec factorLon;  //>! 1/2/dlon/R/cos(lat)
    vec factorLat;  //>! 1/2/dlat/R/cos(lat)
    vec ghts;  //>! \Phi(\theta, t)
    ivec mLon;  //>! m^n_{j}

    TimeLevelIndex oldTimeIdx, halfTimeIdx, newTimeIdx;
public:
    BarotropicModel_Analytic();
    virtual ~BarotropicModel_Analytic();

    virtual void init(int numLon, int numLat);
    
    virtual void input(const std::string &fileName) {}

    virtual void run(TimeManager &timeManager);

    virtual void integrate(const TimeLevelIndex &oldTimeIdx, double dt);
private:
    double calcTotalEnergy(const TimeLevelIndex &timeIdx) const;

    double calcTotalMass(const TimeLevelIndex &timeIdx) const;

    void calcGeopotentialDepthTendency(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &timeIdx);
    
    void calcGeopotentialDepthTendencyL2(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &timeIdx);

    void calcZonalWindTendency(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &timeIdx);
    void calcL1(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx);
    void calcZonalWindTendencyL2(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &timeIdx);

    void calcMeridionalWindTendency(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &timeIdx);

    void calcZonalWindAdvection(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &timeIdx);

    void calcMeridionalWindAdvection(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &timeIdx);

    void calcZonalWindCoriolis(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &timeIdx);

    void calcMeridionalWindCoriolis(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &timeIdx);

    void calcZonalWindPressureGradient(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &timeIdx);
    void calcZonalWindPressureGradientL2(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &timeIdx);

    void calcMeridionalWindPressureGradient(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &timeIdx);
    void semiAnalytic(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &halfTimeIdx, const TimeLevelIndex &tmpTimeIdx, double dt);
    void check_antisym(const TimeLevelIndex &TimeIdx, const TimeLevelIndex &testTimeIdx);
};

}

#endif // __BarotropicModel_Analytic__
