#include "barotropic_model.h"

using namespace barotropic_model;

int main(int argc, const char *argv[])
{
    BarotropicModel_Explicit model;
    RossbyHaurwitzTestCase testCase;

    TimeManager timeManager;
    Time startTime, endTime(68*DAYS);

    timeManager.init(startTime, endTime, 0.1*MINUTES);

    model.init(80, 41);
    testCase.calcInitCond(model);

    model.run(timeManager);

    return 0;
}
