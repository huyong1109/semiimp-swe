#include "BarotropicModel_Analytic.h"

namespace barotropic_model {

BarotropicModel_Analytic::BarotropicModel_Analytic() {
    REPORT_ONLINE;
}

BarotropicModel_Analytic::~BarotropicModel_Analytic() {
    REPORT_OFFLINE;
}

void BarotropicModel_Analytic::init(int numLon, int numLat) {
    // -------------------------------------------------------------------------
    // initialize domain
    domain = new Domain(2);
    domain->setRadius(6.371e6);
    // -------------------------------------------------------------------------
    // initialize mesh
    mesh = new Mesh(*domain);
    mesh->init(numLon, numLat);
    dlon = mesh->getGridInterval(0, FULL, 0);
    dlat = mesh->getGridInterval(1, FULL, 0); // assume equidistance grids
    // -------------------------------------------------------------------------
    // create variables
    u.create("u", "m s-1", "zonal wind speed", *mesh, CENTER, HAS_HALF_LEVEL);
    v.create("v", "m s-1", "meridional wind speed", *mesh, CENTER, HAS_HALF_LEVEL);
    gd.create("gd", "m2 s-2", "geopotential depth", *mesh, CENTER, HAS_HALF_LEVEL);

    ut.create("ut", "(m s-1)*m-2", "transformed zonal wind speed", *mesh, CENTER, HAS_HALF_LEVEL);
    vt.create("vt", "(m s-1)*m-2", "transformed meridional wind speed", *mesh, CENTER, HAS_HALF_LEVEL);
    ght.create("ght", "m-2", "transformed geopotential height", *mesh, CENTER, HAS_HALF_LEVEL);

    dut.create("dut", "m s-2", "zonal wind speed tendency", *mesh, CENTER);
    dvt.create("dvt", "m s-2", "meridional zonal speed tendency", *mesh, CENTER);
    dgd.create("dgd", "m-2 s-1", "geopotential depth tendency", *mesh, CENTER);

    ghu.create("ghu", "m2 s-1", "ut * ght", *mesh, CENTER);
    ghv.create("ghv", "m2 s-1", "vt * ght", *mesh, CENTER);

    fu.create("fu", "* m s-1", "* * u", *mesh, CENTER);
    fv.create("fv", "* m s-1", "* * v", *mesh, CENTER);
    
    L1U.create("lu", "* m s-1", "* * v", *mesh, CENTER, HAS_HALF_LEVEL);
    L1V.create("lv", "* m s-1", "* * v", *mesh, CENTER, HAS_HALF_LEVEL);
    L1H.create("lh", "* m s-1", "* * v", *mesh, CENTER, HAS_HALF_LEVEL);
    L2U.create("lu", "* m s-1", "* * v", *mesh, CENTER, HAS_HALF_LEVEL);
    L2V.create("lv", "* m s-1", "* * v", *mesh, CENTER, HAS_HALF_LEVEL);
    L2H.create("lh", "* m s-1", "* * v", *mesh, CENTER, HAS_HALF_LEVEL);
    // -------------------------------------------------------------------------




    int js = 0, jn = mesh->getNumGrid(1, FULL)-1;

    cosLat.set_size(mesh->getNumGrid(1, FULL));
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        cosLat[j] = mesh->getCosLat(FULL, j);
    }
    cosLat[js] = mesh->getCosLat(HALF,   js)*0.25;
    cosLat[jn] = mesh->getCosLat(HALF, jn-1)*0.25;

    tanLat.set_size(mesh->getNumGrid(1, FULL));
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        tanLat[j] = mesh->getTanLat(FULL, j);
    }
    tanLat[js] = -1/cosLat[js];
    tanLat[jn] =  1/cosLat[jn];

    factorCor.set_size(mesh->getNumGrid(1, FULL));
    for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
        factorCor[j] = 2*OMEGA*mesh->getSinLat(FULL, j);
    }

    factorCur.set_size(mesh->getNumGrid(1, FULL));
    for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
        factorCur[j] = tanLat[j]/domain->getRadius();
    }

    factorLon.set_size(mesh->getNumGrid(1, FULL));
    for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
        factorLon[j] = 1/(2*dlon*domain->getRadius()*cosLat[j]);
    }

    factorLat.set_size(mesh->getNumGrid(1, FULL));
    for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
        factorLat[j] = 1/(2*dlat*domain->getRadius()*cosLat[j]);
    }
    ghts.set_size(mesh->getNumGrid(1, FULL));
    mLon.set_size(mesh->getNumGrid(1, FULL));
    // -------------------------------------------------------------------------
    // set variables in Poles
    for (int i = -1; i < mesh->getNumGrid(0, FULL)+1; ++i) {
        dut(i, js) = 0.0; dut(i, jn) = 0.0;
        dvt(i, js) = 0.0; dvt(i, jn) = 0.0;
        ghu(i, js) = 0.0; ghu(i, jn) = 0.0;
        ghv(i, js) = 0.0; ghv(i, jn) = 0.0;
    }
}

void BarotropicModel_Analytic::run(TimeManager &timeManager) {
    // -------------------------------------------------------------------------
    // initialize IO manager
    io.init(timeManager);
    int fileIdx = io.registerOutputFile(*mesh, "output", IOFrequencyUnit::HOURS, 1);
    io.file(fileIdx).registerOutputField<double, 2, FULL_DIMENSION>(3, &u, &v, &gd);
    // -------------------------------------------------------------------------
    // output initial condition
    io.create(fileIdx);
    io.output<double, 2>(fileIdx, oldTimeIdx, 3, &u, &v, &gd);
    io.close(fileIdx);
    // -------------------------------------------------------------------------
    // main integration loop
    while (!timeManager.isFinished()) {
        integrate(oldTimeIdx, timeManager.getStepSize());
        timeManager.advance();
        oldTimeIdx.shift();
        io.create(fileIdx);
        io.output<double, 2>(fileIdx, oldTimeIdx, 3, &u, &v, &gd);
        io.close(fileIdx);
    }
}

void BarotropicModel_Analytic::integrate(const TimeLevelIndex &oldTimeIdx,
                                                   double dt) {
    // -------------------------------------------------------------------------
    // set time level indices
    halfTimeIdx = oldTimeIdx+0.5;
    newTimeIdx = oldTimeIdx+1;
    // -------------------------------------------------------------------------
    // copy states
    if (firstRun) {
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = -1; i < mesh->getNumGrid(0, FULL)+1; ++i) {
                ght(oldTimeIdx, i, j) = sqrt(gd(oldTimeIdx, i, j));
                ut(oldTimeIdx, i, j) = u(oldTimeIdx, i, j)*ght(oldTimeIdx, i, j);
                vt(oldTimeIdx, i, j) = v(oldTimeIdx, i, j)*ght(oldTimeIdx, i, j);
            }
        }
        firstRun = false;
	// get old total energy and mass
    	double e0 = calcTotalEnergy(oldTimeIdx);
    	double m0 = calcTotalMass(oldTimeIdx);
    	cout << "energy: ";
    	cout << std::fixed << setw(20) << setprecision(2) << e0 << "  ";
    	cout << "mass: ";
    	cout << setw(20) << setprecision(2) << m0 << endl;
    }
        // ---------------------------------------------------------------------
	// Caculate ghts : horizontal avarage of ght
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
	    ghts[j] = 0.0;
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                ghts[j] += ght(oldTimeIdx, i, j);
            }
	    ghts[j] /= mesh->getNumGrid(0, FULL);
	    mLon[j] = (int)(floor(2.0*ghts[j]*dt*factorLon[j]));
	    //cout<<"lat = "<<setw(4)<<j<<setw(20)<<ghts[j] <<setw(4)<<mLon[j]<<endl;
	    ghts[j] = (double)(mLon[j]/(dt*factorLon[j]*2.0));
	    cout<<"lat = "<<setw(4)<<j<<setw(20)<<ghts[j] <<setw(4)<<mLon[j]<<endl;
        }
	//two half semi-implicit step 
	cout<<"First step"<<endl;
	semiAnalytic(oldTimeIdx, halfTimeIdx, oldTimeIdx,  dt);	
	cout<<"check antisymetrical "<<endl;
	check_antisym(oldTimeIdx, halfTimeIdx);	
	cout<<"Second step"<<endl;
	semiAnalytic(oldTimeIdx, newTimeIdx, halfTimeIdx,  dt);	
	cout<<"Third step"<<endl;
	// keep L2 P^{n+1/2}
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                L2H(i, j) = dgd(i, j);
                L2U(i, j) = dut(i, j);
                L2V(i, j) = dvt(i, j);
            }
	}
	
	check_antisym(oldTimeIdx, newTimeIdx);	
	//Caculate \beta  = lqlp/lq2;
	long double beta = 0.0;
	long double lplq = 0.0;
	long double lp2  = 0.0;

        calcGeopotentialDepthTendencyL2(oldTimeIdx, newTimeIdx);
        calcZonalWindTendencyL2(oldTimeIdx, newTimeIdx);
	calcMeridionalWindTendency(oldTimeIdx, newTimeIdx);
	// (H,H) 
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
		lplq += (L1H( i,j)/dt+ dgd(i,j))*(L1H(i,j)/dt + L2H(i,j))*cosLat[j];
		lp2  += (L1H( i,j)/dt+ dgd(i,j))*(L1H(i,j)/dt + dgd(i,j))*cosLat[j];
            }
	}
	// (U,U) 
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
		lplq += (L1U(i,j)/dt+ dut(i,j))*(L1U(i,j)/dt+L2U(i,j))*cosLat[j];
		lp2  += (L1U(i,j)/dt+ dut(i,j))*(L1U(i,j)/dt+dut(i,j))*cosLat[j];
            }
	}
	// (V,V) 
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
		lplq += (L1V(i,j)/dt+ dvt(i,j))*(L1V(i,j)/dt+L2V(i,j))*cosLat[j];
		lp2  += (L1V(i,j)/dt+ dvt(i,j))*(L1V(i,j)/dt+dvt(i,j))*cosLat[j];
            }
	}
	beta = lplq/lp2;
	cout<<"beta = " <<setprecision(15)<<beta<<endl;
//update F^{n+1}
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                gd(newTimeIdx, i, j) = gd(oldTimeIdx, i, j)-beta*(L1H(i,j)+ dt*dgd(i,j));
                ut(newTimeIdx, i, j) = ut(oldTimeIdx, i, j)-beta*(L1U(i,j)+ dt*dut(i,j));
                vt(newTimeIdx, i, j) = vt(oldTimeIdx, i, j)-beta*(L1V(i,j)+ dt*dvt(i,j));
            }
	}

        ut.applyBndCond(newTimeIdx);
        vt.applyBndCond(newTimeIdx);
        gd.applyBndCond(newTimeIdx);

        // ---------------------------------------------------------------------
        // transform geopotential height
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                ght(newTimeIdx, i, j) = sqrt(gd(newTimeIdx, i, j));
            }
        }
        ght.applyBndCond(newTimeIdx);
       // ---------------------------------------------------------------------
       // transform back velocity on half time level
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                u(newTimeIdx, i, j) = ut(newTimeIdx, i, j)/ght(newTimeIdx, i, j);
                v(newTimeIdx, i, j) = vt(newTimeIdx, i, j)/ght(newTimeIdx, i, j);
            }
        }
        u.applyBndCond(newTimeIdx);
        v.applyBndCond(newTimeIdx);

	//
        // get new total energy and mass
        double e1 = calcTotalEnergy(newTimeIdx);
#ifndef NDEBUG
        double m1 = calcTotalMass(newTimeIdx);
        cout << "energy: ";
        cout << std::fixed << setw(20) << setprecision(2) << e1 << "  ";
        cout << "mass: ";
        cout << setw(20) << setprecision(2) << m1 << " " << endl;
#endif
}

double BarotropicModel_Analytic::calcTotalEnergy(const TimeLevelIndex &timeIdx) const {
    double totalEnergy = 0.0;
    for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            totalEnergy += (pow(ut(timeIdx, i, j), 2)+
                            pow(vt(timeIdx, i, j), 2)+
                            pow(gd(timeIdx, i, j), 2))*cosLat[j];
        }
    }
    return totalEnergy;
}

double BarotropicModel_Analytic::calcTotalMass(const TimeLevelIndex &timeIdx) const {
    double totalMass = 0.0;
    for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            totalMass += gd(timeIdx, i, j)*cosLat[j];
        }
    }
    return totalMass;
}
/**
 *  Input: ut, vt, ght
 *  Intermediate: ghu, ghv
 *  Output: dgd
 */
void BarotropicModel_Analytic::calcGeopotentialDepthTendency(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    // calculate intermediate variables
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = -1; i < mesh->getNumGrid(0, FULL)+1; ++i) {
            ghu(i, j) = ut(timeIdx, i, j)*ght(oldtimeIdx, i, j);
            ghv(i, j) = vt(timeIdx, i, j)*ght(oldtimeIdx, i, j)*cosLat[j];
        }
    }
    // normal grids
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            dgd(i, j) = (ghu(i+1, j)-ghu(i-1, j))*factorLon[j]+
                        (ghv(i, j+1)-ghv(i, j-1))*factorLat[j];
        }
    }
    // pole grids
    // last character 's' and 'n' mean 'Sorth Pole' and 'North Pole' respectively
    int js = 0, jn = mesh->getNumGrid(1, FULL)-1;
    double dgds = 0.0, dgdn = 0.0;
    for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
        dgds += ghv(i, js+1);
        dgdn -= ghv(i, jn-1);
    }
    dgds *= factorLat[js]/mesh->getNumGrid(0, FULL);
    dgdn *= factorLat[jn]/mesh->getNumGrid(0, FULL);
    for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
        dgd(i, js) = dgds;
        dgd(i, jn) = dgdn;
    }
#ifndef NDEBUG
    double tmp = 0.0;
    for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            tmp += dgd(i, j)*cosLat[j];
        }
    }
    assert(fabs(tmp) < 1.0e-10);
#endif
}

/**
 *  Input: ut, vt, ght
 *  Intermediate: ghu, ghv
 *  Output: dgd
 */
void BarotropicModel_Analytic::calcGeopotentialDepthTendencyL2(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    // calculate intermediate variables
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = -1; i < mesh->getNumGrid(0, FULL)+1; ++i) {
            ghu(i, j) = ut(timeIdx, i, j)*(ght(oldtimeIdx, i, j)-ghts[j]);
            ghv(i, j) = vt(timeIdx, i, j)*ght(oldtimeIdx, i, j)*cosLat[j];
        }
    }
    // normal grids
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            dgd(i, j) = (ghu(i+1, j)-ghu(i-1, j))*factorLon[j]+
                        (ghv(i, j+1)-ghv(i, j-1))*factorLat[j];
        }
    }
    // pole grids
    // last character 's' and 'n' mean 'Sorth Pole' and 'North Pole' respectively
    int js = 0, jn = mesh->getNumGrid(1, FULL)-1;
    double dgds = 0.0, dgdn = 0.0;
    for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
        dgds += ghv(i, js+1);
        dgdn -= ghv(i, jn-1);
    }
    dgds *= factorLat[js]/mesh->getNumGrid(0, FULL);
    dgdn *= factorLat[jn]/mesh->getNumGrid(0, FULL);
    for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
        dgd(i, js) = dgds;
        dgd(i, jn) = dgdn;
    }
#ifndef NDEBUG
    double tmp = 0.0;
    for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            tmp += dgd(i, j)*cosLat[j];
        }
    }
    assert(fabs(tmp) < 1.0e-10);
#endif
}


void BarotropicModel_Analytic::calcZonalWindTendencyL2(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    calcZonalWindAdvection(oldtimeIdx, timeIdx);
    calcZonalWindCoriolis(oldtimeIdx,timeIdx);
    calcZonalWindPressureGradientL2(oldtimeIdx,timeIdx);
}

void BarotropicModel_Analytic::calcZonalWindTendency(const TimeLevelIndex &oldtimeIdx,const TimeLevelIndex &timeIdx) {
    calcZonalWindAdvection(oldtimeIdx,timeIdx);
    calcZonalWindCoriolis(oldtimeIdx,timeIdx);
    calcZonalWindPressureGradient(oldtimeIdx,timeIdx);
}
void BarotropicModel_Analytic::calcMeridionalWindTendency(const TimeLevelIndex &oldtimeIdx,const TimeLevelIndex &timeIdx) {
    calcMeridionalWindAdvection(oldtimeIdx,timeIdx);
    calcMeridionalWindCoriolis(oldtimeIdx,timeIdx);
    calcMeridionalWindPressureGradient(oldtimeIdx,timeIdx);
}

/**
 *  Input: u, v, ut
 *  Output, s1
 */
void BarotropicModel_Analytic::calcZonalWindAdvection(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = -1; i < mesh->getNumGrid(0, FULL)+1; ++i) {
            fu(i, j) = ut(timeIdx, i, j)*u(oldtimeIdx, i, j);
            fv(i, j) = ut(timeIdx, i, j)*v(oldtimeIdx, i, j)*cosLat[j];
        }
    }
    // normal grids
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            double dx1 = fu(i+1, j)-fu(i-1, j);
            double dy1 = fv(i, j+1)-fv(i, j-1);
            double dx2 = u(oldtimeIdx, i, j)*(ut(timeIdx, i+1, j)-ut(timeIdx, i-1, j));
            double dy2 = v(oldtimeIdx, i, j)*(ut(timeIdx, i, j+1)-ut(timeIdx, i, j-1))*cosLat[j];
            dut(i, j) = 0.5*((dx1+dx2)*factorLon[j]+(dy1+dy2)*factorLat[j]);
        }
    }
}

/**
 *  Input: u, v, vt
 *  Output, dvt
 */
void BarotropicModel_Analytic::calcMeridionalWindAdvection(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = -1; i < mesh->getNumGrid(0, FULL)+1; ++i) {
            fu(i, j) = vt(timeIdx, i, j)*u(oldtimeIdx, i, j);
            fv(i, j) = vt(timeIdx, i, j)*v(oldtimeIdx, i, j)*cosLat[j];
        }
    }
    // normal grids
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            double dx1 = fu(i+1,j)-fu(i-1,j);
            double dy1 = fv(i,j+1)-fv(i,j-1);
            double dx2 = u(oldtimeIdx, i, j)*(vt(timeIdx, i+1, j)-vt(timeIdx, i-1, j));
            double dy2 = v(oldtimeIdx, i, j)*(vt(timeIdx, i, j+1)-vt(timeIdx, i, j-1))*cosLat[j];
            dvt(i, j) = 0.5*((dx1+dx2)*factorLon[j]+(dy1+dy2)*factorLat[j]);
        }
    }
}

/**
 *  Input: u, vt
 *  Output: dut
 */
void BarotropicModel_Analytic::calcZonalWindCoriolis(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            double f = factorCor[j]+u(oldtimeIdx, i, j)*factorCur[j];
            dut(i, j) -= f*vt(timeIdx, i, j);
        }
    }
}

/**
 *  Input: u, ut
 *  Output: dvt
 */
void BarotropicModel_Analytic::calcMeridionalWindCoriolis(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            double f = factorCor[j]+u(oldtimeIdx, i, j)*factorCur[j];
            dvt(i, j) += f*ut(timeIdx, i, j);
        }
    }
}

/*
 *  Input: gd,  ght
 *  Output: dut
 */
void BarotropicModel_Analytic::calcZonalWindPressureGradient(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            dut(i, j) += (gd(timeIdx, i+1, j)-gd(timeIdx, i-1, j))*
                         factorLon[j]*ght(oldtimeIdx, i, j);
        }
    }
}

/*
 *  Input: gd,  ght
 *  Output: dut
 */
void BarotropicModel_Analytic::calcZonalWindPressureGradientL2(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            dut(i, j) += (gd(timeIdx, i+1, j)-gd(timeIdx, i-1, j))*
                         factorLon[j]*(ght(oldtimeIdx, i, j)-ghts[j]);
        }
    }
}
/*
 *  Input: gd,  ght
 *  Output: dvt
 */
void BarotropicModel_Analytic::calcMeridionalWindPressureGradient(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            dvt(i, j) += (gd(timeIdx, i, j+1)-gd(timeIdx, i, j-1))*
                         factorLat[j]*cosLat[j]*ght(oldtimeIdx, i, j);
        }
    }
}

void BarotropicModel_Analytic::calcL1(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    int ilen = mesh->getNumGrid(0, FULL);
    int iw, ie;
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
		ie = (i+mLon[j]+ilen)%ilen; 
		iw = (i-mLon[j]+ilen)%ilen; 
                L1U( i, j) = 0.5*(2.0*ut(oldtimeIdx, i,j)-ut(oldtimeIdx, ie, j)-ut(oldtimeIdx, iw, j))+0.5*(gd(oldtimeIdx, ie, j)-gd(oldtimeIdx, iw, j));
                L1V( i, j) = 0.0;
                L1H( i, j) = 0.5*(2.0*gd(oldtimeIdx, i,j)-gd(oldtimeIdx, ie, j)-gd(oldtimeIdx, iw, j))+0.5*(ut(oldtimeIdx, ie, j)-ut(oldtimeIdx, iw, j));
        }
    }
    // Polar 
    int jn = mesh->getNumGrid(1, FULL)-1;
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                L1U(i,0 ) = 0.0; 
                L1V(i,0 ) = 0.0;
                L1H(i,0 ) = 0.0;
                L1U(i,jn ) = 0.0; 
                L1V(i,jn ) = 0.0;
                L1H(i,jn ) = 0.0;
	}

}
void BarotropicModel_Analytic::semiAnalytic(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &halfTimeIdx, const TimeLevelIndex &tmpTimeIdx, double dt) {
	// Compute L_2F 
        // update geopotential height
        calcGeopotentialDepthTendencyL2(oldTimeIdx, tmpTimeIdx);
        calcZonalWindTendencyL2(oldTimeIdx, tmpTimeIdx);
	calcMeridionalWindTendency(oldTimeIdx, tmpTimeIdx);
	// compute L_1 \hat{P} ^{n+1/2}
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                gd(halfTimeIdx, i, j) = gd(oldTimeIdx, i, j)-dt/2*dgd(i, j);
                ut(halfTimeIdx, i, j) = ut(oldTimeIdx, i, j)-dt/2*dut(i, j);
                vt(halfTimeIdx, i, j) = vt(oldTimeIdx, i, j)-dt/2*dvt(i, j);
            }
	}

        gd.applyBndCond(halfTimeIdx);
        vt.applyBndCond(halfTimeIdx);
	ut.applyBndCond(halfTimeIdx);
	cout<<"Calculate L1"<<endl;
	calcL1(halfTimeIdx, tmpTimeIdx);

	// update P^{n+1/2} =F -tau/2  (L_1 \hat(P)^{n+1/2} + L_2(P)^{n+1/2}  .or. L_2F^n
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                gd(halfTimeIdx, i, j) = gd(oldTimeIdx, i, j)-0.5*(L1H(i,j) + dt*dgd(i, j));
                ut(halfTimeIdx, i, j) = ut(oldTimeIdx, i, j)-0.5*(L1U(i,j) + dt*dut(i, j));
                vt(halfTimeIdx, i, j) = vt(oldTimeIdx, i, j)-0.5*(L1V(i,j) + dt*dvt(i, j));
            }
        }
        gd.applyBndCond(halfTimeIdx);
        vt.applyBndCond(halfTimeIdx);
	ut.applyBndCond(halfTimeIdx);
        
}
void BarotropicModel_Analytic::check_antisym(const TimeLevelIndex &TimeIdx, const TimeLevelIndex &testTimeIdx)
{

    double eTotal = 0.0;
    double eZonalAdv = 0.0;
    double eMeridAdv= 0.0;
    double eCor = 0.0;
    double eZonalCor = 0.0;
    double eMeridCor= 0.0;
    double ePres = 0.0;
    double eZonalPres = 0.0;
    double eMeridPres = 0.0;
    // (T(U) , U) ; (T(V),V)
    calcZonalWindAdvection(TimeIdx, testTimeIdx);
    calcMeridionalWindAdvection(TimeIdx,testTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                eZonalAdv += ut(testTimeIdx, i, j)*dut(i, j)*cosLat[j];
		dut(i,j) = 0.0;
                eMeridAdv += vt(testTimeIdx, i, j)*dvt(i, j)*cosLat[j];
		dvt(i,j) = 0.0;
            }
	}
    // (-f*V, U)+(fU,V) 
    calcZonalWindCoriolis(TimeIdx,testTimeIdx);
    calcMeridionalWindCoriolis(TimeIdx,testTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                eCor += ut(testTimeIdx, i, j)*dut(i, j)*cosLat[j];
        	dut(i,j) = 0.0;
                eCor += vt(testTimeIdx, i, j)*dvt(i, j)*cosLat[j];
        	dvt(i,j) = 0.0;
            }
            }
    // (RU,U) + (R*phi, phi) ; (RV, V) + (R*phi, phi) 
    calcZonalWindPressureGradient(TimeIdx,testTimeIdx);
    calcMeridionalWindPressureGradient(TimeIdx,testTimeIdx);
    calcGeopotentialDepthTendency(TimeIdx, testTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                ePres += ut(testTimeIdx, i, j)*dut(i, j)*cosLat[j];
                ePres += vt(testTimeIdx, i, j)*dvt(i, j)*cosLat[j];
		ePres += gd(testTimeIdx, i, j)*dgd(i, j)*cosLat[j];
		dut(i,j) = 0.0;
		dvt(i,j) = 0.0;
		dgd(i,j) = 0.0;
            }
	}
    calcGeopotentialDepthTendency(TimeIdx, testTimeIdx);
    calcZonalWindTendency(TimeIdx, testTimeIdx);
    calcMeridionalWindTendency(TimeIdx, testTimeIdx);
    for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            eTotal += gd(testTimeIdx, i, j)*dgd(i, j)*cosLat[j];
            eTotal += ut(testTimeIdx, i, j)*dut(i, j)*cosLat[j];
            eTotal += vt(testTimeIdx, i, j)*dvt(i, j)*cosLat[j];
	    dut(i,j) = 0.0;
	    dvt(i,j) = 0.0;
	    dgd(i,j) = 0.0;
        }
    }
    cout << "(AF,F) :";
    cout << std::fixed << setw(20) << setprecision(12) << eTotal<<"  ";
    cout << std::fixed << setw(20) << setprecision(12) << eZonalAdv<<"  ";
    cout << std::fixed << setw(20) << setprecision(12) << eMeridAdv<<"  ";
    cout << std::fixed << setw(20) << setprecision(12) << eCor<<"  ";
    cout << std::fixed << setw(20) << setprecision(12) << ePres<<"  ";
    cout << endl;

    // check for semi-implicit 
    double eZonalL1 = 0.0;
    double eZonalL2 = 0.0;
    double ePresL1 = 0.0;
    double ePresL2 = 0.0;
    
    double tmpPres = 0.0;
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
	    tmpPres = 0.0; 
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                tmpPres += ut(testTimeIdx, i, j)*dut(i, j)*cosLat[j];
                ePresL1 += ut(testTimeIdx, i, j)*dut(i, j)*cosLat[j];
		ePresL1 += gd(testTimeIdx, i, j)*dgd(i, j)*cosLat[j];
		tmpPres += gd(testTimeIdx, i, j)*dgd(i, j)*cosLat[j];
		dut(i,j) = 0.0;
		dgd(i,j) = 0.0;
            }
//	    cout<<"J= "<< j<<setw(10)<<"  "<< tmpPres<<endl;
	}
    calcZonalWindPressureGradientL2(TimeIdx,testTimeIdx);
    calcMeridionalWindPressureGradient(TimeIdx,testTimeIdx);
    calcGeopotentialDepthTendencyL2(TimeIdx, testTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                ePresL2 += ut(testTimeIdx, i, j)*dut(i, j)*cosLat[j];
                ePresL2 += vt(testTimeIdx, i, j)*dvt(i, j)*cosLat[j];
		ePresL2 += gd(testTimeIdx, i, j)*dgd(i, j)*cosLat[j];
		dut(i,j) = 0.0;
		dvt(i,j) = 0.0;
		dgd(i,j) = 0.0;
            }
	}
    cout << "(L_1F,F) :";
    cout << std::fixed << setw(20) << setprecision(12) << ePresL1<<"  ";
    cout << std::fixed << setw(20) << setprecision(12) << ePresL2<<"  ";
    cout << endl;

}

}
