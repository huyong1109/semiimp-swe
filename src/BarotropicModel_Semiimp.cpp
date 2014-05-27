#include "BarotropicModel_Semiimp.h"

namespace barotropic_model {

BarotropicModel_Semiimp::BarotropicModel_Semiimp() {
    REPORT_ONLINE;
}

BarotropicModel_Semiimp::~BarotropicModel_Semiimp() {
    REPORT_OFFLINE;
}

void BarotropicModel_Semiimp::init(int numLon, int numLat) {
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
    u.create("u", "m s-1", "zonal wind speed", *mesh, CENTER);
    v.create("v", "m s-1", "meridional wind speed", *mesh, CENTER);
    gd.create("gd", "m2 s-2", "geopotential depth", *mesh, CENTER, HAS_HALF_LEVEL);

    ut.create("ut", "(m s-1)*m-2", "transformed zonal wind speed", *mesh, CENTER, HAS_HALF_LEVEL);
    vt.create("vt", "(m s-1)*m-2", "transformed meridional wind speed", *mesh, CENTER, HAS_HALF_LEVEL);
    ght.create("ght", "m-2", "transformed geopotential height", *mesh, CENTER);

    dut.create("dut", "m s-2", "zonal wind speed tendency", *mesh, CENTER);
    dvt.create("dvt", "m s-2", "meridional zonal speed tendency", *mesh, CENTER);
    dgd.create("dgd", "m-2 s-1", "geopotential depth tendency", *mesh, CENTER);

    ghu.create("ghu", "m2 s-1", "ut * ght", *mesh, CENTER);
    ghv.create("ghv", "m2 s-1", "vt * ght", *mesh, CENTER);

    fu.create("fu", "* m s-1", "* * u", *mesh, CENTER);
    fv.create("fv", "* m s-1", "* * v", *mesh, CENTER);
    lp.create("lp", "* m s-1", "* * v", *mesh, CENTER);
    LU.create("lp", "* m s-1", "* * v", *mesh, CENTER);
    LV.create("lp", "* m s-1", "* * v", *mesh, CENTER);
    LH.create("lp", "* m s-1", "* * v", *mesh, CENTER);
    //a.create("a", "Tri-diagonal a", *mesh, CENTER);
    //b.create("b", "Tri-diagonal b", *mesh, CENTER);
    //c.create("c", "Tri-diagonal c", *mesh, CENTER);
    //d.create("d", "Tri-diagonal d", *mesh, CENTER);
    //e.create("e", "Tri-diagonal e", *mesh, CENTER);
    //r.create("r", "Tri-diagonal r", *mesh, CENTER);
    //x.create("x", "Tri-diagonal x", *mesh, CENTER);
    // -------------------------------------------------------------------------
    // set coefficients
    // Note: Some coefficients containing cos(lat) will be specialized at Poles.
    for (int i = 0; i<5; i++){
	a[i].set_size(mesh->getNumGrid(0,FULL));
	//tmpa[i].set_size(mesh->getNumGrid(0,FULL)/2);
    }
	b.set_size(mesh->getNumGrid(0,FULL));
	//tmpb.set_size(mesh->getNumGrid(0,FULL)/2);
	x.set_size(mesh->getNumGrid(0,FULL));




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
    // -------------------------------------------------------------------------
    // set variables in Poles
    for (int i = -1; i < mesh->getNumGrid(0, FULL)+1; ++i) {
        dut(i, js) = 0.0; dut(i, jn) = 0.0;
        dvt(i, js) = 0.0; dvt(i, jn) = 0.0;
        ghu(i, js) = 0.0; ghu(i, jn) = 0.0;
        ghv(i, js) = 0.0; ghv(i, jn) = 0.0;
    }
}

void BarotropicModel_Semiimp::run(TimeManager &timeManager) {
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

void BarotropicModel_Semiimp::integrate(const TimeLevelIndex &oldTimeIdx,
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
	    //cout<<"lat = " <<j<< ": " <<ghts[j]<<endl;
        }
	//two half semi-implicit step 
	semiimplicit(oldTimeIdx, halfTimeIdx, oldTimeIdx,  dt);	
	//check_antisym(oldTimeIdx, halfTimeIdx);	
	semiimplicit(oldTimeIdx, newTimeIdx, halfTimeIdx,  dt);	
	
	//check_antisym(oldTimeIdx, newTimeIdx);	
	//Caculate \beta  = lqlp/lq2;
	long double beta = 0.0;
	long double lqlp = 0.0;
	long double lq2 = 0.0;

	// H 
        calcGeopotentialDepthTendency(oldTimeIdx, newTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                LH(i, j) = dgd(i, j);
		lq2 += LH(i,j)*LH(i,j)*cosLat[j];
            }
	}
	
        calcGeopotentialDepthTendencyL1( newTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                lp(i, j) = dgd(i, j);
            }
	}
        calcGeopotentialDepthTendencyL2(oldTimeIdx, halfTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                lp(i, j) += dgd(i, j);
		lqlp += LH(i,j)*lp(i,j)*cosLat[j];
            }
	}
	// U
        calcZonalWindTendency(oldTimeIdx, newTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                LU(i, j) = dut(i, j);
		lq2 += LU(i,j)*LU(i,j)*cosLat[j];
            }
	}
        calcZonalWindTendencyL1( newTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                lp(i, j) = dut(i, j);
            }
	}
        calcZonalWindTendencyL2(oldTimeIdx, halfTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                lp(i, j) += dut(i, j);
		lqlp += LU(i,j)*lp(i,j)*cosLat[j];
            }
	}
	// V
        calcMeridionalWindTendency(oldTimeIdx, newTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                LV(i, j) = dvt(i, j);
		lq2 += LV(i,j)*LV(i,j)*cosLat[j];
            }
	}
        calcMeridionalWindTendency(oldTimeIdx, halfTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                lp(i, j) = dvt(i, j);
		lqlp += LV(i,j)*lp(i,j)*cosLat[j];
            }
	}
	beta = lqlp/lq2;
	cout<<"beta = " <<setprecision(15)<<beta<<endl;
//update F^{n+1}
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                ut(newTimeIdx, i, j) = ut(oldTimeIdx, i, j)-dt*beta*LU(i,j);
                vt(newTimeIdx, i, j) = vt(oldTimeIdx, i, j)-dt*beta*LV(i,j);
                gd(newTimeIdx, i, j) = gd(oldTimeIdx, i, j)-dt*beta*LH(i,j);
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

double BarotropicModel_Semiimp::calcTotalEnergy(const TimeLevelIndex &timeIdx) const {
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

double BarotropicModel_Semiimp::calcTotalMass(const TimeLevelIndex &timeIdx) const {
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
void BarotropicModel_Semiimp::calcGeopotentialDepthTendency(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
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
void BarotropicModel_Semiimp::calcGeopotentialDepthTendencyL1(const TimeLevelIndex &timeIdx) {
    // calculate intermediate variables
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = -1; i < mesh->getNumGrid(0, FULL)+1; ++i) {
            ghu(i, j) = ut(timeIdx, i, j)*ghts[j];
        }
    }
    // normal grids
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            dgd(i, j) = (ghu(i+1, j)-ghu(i-1, j))*factorLon[j];
        }
    }
    // polar grids
	int jn= mesh->getNumGrid(1, FULL)-1;
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            dgd(i, 0) = 0.0;
            dgd(i, jn) = 0.0;
	}
}
/**
 *  Input: ut, vt, ght
 *  Intermediate: ghu, ghv
 *  Output: dgd
 */
void BarotropicModel_Semiimp::calcGeopotentialDepthTendencyL2(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
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

void BarotropicModel_Semiimp::calcZonalWindTendencyL1(const TimeLevelIndex &timeIdx) {
    calcZonalWindPressureGradientL1(timeIdx);
}

void BarotropicModel_Semiimp::calcZonalWindTendencyL2(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    calcZonalWindAdvection(oldtimeIdx, timeIdx);
    calcZonalWindCoriolis(oldtimeIdx,timeIdx);
    calcZonalWindPressureGradientL2(oldtimeIdx,timeIdx);
}

void BarotropicModel_Semiimp::calcZonalWindTendency(const TimeLevelIndex &oldtimeIdx,const TimeLevelIndex &timeIdx) {
    calcZonalWindAdvection(oldtimeIdx,timeIdx);
    calcZonalWindCoriolis(oldtimeIdx,timeIdx);
    calcZonalWindPressureGradient(oldtimeIdx,timeIdx);
}
void BarotropicModel_Semiimp::calcMeridionalWindTendency(const TimeLevelIndex &oldtimeIdx,const TimeLevelIndex &timeIdx) {
    calcMeridionalWindAdvection(oldtimeIdx,timeIdx);
    calcMeridionalWindCoriolis(oldtimeIdx,timeIdx);
    calcMeridionalWindPressureGradient(oldtimeIdx,timeIdx);
}

/**
 *  Input: u, v, ut
 *  Output, s1
 */
void BarotropicModel_Semiimp::calcZonalWindAdvection(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
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
void BarotropicModel_Semiimp::calcMeridionalWindAdvection(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
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
void BarotropicModel_Semiimp::calcZonalWindCoriolis(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
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
void BarotropicModel_Semiimp::calcMeridionalWindCoriolis(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
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
void BarotropicModel_Semiimp::calcZonalWindPressureGradient(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
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
void BarotropicModel_Semiimp::calcZonalWindPressureGradientL1(const TimeLevelIndex &timeIdx) {
    
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            dut(i, j) = (gd(timeIdx, i+1, j)-gd(timeIdx, i-1, j))*
                         factorLon[j]*ghts[j];
        }
    }
}
/*
 *  Input: gd,  ght
 *  Output: dut
 */
void BarotropicModel_Semiimp::calcZonalWindPressureGradientL2(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
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
void BarotropicModel_Semiimp::calcMeridionalWindPressureGradient(const TimeLevelIndex &oldtimeIdx, const TimeLevelIndex &timeIdx) {
    for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            dvt(i, j) += (gd(timeIdx, i, j+1)-gd(timeIdx, i, j-1))*
                         factorLat[j]*cosLat[j]*ght(oldtimeIdx, i, j);
        }
    }
}
vec BarotropicModel_Semiimp::Gauss(vec *a, vec b,  int ie) {
	vec x = b;
	double y;
	int iee = ie-1;
	a[3][0]= a[2][iee];
	a[4][0]= a[0][0];
	a[4][iee-1]= a[2][iee-1];
	a[3][iee-1]= a[0][iee];
        for (int i = 1; i < ie-1; ++i) {
	   y = a[0][i]/a[1][i-1];
	   a[0][i] = 0.0;
	   a[1][i] = a[1][i]-a[2][i-1]*y;
	   a[4][i] = a[4][i] -a[4][i-1]*y;
	   b[i]= b[i]- b[i-1]*y;

	   y = a[3][i-1]/a[1][i-1];
	   a[3][i] = a[3][i] -a[2][i-1]*y;
	   a[1][iee]= a[1][iee] - a[4][i-1]*y;
	   b[iee] = b[iee]-b[i-1]*y;
	
        }
    
	a[1][iee]= a[1][iee] - (a[4][iee-1])*(a[3][iee-1])/a[1][iee-1];
	b[iee]= b[iee] - b[iee-1]*(a[3][iee-1])/a[1][iee-1];

	x[iee] = b[iee]/a[1][iee];
	

	x[iee-1] = (b[iee-1]-a[4][iee-1]*x[iee])/a[1][iee-1];
		
	for (int i = ie-3; i>=0; --i) {
		x[i]= (b[i]-a[2][i]*x[i+1]-a[3][i]*x[iee])/a[1][i];
	}
	return x;
}

vec BarotropicModel_Semiimp::Gaussreorder(vec *a,  vec b,  int ie) {
		    
	vec tmpa[5], tmpb;
	for (int i = 0; i<5; i++){
    	    tmpa[i].set_size(ie/2);
    	}
	    tmpb.set_size(ie/2);
	vec x = b;
	for (int i = 0; i < ie/2; ++i) {
	    for(int j =0; j < 5; ++j){
	//	cout<<" i, j = "<<i<<j<<endl;
	//	 cout<<"tmpa "<<tmpa[j][i]<<endl;
	//	 cout<<"a "<<a[j][i*2]<<endl;
		 tmpa[j][i] = a[j][i*2];
		//cout<<"b "<<endl;
	    }
	   // cout<<"bbb ";
	   // cout<<"b "<<tmpb[i]<<endl;
	    tmpb[i] = b[i*2];
	}
	//cout<< "gauss1"<<endl;
	tmpb = Gauss(tmpa, tmpb, ie/2);
	for (int i = 0; i < ie/2; ++i) {
	    x[i*2] = tmpb[i];
	}
	for (int i = 0; i < ie/2; ++i) {
	    for(int j =0; j < 5; ++j){
		 tmpa[j][i] = a[j][i*2+1];
	    }
	    tmpb[i] = b[i*2+1];
	}
	//cout<< "gauss2"<<endl;
	tmpb = Gauss(tmpa, tmpb, ie/2);
	for (int i = 0; i < ie/2; ++i) {
	    x[i*2+1] = tmpb[i];
	}

	//cout<< "end gauss2"<<endl;
	return x;
}


//vec BarotropicModel_Semiimp::Gaussreorder(vec *a, vec b,  int ie) {
//	//==============reorder 
//	cout<<"\n"<<"fucn ie = "<<ie<<"\n";
//	for (int i = 0; i < ie; ++i) {
//	    cout<<b[i]<<"\t";
//	}
//	vec tmp = b;
//	for (int i = 1; i <= ie/2; ++i) {
//	    b[i] = tmp[i*2];
//	    b[i+ie/2] = tmp[i*2-1];
//	}
//
//	for(int j =0; j < 5; ++j){
//	    tmp = a[j];
//	    for (int i = 1; i <= ie/2; ++i) {
//	        a[j][i] = tmp[i*2];
//	        a[j][i+ie/2] = tmp[i*2-1];
//	    }
//	}
//	cout<<endl<<"-----------reorder-----------------"<<endl;
//	for (int i = 0; i < ie; ++i) {
//	    cout<<b[i]<<"\t";
//	}
//	cout<<endl<<"-----------done-----------------"<<endl;
//
//
//	vec x = b;
//	double y;
//	int iee = ie-1;
//	a[3][0]= a[2][iee];
//	a[4][0]= a[0][0];
//	a[4][iee-1]= a[2][iee-1];
//	a[3][iee-1]= a[0][iee];
//        for (int i = 1; i < ie-1; ++i) {
//	   y = a[0][i]/a[1][i-1];
//	   a[0][i] = 0.0;
//	   a[1][i] = a[1][i]-a[2][i-1]*y;
//	   a[4][i] = a[4][i] -a[4][i-1]*y;
//	   b[i]= b[i]- b[i-1]*y;
//
//	   y = a[3][i-1]/a[1][i-1];
//	   a[3][i] = a[3][i] -a[2][i-1]*y;
//	   a[1][iee]= a[1][iee] - a[4][i-1]*y;
//	   b[iee] = b[iee]-b[i-1]*y;
//	
//        }
//    
//	a[1][iee]= a[1][iee] - (a[4][iee-1])*(a[3][iee-1])/a[1][iee-1];
//	b[iee]= b[iee] - b[iee-1]*(a[3][iee-1])/a[1][iee-1];
//
//	x[iee] = b[iee]/a[1][iee];
//	
//
//	x[iee-1] = (b[iee-1]-a[4][iee-1]*x[iee])/a[1][iee-1];
//		
//	for (int i = ie-3; i>=0; --i) {
//		x[i]= (b[i]-a[2][i]*x[i+1]-a[3][i]*x[iee])/a[1][i];
//	}
//
//    	for (int i = 0; i < ie; ++i) {
//	    tmp[i] = x[i];
//		}
//	for (int i = 1; i <= ie/2; ++i) {
//	    x[i*2] = tmp[i];
//	    x[i*2-1] = tmp[i+ie/2];
//		}
//	return x;
//}
void BarotropicModel_Semiimp::semiimplicit(const TimeLevelIndex &oldTimeIdx, const TimeLevelIndex &halfTimeIdx, const TimeLevelIndex &tmpTimeIdx, double dt) {
	// Compute L_2F 
        // update geopotential height
	
	double tempfactor;

        calcGeopotentialDepthTendencyL2(oldTimeIdx, tmpTimeIdx);
        calcZonalWindTendencyL2(oldTimeIdx, tmpTimeIdx);
	calcMeridionalWindTendency(oldTimeIdx, tmpTimeIdx);
        for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
		tempfactor = dt*factorLon[j]*ghts[j];
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                gd(halfTimeIdx, i, j) = gd(oldTimeIdx, i, j)-dt/2*dgd(i, j);
                ut(halfTimeIdx, i, j) = ut(oldTimeIdx, i, j)-dt/2*dut(i, j);
                vt(halfTimeIdx, i, j) = vt(oldTimeIdx, i, j)-dt/2*dvt(i, j);
            }
	}
        gd.applyBndCond(halfTimeIdx);
        vt.applyBndCond(halfTimeIdx);
//========================  Reorder  ============================
//        for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
//        	tempfactor = dt/2.0*factorLon[j]*ghts[j];
//            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
//        	a[0][i] = -tempfactor*tempfactor;
//        	a[1][i] = 2.0*tempfactor*tempfactor+1.0;
//        	a[2][i] = -tempfactor*tempfactor;
//        	a[3][i] = 0.0;
//        	a[4][i] = 0.0;
//        	x[i]= 0.0;
//        	b[i] = -tempfactor*(gd(halfTimeIdx, i+1, j)-gd(halfTimeIdx, i-1, j))+ut(halfTimeIdx, i, j);
//        	//Test case for gauss and gaussreorder
//        	//b[i] = a[0][i]*((double)((i-2+80)%80)) +a[1][i]*((double)i)+a[2][i]*((double(((i+2+80)%80))));
//            }
//            
//            x= Gaussreorder(a,b, mesh->getNumGrid(0, FULL));
//           // cout<<endl<<"||---------------------------||"<<endl;
//           // cout<<"after gauss"<<endl;
//            for (int i = 0; i<mesh->getNumGrid(0, FULL); ++i)
//            {
//        // update U		
//              ut(halfTimeIdx, i, j) = x[i];
//        	//cout<<x[i]<<"\t";
//            }
//           //cout<<endl<<"||---------------------------||"<<endl;
//        }

//========================  End reorder  ============================
	
//========================  Middle Point  ============================
        for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
		tempfactor = dt/2.0*factorLon[j]*ghts[j];
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
		a[0][i] = -4.0*tempfactor*tempfactor;
		a[1][i] = 8.0*tempfactor*tempfactor+1.0;
		a[2][i] = -4.0*tempfactor*tempfactor;
		a[3][i] = 0.0;
		a[4][i] = 0.0;
		x[i]= 0.0;
		//b for middle point method
		b[i] = -tempfactor*(gd(halfTimeIdx, i+1, j)-gd(halfTimeIdx, i-1, j))+ut(halfTimeIdx, i, j);
            }
	    
	    x= Gauss(a, b, mesh->getNumGrid(0, FULL));
	   // cout<<endl<<"||---------------------------||"<<endl;
	   // cout<<"after gauss"<<endl;
	    for (int i = 0; i<mesh->getNumGrid(0, FULL); ++i)
	    {
	// update U		
              ut(halfTimeIdx, i, j) = x[i];
		//cout<<x[i]<<"\t";
	    }
	   //cout<<endl<<"||---------------------------||"<<endl;
        }

	//======================== End  Middle Point  ============================
        
	ut.applyBndCond(halfTimeIdx);
	// update P =F +L_2(F) + L_1(P) 
	calcGeopotentialDepthTendencyL1(halfTimeIdx);
        for (int j = 1; j < mesh->getNumGrid(1, FULL)-1; ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
                gd(halfTimeIdx, i, j) = gd(halfTimeIdx, i, j)-dt/2*dgd(i,j);
            }
        }
        gd.applyBndCond(halfTimeIdx);
        
}
void BarotropicModel_Semiimp::check_antisym(const TimeLevelIndex &TimeIdx, const TimeLevelIndex &testTimeIdx)
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
    
    calcZonalWindPressureGradientL1(testTimeIdx);
    calcGeopotentialDepthTendencyL1(testTimeIdx);
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
