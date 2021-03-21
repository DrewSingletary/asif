#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sstream>
#include "CyberTimer.hpp"

#include <asif++.h>
// #include "KernelData_70-75kg.h"
#include "KernelData_70-135kg.h"

using namespace std;
using namespace ASIF;

const uint32_t nx = 2;
const uint32_t nu = 1;

const double lb[nu] = {-20};
const double ub[nu] = {20};

// const double m_max = 75.;
const double m_max = 135.;
const double m_min = 70.;
// const double m_mean = (m_min+m_max)/2.;
const double m_mean = 135.;
const interval_t mInt = interval(m_min,m_max);

const double K = 5.7;
const double dK = 0.1;
const interval_t KInt = interval(K-dK,K+dK);

const double F = 23;
const double DF = 2;
const double eps = 1.0e-3;
const interval_t FInt = interval(F-DF,F-DF);

ASIFrobust *asif;

void safetySet(const double *x, double *h, double *Dh)
{
	for(uint32_t i=0; i<SafetySetData.size(); i++)
	{
		h[i] = 1. - SafetySetData[i][0]*x[0] - SafetySetData[i][1]*x[1];
		Dh[i] = -SafetySetData[i][0];
		Dh[i+SafetySetData.size()] = -SafetySetData[i][1];
	}
}

void dynamics(const interval_t *x, interval_t *f, interval_t *g)
{
	f[0] = x[1];
	f[1] = -FInt*x[1]/mInt;

	g[0] = 0.;
	g[1] = KInt/mInt;
}

void dynamicsExact(const double *x, double *f, double *g)
{
	f[0] = x[1];
	f[1] = -F*x[1]/m_mean;

	g[0] = 0.;
	g[1] = K/m_mean;
}

int main()
{
	CyberTimer<500> timer;

	ofstream myfile;
	myfile.open("DI_asif_traj.csv", ofstream::out | ofstream::trunc);
	myfile << "tNow" <<','<< "x" <<','<<  "v" <<',' <<  "uDes" <<','<<  "uAct" <<','<< "relax" <<','<< "rc" << endl;
	myfile << std::fixed;
	myfile << std::setprecision(10);

	ASIFrobust::Options opts;

	opts.relaxCost = 50.0;
	opts.relaxLb = 5.0;

	asif = new ASIFrobust(nx,nu,SafetySetData.size(),safetySet,dynamics,5);

	asif->initialize(lb,ub,opts);

	const double tEnd = 6.0;
	const double dt = 0.01;
	double xNow[2] = {0.0,0.0};
	double uDesNow[1] = {20.0};
	double uActNow[1] = {0.0};
	double tNow = 0.0;
	double relax;
	int32_t rc;

	while(tNow<tEnd)
	{
		tNow += dt;
		timer.tic();
		rc = asif->filter(xNow,uDesNow,uActNow,relax);
		timer.toc();

		// Integrate
		double fCl[nx] = {0.0};
		double f[nx];
		double g[nx];
		dynamicsExact(xNow, f, g);
		for(uint32_t i = 0; i<nx; i++)
		{
			fCl[i]+=f[i];
			for(uint32_t k = 0; k<nu; k++)
			{
				fCl[i]+= g[i+k*nx]*uActNow[k];
			}
		}

		for(uint32_t i = 0; i<nx; i++)
		{
			xNow[i] += dt*fCl[i];
		}

		// Save data
		myfile << tNow << "," ;
		myfile << xNow[0] << "," ;
		myfile << xNow[1] << "," ;
		myfile << uDesNow[0] << "," ;
		myfile << uActNow[0] << "," ;
		myfile << relax << "," ;
		myfile << rc << endl ;
		cout << "t: " << tNow << "/" << tEnd << ", uAct: " << uActNow[0] << ", rc:" << rc << endl;
	}

	cout << "dt average:  " << timer.getAverage() << endl << endl;

	myfile.close();

	cout << "Finished" << endl;
	return 0;
}
