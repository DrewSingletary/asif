#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sstream>

#include <asif++.h>
#include "CyberTimer.hpp"

#define STANDARD // variable alpha
// #define ROBUST // variable p
// #define ROBUST_2 // variable pBounds
// #define ROBUST_3 // variable npSSmax

using namespace std;
using namespace ASIF;

const uint32_t nx = 2;
const uint32_t nu = 1;

const double lb[nu] = {-1.5};
const double ub[nu] = {1.5};
const double xBound[2] = {-M_PI,M_PI};
const double vBound[2] = {-M_PI,M_PI};

double p = 1.0;

#ifdef STANDARD
const double pMin = 1.0;
const double pMax = 1.0;
#endif

#ifdef ROBUST
const double pMin = 0.8;
const double pMax = 1.2;
#endif

#ifdef ROBUST_2
double pMin = 0.8;
double pMax = 1.2;
#endif

#ifdef ROBUST_3
const double pMin = 0.8;
const double pMax = 1.2;
#endif

ASIFrobust *asif;
vector<vector<double>> SafetySetData;

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
	f[1] = sin(x[0]);

	g[0] = 0.;
	g[1] = interval(pMin,pMax);
}

void dynamicsExact(const double *x, double *f, double *g)
{
	f[0] = x[1];
	f[1] = sin(x[0]);

	g[0] = 0.;
	g[1] = p;
}

int main()
{
	CyberTimer<500> timer;

	for(uint32_t i=0; i<3; i++)
	{
		#ifdef ROBUST
		p = pMin + static_cast<double>(i)*(pMax-pMin)/2.;
		#endif
		#ifdef ROBUST_2
		pMin = 1.0 - static_cast<double>(i)*(0.2)/2.;
		pMax = 1.0 + static_cast<double>(i)*(0.2)/2.;
		#endif

		cout << "File " << i+1 << "/3" << endl;
		ofstream myfile;
		#ifdef STANDARD
		myfile.open(string("asif_standard_")+to_string(i)+string(".csv"), ofstream::out | ofstream::trunc);
		#endif
		#ifdef ROBUST
		myfile.open(string("asif_robust_")+to_string(i)+string(".csv"), ofstream::out | ofstream::trunc);
		#endif
		#ifdef ROBUST_2
		myfile.open(string("asif_robust2_")+to_string(i)+string(".csv"), ofstream::out | ofstream::trunc);
		#endif
		#ifdef ROBUST_3
		myfile.open(string("asif_robust3_")+to_string(i)+string(".csv"), ofstream::out | ofstream::trunc);
		#endif
		myfile << "tNow" <<','<< "x" <<','<<  "v" <<',' <<  "uDes" <<','<<  "uAct" <<','<< "relax" <<','<< "rc" << endl;
		myfile << std::fixed;
		myfile << std::setprecision(10);

		ASIFrobust::Options opts;
		

		#ifdef STANDARD
		opts.relaxCost = 1.0;
		opts.relaxLb = .5 + static_cast<double>(i)*(4.5)/2.;
		#else
		opts.relaxCost = 50.0;
		opts.relaxLb = 5.0;
		#endif

		#ifdef ROBUST_3
		const uint32_t npSSmaxTmp = round(2.+static_cast<double>(i)*(static_cast<double>(SafetySetData.size())-2.)/2.);
		asif = new ASIFrobust(nx,nu,SafetySetData.size(),safetySet,dynamics,npSSmaxTmp);
		cout << "npSSmax:  " << npSSmaxTmp << endl;
		#else
		asif = new ASIFrobust(nx,nu,SafetySetData.size(),safetySet,dynamics);
		#endif

		asif->initialize(lb,ub,opts);

		const double tEnd = 6.0;
		const double dt = 0.01;
		double xNow[2] = {0.5,0.0};
		double uDesNow[1] = {0.0};
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
			// cout << "t: " << tNow << "/" << tEnd << ", uAct: " << uActNow[0] << ", rc:" << rc << endl;
		}

		cout << "dt average:  " << timer.getAverage() << endl << endl;

		myfile.close();
	}
	cout << "Finished" << endl;
	return 0;
}

