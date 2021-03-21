#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <asif++.h>
#include "CyberTimer.hpp"


using namespace std;

const uint32_t nx = 2;
const uint32_t nu = 1;
const uint32_t npSS = 4;
const double lb[nu] = {-1.0};
const double ub[nu] = {1.0};
const double xBound[2] = {-1.0,1.0};
const double vBound[2] = {-1.0,1.0};

CyberTimer<1> timer;

ASIF::ASIF *asif;

void safetySet(const double *x, double *h, double *Dh)
{
	if(x[1]>0)
	{
		h[0] = xBound[1] - x[0] - (x[1]*x[1])/2.0; 	Dh[0] = -1.0; Dh[4] = -x[1];
		h[1] = x[0] - xBound[0]; 	                Dh[1] = 1.0; Dh[5] = 0.0; 
	}
	else
	{
		h[0] = -x[0] + xBound[1];                   Dh[0] = -1.0; Dh[4] = 0.0;
		h[1] = x[0] - xBound[0] - (x[1]*x[1])/2.0;  Dh[1] =  1.0; Dh[5] = -x[1];
	}
	h[2] = x[1] - vBound[0];  Dh[2] = 0.0; Dh[6] = 1.0;
	h[3] = -x[1] + vBound[1]; Dh[3] = 0.0; Dh[7] = -1.0;
}

void dynamics(const double *x, double *f, double *g)
{
	const double A[nx*nx] = {0.0, 0.0, 1.0, 0.0};
	const double b[nx*nu] = {0.0, 1.0};

	for(uint32_t i=0;i<nx;i++)
	{		
		f[i] = 0.0;
		for(uint32_t j=0;j<nx;j++)
		{
			f[i] += A[i+(j*nx)]*x[j];
		}
	}

	for(uint32_t i=0;i<nx;i++)
	{		
		for(uint32_t j=0;j<nu;j++)
		{
			g[i+(j*nx)] = b[i+(j*nx)];
		}
	}
}

int main()
{
	ASIF::ASIF::Options opts;
	opts.relaxLb = 6.0;

	asif = new ASIF::ASIF(nx,nu,npSS,safetySet,dynamics);
	asif->initialize(lb,ub);

	const double tEnd = 5.0;
	const double dt = 0.001;
	double xNow[2] = {0.0,0.0};
	double uDesNow[1] = {1.0};
	double uActNow[1] = {0.0};
	double tNow = 0.0;
	double relax = 0.0;
	int32_t rc;

	bool optionsUpdated = false;
	while(tNow<tEnd)
	{
		if(!optionsUpdated && tNow>(tEnd/2))
		{
			optionsUpdated = true;
			asif->updateOptions(opts);
		}

		timer.tic();
		rc = asif->filter(xNow,uDesNow,uActNow,relax);
		timer.toc();

		if(rc!=1)
		{
			cout << "ASIF failed" << endl;
		}	

		// Integrate
		double fCl[nx] = {0.0};
		double f[nx];
		double g[nx];
		dynamics(xNow, f, g);
		for(uint32_t i = 0; i<nx; i++)
		{
			fCl[i]+=f[i];
			for(uint32_t k = 0; k<nu; k++)
			{
				fCl[i]+= uActNow[k]*g[i+k*nx];
			}
		}

		for(uint32_t i = 0; i<nx; i++)
		{
			xNow[i] += dt*fCl[i];
		}
		tNow += dt;

		// Print data
		uint32_t printAcc = 3U;
		cout << std::fixed;
		cout << std::setprecision(3);
		cout << "t: " << tNow ;
		cout << std::scientific;
		cout << std::setprecision(printAcc);
		cout << " | x:";
		cout << std::setw(printAcc+7) << xNow[0];
		cout << " | v:";
		cout << std::setw(printAcc+7) << xNow[1];
		cout <<  " | uDes:";
		cout << std::setw(printAcc+7) << uDesNow[0];
		cout << " | uAct:";
		cout << std::setw(printAcc+7) << uActNow[0];
		cout << " | gamma:";
		cout << std::setw(printAcc+7) << relax;
		cout << " | dt(us):";
		cout << std::fixed;
		cout << std::setprecision(1);
		cout << timer.getLatest()*1e6;
		cout << endl;
	}

	cout << "Finished" << endl;

	return 0;
}

