#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <asif++.h>
#include "customTimer.h"

#define DYNAMICS_WITH_GRADIENT true

using namespace std;

const uint32_t nx = 2;
const uint32_t nu = 1;
const uint32_t npSS = 4;
const uint32_t npBS = 1;
const uint32_t npBTSS = 4;

const double lb[nu] = {-1.0};
const double ub[nu] = {1.0};
const double xBound[2] = {-1.0,1.0};
const double vBound[2] = {-1.0,1.0};

const double A[nx*nx] = {0.0, 0.0, 1.0, 0.0};
const double B[nx*nu] = {0.0, 1.0};

const double K[nx] = {-10.0,-20};
const double P[nx*nx] = {0.500000000000000,0.288675134594813,0.288675134594813,0.5};
const double mPpPt[nx*nx] = {-1.0,-0.577350269189626,-0.577350269189626,1.0};
const double Pv = 0.002;

ASIF::ASIFimplicit *asif;

void safetySet(const double *x, double *h, double *Dh)
{
	h[0] = -x[0] + xBound[1]; Dh[0] = -1.0; Dh[4] =  0.0;
	h[1] =  x[0] - xBound[0]; Dh[1] =  1.0; Dh[5] =  0.0; 
	h[2] =  x[1] - vBound[0]; Dh[2] =  0.0; Dh[6] =  1.0;
	h[3] = -x[1] + vBound[1]; Dh[3] =  0.0; Dh[7] = -1.0;
}

void backupSet(const double *x, double *h, double *Dh)
{
	h[0] = Pv;
	for(uint32_t i = 0; i<nx; i++)
	{
		for(uint32_t j = 0; j<nx; j++)
		{
			h[0]-=P[i+j*nx]*x[i]*x[j];			
		}
	}
	ASIF::matrixVectorMultiply(mPpPt,nx,nx,
							   x,nx,
							   Dh);
}

void dynamics(const double *x, double *f, double *g)
{
	ASIF::matrixVectorMultiply(A,nx,nx,
							   x,nx,
							   f);
	memcpy(g,B,sizeof(B));
}

void backupController(const double *x, double *u, double *Du)
{

	ASIF::matrixVectorMultiply(K,nu,nx,
							   x,nx,
							   u);
	memcpy(Du,K,sizeof(K));
}

void dynamicsGradients(const double *x, double *Df, double *Dg)
{
	memcpy(Df,A,sizeof(A));
	for(uint32_t i=0; i<nx*nu*nx; i++)
		Dg[i]=0.0;
}

void dynamicsWithGradient(const double *x, const double *u, double *f, double *g, double *d_fcl_dx)
{
	dynamics(x,f,g);
	memcpy(d_fcl_dx,A,sizeof(A));
}

int main()
{

	ASIF::ASIFimplicit::Options opts;
	opts.backTrajHorizon = 2.0;
	opts.backTrajDt = 0.01;
	opts.relaxReachLb = 5.0;
	opts.relaxSafeLb = 10.0;

	#ifdef DYNAMICS_WITH_GRADIENT
	asif = new ASIF::ASIFimplicit(nx,nu,npSS,npBS,npBTSS,
											safetySet,backupSet,dynamicsWithGradient,backupController);
	#else
	asif = new ASIF::ASIFimplicit(nx,nu,npSS,npBS,npBTSS,
	 										safetySet,backupSet,dynamics,dynamicsGradients,backupController);
	#endif
	
	asif->initialize(lb,ub,opts);

	const double tEnd = 5.0;
	const double dt = 0.001;
	double xNow[2] = {0.0,0.0};
	double uDesNow[1] = {1.0};
	double uActNow[1] = {0.0};
	double tNow = 0.0;
	double relax[2];
	int32_t rc;

	bool optionsUpdated = false;
	CustomTimer timer;

	opts.backTrajHorizon = 5.0;
	while(tNow<tEnd)
	{
		if(!optionsUpdated && tNow>(tEnd/2))
		{
			optionsUpdated = true;
			asif->updateOptions(opts);
		}
	
		tic(&timer);
		rc = asif->filter(xNow,uDesNow,uActNow,relax);
		toc(&timer);

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
				fCl[i]+= g[i+k*nx]*uActNow[k];
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
		cout << " | gammaSafe:";
		cout << std::setw(printAcc+7) << relax[0];
		cout << " | gammaReach:";
		cout << std::setw(printAcc+7) << relax[1];
		cout << " | hMinIdx:";
		if(rc>0)
			cout << std::setw(printAcc+7) << asif->backTrajCritIdx_[0];
		else
			cout << std::setw(printAcc+7) << 0;
		cout << " | dt(us):";
		cout << std::fixed;
		cout << std::setprecision(1);
		cout << timer.dt*1e6;
		cout << endl;
	}

	cout << "Finished" << endl;

	return 0;
}

