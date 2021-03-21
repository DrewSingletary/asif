#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <asif++.h>

using namespace std;
using namespace ASIF;

const uint32_t nx = 2;
const uint32_t nu = 1;
const uint32_t npSS = 4;
const uint32_t npBS = 1;
const uint32_t npBTSS = 10;

const double lb[nu] = {-1.5};
const double ub[nu] = {1.5};
const double xBound[2] = {-M_PI,M_PI};
const double vBound[2] = {-M_PI,M_PI};

const double K[nx] = {-3.0,-3.0};
const double P[nx*nx] = {1.25,0.25,0.25,0.25};
const double mPpPt[nx*nx] = {-2.5,-0.5,-0.5,-0.5};
const double Pv = 0.05;

ASIFimplicit *asif;

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
	f[0] = x[1];
	f[1] = sin(x[0]);

	g[0] = 0.;
	g[1] = 1.;
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
	Df[0] = 0.;        Df[2] = 1.;
	Df[1] = cos(x[0]); Df[3] = 0.;

	for(uint32_t i=0; i<nx*nu*nx; i++)
		Dg[i]=0.0;
}

int main()
{
	for(uint32_t i=0; i<10; i++)
	{
		cout << "File " << i << "/10" << endl;
		ofstream myfile;
		myfile.open(string("asif_")+to_string(i)+string(".csv"), ofstream::out | ofstream::trunc);
		myfile << "tNow" <<','<< "x" <<','<<  "v" <<',' <<  "uDes" <<','<<  "uAct" <<','<< "gammaSafe" <<','<< "gammaReach" <<','<< "rc" << endl;
		myfile << std::fixed;
		myfile << std::setprecision(10);

		ASIFimplicit::Options opts;
		opts.backTrajHorizon = 5.0;
		opts.backTrajDt = 0.001;
		opts.relaxReachLb = 5.0;
		opts.relaxSafeLb = 10.0;

		asif = new ASIFimplicit(nx,nu,npSS,npBS,npBTSS,
			safetySet,backupSet,dynamics,dynamicsGradients,backupController);

		asif->initialize(lb,ub,opts);

		const double tEnd = 5.0;
		const double dt = 0.001;
		double xNow[2] = {0.1+double(i)*0.29,0.0};
		double uDesNow[1] = {0.0};
		double uActNow[1] = {0.0};
		double tNow = 0.0;
		double relax[2];
		int32_t rc;

		while(tNow<tEnd)
		{

			tNow += dt;
			rc = asif->filter(xNow,uDesNow,uActNow,relax);

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


			// Save data
			myfile << tNow << "," ;
			myfile << xNow[0] << "," ;
			myfile << xNow[1] << "," ;
			myfile << uDesNow[0] << "," ;
			myfile << uActNow[0] << "," ;
			myfile << relax[0] << "," ;
			myfile << relax[1] << "," ;
			myfile << rc << endl ;
		}

		myfile.close();
	}
	cout << "Finished" << endl;
	return 0;
}

