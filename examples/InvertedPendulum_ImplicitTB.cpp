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
const uint32_t npBTSS = 4;
const double uMax = 1.5;

const double lb[nu] = {-uMax};
const double ub[nu] = {uMax};
const double xBound[2] = {-M_PI/2.,M_PI};
const double vBound[2] = {-M_PI/2.,M_PI/2.};

ASIFimplicitTB *asifTB;

void safetySet(const double *x, double *h, double *Dh)
{
	h[0] = -x[0] + xBound[1]; Dh[0] = -1.0; Dh[4] =  0.0;
	h[1] =  x[0] - xBound[0]; Dh[1] =  1.0; Dh[5] =  0.0; 
	h[2] =  x[1] - vBound[0]; Dh[2] =  0.0; Dh[6] =  1.0;
	h[3] = -x[1] + vBound[1]; Dh[3] =  0.0; Dh[7] = -1.0;
}

void backupSetTB(const double *x, double *h, double *Dh, double*DDh)
{
	const double x0 = M_PI/2.;
	// const double xh = 0.1;
	// const double vh = M_PI/4.;

	// const double xh4 = xh*xh*xh*xh;
	// const double vh4 = vh*vh*vh*vh;
	// double tmpx = (x[0]-x0);
	// double tmpv = x[1];

	// h[0] = 1. - tmpx*tmpx*tmpx*tmpx/xh4 - tmpv*tmpv*tmpv*tmpv/vh4;

	// Dh[0] = -4.*tmpx*tmpx*tmpx/xh4;
	// Dh[1] = -4.*tmpv*tmpv*tmpv/vh4;

	// DDh[0] = -12.*tmpx*tmpx/xh4;
	// DDh[1] = 0.;
	// DDh[2] = 0.;
	// DDh[3] = -12.*tmpv*tmpv/vh4;

	h[0] = x[0] - x0 + 0.1;

	Dh[0] = 1.;
	Dh[1] = 0.;

	DDh[0] = 0.;
	DDh[1] = 0.;
	DDh[2] = 0.;
	DDh[3] = 0.;

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

	const double vDes = (M_PI/10.);
	const double K = 10.;

	u[0] = K*(vDes-x[1]);
	Du[0] = 0.;
	Du[1] = -K;
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
	for(uint32_t i=0; i<2; i++)
	{
		ofstream myfiletb;
		myfiletb.open(string("asiftb_")+to_string(i)+string(".csv"), ofstream::out | ofstream::trunc);
		myfiletb << "tNow" <<','<< "x" <<','<<  "v" <<',' <<  "uDes" <<','<<  "uAct" <<','<< "relax" <<','<< "TTS" <<','<< "BTorthoBS" <<',' << "rc"<< endl;
		myfiletb << std::fixed;
		myfiletb << std::setprecision(10);

		ASIFimplicitTB::Options optsTB;
		optsTB.backTrajHorizon = 11.0;
		optsTB.backTrajDt = 0.001;
		optsTB.relaxCost = 10.;
		optsTB.relaxSafeLb = 10.0;
		optsTB.relaxTTS = 30.0;
		optsTB.relaxMinOrtho = 60.0;
		optsTB.backTrajMinOrtho = 0.001;

		asifTB = new ASIFimplicitTB(nx,nu,npSS,npBTSS,
		                            safetySet,backupSetTB,dynamics,dynamicsGradients,backupController);

		asifTB->initialize(lb,ub,optsTB);

		const double tEnd = 5.;
		const double dt = 0.001;
		double xNowTB[2] = {-0.1+static_cast<double>(i)*0.2,0.0};
		// double uDesNow[1] = {-uMax + static_cast<double>(i)*2.*uMax};
		double uDesNow[1] = {0.0};
		double uActNow[1] = {0.0};
		double tNow = 0.0;
		double relaxTB;
		int32_t rc;

		while(tNow<tEnd)
		{
			cout << "File " << i+1 << "/2, t " << tNow << "/" << tEnd << endl;
			tNow += dt;

			rc = asifTB->filter(xNowTB,uDesNow,uActNow,relaxTB);

			// Integrate
			double fCl[nx] = {0.0};
			double f[nx];
			double g[nx];
			dynamics(xNowTB, f, g);
			for(uint32_t i = 0; i<nx; i++)
			{
				fCl[i] += f[i];
				for(uint32_t k = 0; k<nu; k++)
				{
					fCl[i] += g[i+k*nx]*uActNow[k];
				}
			}

			for(uint32_t i = 0; i<nx; i++)
			{
				xNowTB[i] += dt*fCl[i];
			}

			// Save data
			myfiletb << tNow << "," ;
			myfiletb << xNowTB[0] << "," ;
			myfiletb << xNowTB[1] << "," ;
			myfiletb << uDesNow[0] << "," ;
			myfiletb << uActNow[0] << "," ;
			myfiletb << relaxTB << "," ;
			myfiletb << asifTB->TTS_ << "," ;
			myfiletb << asifTB->BTorthoBS_ << "," ;
			myfiletb << rc << endl;

			// myfiletb << asifTB->backTraj_.size() << endl;
			// for(const auto &logSample: asifTB->backTraj_)
			// {
			// 	myfiletb << logSample.first << "," ;
			// 	for(uint32_t i=0; i<2; i++)
			// 		myfiletb << logSample.second[i] << "," ;
			// 	myfiletb << endl;
			// }

			double hTmp;
			double DhTmp[2];
			double DDhTmp[4];
			backupSetTB(xNowTB,&hTmp,DhTmp,DDhTmp);
			
			if(hTmp>0)
				break;
		}

		myfiletb.close();
	}
	cout << "Finished" << endl;
	return 0;
}

