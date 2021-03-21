#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sstream>
#include "CyberTimer.hpp"

#include <asif++.h>
// #include "RealizableKernelData_10Hz.h"
#include "RealizableKernelData_100Hz.h"

using namespace std;
using namespace ASIF;

const uint32_t nx = 2;
const uint32_t nu = 1;

const double lb[nu] = {-20.};
const double ub[nu] = {20.};

// const double xUncertainty[2] = {0.32,0.28};
const double xUncertainty[2] = {0.031,0.028};

// const uint32_t dtPerSample = 100;
const uint32_t dtPerSample = 10;

const double m_max = 75.;
const double m_min = 70.;
const double m_mean = (m_min+m_max)/2.;
// const double m_mean = 140.;

const interval_t mInt = interval(m_min,m_max);

const double K = 5.7;
const double dK = 0.1;
const interval_t KInt = interval(K-dK,K+dK);

const double F = 23;
const double DF = 2;
const double eps = 1.0e-3;
const interval_t FInt = interval(F-DF,F-DF);

ASIFrealizable *asif;

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

	cout << "Loaded safety kernel! " << endl;
	cout << "__" << kernel.vertices.size() << " vertices." << endl;
	cout << "__" << kernel.facets.size() << " facets." << endl;
	cout << "__" << kernel.maxCriticalFacets << " max critical facets." << endl;
	cout << "__" << kernel.maxActiveConstraints << " max active constraints." << endl;

	ofstream myfile;
	myfile.open("DI_asif_realizableSampled.csv", ofstream::out | ofstream::trunc);

	myfile << "tNow" <<','<< "x" <<','<<  "v" <<','<< "xEstim" <<','<<  "vEstim" <<',' <<  "uDes" <<','<<  "uFilter" <<','<<  "uAct" <<','<< "relax" <<',' << "relaxHard" <<','<< "rc" <<',' << "smoothBoundMin" <<','<< "smoothBoundMax" << endl;
	myfile << std::fixed;
	myfile << std::setprecision(10);

	ASIFrealizable::Options opts;

	opts.relaxCost = 100.0;
	opts.relaxOffset = 0.0;
	opts.relaxDes = 10.0;

	asif = new ASIFrealizable(nx,nu,xUncertainty,kernel,dynamics,2);

	timer.tic();
	asif->initialize(lb,ub,opts);
	cout << "Initialization time: " << 1000.*timer.toc() << "ms" << endl;
	timer.reset();

	const double tEnd = 6.0;
	const double dt = 0.001;

	double xNow[2] = {0.0,0.0};
	double xEstim[2] = {0.0};
	double uDesNow[1] = {20.0};
	double uActNow[1] = {0.0};
	double uFilterNow[1] = {0.0};
	double tNow = 0.0;
	double relax[2] = {0.0};
	int32_t rc = 0;
	double smoothBounds[2] = {-20.,20.};
	const double sc = 20.*dtPerSample*0.001;

	uint32_t skipCounter = 0;
	while(tNow<tEnd)
	{
		if(skipCounter==0)
		{
			xEstim[0] = xNow[0];
			xEstim[1] = xNow[1];
			timer.tic();
			rc = asif->filter(xEstim,uDesNow,uFilterNow,relax);
			uActNow[0] = uFilterNow[0];

			if(uActNow[0]>uDesNow[0] && uActNow[0]>smoothBounds[0])
			{
				smoothBounds[0] = uActNow[0];
				if(smoothBounds[0]>smoothBounds[1])
					smoothBounds[1] = smoothBounds[0];
				else
					smoothBounds[1] += sc;
			}
			else if(uActNow[0]<uDesNow[0] && uActNow[0]<smoothBounds[1])
			{
				smoothBounds[1] = uActNow[0];
				if(smoothBounds[0]>smoothBounds[1])
					smoothBounds[0] = smoothBounds[1];
				else
					smoothBounds[0] -= sc;
			}
			else
			{
				smoothBounds[0] -= sc;
				smoothBounds[1] += sc; 
			}

			if(smoothBounds[0]<-20.)
				smoothBounds[0] = -20.;
			if(smoothBounds[1]>20.)
				smoothBounds[1] = 20.;

			if(uActNow[0] > smoothBounds[1])
				uActNow[0] = smoothBounds[1];
			if(uActNow[0] < smoothBounds[0])
				uActNow[0] = smoothBounds[0];

			timer.toc();
		}

		skipCounter++;
		if(skipCounter==dtPerSample)
			skipCounter=0;

		// Save data
		myfile << tNow << "," ;
		myfile << xNow[0] << "," ;
		myfile << xNow[1] << "," ;
		myfile << xEstim[0] << "," ;
		myfile << xEstim[1] << "," ;
		myfile << uDesNow[0] << "," ;
		myfile << uFilterNow[0] << "," ;
		myfile << uActNow[0] << "," ;
		myfile << relax[0] << "," ;
		myfile << relax[1] << "," ;
		myfile << rc << "," ;
		myfile << smoothBounds[0] << "," ;
		myfile << smoothBounds[1] << endl;
		cout << "t: " << tNow << "/" << tEnd << ", uAct: " << uActNow[0] << ", rc:" << rc << endl;

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
			xNow[i] += dt*fCl[i];

		tNow += dt;
	}

	cout << "dt average:  " << timer.getAverage() << endl << endl;

	myfile.close();


	cout << "Finished" << endl;
	return 0;
}

