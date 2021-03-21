#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sstream>
#include <random>

#include <asif++.h>
#include "CyberTimer.hpp"

using namespace std;
using namespace ASIF;

const uint32_t nx = 2;
const uint32_t nu = 1;

const double lb[nu] = {-1.5};
const double ub[nu] = {1.5};
const double xBound[2] = {-M_PI,M_PI};
const double vBound[2] = {-M_PI,M_PI};
const vector<vector<double>> xUncertaintyList {{0.032,0.027},{0.164,0.133},{0.341,0.266}};
uint32_t dtPerSampleList[3] {10, 50, 100};

double p = 1.0;
const double pMin = 0.9;
const double pMax = 1.1;

ASIFrealizable *asif;
ASIFrealizable::kernel_t kernel;

ASIFrobust *asifRobust;
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

int32_t loadRealizableKernel(const string &filename)
{
	ifstream fin(filename);
	if(!fin.is_open()) throw runtime_error("Could not open file: ");

	string line, word;

	getline(fin, line);
	stringstream s(line); 
	if(getline(s, word, ','))
		kernel.vertices.resize(stoi(word),std::vector<double>(nx,0.0));
	else
	{
		throw runtime_error("Could not parse number of vertices");
		return -1;
	}

	getline(fin, line);
	s = stringstream(line); 
	if(getline(s, word, ','))
		kernel.facets.resize(stoi(word),ASIFrealizable::facet_t{std::vector<uint32_t>(nx),
		                                                        std::vector<double>(nx),
		                                                        std::vector<uint32_t>(),
		                                                        std::vector<interval_t>(nx)});
	else
	{
		throw runtime_error("Could not parse number of facets");
		return -1;
	}

	getline(fin, line);
	s = stringstream(line); 
	if(getline(s, word, ','))
	{
		kernel.maxCriticalFacets = stoi(word);
	}
	else
	{
		throw runtime_error("Could not parse max number of critical facets");
		return -1;
	}

	getline(fin, line);
	s = stringstream(line); 
	if(getline(s, word, ','))
	{
		kernel.maxActiveConstraints = stoi(word);
	}
	else
	{
		throw runtime_error("Could not parse max number of active constraints");
		return -1;
	}

	for(uint32_t i=0; i<kernel.vertices.size(); i++)
	{
		if(fin.eof()) throw runtime_error("EOF reached");
		getline(fin, line);
		stringstream ss(line);
		uint32_t j=0;
		while (getline(ss, word, ','))
		{
			kernel.vertices[i][j] = stod(word);
			j++;
		}
	}

	for(uint32_t i=0; i<kernel.facets.size(); i++)
	{
		if(fin.eof()) throw runtime_error("EOF reached");
		getline(fin, line);
		stringstream ss(line);
		uint32_t j=0;
		while (getline(ss, word, ','))
		{
			kernel.facets[i].verticesIdx[j] = stoi(word);
			j++;
		}
		
		if(fin.eof()) throw runtime_error("EOF reached");
		getline(fin, line);
		ss = stringstream(line);
		j=0;
		while (getline(ss, word, ','))
		{
			kernel.facets[i].normal[j] = stod(word);
			j++;
		}
		
		if(fin.eof()) throw runtime_error("EOF reached");
		getline(fin, line);
		ss = stringstream(line);
		j=0;
		kernel.facets[i].activeConstraintsSet.clear();
		while (getline(ss, word, ','))
		{
			kernel.facets[i].activeConstraintsSet.emplace_back(stoi(word));
			j++;
		}
	}

	fin.close();

	cout << "Loaded safety kernel! " << endl;
	cout << "__" << kernel.vertices.size() << " vertices." << endl;
	cout << "__" << kernel.facets.size() << " facets." << endl;
	cout << "__" << kernel.maxCriticalFacets << " max critical facets." << endl;
	cout << "__" << kernel.maxActiveConstraints << " max active constraints." << endl;

	return 1;
}

int32_t loadRobustKernel(const string &filename)
{
	// Load set data
	vector<string> rowString;
	vector<double> rowNumber;

	ifstream fin(filename);

	if(!fin.is_open()) throw runtime_error("Could not open file");

	string line, word;
	while (!fin.eof())
	{
		rowString.clear();
		getline(fin, line);
		stringstream s(line); 
		while (getline(s, word, ','))
		{
			rowString.push_back(word); 
		}
		if(rowString.size()>0)
		{
			rowNumber.resize(rowString.size());
			for(uint32_t i=0; i<rowString.size();i++)
			{
				rowNumber[i] = stod(rowString[i]); 
			}
			SafetySetData.push_back(rowNumber);
		}
	} 
	fin.close();
	cout << "Loaded set with " << SafetySetData.size() << " parts" << endl;

	return 1;
}

int main()
{
	CyberTimer<500> timer;

	if(loadRealizableKernel("examples/IP_realizable_sampled_set_export.csv")<0)
		return -1;

	if(loadRobustKernel("examples/IP_RealizableRobust_set_export.csv")<0)
		return -1;

	for(uint32_t i=3; i<3; i++)
	{

		cout << "File " << i+1 << "/3" << endl;
		ofstream myfile;

		myfile.open(string("asif_realizableSampled_")+to_string(i)+string(".csv"), ofstream::out | ofstream::trunc);

		myfile << "tNow" <<','<< "x" <<','<<  "v" <<','<< "xEstim" <<','<<  "vEstim" <<',' <<  "uDes" <<','<<  "uAct" <<','<< "relax" <<',' << "relaxHard" <<','<< "rc" << endl;
		myfile << std::fixed;
		myfile << std::setprecision(10);

		ASIFrealizable::Options opts;

		opts.relaxCost = 50.0;
		opts.relaxDes = 1.;

		asif = new ASIFrealizable(nx,nu,xUncertaintyList[i].data(),kernel,dynamics);
		
		timer.tic();
		asif->initialize(lb,ub,opts);
		cout << "Initialization time: " << 1000.*timer.toc() << "ms" << endl;
		timer.reset();
		
		const double tEnd = 6.0;
		const double dt = 0.001;

		double xNow[2] = {0.5,0.0};
		double xEstim[2] = {0.0};
		double uDesNow[1] = {0.0};
		double uActNow[1] = {0.0};
		double tNow = 0.0;
		double relax[2] = {0.0};
		int32_t rc = 0;

		uint32_t skipCounter = 0;
		while(tNow<tEnd)
		{
			if(skipCounter==0)
			{
				xEstim[0] = xNow[0];
				xEstim[1] = xNow[1];
				timer.tic();
				rc = asif->filter(xEstim,uDesNow,uActNow,relax);
				timer.toc();
			}

			skipCounter++;
			if(skipCounter==dtPerSampleList[i])
				skipCounter=0;


			// Save data
			myfile << tNow << "," ;
			myfile << xNow[0] << "," ;
			myfile << xNow[1] << "," ;
			myfile << xEstim[0] << "," ;
			myfile << xEstim[1] << "," ;
			myfile << uDesNow[0] << "," ;
			myfile << uActNow[0] << "," ;
			myfile << relax[0] << "," ;
			myfile << relax[1] << "," ;
			myfile << rc << endl ;
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
		delete asif;
	}

	ofstream myfile;
	myfile.open("asif_realizableSampledRobust.csv", ofstream::out | ofstream::trunc);
	myfile << "tNow" <<','<< "x" <<','<<  "v" <<','<< "xEstim" <<','<<  "vEstim" <<',' <<  "uDes" <<','<<  "uAct" <<','<< "relax" <<',' << "rc" << endl;
	myfile << std::fixed;
	myfile << std::setprecision(10);

	ASIFrobust::Options opts;

	opts.relaxCost = 50.0;
	opts.relaxLb = 12.0;

	asifRobust = new ASIFrobust(nx,nu,SafetySetData.size(),safetySet,dynamics);

	timer.tic();
	asifRobust->initialize(lb,ub,opts);
	cout << "Initialization time: " << 1000.*timer.toc() << "ms" << endl;
	timer.reset();

	const double tEnd = 6.0;
	const double dt = 0.001;

	double xNow[2] = {0.5,0.0};
	double xEstim[2] = {0.0};
	double uDesNow[1] = {0.0};
	double uActNow[1] = {0.0};
	double tNow = 0.0;
	double relax = 0;
	int32_t rc = 0;

	uint32_t skipCounter = 0;
	while(tNow<tEnd)
	{
		if(skipCounter==0)
		{
			xEstim[0] = xNow[0];
			xEstim[1] = xNow[1];
			timer.tic();
			rc = asifRobust->filter(xEstim,uDesNow,uActNow,relax);
			timer.toc();
		}

		if(rc<0)
		{
			myfile.close();
			delete asifRobust;
			return -1;
		}

		skipCounter++;
		if(skipCounter==dtPerSampleList[2])
			skipCounter=0;

		// Save data
		myfile << tNow << "," ;
		myfile << xNow[0] << "," ;
		myfile << xNow[1] << "," ;
		myfile << xEstim[0] << "," ;
		myfile << xEstim[1] << "," ;
		myfile << uDesNow[0] << "," ;
		myfile << uActNow[0] << "," ;
		myfile << relax << "," ;
		myfile << rc << endl ;
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

	delete asifRobust;

	cout << "Finished" << endl;
	return 0;
}

