#include "asif_implicit_tb.h"

using namespace std;

namespace ASIF
{
	ASIFimplicitTB::ASIFimplicitTB(const uint32_t nx,
	                               const uint32_t nu,
	                               const uint32_t npSS,
	                               const uint32_t npBTSS,
	                               function<void(const double* /*x*/,
	                                                   double* /*h*/,
	                                                   double* /*Dh*/)> safetySet,
	                               function<void(const double* /*x*/,
	                                                   double* /*h*/,
	                                                   double* /*Dh*/,
	                                                   double* /*DDh*/)> backupSet,
	                               function<void(const double* /*x*/,
	                                                   double* /*f*/,
	                                                   double* /*g*/)> dynamics,
	                               function<void(const double* /*x*/,
	                                                   double* /*Df*/,
	                                                   double* /*Dg*/)> dynamicsGradients,
	                               function<void(const double* /*x*/,
	                                                   double* /*u*/,
	                                                   double* /*Du*/)> backupController,
	                               const QPSOLVER qpSolverType,
	                               const bool diagonalCost):
	ASIFimplicitTB(false,
	             nx,
	             nu,
	             npSS,
	             npBTSS,
	             safetySet,
	             backupSet,
	             dynamics,
	             dynamicsGradients,
	             nullptr,
	             backupController,
	             qpSolverType,
	             diagonalCost)
	{}

	ASIFimplicitTB::ASIFimplicitTB(const uint32_t nx,
	                               const uint32_t nu,
	                               const uint32_t npSS,
	                               const uint32_t npBTSS,
	                               function<void(const double* /*x*/,
	                                                   double* /*h*/,
	                                                   double* /*Dh*/)> safetySet,
	                               function<void(const double* /*x*/,
	                                                   double* /*h*/,
	                                                   double* /*Dh*/,
	                                                   double* /*DDh*/)> backupSet,
	                               function<void(const double* /*x*/,
	                                             const double* /*u*/,
	                                                   double* /*f*/,
	                                                   double* /*g*/,
	                                                   double* /*d_fcl_dx*/)> dynamicsWithGradient,
	                               function<void(const double* /*x*/,
	                                                   double* /*u*/,
	                                                   double* /*Du*/)> backupController,
	                               const QPSOLVER qpSolverType,
	                               const bool diagonalCost):
	ASIFimplicitTB(true,
	             nx,
	             nu,
	             npSS,
	             npBTSS,
	             safetySet,
	             backupSet,
	             nullptr,
	             nullptr,
	             dynamicsWithGradient,
	             backupController,
	             qpSolverType,
	             diagonalCost)
	{
		dynamics_ = [this]
		(const double* x, double* f, double* g)
		->void
		{
			const double uTemp[nu_] = {0.0};
			double DfclTemp[nx_*nx_] = {0.0};
			dynamicsWithGradient_(x,uTemp,f,g,DfclTemp);
		};
	}

	ASIFimplicitTB::ASIFimplicitTB(const bool dynamicsHasGradient,
	                               const uint32_t nx,
	                               const uint32_t nu,
	                               const uint32_t npSS,
	                               const uint32_t npBTSS,
	                               function<void(const double* /*x*/,
	                                                   double* /*h*/,
	                                                   double* /*Dh*/)> safetySet,
	                               function<void(const double* /*x*/,
	                                                   double* /*h*/,
	                                                   double* /*Dh*/,
	                                                   double* /*DDh*/)> backupSet,
	                               function<void(const double* /*x*/,
	                                                   double* /*f*/,
	                                                   double* /*g*/)> dynamics,
	                               function<void(const double* /*x*/,
	                                                   double* /*Df*/,
	                                                   double* /*Dg*/)> dynamicsGradients,
	                               function<void(const double* /*x*/,
	                                             const double* /*u*/,
	                                                   double* /*f*/,
	                                                   double* /*g*/,
	                                                   double* /*d_fcl_dx*/)> dynamicsWithGradient,
	                               function<void(const double* /*x*/,
	                                                   double* /*u*/,
	                                                   double* /*Du*/)> backupController,
	                               const QPSOLVER qpSolverType,
	                               const bool diagonalCost):
	hBackupEnd_(0.),
	hSafetyNow_(0.),
	dynamicsHasGradient_(dynamicsHasGradient),
	nx_(nx),
	nu_(nu),
	nv_(nu+1),
	npSS_(npSS),
	npBTSS_(npBTSS),
	npTC_(npBTSS*npSS+2),
	safetySet_(safetySet),
	backupSet_(backupSet),
	dynamics_(dynamics),
	dynamicsGradients_(dynamicsGradients),
	dynamicsWithGradient_(dynamicsWithGradient),
	backupController_(backupController),
	options_()
	{
		switch(qpSolverType)
		{
			case QPSOLVER::OSQP:
			{
				QPsolver_ = new QPWrapperOsqp(nv_,npTC_,diagonalCost);
				break;
			}
			default :
			{
				QPsolver_ = nullptr;
				break;
			}
		}

		H_ = new double[nv_*nv_]{0.0};
		c_ = new double[nv_]{0.0};
		A_ = new double[npTC_*nv_]{0.0};
		b_ = new double[npTC_]{0.0};
		lb_ = new double[nv_]{0.0};
		ub_ = new double[nv_]{0.0};
	}

	ASIFimplicitTB::~ASIFimplicitTB(void)
	{
		delete[] H_;
		delete[] c_;
		delete[] A_;
		delete[] b_;
		delete[] lb_;
		delete[] ub_;

		if(QPsolver_!=nullptr)
			delete QPsolver_;
	}

	int32_t ASIFimplicitTB::initialize(const double lb[],
	                                   const double ub[])
	{
		#ifdef ASIF_DEBUG
		cout << "Initializing ASIFimplicitTB" << endl;
		#endif

		// Initialize backup trajectory
		npBT_ = round(options_.backTrajHorizon*(1.0 + options_.backTrajExtend)/options_.backTrajDt)+1;
		if(npBT_<npBTSS_)
		{
			npBT_ = npBTSS_;
			options_.backTrajDt = options_.backTrajHorizon*(1.0 + options_.backTrajExtend)/static_cast<double>(npBT_-1);
		}
		backTraj_ = vector<pair<double,state_t>>(npBT_,pair<double,state_t>(0.0,state_t(nx_+nx_*nx_)));

		// Initialize cost hessian
		double Hidentity[nu_*nu_];
		for(uint32_t j=0; j<nu_;j++)
		{
			for(uint32_t i=0; i<nu_;i++)
			{
				if(i==j)
					Hidentity[i+(j*nu_)] = 1.0;
				else
					Hidentity[i+(j*nu_)] = 0.0;
			}
		}
		updateH(Hidentity);
		H_[(nv_-1)+((nv_-1)*nv_)] = options_.relaxCost;

		// Initialize bounds
		memcpy(lb_,lb,nu_*sizeof(double));
		lb_[nv_-1] = options_.relaxSafeLb;
		memcpy(ub_,ub,nu_*sizeof(double));
		ub_[nv_-1] = options_.inf;

		// Initialize constraints
		updateConstraintsTrivial();

		// Initialize linear cost
		c_[nv_-1] = -2.0*options_.relaxCost*options_.relaxSafeLb;
		double uDes0[nu_] = {0.0};
		updateCost(uDes0);

		QPsolver_->initialize(H_, c_, A_, b_, lb_, ub_);
		QPsolver_->solve();

		#ifdef ASIF_DEBUG
		double solutionFull[nv_];
		QPsolver_->getSolution(solutionFull);
		#endif

		return 1;
	}


	int32_t ASIFimplicitTB::initialize(const double lb[],
	                                   const double ub[],
	                                   const Options &options)
	{
		options_ = options;
		return initialize(lb,ub);
	}

	int32_t ASIFimplicitTB::filter(const double x[],
	                               const double uDes[],
	                                     double uAct[])
	{
		double relax;
		updateCost(uDes);
		return filter(x,nullptr,c_,uAct,relax);
	}

	int32_t ASIFimplicitTB::filter(const double x[],
	                               const double uDes[],
	                                     double uAct[],
	                                     double &relax)
	{
		updateCost(uDes);
		return filter(x,nullptr,c_,uAct,relax);
	}

	int32_t ASIFimplicitTB::filter(const double x[],
	                               const double H[],
	                               const double c[],
	                                     double uAct[])
	{
		double relax;
		return filter(x,H,c,uAct,relax);
	}

	int32_t ASIFimplicitTB::filter(const double x[],
	                               const double H[],
	                               const double c[],
	                                     double uAct[],
	                                     double &relax)
	{
		#ifdef ASIF_DEBUG
		cout << "Filtering"<< endl;
		#endif

		double h[1];
		double Dh[nx_];
		double DDh[nx_*nx_];

		backupSet_(backTraj_.back().second.data(),h,Dh,DDh);
		hBackupEnd_ = h[0];

		backupSet_(x,h,Dh,DDh);
		
		double hSafetyNow[npSS_];
		double DhSafetyNow[npSS_*nx_];
		safetySet_(x,hSafetyNow,DhSafetyNow);
		hSafetyNow_ = *std::min_element(hSafetyNow,hSafetyNow+npSS_);

		if(h[0]>=0)
		{
			updateConstraintsTrivial();
			if(H==nullptr)
				QPsolver_->updateCost(nullptr,c);
			else
			{
				updateH(H);
				QPsolver_->updateCost(H_,c);
			}

			QPsolver_->updateA(A_);
			QPsolver_->updateb(b_);		
			int32_t rtCode = QPsolver_->solve();

			if(static_cast<QPWrapperAbstract::SOLVER_STATUS>(rtCode)==QPWrapperAbstract::SOLVER_STATUS::FEASIBLE)
			{
				double solutionFull[nv_];
				QPsolver_->getSolution(solutionFull);
				memcpy(uAct,solutionFull,nu_*sizeof(double));
				inputSaturate(uAct);
				relax = solutionFull[nu_];
				return 2;
			}
			else
			{			
				// Compute backup control action
				double Du[nu_*nx_];
				backupController_(x,uAct,Du);
				inputSaturate(uAct);
				return -1;
			}	
		}
		else
		{
			if(updateConstraints(x)==1)
			{
				backupSet_(backTraj_.back().second.data(),&hBackupEnd_,Dh,DDh);			

				if(H==nullptr)
					QPsolver_->updateCost(nullptr,c);
				else
				{
					updateH(H);
					QPsolver_->updateCost(H_,c);
				}

				QPsolver_->updateA(A_);
				QPsolver_->updateb(b_);		
				int32_t rtCode = QPsolver_->solve();

				if(rtCode==1)
				{
					double solutionFull[nv_];
					QPsolver_->getSolution(solutionFull);
					memcpy(uAct,solutionFull,nu_*sizeof(double));
					inputSaturate(uAct);
					relax = solutionFull[nu_];
					return 1;
				}
				else
				{			
					// Compute backup control action
					double Du[nu_*nx_];
					backupController_(x,uAct,Du);
					inputSaturate(uAct);
					return rtCode;
				}	
			}
			else
			{
				// Compute backup control action
				double Du[nu_*nx_];
				backupController_(x,uAct,Du);
				inputSaturate(uAct);
				return -3;
			}		
		}
	}

	int32_t ASIFimplicitTB::updateOptions(const Options &options)
	{
		#ifdef ASIF_DEBUG
		cout << "Updating options"<< endl;
		#endif

		options_ = options;
		return updateOptions();
	}

	int32_t ASIFimplicitTB::updateOptions(void)
	{
		npBT_ = round(options_.backTrajHorizon/options_.backTrajDt)+1;
		if(npBT_<npBTSS_)
		{
			npBT_ = npBTSS_;
			options_.backTrajDt = options_.backTrajHorizon/static_cast<double>(npBT_-1);
		}

		backTraj_ = vector<pair<double,state_t>>(npBT_,pair<double,state_t>(0.0,state_t(nx_+nx_*nx_)));
		H_[(nv_-1)+((nv_-1)*nv_)] = options_.relaxCost;
		c_[nv_-1] = -2.0*options_.relaxCost*options_.relaxSafeLb;
		lb_[nv_-1] = options_.relaxSafeLb;
		QPsolver_->updateBounds(lb_,nullptr);	
		QPsolver_->updateCost(H_,c_);

		if(options_.satSharpness>2)
		{
			options_.satSharpness = 2;
			return 2;
		}
		else if(options_.satSharpness<0.01)
		{
			options_.satSharpness = 0.01;
			return 3;
		}
		else
		{
			return 1;
		}
	}

	int32_t ASIFimplicitTB::updateConstraints(const double x[])
	{
		#ifdef ASIF_DEBUG
			cout << "Updating constraints: ";
			ASIFtimer timer;
			tic(&timer);
		#endif

		//Compute dynamics
		double f[nx_];
		double g[nx_*nu_];
		dynamics_(x,f,g);

		//Compute backup trajectory
		state_t x0BackupTraj(nx_+nx_*nx_);
		for(uint32_t i=0;i<nx_;i++)
		{
			x0BackupTraj[i]=x[i];
		}
		for(uint32_t i=nx_;i<(nx_+nx_*nx_);i+=nx_+1)
		{
			x0BackupTraj[i]=1.0;
		}

		#ifdef USE_ODEINT
		auto rhsBind = bind(&ASIFimplicitTB::ODErhs, this,
		                         placeholders::_1,
		                         placeholders::_2,
		                         placeholders::_3);	
		auto stepper = make_dense_output(options_.backTrajAbsTol,
		                                 options_.backTrajRelTol,
		                                 stepper_t());
		auto itBegin = make_n_step_iterator_begin(stepper, rhsBind, x0BackupTraj, 0.0, options_.backTrajDt, npBT_-1);
		auto itEnd = make_n_step_iterator_end(stepper, rhsBind, x0BackupTraj);
		#endif

		// Integrate trajectory forward in time
		double *btX;
		double hFull[npBT_*npSS_];
		double DhFull[npBT_*npSS_*nx_];
		double hFullMin[npBT_];
		vector<uint32_t> idxOfMinH(npBT_);

		#ifdef USE_ODEINT
			uint32_t iBT = 0;
			for(auto it = itBegin; it!=itEnd; it++)
			{
				backTraj_[iBT].first = options_.backTrajDt*static_cast<double>(iBT);
				copy(it->begin(), it->end(), backTraj_[iBT].second.begin());
				btX = backTraj_[iBT].second.data();
				safetySet_(btX,
				           &(hFull[iBT*npSS_]),
				           &(DhFull[iBT*npSS_*nx_]));
				hFullMin[iBT] = *min_element(&(hFull[iBT*npSS_]),&(hFull[iBT*npSS_+npSS_]));
				idxOfMinH[iBT] = iBT;
				iBT++;
			}
		#else
			backTraj_[0].first = 0.0;
			backTraj_[0].second = x0BackupTraj;
			btX = backTraj_[0].second.data();
			safetySet_(btX,
			           &(hFull[0]),
			           &(DhFull[0]));
			hFullMin[0] = *min_element(&(hFull[0]),&(hFull[0+npSS_]));
			idxOfMinH[0] = 0;
			for(uint32_t i = 1; i<npBT_; i++)
			{
				backTraj_[i].first = backTraj_[i-1].first + options_.backTrajDt;
				ODErhs(backTraj_[i-1].second,backTraj_[i].second,0.0);
				transform(backTraj_[i].second.begin(), backTraj_[i].second.end(), backTraj_[i].second.begin(),
				          bind(multiplies<double>(), placeholders::_1, options_.backTrajDt));
				transform (backTraj_[i].second.begin(), backTraj_[i].second.end(), backTraj_[i-1].second.begin(), backTraj_[i].second.begin(),
					        plus<double>());
				btX = backTraj_[i].second.data();
				safetySet_(btX,
				           &(hFull[i*npSS_]),
				           &(DhFull[i*npSS_*nx_]));
				hFullMin[i] = *min_element(&(hFull[i*npSS_]),&(hFull[i*npSS_+npSS_]));
				idxOfMinH[i] = i;
			}			
		#endif

		// Find time-to-safety
		double hBSnm1 = -1.0;
		double hBS[1];
		double DhBS[nx_];
		double DDhBS[nx_*nx_];
		double hReach;
		bool BSHit = false;
		uint32_t idxHit;
		double cosTilde[1];
		double fClBS[nx_];
		double DfClBS[nx_*nx_];
		double den1;
		double den2; 
		double den;
		
		for(uint32_t i = 1; i<npBT_; i++)
		{
			if(hBSnm1<0.0)
			{
				backupSet_(backTraj_[i].second.data(),hBS,DhBS,DDhBS);
				if(hBS[0]>=0.0)
				{
					idxHit = i;
					BSHit = true;
					btX = backTraj_[i].second.data();
					backupCLdynamics(btX,fClBS,DfClBS);
					matrixMultiply(DhBS,1,nx_,
					               fClBS,nx_,1,
					               cosTilde);
					den1 = vectorNorm(DhBS,nx_);
					den2 = vectorNorm(fClBS,nx_);
					den = den1*den2;
					BTorthoBS_ = cosTilde[0]/den;
					if(BTorthoBS_>2.0*options_.backTrajMinOrtho)
						break;				
				}
			}
			hBSnm1 = hBS[0];
		}
		if(!BSHit)
		{
			#ifdef ASIF_DEBUG
			cout << "Failed to determine time to safety." << endl;
			#endif
			BTorthoBS_ = 0;
			return -1;
		}

		// Sort values of safety set along trajectory
		sort(idxOfMinH.begin(), idxOfMinH.begin() + idxHit + 1, [&hFullMin](uint32_t a, uint32_t b) -> bool {return hFullMin[a]<hFullMin[b];});

		// Compute contigent cone
		double h[npTC_] = {0.0};
		double Dh[npTC_*nx_] = {0.0};
		double *DhSS;
		double DhSSDx[npSS_*nx_];
		uint32_t idxMinCurr;	
		double *btDX;
		double DhBSDx[nx_];
		if(npBTSS_>idxHit)
			backTrajCritIdx_.assign(idxOfMinH.begin(),idxOfMinH.begin() + idxHit + 1);
		else
			backTrajCritIdx_.assign(idxOfMinH.begin(),idxOfMinH.begin()+ npBTSS_);

		for(uint32_t idx = 0; idx<npBTSS_; idx++)
		{
			if(idx>idxHit)
			{
				for(uint32_t i = 0; i<npSS_; i++)
				{
					h[idx*npSS_+i] = 1.0;
					for(uint32_t j = 0; j<nx_; j++)
					{
						Dh[(idx*npSS_+i)+(j*npTC_)] = 0.0;
					}
				}
			}
			else
			{
				idxMinCurr = idxOfMinH[idx];
				memcpy(&(h[idx*npSS_]),&(hFull[idxMinCurr*npSS_]),npSS_*sizeof(double));
				DhSS = &(DhFull[idxMinCurr*npSS_*nx_]);

				btDX = backTraj_[idxMinCurr].second.data()+nx_;
				matrixMultiply(DhSS,npSS_,nx_,
				               btDX,nx_,nx_,
				               DhSSDx);

				for(uint32_t i = 0; i<npSS_; i++)
				{
					for(uint32_t j = 0; j<nx_; j++)
					{
						Dh[(idx*npSS_+i)+(j*npTC_)] = DhSSDx[i+(j*npSS_)];
					}
				}				
			}
		}

		TTS_ = backTraj_[idxHit].first;
		hReach = options_.backTrajHorizon - backTraj_[idxHit].first;
		btDX = btX+nx_;
		matrixMultiply(DhBS,1,nx_,
		               btDX,nx_,nx_,
		               DhBSDx);

		h[npBTSS_*npSS_] = hReach;	
		for(uint32_t i = 0; i<nx_; i++)
		{
			Dh[(npBTSS_*npSS_)+(i*npTC_)] = DhBSDx[i]/cosTilde[0];
		}

		double denSquared = den*den;
		h[npBTSS_*npSS_+1] = BTorthoBS_ - options_.backTrajMinOrtho;
		double DxHit[nx_*nx_];
		matrixMultiply(fClBS,nx_,1,
		               DhBSDx,1,nx_,
		               DxHit);
		for(uint32_t i = 0; i<nx_*nx_; i++)
		{
			DxHit[i] = btDX[i] - DxHit[i];
		}
		double Dnum[nx_] = {0.0};
		double Dden1[nx_] = {0.0};
		double Dden2[nx_] = {0.0};
		double Dden[nx_];
		for(uint32_t i = 0; i<nx_; i++)
		{
			for(uint32_t k = 0; k<nx_; k++)
			{
				double temp1 = 0.0;
				double temp2 = 0.0; 
				for(uint32_t l = 0; l<nx_; l++)
				{
					temp1 += DDhBS[k+l*nx_]*DxHit[l+i*nx_];
					temp2 += DfClBS[k+l*nx_]*DxHit[l+i*nx_];					
				}
				double temp3 = DhBS[k]*temp2;
				double temp4 = temp1*fClBS[k];
				Dden1[i] += temp3;
				Dden2[i] += temp4;
				Dnum[i]  += temp3 + temp4;	
			}
		}
		for(uint32_t i = 0; i<nx_; i++)
		{
			Dden[i] = den2*Dden1[i]/den1 + den1*Dden2[i]/den2;		
		}
		for(uint32_t i = 0; i<nx_; i++)
		{
			DhBSDx[i] = (Dnum[i]*den - cosTilde[0]*Dden[i])/denSquared;	
			Dh[(npBTSS_*npSS_+1)+(i*npTC_)] = DhBSDx[i];	
		}		      

		// Compute Lfh
		double Lfh[npTC_];
		matrixVectorMultiply(Dh,npTC_,nx_,
		                     f,nx_,
		                     Lfh);

		// Compute Lgh
		double Lgh[npTC_*nu_];
		matrixMultiply(Dh,npTC_,nx_,
		               g,nx_,nu_,
		               Lgh);	

		// Fillup A
		for(uint32_t i=0;i<npTC_;i++)
		{
			for(uint32_t j=0;j<nu_;j++)
			{
				A_[i+(j*npTC_)] = Lgh[i+(j*npTC_)];
			}
		}
		for(uint32_t i=0;i<npBTSS_*npSS_;i++)
		{
			A_[i+(nu_*npTC_)] = h[i];		
		}

		// Fillup b
		for(uint32_t i=0;i<npTC_;i++)
		{
			b_[i] = -Lfh[i];
		}
		b_[npBTSS_*npSS_] -= options_.relaxTTS*h[npBTSS_*npSS_];
		b_[npBTSS_*npSS_+1] -= options_.relaxMinOrtho*(h[npBTSS_*npSS_+1]);

		#ifdef ASIF_DEBUG
			toc(&timer);
			cout << timer.dt*1e6 << "us" << endl;

			cout << "x0: " << endl;
			printVector(backTraj_[0].second.data(),nx_,6);

			cout << "xend: " << endl;
			printVector(backTraj_[npBT_-1].second.data(),nx_,6);

			cout << "h:" << endl;
			printVector(h,npTC_);

			cout << "Dh:" << endl;
			printMatrix(Dh,npTC_,nx_);

			cout << "f:" << endl;
			printVector(f,nx_);

			cout << "g:" << endl;
			printMatrix(g,nx_,nu_);

			cout << "Lfh:" << endl;
			printVector(Lfh,npTC_);

			cout << "Lgh:" << endl;
			printMatrix(Lgh,npTC_,nu_);

			cout << "A:" << endl;
			printMatrix(A_,npTC_,nv_);

			cout << "b:" << endl;
			printVector(b_,npTC_);
		#endif



		return 1;		
	}

	int32_t ASIFimplicitTB::updateConstraintsTrivial(void)
	{
		// Fillup A
		for(uint32_t i=0;i<npTC_*nv_;i++)
		{
			A_[i] = 0.0;
		}

		// Fillup b
		for(uint32_t i=0;i<npTC_;i++)
		{
			b_[i] = -options_.inf;
		}

		TTS_ = 0.0;
		BTorthoBS_ = 1.0;
		return 1;
	}

	int32_t ASIFimplicitTB::updateCost(const double uDes[])
	{
		#ifdef ASIF_DEBUG
		cout << "Updating cost"<< endl;
		#endif

		for(uint32_t i=0; i<nu_;i++)
		{
			c_[i] = -2.0*uDes[i];
		}
		return 1;	
	}

	int32_t ASIFimplicitTB::updateH(const double H[])
	{
		#ifdef ASIF_DEBUG
		cout << "Building H"<< endl;
		#endif

		for(uint32_t j=0; j<nu_;j++)
		{
			for(uint32_t i=0; i<nu_;i++)
			{
				H_[i+(j*nv_)] = H[i+(j*nu_)];
			}
		}
		return 1;	
	}

	void ASIFimplicitTB::inputSaturateSoft(const double u[],
	                                             double uSat[],
	                                             double DuSat[])
	{
		double uc, mi, ma, range, middle, bevelL, bevelXc, bevelYc, bevelStart, bevelStop;
		const double &r = options_.satSharpness;
		const double alpha = M_PI/8;
		const double beta = M_PI/4;
		for(uint32_t i = 0; i<nu_; i++)
		{
			mi = lb_[i];
			ma = ub_[i];
			range = ma - mi;
			middle = (ma+mi)/2;
			uc = 2*(u[i]-middle)/range;

			bevelL = r*tan(alpha);
			bevelStart = 1-cos(beta)*bevelL;
			bevelStop = 1+bevelL;
			bevelXc = bevelStop;
			bevelYc = 1 - r;

			if(uc>=bevelStop)
			{
				uSat[i] = ma;
				DuSat[i] = 0;
			}
			else if(uc<=-bevelStop)
			{
				uSat[i] = mi;
				DuSat[i] = 0;
			}
			else if(uc<=bevelStart && uc>=-bevelStart)
			{
				uSat[i] = u[i];
				DuSat[i] = 1;
			}
			else if(uc>bevelStart)
			{
				uSat[i] = sqrt(r*r - (uc-bevelXc)*(uc-bevelXc)) + bevelYc;
				DuSat[i] = (bevelXc-uc)/sqrt(r*r - (uc-bevelXc)*(uc-bevelXc));
				uSat[i] = 0.5*uSat[i]*range + middle;
			}
			else if(uc<-bevelStart)
			{
				uSat[i] = -sqrt(r*r - (uc+bevelXc)*(uc+bevelXc)) - bevelYc;
				DuSat[i] = (bevelXc+uc)/sqrt(r*r - (uc+bevelXc)*(uc+bevelXc));
				uSat[i] = 0.5*uSat[i]*range + middle;
			}
			else
			{
				DuSat[i] = 1;
				uSat[i] = u[i];			
			}
		}
	}

	void ASIFimplicitTB::inputSaturate(double u[])
	{
		for(uint32_t i = 0; i<nu_; i++)
		{
			if(u[i]>ub_[i])
				u[i] = ub_[i];
			else if(u[i]<lb_[i])
				u[i] = lb_[i];
		}
	}


	void ASIFimplicitTB::backupCLdynamics(const double x[],
	                                            double fCL[],
	                                            double DfCL[])
	{
		double f[nx_];
		double g[nx_*nu_];
		double u[nu_];
		double Du[nu_*nx_];
		double uSat[nu_];
		double DuSat[nu_];

		// Compute control action
		backupController_(x,u,Du);

		// Saturate inputs
		inputSaturateSoft(u,uSat,DuSat);

		// Compute gradient of closed loop dynamics
		uint32_t idx;
		if(dynamicsHasGradient_)
		{
			double d_fcl_dx[nx_*nx_];
			dynamicsWithGradient_(x,uSat,f,g,d_fcl_dx);

			for(uint32_t i = 0; i<nx_; i++)
			{
				for(uint32_t j = 0; j<nx_; j++)
				{
					idx = i+(j*nx_);
					DfCL[idx] = d_fcl_dx[idx];
					for(uint32_t k = 0; k<nu_; k++)
					{
						DfCL[idx]+= g[i+(k*nx_)]*DuSat[k]*Du[k+(j*nu_)];
					}
				}			
			}
		}
		else
		{
			double Df[nx_*nx_];
			double Dg[nx_*nu_*nx_];
			
			dynamics_(x,f,g);
			dynamicsGradients_(x,Df,Dg);
			for(uint32_t i = 0; i<nx_; i++)
			{
				for(uint32_t j = 0; j<nx_; j++)
				{
					idx = i+(j*nx_);
					DfCL[idx] = Df[idx];
					for(uint32_t k = 0; k<nu_; k++)
					{
						DfCL[idx]+= Dg[i+(k*nx_)+(j*nx_*nu_)]*uSat[k] + g[i+(k*nx_)]*DuSat[k]*Du[k+(j*nu_)];
					}
				}			
			}
		}

		// Compute closed loop dynamics
		matrixVectorMultiply(g,nx_,nu_,
		                     uSat,nu_,
		                     fCL);
		for(uint32_t i = 0; i<nx_; i++)
			fCL[i]+=f[i];
	}

	void ASIFimplicitTB::ODErhs(const vector<double> &x,
	                                  vector<double> &xDot,
	                            const double t)
	{
		xDot.resize(nx_+nx_*nx_);
		double DfCL[nx_*nx_];
		backupCLdynamics(x.data(),xDot.data(),DfCL);
		matrixMultiply(DfCL,nx_,nx_,
		               x.data()+nx_,nx_,nx_,
		               xDot.data()+nx_);
	}

	string ASIFimplicitTB::filterErrorMsgString(const int32_t rc)
	{
		switch(rc)
		{
			case 1:
				return "Success";

			case 2:
				return "Inside backup set";

			case -1:
				return "Inside backup set QP failed";

			case -2:
				return "QP failed";

			case -3:
				return "Backup set not reached";				

			default:	
				return "Unkown";
		}
	}	
} //end ASIF namespace