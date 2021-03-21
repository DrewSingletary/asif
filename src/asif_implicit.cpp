#include "asif_implicit.h"
#include <iostream>

using namespace std;

namespace ASIF
{
	ASIFimplicit::ASIFimplicit(const uint32_t nx,
	                           const uint32_t nu,
	                           const uint32_t npSS,
	                           const uint32_t npBS,
	                           const uint32_t npBTSS,
	                           function<void(const double* /*x*/,
	                                               double* /*h*/,
	                                               double* /*Dh*/)> safetySet,
	                           function<void(const double* /*x*/,
	                                               double* /*h*/,
	                                               double* /*Dh*/)> backupSet,
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
	ASIFimplicit(false,
	             nx,
	             nu,
	             npSS,
	             npBS,
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

	ASIFimplicit::ASIFimplicit(const uint32_t nx,
	                           const uint32_t nu,
	                           const uint32_t npSS,
	                           const uint32_t npBS,
	                           const uint32_t npBTSS,
	                           function<void(const double* /*x*/,
	                                               double* /*h*/,
	                                               double* /*Dh*/)> safetySet,
	                           function<void(const double* /*x*/,
	                                               double* /*h*/,
	                                               double* /*Dh*/)> backupSet,
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
	ASIFimplicit(true,
	             nx,
	             nu,
	             npSS,
	             npBS,
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

	ASIFimplicit::ASIFimplicit(const bool dynamicsHasGradient,
	                           const uint32_t nx,
	                           const uint32_t nu,
	                           const uint32_t npSS,
	                           const uint32_t npBS,
	                           const uint32_t npBTSS,
	                           function<void(const double* /*x*/,
	                                               double* /*h*/,
	                                               double* /*Dh*/)> safetySet,
	                           function<void(const double* /*x*/,
	                                               double* /*h*/,
	                                               double* /*Dh*/)> backupSet,
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
	nv_(nu+2),
	npSS_(npSS),
	npBS_(npBS),
	npBTSS_(npBTSS),
	npTC_(npBTSS*npSS+npBS),
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

	ASIFimplicit::~ASIFimplicit(void)
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

	int32_t ASIFimplicit::initialize(const double lb[],
	                                 const double ub[])
	{

		// Initialize Lfh and Lgh outputs (for learning)
		for (int i = 0; i < (int) npTC_; i++)
		{
			Lfh_out_.push_back(0.0);
			std::vector<double> tmp;
			std::vector<double> tmp2;
				for (int j = 0; j < (int) nu_; j++)
				{
					tmp.push_back(0.0);
				}
				for (int j = 0; j < (int) nx_; j++)
				{
					tmp2.push_back(0.0);
				}
			Lgh_out_.push_back(tmp);
			Dh_out_.push_back(tmp2);
		}


		#ifdef ASIF_DEBUG
		cout << "Initializing ASIFimplicit" << endl;
		#endif

		if (options_.x0 == nullptr)
			x0_ = new double[nx_]();
		else
			x0_ = options_.x0;

		for (int i = 0; i < (int) nx_*(int)npSS_; i++)
			Dh_index_.push_back(0.0);
		for (int i = 0; i < (int)npSS_; i++)
			h_index_.push_back(0.0);

		// Initialize backup trajectory
		npBT_ = round(options_.backTrajHorizon/options_.backTrajDt)+1;
		if(npBT_<npBTSS_)
		{
			npBT_ = npBTSS_;
			options_.backTrajDt = options_.backTrajHorizon/static_cast<double>(npBT_-1);
		}
		backTraj_ = vector<pair<double,state_t>>(npBT_,pair<double,state_t>(0.0,state_t(nx_+nx_*nx_)));

		if (options_.n_debug != -1){
			if (options_.n_debug > -1 && options_.n_debug < (int)npBT_ -1)
				index_debug_ = options_.n_debug;
			else
				options_.n_debug = -1;
		}
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
		H_[(nv_-2)+((nv_-2)*nv_)] = options_.relaxCost;
		H_[(nv_-1)+((nv_-1)*nv_)] = options_.relaxCost;

		// Initialize bounds
		memcpy(lb_,lb,nu_*sizeof(double));
		lb_[nv_-2] = options_.relaxSafeLb;
		lb_[nv_-1] = options_.relaxReachLb;
		memcpy(ub_,ub,nu_*sizeof(double));
		ub_[nv_-2] = options_.inf;
		ub_[nv_-1] = options_.inf;

		// Initialize constraints
		updateConstraints(x0_);

		// Initialize linear cost
		c_[nv_-2] = -2.0*options_.relaxCost*options_.relaxSafeLb;
		c_[nv_-1] = -2.0*options_.relaxCost*options_.relaxReachLb;
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


	int32_t ASIFimplicit::initialize(const double lb[],
	                                 const double ub[],
	                                 const Options &options)
	{
		options_ = options;
		return initialize(lb,ub);
	}

	int32_t ASIFimplicit::filter(const double x[],
	                             const double uDes[],
	                                   double uAct[])
	{
		double relax[2];
		updateCost(uDes);
		return filter(x,nullptr,c_,uAct,relax);
	}

	int32_t ASIFimplicit::filter(const double x[],
	                             const double uDes[],
	                                   double uAct[],
	                                   double relax[2])
	{
		updateCost(uDes);
		return filter(x,nullptr,c_,uAct,relax);
	}

	int32_t ASIFimplicit::filter(const double x[],
	                             const double H[],
	                             const double c[],
	                                   double uAct[])
	{
		double relax[2];
		return filter(x,H,c,uAct,relax);
	}

	int32_t ASIFimplicit::filter(const double x[],
	                             const double H[],
	                             const double c[],
	                                   double uAct[],
	                                   double relax[2])
	{
		#ifdef ASIF_DEBUG
		cout << "Filtering"<< endl;
		#endif

		double h[npBS_];
		double Dh[nx_*npBS_];
		backupSet_(backTraj_.back().second.data(),h,Dh);
		hBackupEnd_ = *std::min_element(h,h+npBS_);

		double hSafetyNow[npSS_];
		double DhSafetyNow[npSS_*nx_];
		safetySet_(x,hSafetyNow,DhSafetyNow);
		hSafetyNow_ = *std::min_element(hSafetyNow,hSafetyNow+npSS_);

		updateConstraints(x);
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
			relax[0] = solutionFull[nu_];
			relax[1] = solutionFull[nu_+1];
			return 1;
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

	int32_t ASIFimplicit::updateOptions(const Options &options)
	{
		#ifdef ASIF_DEBUG
		cout << "Updating options"<< endl;
		#endif

		options_ = options;
		return updateOptions();
	}

	int32_t ASIFimplicit::updateOptions(void)
	{
		npBT_ = round(options_.backTrajHorizon/options_.backTrajDt)+1;
		if(npBT_<npBTSS_)
		{
			npBT_ = npBTSS_;
			options_.backTrajDt = options_.backTrajHorizon/static_cast<double>(npBT_-1);
		}

		backTraj_ = vector<pair<double,state_t>>(npBT_,pair<double,state_t>(0.0,state_t(nx_+nx_*nx_)));
		H_[(nv_-2)+((nv_-2)*nv_)] = options_.relaxCost;
		H_[(nv_-1)+((nv_-1)*nv_)] = options_.relaxCost;
		c_[nv_-2] = -2.0*options_.relaxCost*options_.relaxSafeLb;
		c_[nv_-1] = -2.0*options_.relaxCost*options_.relaxReachLb;
		lb_[nv_-2] = options_.relaxSafeLb;
		lb_[nv_-1] = options_.relaxReachLb;
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

	int32_t ASIFimplicit::updateConstraints(const double x[])
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
		auto rhsBind = bind(&ASIFimplicit::ODErhs, this,
		                         placeholders::_1,
		                         placeholders::_2,
		                         placeholders::_3);
		auto stepper = make_dense_output(options_.backTrajAbsTol,
		                                 options_.backTrajRelTol,
		                                 stepper_t());
		auto itBegin = make_n_step_iterator_begin(stepper, rhsBind, x0BackupTraj, 0.0, options_.backTrajDt, npBT_-1);
		auto itEnd = make_n_step_iterator_end(stepper, rhsBind, x0BackupTraj);
		#endif

		// Select critical points along trajectory
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

		sort(idxOfMinH.begin(), idxOfMinH.end(), [&hFullMin](int a, int b) -> bool {return hFullMin[a]<hFullMin[b];});

		// Compute contigent cone
		double h[npTC_] = {0.0};
		double Dh[npTC_*nx_] = {0.0};
		double *DhSS;
		double DhSSDx[npSS_*nx_];
		double DhBS[npBS_*nx_];
		double DhBSDx[npBS_*nx_];
		uint32_t idxMinCurr;
		double *btDX;

		// get index_debug value of h and Dh
		if (options_.n_debug != -1) {
			double *DhSS_first_;
			double *btDX_first_;
			DhSS_first_ = &(DhFull[index_debug_*npSS_*nx_]);
			btDX_first_ = backTraj_[index_debug_].second.data()+nx_;
			matrixMultiply(DhSS_first_,npSS_,nx_,
										 btDX_first_,nx_,nx_,
										 Dh_index_.data());
			btX = backTraj_[index_debug_].second.data();
			double hLearn[npSS_];
			double DhLearn[npSS_*nx_];
			safetySet_(btX,
								 hLearn,
								 DhLearn);
			memcpy(h_index_.data(),hLearn,npSS_*sizeof(double));
		}
		////////////////////////////////////////////////////
		backTrajCritIdx_.assign(idxOfMinH.begin(),idxOfMinH.begin()+ npBTSS_);
		for(uint32_t idx = 0; idx<npBTSS_; idx++)
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
					if (idx == 0 && options_.n_debug == -1) {
						Dh_index_[j*npSS_+i] = DhSSDx[i+(j*npSS_)];
						h_index_[i] = h[i+idx*npSS_];
						index_debug_ = idxMinCurr;
					}
				}
			}
		}

		btX = backTraj_.back().second.data();
		btDX = btX+nx_;
		backupSet_(btX,&(h[npBTSS_*npSS_]),DhBS);
		matrixMultiply(DhBS,npBS_,nx_,
		               btDX,nx_,nx_,
		               DhBSDx);
		for(uint32_t i = 0; i<npBS_; i++)
		{
			for(uint32_t j = 0; j<nx_; j++)
			{
				Dh[(npBTSS_*npSS_+i)+(j*npTC_)] = DhBSDx[i+(j*npBS_)];
			}
		}

		for (int i = 0; i < (int) npTC_; i++)
		{
			for (int j = 0; j < (int) nx_; j++){
				Dh_out_[i][j] = Dh[nx_*i+j];
			}
		}

		// Compute Lfh
		double Lfh[npTC_];
		matrixVectorMultiply(Dh,npTC_,nx_,
							  f,nx_,
							  Lfh);

		for (int i = 0; i < (int) npTC_; i++)
			Lfh_out_[i] = Lfh[i];


		// Compute Lgh
		double Lgh[npTC_*nu_];
		matrixMultiply(Dh,npTC_,nx_,
					   g,nx_,nu_,
					   Lgh);
		for (int i = 0; i < (int) npTC_; i++)
		{
			for (int j = 0; j < (int) nu_; j++){
				Lgh_out_[i][j] = Lgh[nu_*i+j];
			}
		}

		// Update Lfh and Lgh with learning if necessary
		if (options_.use_learning) {
			update_weights(&learning_data_,x,nx_,Dh_index_.data(),Lfh,Lgh,nu_);
		}

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
		for(uint32_t i=npBTSS_*npSS_;i<npTC_;i++)
		{
			A_[i+((nu_+1)*npTC_)] = h[i];
		}

		// Fillup b
		for(uint32_t i=0;i<npTC_;i++)
		{
			b_[i] = -Lfh[i];
		}

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

	int32_t ASIFimplicit::updateCost(const double uDes[])
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

	int32_t ASIFimplicit::updateH(const double H[])
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

	void ASIFimplicit::inputSaturateSoft(const double u[],
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

	void ASIFimplicit::inputSaturate(double u[])
	{
		for(uint32_t i = 0; i<nu_; i++)
		{
			if(u[i]>ub_[i])
				u[i] = ub_[i];
			else if(u[i]<lb_[i])
				u[i] = lb_[i];
		}
	}


	void ASIFimplicit::backupCLdynamics(const double x[],
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

	void ASIFimplicit::ODErhs(const vector<double> &x,
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

	string ASIFimplicit::filterErrorMsgString(const int32_t rc)
	{
		switch(rc)
		{
			case 1:
				return "Success";

			case -1:
				return "QP failed";

			default:
				return "Unkown";
		}
	}
} //end ASIF namespace
