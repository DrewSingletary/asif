#include "asif_robust.h"

namespace ASIF
{
	ASIFrobust::ASIFrobust(const uint32_t nx,
	                       const uint32_t nu,
	                       const uint32_t npSS,
	                       std::function<void(const double* /*x*/,
	                                                double* /*h*/,
	                                                double* /*Dh*/)> safetySet,
	                       std::function<void(const interval_t* /*x*/,
	                                                interval_t* /*f*/,
	                                                interval_t* /*g*/)> dynamics,
	                       const uint32_t npSSmax,
	                       const QPSOLVER qpSolverType,
	                       const bool diagonalCost):
	nx_(nx),
	nu_(nu),
	npSS_(npSS),
	npSSmax_(npSSmax>npSS?npSS:npSSmax),
	nv_(nu+1+npSSmax_*2*(nu+1)),
	nc_(npSSmax_*(1+(nu+1))),
	safetySet_(safetySet),
	dynamics_(dynamics),
	options_()
	{
		H_ = new double[nv_*nv_]{0.0};
		c_ = new double[nv_]{0.0};
		A_ = new double[nc_*nv_]{0.0};
		b_ = new double[nc_]{0.0};
		lb_ = new double[nv_]{0.0};
		ub_ = new double[nv_]{0.0};

		switch(qpSolverType)
		{
			case QPSOLVER::OSQP:
			{
				QPsolver_ = new QPWrapperOsqp(nv_,nc_,diagonalCost);
				break;
			}
			default :
			{
				QPsolver_ = nullptr;
				break;
			}
		}
	}

	ASIFrobust::~ASIFrobust(void)
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


	int32_t ASIFrobust::initialize(const double lb[],
	                               const double ub[])
	{
		#ifdef ASIF_DEBUG
		std::cout << "Initializing ASIF" << std::endl;

		std::cout << "nx:" << nx_ << std::endl;
		std::cout << "nu:" << nu_ << std::endl;
		std::cout << "npSS:" << npSS_ << std::endl;
		std::cout << "npSSmax:" << npSSmax_ << std::endl;
		std::cout << "nv:" << nv_ << std::endl;
		std::cout << "nc:" << nc_ << std::endl;
		#endif

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
		H_[nu_+(nu_*nv_)] = options_.relaxCost;

		// Initialize bounds
		memcpy(lb_,lb,nu_*sizeof(double));
		lb_[nu_] = options_.relaxLb;
		memcpy(ub_,ub,nu_*sizeof(double));
		ub_[nu_] = options_.inf;
		for(uint32_t i=nu_+1; i<nv_; i++)
		{
			lb_[i] = 0.0;
			ub_[i] = options_.inf;
		}

		// Initialize contraint matrices
		uint32_t iColLambda = nu_+1;
		for(uint32_t iRow = 0; iRow<nc_; iRow+=(nu_+2))
		{
			// Fillup A
			for(uint32_t i=0;i<nu_;i++)
			{
				for(uint32_t j=0;j<nu_;j++)
				{
					A_[(iRow+1+i)+(j*nc_)] = -1.0;
				}
			}

			// Fillup A
			for(uint32_t i=0;i<(nu_+1);i++)
			{
				for(uint32_t j=0;j<(nu_+1);j++)
				{
					if(i==j)
					{
						A_[(iRow+1+i)+((iColLambda+j)*nc_)] = 1.0;
						A_[(iRow+1+i)+((iColLambda+nu_+1+j)*nc_)] = -1.0;
					}
				}
			}

			// Fillup b
			b_[iRow+nu_+1] = 1.0;		

			iColLambda+=2*(nu_+1);
		}

		// Initialize constraints
		double x0[nx_] = {0.0};
		updateConstraints(x0);

		// Initialize linear cost
		double uDes0[nu_] = {0.0};
		updateCost(uDes0);
		c_[nu_] = -2.0*options_.relaxCost*options_.relaxLb;

		// Initialize contraint matrix structure
		bool be[nc_];
		std::fill(be,be+nc_,true);
		for(uint32_t i=0; i<nc_; i+=(nu_+2))
			be[i] = false;

		QPsolver_->initialize(H_, c_, A_, b_, lb_, ub_, be);
		QPsolver_->solve();

		#ifdef ASIF_DEBUG
		double solutionFull[nv_];
		QPsolver_->getSolution(solutionFull);

		std::cout << "H:" << std::endl;
		printMatrix(H_,nv_,nv_);

		std::cout << "c:" << std::endl;
		printVector(c_,nv_);

		std::cout << "A:" << std::endl;
		printMatrix(A_,nc_,nv_);

		std::cout << "b:" << std::endl;
		printVector(b_,nc_);

		std::cout << "lb:" << std::endl;
		printVector(lb_,nc_);

		std::cout << "ub:" << std::endl;
		printVector(ub_,nc_);

		std::cout << "be:" << std::endl;
		printVector(be,nc_);

		#endif

		return 1;
	}

	int32_t ASIFrobust::initialize(const double lb[],
	                               const double ub[],
	                               const Options &options)
	{
		options_ = options;
		return initialize(lb,ub);
	}

	int32_t ASIFrobust::filter(const double x[],
	                           const double uDes[],
	                                 double uAct[])
	{
		double relax;
		updateCost(uDes);
		return filter(x,nullptr,c_,uAct,relax);
	}

	int32_t ASIFrobust::filter(const double x[],
	                           const double uDes[],
	                                 double uAct[],
	                                 double &relax)
	{
		updateCost(uDes);
		return filter(x,nullptr,c_,uAct,relax);
	}

	int32_t ASIFrobust::filter(const double x[],
	                           const double H[],
	                           const double c[],
	                                 double uAct[])
	{
		double relax;
		return filter(x,H,c,uAct,relax);
	}

	int32_t ASIFrobust::filter(const double x[],
	                           const double H[],
	                           const double c[],
	                                 double uAct[],
	                                 double &relax)
	{
		#ifdef ASIF_DEBUG
		std::cout << "Filtering"<< std::endl;
		#endif

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
			relax = solutionFull[nu_];
			return 1;
		}
		else
			return -1;
	}

	
	int32_t ASIFrobust::updateOptions(const Options &options)
	{
		#ifdef ASIF_DEBUG
		std::cout << "Updating options"<< std::endl;
		#endif

		options_ = options;
		return updateOptions();
	}

	int32_t ASIFrobust::updateOptions(void)
	{
		H_[nu_+(nu_*nv_)] = options_.relaxCost;
		c_[nu_] = -2.0*options_.relaxCost*options_.relaxLb;
		lb_[nu_] = options_.relaxLb;
		QPsolver_->updateBounds(lb_,nullptr);	
		QPsolver_->updateCost(H_,c_);
		return 1;
	}

	int32_t ASIFrobust::updateConstraints(const double x[])
	{
		#ifdef ASIF_DEBUG
		std::cout << "Updating constraints"<< std::endl;
		#endif

		// Convert states to interval values
		interval_t xInt[nx_];
		for(uint32_t i=0; i<nx_; i++)
			xInt[i] = interval(x[i]);

		double hFull[npSS_];
		double DhFull[npSS_*nx_];
		safetySet_(x,hFull,DhFull);

		interval_t f[nx_];
		interval_t g[nx_*nu_];
		dynamics_(xInt,f,g);

		double *h;
		double *Dh;
		if(npSSmax_<npSS_)
		{
			h = new double[npSSmax_];
			Dh = new double[npSSmax_*nx_];

			std::vector<uint32_t> idxOfMinH(npSS_);
			for(uint32_t i=0; i<npSS_;i++)
				idxOfMinH[i] = i;

			sort(idxOfMinH.begin(), idxOfMinH.end(), [&hFull](int a, int b) -> bool {return hFull[a]<hFull[b];});
			for(uint32_t i=0; i<npSSmax_; i++)
			{
				h[i] = hFull[idxOfMinH[i]];
				for(uint32_t j=0; j<nx_; j++)
				{
					Dh[i+j*npSSmax_] = DhFull[idxOfMinH[i]+j*npSS_];
				}
			}
		}
		else
		{
			h = hFull;
			Dh = DhFull;
		}

		// Convert Dh to interval values
		interval_t DhInt[npSSmax_*nx_];
		for(uint32_t i=0; i<npSSmax_*nx_; i++)
			DhInt[i] = interval(Dh[i]);

		interval_t Lfh[npSSmax_];
		interval_t Lgh[npSSmax_*nu_];

		// Compute Lfh
		matrixVectorMultiply(DhInt,npSSmax_,nx_,
		                     f,nx_,
		                     Lfh);

		// Compute Lgh
		matrixMultiply(DhInt,npSSmax_,nx_,
		               g,nx_,nu_,
		               Lgh);		

		// Fillup contraint matrices
		uint32_t iColLambda = nu_+1;
		uint32_t inPSS = 0;
		for(uint32_t iRow = 0; iRow<nc_; iRow+=(nu_+2))
		{
			// Fillup A
			A_[iRow+(nu_*nc_)] = h[inPSS];
			for(uint32_t j=0;j<nu_;j++)
			{
				const interval tmp = Lgh[inPSS+(j*npSSmax_)].convert();
				A_[iRow+((iColLambda+j)*nc_)] = tmp.left();
				A_[iRow+((iColLambda+(nu_+1)+j)*nc_)] = -tmp.right();
			}
			const interval tmp = Lfh[inPSS].convert();
			A_[iRow+((iColLambda+nu_)*nc_)] = tmp.left();
			A_[iRow+((iColLambda+(nu_+1)+nu_)*nc_)] = -tmp.right();

			iColLambda+=2*(nu_+1);
			inPSS++;
		}

		if(npSSmax_<npSS_)
		{
			delete[] h;
			delete[] Dh;
		}
		 
		return 1;
	}

	int32_t ASIFrobust::updateCost(const double uDes[])
	{
		#ifdef ASIF_DEBUG
		std::cout << "Updating cost"<< std::endl;
		#endif

		for(uint32_t i=0; i<nu_;i++)
		{
			c_[i] = -2.0*uDes[i];
		}
		return 1;	
	}

	int32_t ASIFrobust::updateH(const double H[])
	{
		#ifdef ASIF_DEBUG
		std::cout << "Building H"<< std::endl;
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

	void ASIFrobust::inputSaturate(double u[])
	{
		for(uint32_t i = 0; i<nu_; i++)
		{
			if(u[i]>ub_[i])
				u[i] = ub_[i];
			else if(u[i]<lb_[i])
				u[i] = lb_[i];
		}
	}
} //end ASIF namespace