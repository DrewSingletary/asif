#include "asif_realizable.h"

namespace ASIF
{
	ASIFrealizable::ASIFrealizable(const uint32_t nx,
	                               const uint32_t nu,
	                               const double uncertaintyBounds[],
											 const kernel_t &kernel,
	                               std::function<void(const interval_t* /*x*/,
	                                                        interval_t* /*f*/,
	                                                        interval_t* /*g*/)> dynamics,
	                               const uint32_t npSSmax,
	                               const QPSOLVER qpSolverType,
	                               const bool diagonalCost):
	nx_(nx),
	nu_(nu),
	kernel_(kernel),
	nFacets_(kernel.facets.size()),
	npSS_(kernel.maxCriticalFacets*kernel.maxActiveConstraints),
	npSSmax_((npSSmax>nFacets_) ? nFacets_ : npSSmax),
	nv_((npSSmax_>0) ? (nu+npSS_*2*(nu+1)+1) : (nu+npSS_*2*(nu+1))),
	nc_(npSS_*(nu+2) + npSSmax_),
	dynamics_(dynamics),
	options_(),
	criticalFacets_(kernel.maxCriticalFacets),
	nCriticalFacets_(0),
	criticalBarrierFacets_(),
	hBarrier_(),
	DhBarrier_()
	{
		H_ = new double[nv_*nv_]{0.0};
		c_ = new double[nv_]{0.0};
		A_ = new double[nc_*nv_]{0.0};
		b_ = new double[nc_]{0.0};
		lb_ = new double[nv_]{0.0};
		ub_ = new double[nv_]{0.0};
		be_ = new bool[nc_]{false};

		A_facet_ = new double[(2*nx_+1)*nx_]{0.0};
		b_facet_ = new double[2*nx_+1]{0.0};

		uncertaintyBounds_ = new double[nx_]{0.0};

		std::copy(uncertaintyBounds,uncertaintyBounds+nx_,uncertaintyBounds_);

		switch(qpSolverType)
		{
			case QPSOLVER::OSQP:
			{
				QPsolver_ = new QPWrapperOsqp(nv_,nc_,diagonalCost);
				facetSolver_ = new QPWrapperOsqp(nx_,2*nx_+1,true);
				break;
			}
			default :
			{
				QPsolver_ = nullptr;
				facetSolver_ = nullptr;
				break;
			}
		}

		if(npSSmax_>0)
		{
			criticalBarrierFacets_.resize(npSSmax_);
			hBarrier_.resize(npSSmax_);
			DhBarrier_.resize(npSSmax_*nx_);
		}
		if(npSSmax==nFacets_)
		{
			for(uint32_t i=0; i<nFacets_; i++)
				criticalBarrierFacets_[i] = i;
		}

	}

	ASIFrealizable::~ASIFrealizable(void)
	{
		delete[] H_;
		delete[] c_;
		delete[] A_;
		delete[] b_;
		delete[] lb_;
		delete[] ub_;
		delete[] be_;

		delete[] A_facet_;
		delete[] b_facet_;

		delete[] uncertaintyBounds_;
		
		if(QPsolver_!=nullptr)
			delete QPsolver_;

		if(facetSolver_!=nullptr)
			delete facetSolver_;
	}

	int32_t ASIFrealizable::initialize(const double lb[],
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

		// Initialize facetSolver
		double H_facet[nx_*nx_] = {0.0};
		for(uint32_t j=0; j<nx_;j++)
		{
			for(uint32_t i=0; i<nx_;i++)
			{
				if(i==j)
					H_facet[i+(j*nx_)] = 1.0;
			}
		}
		double c_facet[nx_] = {0.0};
		double lb_facet[nx_] = {0.0};
		double ub_facet[nx_];
		std::fill(ub_facet,ub_facet+nx_,1.0);
		bool be_facet[2*nx_+1] = {0};
		be_facet[0] = true;
		for(uint32_t i=0; i<nx_; i++)
		{
			A_facet_[i*(2*nx_+1)] = 1.0;
		}
		b_facet_[0] = 1.0;

		facetSolver_->initialize(H_facet, c_facet, A_facet_, b_facet_, lb_facet, ub_facet,be_facet);
		facetSolver_->solve();

		// Initialize xFaceInt
		const std::vector<std::vector<double>> &vertices = kernel_.vertices;
		for(uint32_t i=0; i<nFacets_; i++)
		{
			std::vector<interval_t> &xFaceInt = kernel_.facets[i].xFaceInt;
			xFaceInt.resize(nx_);
			const std::vector<uint32_t> &verticesIdx = kernel_.facets[i].verticesIdx;

			for(uint32_t j=0; j<nx_; j++)
				xFaceInt[j] = vertices[verticesIdx[0]][j];

			const uint32_t nLambdas = nx_-1;

			for(uint32_t j=1; j<=nLambdas; j++)
			{
				interval_t lamdaCurr = interval_t(0.,1.);
				const std::vector<double> &vertex = vertices[verticesIdx[j]];
				for(uint32_t k=0; k<nx_; k++)
					xFaceInt[k] = lamdaCurr*xFaceInt[k] + (1.-lamdaCurr)*vertex[k];
			}
		}

		// Initialize Bounding Box
		for(uint32_t i=0; i<nFacets_; i++)
		{
			std::vector<std::pair<double,double>> &boundingBox = kernel_.facets[i].boundingBox;
			boundingBox.resize(nx_);
			const std::vector<uint32_t> &verticesIdx = kernel_.facets[i].verticesIdx;
			double verticesComponents[nx_];
			for(uint32_t j=0; j<nx_; j++)
			{
				for(uint32_t k=0; k<nx_; k++)
					verticesComponents[k] = vertices[verticesIdx[k]][j];

				boundingBox[j].first = *std::min_element(verticesComponents,verticesComponents+nx_);
				boundingBox[j].second = *std::max_element(verticesComponents,verticesComponents+nx_);
			}
		}

		// Initialize cost hessian
		double Hidentity[nu_*nu_] = {0.0};
		for(uint32_t j=0; j<nu_;j++)
		{
			for(uint32_t i=0; i<nu_;i++)
			{
				if(i==j)
					Hidentity[i+(j*nu_)] = 1.0;
			}
		}
		updateH(Hidentity);
		if(npSSmax_>0)
			H_[nv_*nv_-1] = options_.relaxCost;

		// Initialize bounds
		memcpy(lb_,lb,nu_*sizeof(double));
		memcpy(ub_,ub,nu_*sizeof(double));
		for(uint32_t i=nu_; i<nv_; i++)
		{
			lb_[i] = 0.0;
			ub_[i] = options_.inf;
		}

		// Initialize contraint matrices
		uint32_t iColLambda = nu_;
		for(uint32_t iRow = 0; iRow<(npSS_*(nu_+2)); iRow+=(nu_+2))
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

		// Initialize contraint matrix structure
		std::fill(be_,be_+(npSS_*(nu_+2)),true);
		for(uint32_t i=0; i<(npSS_*(nu_+2)); i+=(nu_+2))
			be_[i] = false;

		// Initialize constraints
		double x0[nx_] = {0.0};
		updateConstraints(x0);

		// Initialize linear cost
		double uDes0[nu_] = {0.0};
		updateCost(uDes0);



		QPsolver_->initialize(H_, c_, A_, b_, lb_, ub_, be_);
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

	int32_t ASIFrealizable::initialize(const double lb[],
	                                   const double ub[],
	                                   const Options &options)
	{
		options_ = options;
		return initialize(lb,ub);
	}

	int32_t ASIFrealizable::filter(const double x[],
	                               const double uDes[],
	                                     double uAct[])
	{
		double relax[2];
		updateCost(uDes);
		return filter(x,nullptr,c_,uAct,relax);
	}

	int32_t ASIFrealizable::filter(const double x[],
	                               const double uDes[],
	                                     double uAct[],
	                                     double relax[2])
	{
		updateCost(uDes);
		return filter(x,nullptr,c_,uAct,relax);
	}

	int32_t ASIFrealizable::filter(const double x[],
	                               const double H[],
	                               const double c[],
	                                     double uAct[])
	{
		double relax[2];
		return filter(x,H,c,uAct,relax);
	}

	int32_t ASIFrealizable::filter(const double x[],
	                               const double H[],
	                               const double c[],
	                                     double uAct[],
	                                     double relax[2])
	{
		#ifdef ASIF_DEBUG
		std::cout << "Filtering"<< std::endl;
		#endif

		int32_t rtCode = updateConstraints(x);
		if(rtCode<0)
			return -2;

		if(H==nullptr)
			QPsolver_->updateCost(nullptr,c);
		else
		{
			updateH(H);
			QPsolver_->updateCost(H_,c);
		}

		QPsolver_->updateA(A_);
		QPsolver_->updateb(b_);		
		rtCode = QPsolver_->solve();

		if(static_cast<QPWrapperAbstract::SOLVER_STATUS>(rtCode)==QPWrapperAbstract::SOLVER_STATUS::FEASIBLE)
		{
			double solutionFull[nv_];
			QPsolver_->getSolution(solutionFull);
			memcpy(uAct,solutionFull,nu_*sizeof(double));
			inputSaturate(uAct);
			relax[0] = solutionFull[nu_];
			relax[1] = solutionFull[nv_-1];
			return 1;
		}
		else
			return -1;
	}

	
	int32_t ASIFrealizable::updateOptions(const Options &options)
	{
		#ifdef ASIF_DEBUG
		std::cout << "Updating options"<< std::endl;
		#endif

		options_ = options;
		return updateOptions();
	}

	int32_t ASIFrealizable::updateOptions(void)
	{
		if(npSSmax_>0)
			H_[nv_*nv_-1] = options_.relaxCost;

		QPsolver_->updateBounds(lb_,nullptr);	
		QPsolver_->updateCost(H_,c_);
		return 1;
	}

	int32_t ASIFrealizable::updateConstraints(const double x[])
	{
		#ifdef ASIF_DEBUG
		std::cout << "Updating constraints"<< std::endl;
		#endif

		// Find critical facets
		nCriticalFacets_ = 0;
		for(uint32_t i=0; i<nx_; i++)
		{
			b_facet_[1+i] = x[i] - uncertaintyBounds_[i];
			b_facet_[1+i+nx_] = -x[i] - uncertaintyBounds_[i];
		}
		facetSolver_->updateb(b_facet_);

		//Evaluate safety set
		double hFull[nFacets_];
		for(uint32_t i=0; i<nFacets_; i++)
		{
			const std::vector<double> &normal = kernel_.facets[i].normal;
			hFull[i] = 1.;
			for(uint32_t j=0; j<nx_; j++)
				hFull[i] -= normal[j]*x[j];
		}

		const std::vector<std::vector<double>> &vertices = kernel_.vertices;
		for(uint32_t i=0; i<nFacets_; i++)
		{
			const std::vector<uint32_t> &verticesIdx = kernel_.facets[i].verticesIdx;
			const std::vector<std::pair<double,double>> &boundingBox = kernel_.facets[i].boundingBox;
			bool potentialFacet = true;

			for(uint32_t j=0; j<nx_; j++)
			{
				if(x[j] < (boundingBox[j].first - uncertaintyBounds_[j])  ||
					x[j] > (boundingBox[j].second + uncertaintyBounds_[j]))
				{
					potentialFacet = false;
					break;
				}
			}

			if(potentialFacet)
			{
				for(uint32_t j=0; j<nx_; j++)
				{
					interval_t lamdaCurr = interval_t(0.,1.);
					const std::vector<double> &vertex = vertices[verticesIdx[j]];

					for(uint32_t k=0; k<nx_; k++)
					{
						A_facet_[1+k+(j*(2*nx_+1))] = vertex[k];
						A_facet_[1+k+nx_+(j*(2*nx_+1))] = -vertex[k];
					}
				}

				facetSolver_->updateA(A_facet_);
				int32_t rtCode = facetSolver_->solve();
				if(static_cast<QPWrapperAbstract::SOLVER_STATUS>(rtCode)==QPWrapperAbstract::SOLVER_STATUS::FEASIBLE)
				{
					criticalFacets_[nCriticalFacets_] = i;
					nCriticalFacets_++;
					if(nCriticalFacets_>=kernel_.maxCriticalFacets)
						break;
				}
			}
		}

		// Generate critical constraints
		interval_t Lfh[npSS_];
		interval_t Lgh[npSS_*nu_];
		std::fill(Lfh,Lfh+npSS_,0.);
		std::fill(Lgh,Lgh+(npSS_*nu_),0.);
		if(nCriticalFacets_>0)
		{
			const uint32_t DhNrows = nCriticalFacets_*kernel_.maxActiveConstraints;
			double Dh[DhNrows*nx_] = {0.0};
			uint32_t criticalFacetMap[DhNrows];

			// Compute Dh
			uint32_t activeConstraintsTotal = 0;
			for(uint32_t i=0; i<nCriticalFacets_; i++)
			{
				const uint32_t &iFacet = criticalFacets_[i];
				const std::vector<uint32_t>  &activeConstraintsSet = kernel_.facets[iFacet].activeConstraintsSet;

				for(uint32_t j=0; j<activeConstraintsSet.size(); j++)
				{
					criticalFacetMap[activeConstraintsTotal] = i;
					const std::vector<double> &normal = kernel_.facets[activeConstraintsSet[j]].normal;
					for(uint32_t k=0; k<nx_; k++)
					{
						Dh[activeConstraintsTotal+k*DhNrows] = -normal[k];
					}
					activeConstraintsTotal++;
				}
			}

			// Convert Dh to interval values
			interval_t DhInt[activeConstraintsTotal*nx_];
			for(uint32_t i=0; i<activeConstraintsTotal; i++)
			{
				for(uint32_t j=0; j<nx_; j++)
					DhInt[i+j*activeConstraintsTotal] = interval(Dh[i+j*DhNrows]);
			}

			// Compute Lfh and Lgh
			for(uint32_t i=0; i<activeConstraintsTotal; i++)
			{
				const uint32_t &iFacet = criticalFacets_[criticalFacetMap[i]];
				const std::vector<interval_t> &xFaceInt = kernel_.facets[iFacet].xFaceInt;
				interval_t f[nx_];
				interval_t g[nx_*nu_];
				dynamics_(xFaceInt.data(),f,g);

				Lfh[i] = 0.;
				for(uint32_t j=0; j<nx_; j++)
					Lfh[i] = Lfh[i] + f[j]*DhInt[i+(j*activeConstraintsTotal)];

				for(uint32_t j=0; j<nu_; j++)
				{
					Lgh[i+j*npSS_] = 0.;
					for(uint32_t k=0; k<nx_; k++)
						Lgh[i+j*npSS_] = Lgh[i+j*npSS_]+ g[k+j*nx_]*DhInt[i+(k*activeConstraintsTotal)];
				}
			}
		}

		// Fillup contraint matrices
		uint32_t iColLambda = nu_;
		uint32_t inPSS = 0;
		for(uint32_t iRow = 0; iRow<(npSS_*(nu_+2)); iRow+=(nu_+2))
		{
			// Fillup A
			for(uint32_t j=0;j<nu_;j++)
			{
				const interval tmp = Lgh[inPSS+(j*npSS_)].convert();
				A_[iRow+((iColLambda+j)*nc_)] = tmp.left();
				A_[iRow+((iColLambda+(nu_+1)+j)*nc_)] = -tmp.right();
			}
			const interval tmp = Lfh[inPSS].convert();
			A_[iRow+((iColLambda+nu_)*nc_)] = tmp.left();
			A_[iRow+((iColLambda+(nu_+1)+nu_)*nc_)] = -tmp.right();

			iColLambda+=2*(nu_+1);
			inPSS++;
		}

		// Add barrier constraints
		if(npSSmax_>0)
		{
			// Convert states to interval values
			interval_t xInt[nx_];
			for(uint32_t i=0; i<nx_; i++)
				xInt[i] = interval(x[i]);

			double DhFull[nFacets_*nx_];
			double f[nx_];
			double g[nx_*nu_];
			interval_t fInt[nx_];
			interval_t gInt[nx_*nu_];

			for(uint32_t i=0; i<nFacets_; i++)
			{
				const std::vector<double> &normal = kernel_.facets[i].normal;
				for(uint32_t j=0; j<nx_; j++)
					DhFull[i+j*nFacets_] = -normal[j];
			}
			dynamics_(xInt,fInt,gInt);
			for(uint32_t i=0; i<nx_; i++)
			{
				f[i] = fInt[i].convert().mid();
				for(uint32_t j=0; j<nu_; j++)
					g[i+j*nx_] = gInt[i+j*nx_].convert().mid();
			}

			double *hBarrier;
			double *DhBarrier;

			if(npSSmax_<nFacets_)
			{
				hBarrier = hBarrier_.data();
				DhBarrier = DhBarrier_.data();

				std::vector<uint32_t> idxOfMinH(nFacets_);
				for(uint32_t i=0; i<nFacets_;i++)
					idxOfMinH[i] = i;

				sort(idxOfMinH.begin(), idxOfMinH.end(), [&hFull](int a, int b) -> bool {return hFull[a]<hFull[b];});
				for(uint32_t i=0; i<npSSmax_; i++)
				{
					criticalBarrierFacets_[i] = idxOfMinH[i];
					hBarrier[i] = hFull[idxOfMinH[i]];
					for(uint32_t j=0; j<nx_; j++)
					{
						DhBarrier[i+j*npSSmax_] = DhFull[idxOfMinH[i]+j*nFacets_];
					}
				}
			}
			else
			{
				std::copy(hFull,hFull+nFacets_,hBarrier_.begin());
				hBarrier = hFull;
				DhBarrier = DhFull;
			}

			double Lfh[npSSmax_];
			double Lgh[npSSmax_*nu_];

			// Compute Lfh
			matrixVectorMultiply(DhBarrier,npSSmax_,nx_,
			                     f,nx_,
			                     Lfh);

			// Compute Lgh
			matrixMultiply(DhBarrier,npSSmax_,nx_,
			               g,nx_,nu_,
			               Lgh);		

			// Fillup A and b
			for(uint32_t i=0;i<npSSmax_;i++)
			{
				for(uint32_t j=0;j<nu_;j++)
				{
					A_[(npSS_*(nu_+2))+i+(j*nc_)] = Lgh[i+(j*npSSmax_)];
				}
				A_[(npSS_*(nu_+2))+i+((nv_-1)*nc_)] = 1.0;
				b_[(npSS_*(nu_+2))+i] = -Lfh[i] -options_.relaxDes*(hBarrier[i] - options_.relaxOffset);
			}
		}

		if(nCriticalFacets_==0 && std::any_of(hFull, hFull+nFacets_, [](double v){return v<0.;}))
			return -1;
		else
			return 1;
		
		#ifdef ASIF_DEBUG
		std::cout << "x: " << std::endl;
		for(uint32_t i=0; i<nx_;i++)
			std::cout << x[i] << std::endl;
		std::cout << std::endl;

		std::cout << "nCriticalFacets:" << nCriticalFacets_ << std::endl;
		std::cout << "Critical Facets:" << std::endl;
		for(uint32_t i=0; i<nCriticalFacets_;i++)
			std::cout << criticalFacets_[i] << std::endl;
		std::cout << std::endl;

		std::cout << "criticalBarrierFacets and value: " << std::endl;
		for(uint32_t i=0; i<npSSmax_;i++)
			std::cout << criticalBarrierFacets_[i] << ", " << hBarrier_[i] << std::endl;
		std::cout << std::endl;

		std::cout << "nCriticalFacets: " << nCriticalFacets_ << std::endl;
		std::cout << "H:" << std::endl;
		for(uint32_t i=0; i<nv_; i++)
		{
			for(uint32_t j=0; j<nv_; j++)	
			{
				std::cout << H_[i+j*nv_] << ", ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;


		std::cout << "c:" << std::endl;
		for(uint32_t i=0; i<nv_; i++)
		{
				std::cout << c_[i] << ", ";
		}
		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "lb:" << std::endl;
		for(uint32_t i=0; i<nv_; i++)
		{
			std::cout << lb_[i] << ", ";
		}
		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "ub:" << std::endl;
		for(uint32_t i=0; i<nv_; i++)
		{
			std::cout << ub_[i] << ", ";
		}
		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "A, b " << std::endl;
		for(uint32_t i=0; i<nc_; i++)
		{
			for(uint32_t j=0; j<nv_; j++)	
			{
				std::cout << A_[i+j*nc_] << ", ";
			}
			if(be_[i])
				std::cout << "= ";
			else
				std::cout << "> ";
			std::cout << b_[i] << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
		#endif
	}

	int32_t ASIFrealizable::updateCost(const double uDes[])
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

	int32_t ASIFrealizable::updateH(const double H[])
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

	void ASIFrealizable::inputSaturate(double u[])
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