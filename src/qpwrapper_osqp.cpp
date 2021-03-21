#include "qpwrapper_osqp.h"

namespace ASIF
{
	QPWrapperOsqp::QPWrapperOsqp(const uint32_t nv,
	                             const uint32_t nc,
	                             const bool diagonalCost):
	QPWrapperAbstract(nv,nc,diagonalCost)
	{
		settings_ = (OSQPSettings *)c_malloc(sizeof(OSQPSettings));
		data_ = (OSQPData *)c_malloc(sizeof(OSQPData));
		if(diagonalCost)
		{
			H_x_ = new c_float[nv_];
			H_r_ = new c_int[nv_];
		}
		else
		{
			H_x_ = new c_float[nv_*nv_];
			H_r_ = new c_int[nv_*nv_];
		}
		H_c_ = new c_int[nv_+1];
		A_x_ = new c_float[nv_*nc_ + nv_];
		A_r_ = new c_int[nv_*nc_ + nv_];
		A_c_ = new c_int[nv_+1];
		c_x_ = new c_float[nv_];
		lb_x_ = new c_float[nc_ + nv_];
		ub_x_ = new c_float[nc_ + nv_];
		Ax_new_idx_ = new c_int[nv_*nc_];
	}

	QPWrapperOsqp::~QPWrapperOsqp(void)
	{
		osqp_cleanup(work_);
		if (data_)
		{
			c_free(data_->A);
			c_free(data_->P);
			c_free(data_);
		}
		if (settings_)  c_free(settings_);

		delete[] H_x_;
		delete[] H_r_;
		delete[] H_c_;
		delete[] A_x_;
		delete[] A_r_;
		delete[] A_c_;
		delete[] c_x_;
		delete[] lb_x_;
		delete[] ub_x_;
		delete[] Ax_new_idx_;
	}

	int32_t QPWrapperOsqp::initialize(const double H[],
	                                  const double c[],
	                                  const double A[],
	                                  const double b[],
	                                  const double lb[],
	                                  const double ub[],
	                                  const bool   be[])
	{
		#ifdef ASIF_DEBUG
		std::cout << "Initializing" << std::endl;
		#endif

		// Setting default settings
		osqp_set_default_settings(settings_);
		settings_->max_iter = 2000;
		
		// Fill up equality indicator vector
		if(be!=nullptr)
		{
			for(uint32_t i=0; i<nc_; i++)
			{
				be_[i] = be[i];
			}
		}

		// Fill up upper bound
		for(uint32_t i=0; i<nc_;i++)
		{
			ub_x_[i] = OSQP_INFTY;
		}

		// Fill up hessian matrix
		buildH(H,false);

		// Fill up linear cost
		buildc(c);

		// Fill up constraint matrix
		buildA(A);

		// Fill up bounds
		buildlb(b,lb);
		buildub(ub);

		// Prepare A update idx
		uint32_t valueCount = 0;
		for(uint32_t j=0; j<nv_;j++)
		{
			for(uint32_t i=0; i<nc_;i++)
			{
				Ax_new_idx_[i+j*nc_] = valueCount;
				valueCount++;
			}
			valueCount++;
		}

		// Fill up osqp data
		data_->n = nv_;
		data_->m = nc_ + nv_;
		data_->P = csc_matrix(data_->n, data_->n, H_nnz_, H_x_, H_r_, H_c_);
		data_->q = c_x_;
		data_->A = csc_matrix(data_->m, data_->n, A_nnz_, A_x_, A_r_, A_c_);
		data_->l = lb_x_;
		data_->u = ub_x_;
		

		return osqp_setup(&work_, data_, settings_);
	}

	int32_t QPWrapperOsqp::updateCost(const double H[],
	                                  const double c[])
	{

		if(H!=nullptr)
		{
			#ifdef ASIF_DEBUG
			std::cout << __func__ << ": ";
			std::cout << "Updating H" << std::endl;
			#endif

			buildH(H);
			if(diagonalCost_)
			{
				osqp_update_P(work_,H_x_,OSQP_NULL,nv_);
			}
			else
			{
				double Htemp[(nv_*(nv_+1))/2];
				uint32_t HtempIdx = 0;
				for(uint32_t j=0; j<nv_;j++)
				{
					for(uint32_t i=0; i<=j;i++)
					{
						Htemp[HtempIdx] = H_x_[j*nv_+i];
						HtempIdx++;
					}
				}
				osqp_update_P(work_,Htemp,OSQP_NULL,(nv_*(nv_+1))/2);			
			}

		}

		if(c!=nullptr)
		{
			#ifdef ASIF_DEBUG
			std::cout << __func__ << ": ";
			std::cout << "Updating c" << std::endl;
			buildc(c);
			#endif
			osqp_update_lin_cost(work_,c);
		}
		
		return 1;
	}

	int32_t QPWrapperOsqp::updateA(const double A[])
	{	
		#ifdef ASIF_DEBUG
		std::cout << __func__ << ": ";
		std::cout << "Updating A" << std::endl;
		buildA(A);
		#endif
		osqp_update_A(work_, A, Ax_new_idx_, nv_*nc_);
		return 1;
	}

	int32_t QPWrapperOsqp::updateb(const double b[])
	{
		#ifdef ASIF_DEBUG
		std::cout << __func__ << ": ";
		std::cout << "Updating b" << std::endl;
		#endif
		buildlb(b);
		osqp_update_lower_bound(work_,lb_x_);
		return 1;
	}

	int32_t QPWrapperOsqp::updateBounds(const double lb[],
		                                 const double ub[])
	{
		if(lb!=nullptr)
		{
			#ifdef ASIF_DEBUG
			std::cout << __func__ << ": ";
			std::cout << "Updating lb" << std::endl;
			#endif
			buildlb(nullptr,lb);
			osqp_update_lower_bound(work_,lb_x_);
		}

		if(ub!=nullptr)
		{		
			#ifdef ASIF_DEBUG
			std::cout << __func__ << ": ";
			std::cout << "Updating ub" << std::endl;
			#endif
			buildub(ub);
			osqp_update_upper_bound(work_,ub_x_);
		}
		return 1;
	}

	int32_t QPWrapperOsqp::solve(void)
	{

		#ifdef ASIF_DEBUG
		std::cout << "Solving... "<< std::endl;
		#endif
		osqp_solve(work_);

		if(work_->info->status_val==OSQP_SOLVED || work_->info->status_val==OSQP_SOLVED_INACCURATE)
		{
			#ifdef ASIF_DEBUG
			std::cout << "...Successfull" << std::endl;
			#endif	
			return (int32_t)QPWrapperAbstract::SOLVER_STATUS::FEASIBLE;
		}
		else
		{	
			#ifdef ASIF_DEBUG
			std::cout << "...Failed" << std::endl;
			#endif
			return (int32_t)work_->info->status_val;
		}
	}

	int32_t QPWrapperOsqp::getSolution(double sol[])
	{
		#ifdef ASIF_DEBUG
		std::cout << "Getting solution : ";
		#endif

		for(uint32_t i = 0; i<nv_;i++)
		{
			sol[i] = static_cast<double>(work_->solution->x[i]);
		}

		#ifdef ASIF_DEBUG
		for(uint32_t i = 0; i<nv_;i++)
		{
			std::cout << sol[i] << " ";
		}
		std::cout << std::endl;
		#endif

		return 1;
	}

	void QPWrapperOsqp::buildH(const double H[],
	                           const bool valuesOnly)
	{

		if(diagonalCost_)
		{
			for(uint32_t i=0; i<nv_;i++)
			{
				H_x_[i] = 2.0*H[i+(i*nv_)];
			}

			if(!valuesOnly)
			{
				H_nnz_ = nv_;
				for(uint32_t i=0; i<nv_;i++)
				{
					H_r_[i] = i;
					H_c_[i] = i;
				}
				H_c_[nv_] = nv_;
			}
		}
		else
		{
			for(uint32_t j=0; j<nv_;j++)
			{
				for(uint32_t i=0; i<nv_;i++)
				{
					H_x_[i+(j*nv_)] = 2.0*H[i+(j*nv_)];
				}	
			}

			if(!valuesOnly)
			{
				H_nnz_ = nv_*nv_;
				for(uint32_t j=0; j<nv_;j++)
				{
					for(uint32_t i=0; i<nv_;i++)
					{
						H_r_[i+(j*nv_)] = i;
					}	
					H_c_[j] = j*nv_;
				}
				H_c_[nv_] = nv_*nv_;
			}
		}
	}

	void QPWrapperOsqp::buildc(const double c[])
	{
		for(uint32_t i=0; i<nv_;i++)
		{
			c_x_[i] = c[i];
		}
	}

	void QPWrapperOsqp::buildA(const double A[])
	{
		A_nnz_ = nv_*nc_ + nv_;
		uint32_t valueCount = 0;
		for(uint32_t j=0; j<nv_;j++)
		{
			A_c_[j] = valueCount;
			for(uint32_t i=0; i<nc_;i++)
			{
				A_x_[valueCount] = A[i+j*nc_];
				A_r_[valueCount] = i;
				valueCount++;
			}
			for(uint32_t i=0; i<nv_;i++)
			{
				if(i==j)
				{
					A_x_[valueCount] = 1.0;		
					A_r_[valueCount] = i+nc_;
					valueCount++;
				}
			}	
		}
		A_c_[nv_] = valueCount;
	}


	void QPWrapperOsqp::buildlb(const double b[],
			        			       const double lb[])
	{
		if(b!=nullptr)
		{
			buildlb(b);
		}
		for(uint32_t i=nc_; i<(nc_+nv_);i++)
		{
			lb_x_[i] = lb[i-nc_];
		}
	}

	void QPWrapperOsqp::buildlb(const double b[])
	{
		for(uint32_t i=0; i<nc_;i++)
		{
			lb_x_[i] = b[i];
			if(be_[i])
				ub_x_[i] = b[i];
		}

	}

	void QPWrapperOsqp::buildub(const double ub[])
	{
		for(uint32_t i=nc_; i<(nc_+nv_);i++)
		{
			ub_x_[i] = ub[i-nc_];
		}
	}
	
} //end ASIF namespace
