#ifndef _QPWRAPPER_OSQP_H_
#define _QPWRAPPER_OSQP_H_

#include <osqp.h>
#include "qpwrapper_abstract.h"

namespace ASIF
{
	class QPWrapperOsqp : public QPWrapperAbstract
	{
	public:
		QPWrapperOsqp(const uint32_t nv,
					     const uint32_t nc,
					     const bool diagonalCost);
		virtual ~QPWrapperOsqp(void);

		virtual int32_t initialize(const double H[],
			                        const double c[],
			                        const double A[],
			                        const double b[],
			                        const double lb[],
			                        const double ub[],
			                        const bool   be[] = nullptr);
		virtual int32_t updateCost(const double H[], const double c[]);
		virtual int32_t updateA(const double A[]);
		virtual int32_t updateb(const double b[]);
		virtual int32_t updateBounds(const double lb[], const double ub[]);
		virtual int32_t solve(void);
		virtual int32_t getSolution(double sol[]);

	protected:
		void buildH(const double H[],
			         const bool valuesOnly = true);
		void buildc(const double c[]);
		void buildA(const double A[]);
		void buildlb(const double b[],
			          const double lb[]);
		void buildlb(const double b[]);
		void buildub(const double ub[]);

		OSQPSettings *settings_;
		OSQPWorkspace *work_;
		OSQPData *data_;

		c_float *H_x_;
		c_int *H_r_;
		c_int *H_c_;
		c_int H_nnz_;

		c_float *c_x_;

		c_float *A_x_;
		c_int *A_r_;
		c_int *A_c_;
		c_int A_nnz_;

		c_float *lb_x_;
		c_float *ub_x_;

		c_int *Ax_new_idx_;
	};
} //end ASIF namespace
#endif
