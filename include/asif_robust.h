#ifndef _ASIF_ROBUST_H_
#define _ASIF_ROBUST_H_

#include "asif_utils.h"
#include "aa.h"

typedef AAF interval_t;

namespace ASIF
{
	class ASIFrobust
	{
	public:
		typedef struct
		{
			double relaxLb = 5.0;
			double relaxCost = 50.0;
			double inf = 1e20;
		} Options;
	public:
		ASIFrobust(const uint32_t nx,
		           const uint32_t nu,
		           const uint32_t npSS,
		           std::function<void(const double* /*x*/,
		                                    double* /*h*/,
		                                    double* /*Dh*/)> safetySet,
		           std::function<void(const interval_t* /*x*/,
		                                    interval_t* /*f*/,
		                                    interval_t* /*g*/)> dynamics,
		           const uint32_t npSSmax = -1,
		           const QPSOLVER qpSolverType = QPSOLVER::OSQP,
		           const bool diagonalCost = true);
		~ASIFrobust(void);

		int32_t initialize(const double lb[],
		                   const double ub[]);
		int32_t initialize(const double lb[],
		                   const double ub[],
		                   const Options &options);

		int32_t filter(const double x[],
		               const double uDes[],
		                     double uAct[]);
		int32_t filter(const double x[],
		               const double uDes[],
		                     double uAct[],
		                     double &relax);
		int32_t filter(const double x[],
		               const double H[],
		               const double c[],
		                     double uAct[]);
		int32_t filter(const double x[],
		               const double H[],
		               const double c[],
		                     double uAct[],
		                     double &relax);

		int32_t updateOptions(void);
		int32_t updateOptions(const Options &options);

	protected:
		int32_t updateConstraints(const double x[]);
		int32_t updateCost(const double uDes[]);
		int32_t updateH(const double H[]);
		void inputSaturate(double u[]);
		const uint32_t nx_;
		const uint32_t nu_;
		const uint32_t npSS_;
		const uint32_t npSSmax_;
		const uint32_t nv_;
		const uint32_t nc_;
		std::function<void(const double* /*x*/,
		                         double* /*h*/,
		                         double* /*Dh*/)> safetySet_;
		std::function<void(const interval_t* /*x*/,
		                         interval_t* /*f*/,
		                         interval_t* /*g*/)> dynamics_;

		Options options_;
		QPWrapperAbstract *QPsolver_;

		double *H_;
		double *c_;
		double *A_;
		double *b_;
		double *lb_;
		double *ub_;

	}; //end ASIFrobust class
} //end ASIF namespace
#endif
