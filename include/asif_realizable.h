#ifndef _ASIF_REALIZABLE_H_
#define _ASIF_REALIZABLE_H_

#include "asif_utils.h"
#include "aa.h"

typedef AAF interval_t;

namespace ASIF
{
	class ASIFrealizable
	{
	public:
		typedef struct
		{
			double relaxDes = 5.0;
			double relaxOffset = 5.0;
			double relaxCost = 50.0;
			double inf = 1e20;
		} Options;

		typedef struct
		{
			std::vector<uint32_t>                 verticesIdx;
			std::vector<double>                   normal;
			std::vector<uint32_t>                 activeConstraintsSet;
			std::vector<interval_t>               xFaceInt;
			std::vector<std::pair<double,double>> boundingBox;
		} facet_t;

		typedef struct
		{
			std::vector<std::vector<double>> vertices;
			std::vector<facet_t> facets;
			uint32_t maxCriticalFacets;
			uint32_t maxActiveConstraints;
		} kernel_t;

	public:
		ASIFrealizable(const uint32_t nx,
		               const uint32_t nu,
		               const double uncertaintyBounds[],
		               const kernel_t &kernel,
		               std::function<void(const interval_t* /*x*/,
		                                        interval_t* /*f*/,
		                                        interval_t* /*g*/)> dynamics,
		               const uint32_t npSSmax = 0,
		               const QPSOLVER qpSolverType = QPSOLVER::OSQP,
		               const bool diagonalCost = true);
		~ASIFrealizable(void);

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
		                     double relax[2]);
		int32_t filter(const double x[],
		               const double H[],
		               const double c[],
		                     double uAct[]);
		int32_t filter(const double x[],
		               const double H[],
		               const double c[],
		                     double uAct[],
		                     double relax[2]);

		int32_t updateOptions(void);
		int32_t updateOptions(const Options &options);

	protected:
		int32_t updateConstraints(const double x[]);
		int32_t updateCost(const double uDes[]);
		int32_t updateH(const double H[]);
		void inputSaturate(double u[]);



		const uint32_t nx_;
		const uint32_t nu_;
		double *uncertaintyBounds_;
		kernel_t kernel_;
		const uint32_t nFacets_;
		const uint32_t npSS_;
		const uint32_t npSSmax_;
		const uint32_t nv_;
		const uint32_t nc_;
		
		std::function<void(const interval_t* /*x*/,
		                         interval_t* /*f*/,
		                         interval_t* /*g*/)> dynamics_;

		Options options_;
		QPWrapperAbstract *QPsolver_;
		QPWrapperAbstract *facetSolver_;

		double *H_;
		double *c_;
		double *A_;
		double *b_;
		double *lb_;
		double *ub_;
		bool *be_;

		double *A_facet_;
		double *b_facet_;

	public:
		std::vector<uint32_t> criticalFacets_;
		uint32_t nCriticalFacets_;
		std::vector<uint32_t> criticalBarrierFacets_;
		std::vector<double> hBarrier_;
		std::vector<double> DhBarrier_;
	}; //end ASIFrobust class
} //end ASIF namespace
#endif
