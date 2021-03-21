#ifndef _ASIF_IMPLICIT_ROBUST_H_
#define _ASIF_IMPLICIT_ROBUST_H_

#include "asif_utils.h"
#include "asif_learning_utils.h"
#include "aa.h"

typedef AAF interval_t;
#ifdef USE_ODEINT
	#include <boost/numeric/odeint.hpp>
	#include <boost/numeric/odeint/iterator/n_step_iterator.hpp>
	using namespace boost::numeric::odeint;
#endif

using namespace std;

namespace ASIF
{
	class ASIFimplicitRB
	{
	public:
		typedef struct
		{
			double *x0 = nullptr;
			double *x_unc = nullptr;
			int n_debug = -1;
			double relaxCost = 50.0;
			double relaxReachLb = 5.0;
			double relaxSafeLb = 5.0;
			double backTrajHorizon = 1.0;
			double backContDt = 0.01;
			double backTrajDt = 0.01;
			double backTrajAbsTol = 1.0e-6;
			double backTrajRelTol = 1.0e-6;
			double satSharpness = 0.1;
			double inf = 1e20;
			bool use_learning = false;
		} Options;

		typedef vector<double> state_t;

		#ifdef USE_ODEINT
			typedef runge_kutta_dopri5<state_t> stepper_t;
		#endif

	public:
		ASIFimplicitRB(const uint32_t nx,
		             const uint32_t nu,
		             const uint32_t npSS,
		             const uint32_t npBS,
		             const uint32_t npBTSS,
								 function<void(const double* /*x*/,
		                                 double* /*h*/,
		                                 double* /*Dh*/)> safetySet,
		             function<void(const interval_t* /*x*/,
		                                 interval_t* /*h*/,
		                                 interval_t* /*Dh*/)> safetySet_int,
		             function<void(const double* /*x*/,
		                                 double* /*h*/,
		                                 double* /*Dh*/)> backupSet,
								 function<void(const interval_t* /*x*/,
		                                 interval_t* /*h*/,
		                                 interval_t* /*Dh*/)> backupSet_int,
								 function<void(const double* /*x*/,
		                                 double* /*f*/,
		                                 double* /*g*/)> dynamics,
								 function<void(const interval_t* /*x*/,
		                                 interval_t* /*f*/,
		                                 interval_t* /*g*/)> dynamics_int,
								 function<void(const double* /*x*/,
								 										 double* /*Df*/,
								 										 double* /*Dg*/)> dynamicsGradients,
								 function<void(const interval_t* /*x*/,
		                                 interval_t* /*Df*/,
		                                 interval_t* /*Dg*/)> dynamicsGradients_int,
		             function<void(const double* /*x*/,
		                                 double* /*u*/,
		                                 double* /*Du*/)> backupController,
								 const QPSOLVER qpSolverType = QPSOLVER::OSQP,
		             const bool diagonalCost = true);

		ASIFimplicitRB(const uint32_t nx,
		             const uint32_t nu,
		             const uint32_t npSS,
		             const uint32_t npBS,
		             const uint32_t npBTSS,
								 function<void(const double* /*x*/,
								 										 double* /*h*/,
								 										 double* /*Dh*/)> safetySet,
								 function<void(const interval_t* /*x*/,
		                                 interval_t* /*h*/,
		                                 interval_t* /*Dh*/)> safetySet_int,
								 function<void(const double* /*x*/,
								 										 double* /*h*/,
								 										 double* /*Dh*/)> backupSet,
								 function<void(const interval_t* /*x*/,
		                                 interval_t* /*h*/,
		                                 interval_t* /*Dh*/)> backupSet_int,
								 function<void(const double* /*x*/,
		                           const double* /*u*/,
		                                 double* /*f*/,
		                                 double* /*g*/,
		                                 double* /*d_fcl_dx*/)> dynamicsWithGradient,
		             function<void(const interval_t* /*x*/,
		                           const double* /*u*/,
		                                 interval_t* /*f*/,
		                                 interval_t* /*g*/,
		                                 interval_t* /*d_fcl_dx*/)> dynamicsWithGradient_int,
		             function<void(const double* /*x*/,
		                                 double* /*u*/,
		                                 double* /*Du*/)> backupController,
		             const QPSOLVER qpSolverType = QPSOLVER::OSQP,
		             const bool diagonalCost = true);

		~ASIFimplicitRB(void);

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

		int32_t updateOptions(const Options &options);

		vector<pair<double,state_t>> backTraj_;
		vector<uint32_t> backTrajCritIdx_;
		double hBackupEnd_;
		double hSafetyNow_;
		int index_debug_;
		std::vector<double> Dh_index_;
		std::vector<double> h_index_;
		LearningData learning_data_;
	protected:
		ASIFimplicitRB(const bool dynamicsHasGradient,
		             const uint32_t nx,
		             const uint32_t nu,
		             const uint32_t npSS,
		             const uint32_t npBS,
		             const uint32_t npBTSS,
								 function<void(const double* /*x*/,
								 									 	 double* /*h*/,
								 										 double* /*Dh*/)> safetySet,
								 function<void(const interval_t* /*x*/,
		                                 interval_t* /*h*/,
		                                 interval_t* /*Dh*/)> safetySet_int,
								 function<void(const double* /*x*/,
		                                 double* /*h*/,
		                                 double* /*Dh*/)> backupSet,
		             function<void(const interval_t* /*x*/,
		                                 interval_t* /*h*/,
		                                 interval_t* /*Dh*/)> backupSet_int,
								 function<void(const double* /*x*/,
		                                 double* /*f*/,
		                                 double* /*g*/)> dynamics,
		             function<void(const interval_t* /*x*/,
		                                 interval_t* /*f*/,
		                                 interval_t* /*g*/)> dynamics_int,
								 function<void(const double* /*x*/,
		                                 double* /*Df*/,
		                                 double* /*Dg*/)> dynamicsGradients,
		             function<void(const interval_t* /*x*/,
		                                 interval_t* /*Df*/,
		                                 interval_t* /*Dg*/)> dynamicsGradients_int,
								 function<void(const double* /*x*/,
		                           const double* /*u*/,
		                                 double* /*f*/,
		                                 double* /*g*/,
		                                 double* /*d_fcl_dx*/)> dynamicsWithGradient,
		             function<void(const interval_t* /*x*/,
		                           const double* /*u*/,
		                                 interval_t* /*f*/,
		                                 interval_t* /*g*/,
		                                 interval_t* /*d_fcl_dx*/)> dynamicsWithGradient_int,
		             function<void(const double* /*x*/,
		                                 double* /*u*/,
		                                 double* /*Du*/)> backupController,
		             const QPSOLVER qpSolverType = QPSOLVER::OSQP,
		             const bool diagonalCost = true);

		int32_t updateOptions(void);
		int32_t updateConstraints(const double x[]);
		int32_t updateCost(const double uDes[]);
		int32_t updateH(const double H[]);
		void inputSaturateSoft(const double u[],
		                             double uSat[],
		                             double DuSat[]);
		void inputSaturate(double u[]);
		void backupCLdynamics(const double x[],
		                            double f[],
		                            double Df[],
																double t);
		void ODErhs(const state_t &x,
		                  state_t &xDot,
		            const double t);

		const bool dynamicsHasGradient_;
		const uint32_t nx_;
		const uint32_t nu_;
		const uint32_t nv_;
		const uint32_t npSS_;
		const uint32_t npBS_;
		const uint32_t npBTSS_;
		const uint32_t npTC_;
		function<void(const double* /*x*/,
		                    double* /*h*/,
		                    double* /*Dh*/)> safetySet_;
		function<void(const interval_t* /*x*/,
												interval_t* /*h*/,
												interval_t* /*Dh*/)> safetySet_int_;
		function<void(const double* /*x*/,
												double* /*h*/,
												double* /*Dh*/)> backupSet_;
		function<void(const interval_t* /*x*/,
		                    interval_t* /*h*/,
		                    interval_t* /*Dh*/)> backupSet_int_;
		function<void(const double* /*x*/,
												double* /*f*/,
												double* /*g*/)> dynamics_;
		function<void(const interval_t* /*x*/,
		                    interval_t* /*f*/,
		                    interval_t* /*g*/)> dynamics_int_;
		function<void(const double* /*x*/,
												double* /*Df*/,
												double* /*Dg*/)> dynamicsGradients_;
		function<void(const interval_t* /*x*/,
		                    interval_t* /*Df*/,
		                    interval_t* /*Dg*/)> dynamicsGradients_int_;
		function<void(const double* /*x*/,
									const double* /*u*/,
												double* /*f*/,
												double* /*g*/,
												double* /*d_fcl_dx*/)> dynamicsWithGradient_;
		function<void(const interval_t* /*x*/,
		              const double* /*u*/,
		                    interval_t* /*f*/,
		                    interval_t* /*g*/,
		                    interval_t* /*d_fcl_dx*/)> dynamicsWithGradient_int_;
		function<void(const double* /*x*/,
		                    double* /*u*/,
		                    double* /*Du*/)> backupController_;

		Options options_;
		double* x0_;
		double* x_unc_;
		double* u_zoh_;
		double* Du_zoh_;
		double t_last_zoh_ = -1;
		QPWrapperAbstract *QPsolver_;

		uint32_t npBT_;

		double *H_;
		double *c_;
		double *A_;
		double *b_;
		double *lb_;
		double *ub_;

	public:
	static string filterErrorMsgString(const int32_t rc);

}; //end ASIFimplicitRB class
} //end ASIF namespace
#endif
