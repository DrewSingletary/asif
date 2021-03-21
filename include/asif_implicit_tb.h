#ifndef _ASIF_IMPLICIT_TB_H_
#define _ASIF_IMPLICIT_TB_H_

#include "asif_utils.h"

#ifdef USE_ODEINT
	#include <boost/numeric/odeint.hpp>
	#include <boost/numeric/odeint/iterator/n_step_iterator.hpp>
	using namespace boost::numeric::odeint;
#endif

using namespace std;

namespace ASIF
{
	class ASIFimplicitTB
	{
	public:
		typedef struct
		{
			double relaxCost = 50.0;
			double relaxSafeLb = 5.0;
			double relaxTTS = 5.0;
			double relaxMinOrtho = 5.0;
			double backTrajHorizon = 1.0;
			double backTrajExtend = 0.05;
			double backTrajDt = 0.01;
			double backTrajMinOrtho = 0.01;
			double backTrajAbsTol = 1.0e-6;
			double backTrajRelTol = 1.0e-6;
			double satSharpness = 0.1;
			double inf = 1e20;
		} Options;

		typedef vector<double> state_t;

		#ifdef USE_ODEINT
			typedef runge_kutta_dopri5<state_t> stepper_t;
		#endif

	public:
		ASIFimplicitTB(const uint32_t nx,
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
		               const QPSOLVER qpSolverType = QPSOLVER::OSQP,
		               const bool diagonalCost = true);
		ASIFimplicitTB(const uint32_t nx,
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
		               const QPSOLVER qpSolverType = QPSOLVER::OSQP,
		               const bool diagonalCost = true);
		~ASIFimplicitTB(void);

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
		
		int32_t updateOptions(const Options &options);

		vector<pair<double,state_t>> backTraj_;
		vector<uint32_t> backTrajCritIdx_;
		double TTS_;
		double BTorthoBS_;
		double hBackupEnd_;
		double hSafetyNow_;
	protected:
		ASIFimplicitTB(const bool dynamicsHasGradient,
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
		               const QPSOLVER qpSolverType = QPSOLVER::OSQP,
		               const bool diagonalCost = true);

		int32_t updateOptions(void);
		int32_t updateConstraints(const double x[]);
		int32_t updateConstraintsTrivial(void);
		int32_t updateCost(const double uDes[]);
		int32_t updateH(const double H[]);
		void inputSaturateSoft(const double u[],
		                             double uSat[],
		                             double DuSat[]);
		void inputSaturate(double u[]);
		void backupCLdynamics(const double x[],
		                            double f[],
		                            double Df[]);
		void ODErhs(const state_t &x,
		                  state_t &xDot,
		            const double t);

		const bool dynamicsHasGradient_;
		const uint32_t nx_;
		const uint32_t nu_;
		const uint32_t nv_;
		const uint32_t npSS_;
		const uint32_t npBTSS_;
		const uint32_t npTC_;
		function<void(const double* /*x*/,
		                    double* /*h*/,
		                    double* /*Dh*/)> safetySet_;
		function<void(const double* /*x*/,
		                    double* /*h*/,
		                    double* /*Dh*/,
								  double* /*DDh*/)> backupSet_;
		function<void(const double* /*x*/,
		                    double* /*f*/,
		                    double* /*g*/)> dynamics_;
		function<void(const double* /*x*/,
		                    double* /*Df*/,
		                    double* /*Dg*/)> dynamicsGradients_;
		function<void(const double* /*x*/,
		              const double* /*u*/,
		                    double* /*f*/,
		                    double* /*g*/,
		                    double* /*d_fcl_dx*/)> dynamicsWithGradient_;		
		function<void(const double* /*x*/,
		                    double* /*u*/,
		                    double* /*Du*/)> backupController_;

		Options options_;
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

	}; //end ASIFimplicit class
} //end ASIF namespace
#endif
