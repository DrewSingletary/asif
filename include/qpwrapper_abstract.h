#ifndef _QPWRAPPER_ABSTRACT_H_
#define _QPWRAPPER_ABSTRACT_H_

#include <stdint.h>
#ifdef ASIF_DEBUG
	#include <iostream>
#endif

namespace ASIF
{
	////////////////////////////////////////////
	//Solves  min  (xT.H.x+cT.x)              //
	//        s.t. A.x>=b                     //
	//             lb<=x<=ub                  //
	////////////////////////////////////////////
	class QPWrapperAbstract
	{
	public:
		enum class SOLVER_STATUS : int32_t
		{
			INFEASIBLE = 0,
			FEASIBLE = 1
		};

	public:
		QPWrapperAbstract(const uint32_t nv,
		                  const uint32_t nc,
		                  const bool diagonalCost);
		virtual ~QPWrapperAbstract(void);

		virtual int32_t initialize(const double H[],
		                           const double c[],
		                           const double A[],
		                           const double b[],
		                           const double lb[],
		                           const double ub[],
		                           const bool   be[] = nullptr) = 0;
		virtual int32_t updateCost(const double H[], const double c[]) = 0;
		virtual int32_t updateA(const double A[]) = 0;
		virtual int32_t updateb(const double b[]) = 0;
		virtual int32_t updateBounds(const double lb[], const double ub[]) = 0;
		virtual int32_t solve(void) = 0;
		virtual int32_t getSolution(double sol[]) = 0;

	protected:

		const uint32_t nv_;
		const uint32_t nc_;
		const bool diagonalCost_;
		bool *be_;
	};
} //end ASIF namespace
#endif
