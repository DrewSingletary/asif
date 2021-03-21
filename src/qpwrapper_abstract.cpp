#include "qpwrapper_abstract.h"

namespace ASIF
{
	QPWrapperAbstract::QPWrapperAbstract(const uint32_t nv,
	                                     const uint32_t nc,
	                                     const bool diagonalCost):
	nv_(nv),
	nc_(nc),
	diagonalCost_(diagonalCost)
	{
		be_ = new bool[nc];
		for(uint32_t i=0; i<nc_; i++)
		{
			be_[i] = false;
		}
	}

	QPWrapperAbstract::~QPWrapperAbstract(void)
	{
		delete[] be_;
	}

} //end ASIF namespace