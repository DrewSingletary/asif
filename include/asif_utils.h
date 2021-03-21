#ifndef _ASIF_UTILS_H_
#define _ASIF_UTILS_H_

#include <assert.h>
#include <stdint.h>
#include <cstring>
#include <vector>
#include <functional>
#include <cmath>
#include "qpwrappers.h"
#include <algorithm>

#ifdef ASIF_DEBUG
#include <iostream>
#include <iomanip>
#include <stdio.h>
#endif

namespace ASIF
{

	template<class T> inline void matrixMultiply(const T A[],
												 const uint32_t nlA,
												 const uint32_t ncA,
												 const T B[],
												 const uint32_t nlB,
												 const uint32_t ncB,
													   T AB[])
	{
		assert(ncA==nlB);
		uint32_t idxAB = 0;
		for(uint32_t i = 0; i<nlA; i++)
		{
			for(uint32_t j = 0; j<ncB; j++)
			{
				idxAB = i + j*nlA;
				AB[idxAB] = 0.0;
				for(uint32_t k = 0; k<ncA; k++)
				{
					AB[idxAB] = AB[idxAB] + A[i+k*nlA]*B[k+j*ncA];
				}				
			}
		}
	}

	template<class T> inline void matrixVectorMultiply(const T A[],
													   const uint32_t nlA,
													   const uint32_t ncA,
													   const T b[],
													   const uint32_t nlb,
															 T Ab[])
	{
		assert(ncA==nlb);
		for(uint32_t i = 0; i<nlA; i++)
		{
			Ab[i] = 0.0;
			for(uint32_t k = 0; k<ncA; k++)
			{
				Ab[i] = Ab[i] + A[i+k*nlA]*b[k];
			}				
		}
	}

	template<class T> inline T vectorNorm(const T v[],
													  const uint32_t len)
	{
		double tmp = 0;
		for(uint32_t i = 0; i<len; i++)
		{
			tmp += v[i]*v[i];
		}
		return sqrt(tmp);
	}

	#ifdef ASIF_DEBUG
	template<class T> inline void printMatrix(const T M[],
											  const uint32_t nlM,
											  const uint32_t ncM,
											  const uint32_t decimals = 2)
	{
		std::cout << std::scientific;
		std::cout << std::setprecision(decimals);
		for(uint32_t i=0; i<((decimals+7)*ncM + ncM + 1); i++)
		{
			std::cout << "_" ;
		}
		std::cout << std::endl;
	
		for(uint32_t i = 0; i<nlM; i++)
		{
			std::cout << '|' ;
			for(uint32_t j = 0; j<ncM; j++)
			{
				std::cout << std::setw(decimals+7) << M[i+j*nlM] << "|" ;
			}
			std::cout << std::endl ;
		}

		for(uint32_t i=0; i<((decimals+7)*ncM + ncM + 1); i++)
		{
			std::cout << "-" ;
		}
		std::cout << std::endl;
	}

	template<class T> inline void printVector(const T v[],
											  const uint32_t nlv,
											  const uint32_t decimals = 2)
	{
		std::cout << std::scientific;
		std::cout << std::setprecision(decimals);;
		for(uint32_t i=0; i<((decimals+7) + 2); i++)
		{
			std::cout << "_" ;
		}
		std::cout << std::endl;

		for(uint32_t i = 0; i<nlv; i++)
		{
			std::cout << '|' ;
			std::cout << std::setw(decimals+7) << v[i];
			std::cout << '|' << std::endl ;
		}

		for(uint32_t i=0; i<((decimals+7) + 2); i++)
		{
			std::cout << "-" ;
		}
		std::cout << std::endl;
	}

	inline void tic(ASIFtimer *t)
	{
		clock_gettime(CLOCK_MONOTONIC, &t->tic);
	}

	inline void toc(ASIFtimer *t)
	{
		struct timespec temp;

		clock_gettime(CLOCK_MONOTONIC, &t->toc);

		if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0) {
			temp.tv_sec  = t->toc.tv_sec - t->tic.tv_sec - 1;
			temp.tv_nsec = 1e9 + t->toc.tv_nsec - t->tic.tv_nsec;
		} else {
			temp.tv_sec  = t->toc.tv_sec - t->tic.tv_sec;
			temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
		}
		t->dt = ((double)temp.tv_sec + (double)temp.tv_nsec / 1e9);
	}
	#endif
} //end ASIF namespace

#endif
