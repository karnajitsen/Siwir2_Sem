#pragma once
#include<iostream>
#include <assert.h>
#include <cmath>
#include <stdlib.h>
#include <memory>
#include <string>
#include<malloc.h>
#include <omp.h>
#include<string.h>
#define LD 64
#define ALLIGNMENT 64
typedef double Real;
//#define offset 8
using namespace std;
//#define M_PI 3.14
class Grid
{

    //__declspec(align(128))
     Real * __restrict data = NULL;
     size_t sizeX, sizeY, ld;
     Real hx, hy;

public:
    explicit Grid()
    {
        data = (Real*)memalign(ALLIGNMENT, 0);
        sizeX = 0;
        sizeY = 0;
        hx = 0.0;
        hy = 0.0;
//alpha = 0.0;
//center = 0.0;
    }

    explicit Grid(const size_t x, const size_t y, const Real& _hx, const Real& _hy , bool bndrYN)
    {
        sizeX = x;
        sizeY = y;
        hx = _hx;
        hy = _hy;
	
	//alpha = 1/hx/hx;
	//center = 4.0 * alpha; 
        ld = x + LD ;
	// size = ld*y*sizeof(Real);
	
         data = (Real*) memalign(ALLIGNMENT, ld*y*sizeof(Real));
         if (bndrYN)
        {
            //Real l = - 1.0 + (sizeX - 1.0)*hx;
	#pragma omp parallel
			{
	#pragma omp for
				for (int j = 0.0; (size_t)j < sizeY; j++)
				{
					Real k = -1.0 + j*hx;
					data[j] = gxy(k, -1.0);
					data[j*ld] = gxy(-1.0, k);
					data[(j + sizeY - 1)] = gxy(k + 1.0, -1.0);
					data[j * ld + (sizeX - 1)] = gxy(1.0, k);				

				}
			}
        }

        
    }
    ~Grid()
    {
        free(data);
    }

    inline Real gxy(const Real x, const Real y)
    {
        if (y == 0.0 && x == 0.0)
	return 0.0;
	//Real r = sqrt(sqrt(x*x + y*y));

		//Real theta = atan2(y, x);
        return (sqrt(sqrt(x*x + y*y)))*sin(M_PI - abs(atan2(y, x)) * 0.5);
    }

    inline void reset()
    {
	memset(data,0.0,ld*sizeY*sizeof(Real));
    }

    inline Real& operator()(const size_t x, const size_t y)
    {
        //assert(x < sizeX);
        //assert(y < sizeY);
        return data[y*ld + x];
    }

    inline Real& operator()(const size_t x, const size_t y) const
    {
        //assert(x < sizeX);
        //assert(y < sizeY);
        return data[y*ld + x];
    }


    inline size_t getXsize() const
    {
         return sizeX;
    }

    inline size_t getYsize() const
    {
        return sizeY;
    }

    inline Real getHx() const
    {
        return hx;
    }

    inline Real getHy() const
    {
        return hy;
    }

};
