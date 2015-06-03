#pragma once
#include<iostream>
#include <assert.h>
#include <cmath>
#include <stdlib.h>
#include <memory>
#include <string>
#include<malloc.h>
#include <omp.h>
#define LD 16
#define ALLIGNMENT 32
using namespace std;
//#define M_PI 3.14
class Grid
{

    //__declspec(align(128))
	double * __restrict data = NULL;
    size_t sizeX, sizeY, ld;
    double hx, hy;

public:
    explicit Grid()
    {
        data = (double*)memalign(ALLIGNMENT, 0);
        sizeX = 0;
        sizeY = 0;
        hx = 0.0;
        hy = 0.0;
    }

    explicit Grid(const size_t x, const size_t y, const double& _hx, const double& _hy , bool bndrYN)
    {
        sizeX = x;
        sizeY = y;
        hx = _hx;
        hy = _hy;
        ld = x + LD;
        //totLength = (x - 2)*(y - 2);
        data = (double*) memalign(ALLIGNMENT, ld*y*sizeof(double));
        //data = (double*) _aligned_malloc(ld*y*sizeof(double), ALLIGNMENT);
        if (bndrYN)
        {
            //double l = - 1.0 + (sizeX - 1.0)*hx;
	#pragma omp parallel
			{
	#pragma omp for
				for (int j = 0.0; (size_t)j < sizeY; j++)
				{
					double k = -1.0 + j*hx;
					data[j] = gxy(k, -1.0);
					data[j*ld] = gxy(-1.0, k);
					data[j + sizeY - 1] = gxy(k + 1.0, -1.0);
					data[j * ld + (sizeX - 1)] = gxy(1.0, k);				

				}
			}
        }

          //data++;
    }
    ~Grid()
    {
        //--data;
        free(data);
    }

    inline double gxy(const double x, const double y)
    {
		if (y == 0.0)
			return 0.0;
		double r = sqrt(sqrt(x*x + y*y));

		double theta = atan2(y, x);
		//std::cout << "Arctan = " << 360 - abs(180 * theta / 3.14159) << '\n';
		return r*sin(M_PI - abs(theta) * 0.5);
    }

    inline void reset()
    {
		size_t x = sizeX;
		size_t y = sizeY ;
#pragma omp parallel
		{
	#pragma omp for 
			for (size_t i = 0; i < y; i++)
			{
				//#pragma omp parallel for
				for (size_t j = 0; j < x; j+=8)
				{
					data[i*ld + j] = 0.0;
					data[i*ld + j+1] = 0.0;
					data[i*ld + j+2] = 0.0;
					data[i*ld + j+3] = 0.0;
					data[i*ld + j + 4] = 0.0;
					data[i*ld + j + 5] = 0.0;
					data[i*ld + j + 6] = 0.0;
					data[i*ld + j + 7] = 0.0;
				}
			}
		}
    }

    inline double& operator()(const size_t x, const size_t y)
    {
        //assert(x < sizeX);
        //assert(y < sizeY);
        return data[y*ld + x];
    }

    inline double& operator()(const size_t x, const size_t y) const
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

    inline double getHx() const
    {
        return hx;
    }

    inline double getHy() const
    {
        return hy;
    }

};
