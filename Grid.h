#pragma once
#include<iostream>
#include <assert.h>
#include <cmath>
#include <stdlib.h>
#include <memory>
#include <string>
#include<malloc.h>
#define LD 16
#define ALLIGNMENT 32
using namespace std;
//#define M_PI 3.14
class Grid
{

    //__declspec(align(128))
    double * __restrict data = NULL;
    size_t sizeX, sizeY, ld, totLength;
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

    explicit Grid(const size_t x, const size_t y, const double& _hx, const double& _hy , bool bndrYN , bool isNeumann)
    {
        sizeX = x;
        sizeY = y;
        hx = _hx;
        hy = _hy;
        ld = x + LD;
        totLength = (x - 2)*(y - 2);
        data = (double*) memalign(ALLIGNMENT, ld*y*sizeof(double));
        //data = (double*) _aligned_malloc(ld*y*sizeof(double), ALLIGNMENT);
        if (bndrYN && !isNeumann)
        {
            double l = (sizeX - 1.0)*hx;
            for (int j = 0.0; (size_t)j < sizeX; j++)
            {
                double k = j*hx;
                data[j] = gxy1(k, 0.0);
                data[j*ld] = gxy1(0.0, k);
                data[j + ld * (sizeX - 1)] = gxy1(k, l);
                data[j * ld + (sizeX - 1)] = gxy1(l, k);
            }
        }

        if (bndrYN && isNeumann)
        {
            for (int j = 0.0; (size_t)j < sizeX; j++)
            {
                double k = j*hx;
                data[j] = gxy2(k);
                data[j + ld * (sizeX - 1)] = gxy2(k);
            }
        }
        //data++;
    }
    ~Grid()
    {
        //--data;
        free(data);
    }

    inline double gxy1(const double x, const double y)
    {
        return sin(M_PI * x) * sinh(M_PI * y);
    }

    inline double gxy2(const double x)
    {
        return (x * (1.0 - x));
    }

    inline void reset()
    {
        for (size_t i = 1; i < sizeX - 1; i++)
        {
            for (size_t j = 1; j < sizeX - 1; j++)
            {
                data[i*ld + j] = 0.0;
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

    inline Grid * operator+=(const Grid * rhs)
    {
        return this;
    }

    inline Grid* operator-=(const Grid* rhs)
    {
        return this;
    }

    inline size_t getXsize() const
    {
        return sizeX;
    }

    inline size_t getYsize() const
    {
        return sizeY;
    }

    inline size_t getTotLength() const
    {
        return totLength;
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
