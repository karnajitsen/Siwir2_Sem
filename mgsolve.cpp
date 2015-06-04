#include<iostream>
#include<fstream>
#include <sys/time.h>
#include <omp.h>
#include <immintrin.h>
#include "Grid.h"
#define XDOMLOW -1.0
#define XDOMHIGH 1.0
#define YDOMLOW -1.0
#define YDOMHIGH 1.0
#define TOLERR 0.0000918
#define V1 2
#define V2 1

Grid ** __restrict xGrids = nullptr;
Grid ** __restrict fGrids = nullptr;
Grid * __restrict sGrid = nullptr;
Grid * __restrict iGrid = nullptr;

void init(double hsize, const size_t level)
{

	size_t je = level;
    size_t xdim = pow(2, je) + 1;
    size_t ydim = (xdim/2)+1;
    bool flag = true;
    xGrids = (Grid**) memalign(ALLIGNMENT, level*sizeof(Grid*));
    fGrids = (Grid**) memalign(ALLIGNMENT, level*sizeof(Grid*));
    iGrid = new Grid(xdim, ydim, hsize, hsize, flag);
    for (size_t i = 0; i < level; i++)
    {
        xGrids[i] = new Grid(xdim, ydim, hsize, hsize, flag);
		fGrids[i] = new Grid(xdim, ydim, hsize, hsize, false);
        xdim = pow(2, --je) + 1;
        ydim = (xdim * 0.5) + 1;
        hsize *= 2.0;
        flag = false;
    }


}

inline void smooth(Grid* __restrict xgrd, const  Grid* __restrict fgrd, const size_t iter)
{
	size_t dimX = (*xgrd).getXsize()-1;
	size_t dimY = (*xgrd).getYsize()-1;
	double hx = (*xgrd).getHx();
    //double hy = (*xgrd).getHy();
    size_t midX = dimX >> 1;
	size_t j,k;	
	for (size_t i = 0; i < iter; i++)
	{
		
#pragma omp parallel private(j,k) shared(dimY,dimX,midX,hx) 
		{
			
#pragma omp for nowait
			for ( j = 1; j < dimY; j++)
			{
                size_t l = ((j + 1) & 0x1) + 1;
                    for ( k = l; k < dimX; k += 2)

					{
			     (*xgrd)(k, j) = (hx*hx*(*fgrd)(k, j) + (*xgrd)(k + 1, j) + (*xgrd)(k - 1, j) + (*xgrd)(k, j + 1) + (*xgrd)(k, j - 1)) * 0.25;
	 	}				
					
			}
#pragma omp for nowait
			for ( k = 1; k < midX; k += 2)
			{
               (*xgrd)(k, dimY) = (hx*hx*(*fgrd)(k, dimY) + (*xgrd)(k + 1, dimY) + (*xgrd)(k - 1, dimY) + 2.0 * (*xgrd)(k, dimY - 1)) * 0.25;

			}

		}

		
#pragma omp parallel private(j,k) shared(dimY,dimX,midX,hx)
		{
#pragma omp for nowait
			for ( j = 1; j < dimY; j++)
			{
				size_t l = (j & 0x1) + 1;
				for ( k = l; k < dimX; k += 2)
				{

                 (*xgrd)(k, j) = (hx*hx*(*fgrd)(k, j) + (*xgrd)(k + 1, j) + (*xgrd)(k - 1, j) + (*xgrd)(k, j + 1) + (*xgrd)(k, j - 1)) * 0.25;

				}


			}
#pragma omp for nowait
			for ( k = 2; k < midX; k += 2)
			{
                (*xgrd)(k, dimY) = (hx*hx*(*fgrd)(k, dimY) + (*xgrd)(k + 1, dimY) + (*xgrd)(k - 1, dimY) + 2.0 * (*xgrd)(k, dimY - 1)) * 0.25;
			}
		
		}
	}
	
}

inline void restriction(const  Grid * __restrict xgrd, const Grid *  __restrict fgrd, Grid* __restrict rgrid)
{
    size_t xlen = (*xgrd).getXsize()-1;
    size_t ylen = (*xgrd).getYsize()-1;
    double hx = (*xgrd).getHx();
    double	alpha = 1.0 / hx / hx;
     double	center = (4.0 * alpha);
    size_t midX = xlen >> 1;

    Grid tmpgrd(xlen+1, ylen+1, hx, hx, false);

#pragma omp parallel
	{

#pragma omp for nowait
		for (size_t i = 1; i < ylen; i++)
		{
			
				for (size_t j = 1; j < xlen; j+=4)
                {

	                
			tmpgrd(j, i) = (*fgrd)(j, i) + alpha*((*xgrd)(j + 1, i) + (*xgrd)(j - 1, i) + (*xgrd)(j, i + 1) + (*xgrd)(j, i - 1)) - (*xgrd)(j, i) * center;
			 tmpgrd(j+1, i) = (*fgrd)(j+1, i) + alpha*((*xgrd)(j + 2, i) + (*xgrd)(j , i) + (*xgrd)(j+1, i + 1) + (*xgrd)(j+1, i - 1)) - (*xgrd)(j+1, i) * center;		
 tmpgrd(j+2, i) = (*fgrd)(j+2, i) + alpha*((*xgrd)(j + 3, i) + (*xgrd)(j+1 , i) + (*xgrd)(j+2, i + 1) + (*xgrd)(j+2, i - 1)) - (*xgrd)(j+2, i) * center;
 tmpgrd(j+3, i) = (*fgrd)(j+3, i) + alpha*((*xgrd)(j + 4, i) + (*xgrd)(j+2 , i) + (*xgrd)(j+3, i + 1) + (*xgrd)(j+3, i - 1)) - (*xgrd)(j+3, i) * center;

}	
				}
		

			
#pragma omp for	nowait	
        for (size_t j = 1; j < midX; j+=2)
				{
                   tmpgrd(j, ylen) = (*fgrd)(j, ylen) + alpha*((*xgrd)(j + 1, ylen) + (*xgrd)(j - 1, ylen) + 2.0 * (*xgrd)(j, ylen - 1)) - (*xgrd)(j, ylen) * center;
		   tmpgrd(j+1, ylen) = (*fgrd)(j+1, ylen) + alpha*((*xgrd)(j + 2, ylen) + (*xgrd)(j , ylen) + 2.0 * (*xgrd)(j+1, ylen - 1)) - (*xgrd)(j+1, ylen) * center;
				}
		
				
}
    size_t rxlen = (*rgrid).getXsize() - 1;
    size_t rylen = (*rgrid).getYsize() -1;
    midX = rxlen >> 1;

#pragma omp parallel
	{

#pragma omp for nowait
		for (size_t i = 1; i < rylen+1; i++)
		{
			//size_t k = i << 1;
				for (size_t j = 1; j < rxlen; j++)
                {
	              //size_t l = j << 1;
			if(i<rylen)
                        (*rgrid)(j, i) = (tmpgrd(2 * j - 1, 2 * i - 1) + tmpgrd(2 * j - 1, 2 * i + 1) +
                        tmpgrd(2 * j + 1, 2 * i - 1) + tmpgrd(2 * j + 1, 2 * i + 1)) * 0.0625 +
                        0.125 *(tmpgrd(2 * j, 2 * i - 1) + tmpgrd(2 * j, 2 * i + 1) +
                        tmpgrd(2 * j - 1, 2 * i) + tmpgrd(2 * j + 1, 2 * i)) + 0.25 * tmpgrd(2 * j, 2 * i);
            		else if(j<midX)
			(*rgrid)(j, rylen) = (tmpgrd(2 * j - 1, 2 * rylen - 1) +
                        tmpgrd(2 * j + 1, 2 * rylen - 1)) * 0.125 +
                        0.125 *(tmpgrd(2 * j, 2 * rylen - 1) * 2.0 + tmpgrd(2 * j - 1, 2 * rylen) + tmpgrd(2 * j + 1, 2 * rylen)) + 0.25 * tmpgrd(2 * j, 2 * rylen);

			}
			}
}
 // delete tmpgrd;
   }

inline void interpolate(Grid * __restrict srcgrd, Grid * __restrict tgtgrd)
{
      size_t txlen = (*tgtgrd).getXsize()-1;
	size_t tylen = (*tgtgrd).getYsize()-1;
    __m256d b;
    b = _mm256_setzero_pd ();

#pragma omp parallel
	{
#pragma omp for nowait
		for (size_t i = 1; i < tylen; i += 2)
		{
			size_t l = i * 0.5;
			for (size_t j = 1; j < txlen; j += 2)
			{
				size_t k = j * 0.5;
				(*tgtgrd)(j, i) += 0.25*((*srcgrd)(k, l) + (*srcgrd)(k + 1, l) + (*srcgrd)(k, l + 1)
					+ (*srcgrd)(k + 1, l + 1));
				(*tgtgrd)(j + 1, i) += 0.5*((*srcgrd)(k + 1, l) + (*srcgrd)(k + 1, l + 1));
				(*tgtgrd)(j, i + 1) += 0.5*((*srcgrd)(k, l + 1) + (*srcgrd)(k + 1, l + 1));
				(*tgtgrd)(j + 1, i + 1) += (*srcgrd)(k + 1, l + 1);

			}
		}

#pragma omp for nowait
        for (size_t j = txlen >> 1; j <= txlen; j+=8)
		{
            _mm256_store_pd(&(*tgtgrd)(j, tylen),b);
            _mm256_store_pd(&(*tgtgrd)(j+4, tylen),b);

		}        


	}

}

inline double errorNorm(const Grid* __restrict  xgrd, const Grid * __restrict sgrd)
{

    size_t dimX = (*xgrd).getXsize();
    size_t dimY = (*xgrd).getYsize();
    double r1 = 0.0,r2=0.0, r3=0.0,r4 = 0.0, sum = 0.0;
   //size_t j,k;
#pragma omp parallel reduction(+: sum) private(r1,r2,r3,r4)
	{
#pragma omp for nowait
		for (size_t  j = 0; j < dimY ; j++)
		{
            for (size_t  k = 0; k < dimX; k+=4)
			{
	r1 = (*sgrd)(k, j) - (*xgrd)(k,j);
	 r2 = (*sgrd)(k+1, j) - (*xgrd)(k+1,j);
 r3 = (*sgrd)(k+2, j) - (*xgrd)(k+2,j);
 r4 = (*sgrd)(k+3, j) - (*xgrd)(k+3,j);

                      if(j<dimY-1)
  sum += 2.0 *(r1*r1 + r2*r2 + r3*r3 + r4*r4);
else
sum += (r1*r1 + r2*r2 + r3*r3 + r4*r4);
			}
		}
	}
	
return sqrt(sum / dimX / ((dimY << 1) - 1.0));
}

void solvemg(size_t level)
{
    size_t xdim = (1 << level) + 1, ydim;
    double newnorm = 1.0;
    double hsize = (XDOMHIGH - XDOMLOW) / (xdim - 1.0);
	size_t i = 0;    
    ydim = (xdim >> 1) + 1;
    sGrid = new Grid(xdim, ydim, hsize, hsize, true);
    //size_t j,k;
#pragma omp parallel //private(j,k) //shared(ydim,xdim,hsize)
	{
#pragma omp for nowait
        for ( size_t k = 0; k < ydim; k++)
		{
            for (size_t j = 0; j < xdim; j++)
			{
				(*sGrid)(j, k) = (*sGrid).gxy(-1.0 + j*hsize, -1.0 + k*hsize);
			}
		}
	}

	for ( i = 1; newnorm > TOLERR; i++)
    {
        for (size_t jl = 0; jl < level - 1; jl++)
        {
            smooth(xGrids[jl], fGrids[jl], V1);
			restriction(xGrids[jl], fGrids[jl], fGrids[jl + 1]);
        }

        for (size_t j = level - 1; j > 0; j--)
        {
            smooth(xGrids[j], fGrids[j], V2);
            interpolate(xGrids[j], xGrids[j - 1]);
            (*xGrids[j]).reset();
            (*fGrids[j]).reset();
        }

		newnorm = errorNorm(xGrids[0], sGrid);
      }
    
    std::cout << "Dirichlet:: Error L2 Norm for h as 1/" << xdim - 1 << " after " << i << " V-Cycle = " << newnorm << "\n\n";
}

int main(int argc, char** argv)
{

    if (argc != 2)
    {
        std::cout << "Invalid number of argument";
        exit(0);
    }	

	
	size_t level = atoi(argv[1]);
    double hsize = (XDOMHIGH - XDOMLOW) / (pow(2, level));
    init(hsize, level);

    string your_alias = "Group_Karnajit_Ramyar";


    std::cout<<"Your Alias: "<<your_alias<<std::endl;
    struct timeval t0, t;
    gettimeofday(&t0, NULL);
     solvemg(level);
    gettimeofday(&t, NULL);
    std::cout << "Wall clock time of MG execution: " <<  ((int64_t)(t.tv_sec - t0.tv_sec) * (int64_t)1000000 +  (int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3  << " ms" << std::endl;


    size_t xdim = (*xGrids[0]).getXsize();
	size_t ydim = (*xGrids[0]).getYsize();

    std::string fname1 = std::string("solution.dat") ;
    std::ofstream	fOut1(fname1);
    std::string fnames1 = std::string("init.dat");
    std::ofstream	fOutsolt1(fnames1);
	std::cout << "\n\nWriting solution to the file...\n\n";

//#pragma omp parallel for
		for (size_t y = 0; y < ydim; ++y) {
			for (size_t x = 0; x < xdim; ++x) {

				fOut1 << x*hsize - 1.0 << "\t" << y*hsize - 1.0 << "\t" << (*xGrids[0])(x, y) << std::endl;
                fOutsolt1 << x*hsize - 1.0 << "\t" << y*hsize - 1.0 << "\t" << (*iGrid)(x, y) << std::endl;
			}
            //fOut1 << std::endl;
              //       fOutsolt1 << std::endl;

		}

         xdim = (*xGrids[0]).getXsize();
		 ydim = (*xGrids[0]).getYsize();

		for (int y = ydim - 2; y >= 0; y = y-1) 
		{
			for (size_t x = 0; x < xdim; ++x) 
			{
				fOut1 << x*hsize - 1.0 << "\t" << (ydim-1-y)*hsize << "\t" << (*xGrids[0])(x, y) << std::endl;
                fOutsolt1 << x*hsize - 1.0 << "\t" << (ydim - 1 - y)*hsize << "\t" << (*iGrid)(x, y) << std::endl;
			}
            //fOut1 << std::endl;
              //       fOutsolt1 << std::endl;
        }
	
		fOut1.close();
		fOutsolt1.close();
	
    std::cout << "\n\n =============== Dirichlet Boundary Value Problem 1 ends here ===================\n\n";

    delete xGrids;
    delete fGrids;
    delete sGrid;
    delete iGrid;

    return 0;
}
