#include "Grid.h"
#include<iostream>
#include<fstream>
#include <sys/time.h>
#include <omp.h>
#define XDOMLOW -1.0
#define XDOMHIGH 1.0
#define YDOMLOW -1.0
#define YDOMHIGH 1.0
//#define TOLERR 0.0000918
#define TOLERR 0.002
#define V1 2
#define V2 1

Grid ** __restrict xGrids = nullptr;
Grid ** __restrict fGrids = nullptr;
Grid * __restrict sGrid = nullptr;
//bool isNeumann = false;

inline void init(double hsize, const size_t level)
{
	//std::cout << "Init";
	size_t je = level;
    size_t xdim = pow(2, je) + 1;
    size_t ydim = (xdim/2)+1;
    bool flag = true;
    xGrids = (Grid**) memalign(ALLIGNMENT, level*sizeof(Grid*));
    fGrids = (Grid**) memalign(ALLIGNMENT, level*sizeof(Grid*));
//#pragma omp parallel for
    for (size_t i = 0; i < level; i++)
    {
        xGrids[i] = new Grid(xdim, ydim, hsize, hsize, flag);
		fGrids[i] = new Grid(xdim, ydim, hsize, hsize, false);
        xdim = pow(2, --je) + 1;
		ydim = (xdim / 2) + 1;
        hsize *= 2.0;
        flag = false;
    }
   
}

inline void smooth(Grid* __restrict xgrd, const  Grid* __restrict fgrd, const size_t iter)
{
	size_t dimX = (*xgrd).getXsize()-1;
	size_t dimY = (*xgrd).getYsize()-1;
	double hx = (*xgrd).getHx();
	double hy = (*xgrd).getHy();
	//size_t midY = (dimY - 1) / 2;
	size_t midX = (dimX - 1) / 2;
	//std::cout << "Smooth";
	//size_t j = 0;
	for (size_t i = 0; i < iter; i++)
	{
		
#pragma omp parallel //firstprivate(dimY,dimX,midX,midY,hx,hy)
		{
			
#pragma omp for
			for (size_t j = 1; j < dimY; j++)
			{
				size_t l = ((j+1) & 0x1) + 1;
				for (size_t k = l; k < dimX ; k += 2)
				{

					(*xgrd)(k, j) = (hx*hy*(*fgrd)(k, j) + (*xgrd)(k + 1, j) + (*xgrd)(k - 1, j) + (*xgrd)(k, j + 1)
						+ (*xgrd)(k, j - 1)) * 0.25;
				}

			}
#pragma omp for
			for (size_t k = 1; k < midX ; k += 2)
			{
				(*xgrd)(k, dimY) = (hx*hy*(*fgrd)(k, dimY) + (*xgrd)(k + 1, dimY) + (*xgrd)(k - 1, dimY) + 2.0 * (*xgrd)(k, dimY - 1)) * 0.25;
			}
		}
#pragma omp parallel //firstprivate(dimY,dimX,midX,midY,hx,hy)
		{
#pragma omp for
			for (size_t j = 1; j < dimY; j++)
			{
				size_t l = (j & 0x1) + 1;
					for (size_t k = l; k < dimX; k += 2)
					{

						(*xgrd)(k, j) = (hx*hy*(*fgrd)(k, j) + (*xgrd)(k + 1, j) + (*xgrd)(k - 1, j) + (*xgrd)(k, j + 1)
							+ (*xgrd)(k, j - 1)) * 0.25;
					}			
				
			}
#pragma omp for
			for (size_t k = 2; k < midX; k += 2)
			{
				(*xgrd)(k, dimY) = (hx*hy*(*fgrd)(k, dimY) + (*xgrd)(k + 1, dimY) + (*xgrd)(k - 1, dimY) + 2.0 * (*xgrd)(k, dimY - 1)) * 0.25;
			}

		}
	}
	
}

inline void restriction(const  Grid * __restrict xgrd, const Grid *  __restrict fgrd, Grid* __restrict rgrid)
{
    size_t xlen = (*xgrd).getXsize()-1;
    size_t ylen = (*xgrd).getYsize()-1;
    double hx = (*xgrd).getHx();
    double hy = (*xgrd).getHy();
    double	alpha = 1.0 / hx / hx;
    double	beta = 1.0 / hy / hy;
    double	center = (2.0 * alpha) + (2.0 * beta);
	
	//std::cout << "***************:: restriction= ";
	Grid tmpgrd(xlen+1, ylen+1, hx, hy, false);
	//size_t i = 1;
#pragma omp parallel //firstprivate(xlen,ylen,midX,midY,alpha,beta,center)
	{
		//size_t i = 0;
#pragma omp for
		for (size_t i = 1; i < ylen; i++)
		{
			for (size_t j = 1; j < xlen; j++)
			{
					tmpgrd(j, i) = (*fgrd)(j, i) + alpha*((*xgrd)(j + 1, i) + (*xgrd)(j - 1, i)) + beta * ((*xgrd)(j, i + 1)
					+ (*xgrd)(j, i - 1)) - (*xgrd)(j, i) * center;
			}


		}

#pragma omp for
		for (size_t j = 1; j < xlen / 2; j++)
		{
			tmpgrd(j, ylen) = (*fgrd)(j, ylen) + alpha*((*xgrd)(j + 1, ylen) + (*xgrd)(j - 1, ylen)) + beta * (2.0 * (*xgrd)(j, ylen - 1)) - (*xgrd)(j, ylen) * center;
		}
	}

    size_t rxlen = (*rgrid).getXsize() - 1;
    size_t rylen = (*rgrid).getYsize() - 1;

#pragma omp parallel //firstprivate(rxlen,rylen,midX,midY)
	{
		//size_t i = 0;
#pragma omp for
		for (size_t i = 1; i < rylen; i++)
		{
			for (size_t j = 1; j < rxlen; j++)
			{
				//if ((i == midY && j < midX) || i != midY)
					(*rgrid)(j, i) = (tmpgrd(2 * j - 1, 2 * i - 1) + tmpgrd(2 * j - 1, 2 * i + 1) +
					tmpgrd(2 * j + 1, 2 * i - 1) + tmpgrd(2 * j + 1, 2 * i + 1)) * 0.0625 +
					0.125 *(tmpgrd(2 * j, 2 * i - 1) + tmpgrd(2 * j, 2 * i + 1) +
					tmpgrd(2 * j - 1, 2 * i) + tmpgrd(2 * j + 1, 2 * i)) + 0.25 * tmpgrd(2 * j, 2 * i);
			}

		}
#pragma omp for
		for (size_t j = 1; j < rxlen - 1; j++)
		{
			(*rgrid)(j, rylen) = (tmpgrd(2 * j - 1, 2 * rylen - 1) +
				tmpgrd(2 * j + 1, 2 * rylen - 1)) * 0.125 +
				0.125 *(tmpgrd(2 * j, 2 * rylen - 1) * 2.0 + tmpgrd(2 * j - 1, 2 * rylen) + tmpgrd(2 * j + 1, 2 * rylen)) + 0.25 * tmpgrd(2 * j, 2 * rylen);
		}

	}
   }

inline void interpolate(Grid * __restrict srcgrd, Grid * __restrict tgtgrd)
{
      size_t txlen = (*tgtgrd).getXsize()-1;
	size_t tylen = (*tgtgrd).getYsize()-1;
	//std::cout << "Interpolate";
#pragma omp parallel  //firstprivate(txlen,tylen) 
	{
#pragma omp for
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
	}

}

inline double errorNorm(const Grid* __restrict  xgrd, const Grid * __restrict sgrd)
{

    size_t dimX = (*xgrd).getXsize();
    size_t dimY = (*xgrd).getYsize();
	double r = 0.0, sum = 0.0;
	//size_t j;
#pragma omp parallel reduction(+: sum) firstprivate(dimY, dimX,r)
	{
#pragma omp for
		for (size_t j = 0; j < dimY; j++)
		{
			for (size_t k = 0; k < dimX; k++)
			{
				r = (*sgrd)(k, j) - (*xgrd)(k, j);
				 
				if (j!= dimY-1)
					sum += 2.0*r*r;
				else
					sum += r*r;
			}

		}
	}

    return sqrt(sum / dimX / (2*dimY- 1));
}

void mgsolve(size_t level)
{
    size_t xdim = pow(2, level) + 1, ydim;
    double newnorm = 1.0;
    double hsize = (XDOMHIGH - XDOMLOW) / (xdim - 1.0);
	size_t i = 0;
    init(hsize, level);
	ydim = (xdim / 2) + 1;
	sGrid = new Grid(xdim, ydim, hsize, hsize, true);
	//size_t k = 0;
	//std::cout << "solution grid";
#pragma omp parallel //firstprivate(gdim,hsize)
	{
#pragma omp for
		for (size_t k = 0; k < ydim; k++)
		{
			//#pragma omp parallel for
			for (size_t j = 0; j < xdim; j++)
			{
				(*sGrid)(j, k) = (*sGrid).gxy(-1.0 + j*hsize, -1.0 + k*hsize);
			}
		}
	}
	
	//std::cout << "solution grid2";
	for ( i = 1; newnorm > TOLERR*4; i++)
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
		std::cout << "Dirichlet:: Error L2 Norm for h as 1/ = " << newnorm << "\n\n";	

    }
   // vcycle = i;   
    
    std::cout << "Dirichlet:: Error L2 Norm for h as 1/" << xdim - 1 << " after " << i << " V-Cycle = " << newnorm << "\n\n";
}

int main(int argc, char** argv)
{

    //std::cout << "1";
    if (argc != 2)
    {
        std::cout << "Invalid number of argument";
        exit(0);
    }

	//omp_set_num_threads(16);

	
	size_t level = atoi(argv[1]);
    //size_t vcycle = 0;

    timeval start, end;

    //std::cout << "Dirichlet:: Level = " << level << "\n\n";

    std::cout << "\n\n =============== Output for Dirichlet Boundary Value Problem 1 ===================\n\n";
	std::cout << "333";
    gettimeofday(&start, 0);
	std::cout << "22222";
    mgsolve(level);

    gettimeofday(&end, 0);
    double elapsed = 0.000001 * ((double)((end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec));
    std::cout << "Dirichlet:: Time spend for Multigrid Solver = " << elapsed << '\n';

    double hsize = (*xGrids[0]).getHx();
    size_t xdim = (*xGrids[0]).getXsize();
	size_t ydim = (*xGrids[0]).getYsize();

    std::string fname1 = std::string("data/Dirichlet/solution_h_") + std::string(to_string(xdim - 1)) + std::string(".txt");
    std::ofstream	fOut1(fname1);
    std::string fnames1 = std::string("data/Dirichlet/exactsolution_h_") + std::string(to_string(xdim - 1)) + std::string(".txt");
    std::ofstream	fOutsolt1(fnames1);
	std::cout << "\n\nWriting solution to the file...\n\n";

//#pragma omp parallel for
		for (size_t y = 0; y < ydim; ++y) {
			for (size_t x = 0; x < xdim; ++x) {

				fOut1 << x*hsize - 1.0 << "\t" << y*hsize - 1.0 << "\t" << (*xGrids[0])(x, y) << std::endl;
				fOutsolt1 << x*hsize - 1.0 << "\t" << y*hsize - 1.0 << "\t" << (*sGrid)(x, y) << std::endl;
			}
			fOut1 << std::endl;
			fOutsolt1 << std::endl;
		}


		std::cout << "\n\nWriting solution to the file..,,,,.\n\n" << ydim ;
		for (size_t y = ydim - 2; y >= 0; y = y-1) 
		{
//#pragma omp parallel for
			for (size_t x = 0; x < xdim; ++x) 
			{
				fOut1 << x*hsize - 1.0 << "\t" << (ydim-1-y)*hsize << "\t" << (*xGrids[0])(x, y) << std::endl;
				fOutsolt1 << x*hsize - 1.0 << "\t" << (ydim - 1 - y)*hsize << "\t" << (*sGrid)(x, y) << std::endl;
			}
			fOut1 << std::endl;
			fOutsolt1 << std::endl;
		}
	
		fOut1.close();
		fOutsolt1.close();
	
    std::cout << "\n\n =============== Dirichlet Boundary Value Problem 1 ends here ===================\n\n";

    return 0;
}
