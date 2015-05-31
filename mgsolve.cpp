#include "Grid.h"
#include<iostream>
#include<fstream>
#include <sys/time.h>
#include <omp.h>
#define XDOMLOW -1.0
#define XDOMHIGH 1.0
#define YDOMLOW -1.0
#define YDOMHIGH 1.0
#define TOLERR 0.0000918
#define V1 2
#define V2 1

Grid ** xGrids = nullptr;
Grid ** fGrids = nullptr;
Grid *sGrid = nullptr;
//bool isNeumann = false;

inline void init(double hsize, const size_t level)
{
    size_t je = level;
    size_t ydim = pow(2, je) + 1;
    size_t xdim = ydim;
    bool flag = true;
    xGrids = (Grid**) memalign(ALLIGNMENT, level*sizeof(Grid*));
    fGrids = (Grid**) memalign(ALLIGNMENT, level*sizeof(Grid*));
//#pragma omp parallel for
    for (size_t i = 0; i < level; i++)
    {
        xGrids[i] = new Grid(xdim, ydim, hsize, hsize, flag);
        fGrids[i] = new Grid(xdim, ydim, hsize, hsize, false);
        ydim = pow(2, --je) + 1;
        xdim = ydim;
        hsize *= 2.0;
        flag = false;
    }

	/*cout << "====After initialization=== \n\n";
	for (size_t j = 0; j < (*xGrids[0]).getYsize(); j++)
	{
		for (size_t k = 0; k < (*xGrids[0]).getXsize(); k++)
	{
		cout << (*xGrids[0])(k, j) << " ";
	}
	cout << '\n';
	}*/
    
}


inline void smooth(Grid* xgrd, const Grid* fgrd, const size_t iter)
{
    size_t dimX = (*xgrd).getXsize();
    size_t dimY = (*xgrd).getYsize();
    double hx = (*xgrd).getHx();
    double hy = (*xgrd).getHy();
   // double	alpha = 1.0;
   // double	beta = 1.0;
    //double	center = 1.0 / (2.0 * alpha + 2.0 * beta);
	size_t midY = (dimY - 1) / 2;
	size_t midX = (dimX - 1) / 2;

	//timeval start, end;
	//gettimeofday(&start, 0);

    for (size_t i = 0; i < iter; i++)
    {
//#pragma omp parallel num_threads(4)
	//	{
		
#pragma omp parallel //firstprivate(dimY,dimX,midX,midY,hx,hy)
		{
		#pragma omp for
			for ( size_t j = 1; j < dimY - 1; j++)
			{
				/*tid1 = omp_get_num_threads();
				tid = omp_get_num_threads();
				std::cout << "inside smooth for " << tid1 << " " << tid << std::endl;
				*/
				size_t l = ((j + 1) & 0x1) + 1;
				for (size_t k = l; k < dimX - 1; k += 2)
				{
					if ((j == midY && k < midX) || j!=midY)
					{
						(*xgrd)(k, j) = (hx*hy*(*fgrd)(k, j) +  ((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) +  ((*xgrd)(k, j + 1)
							+ (*xgrd)(k, j - 1))) * 0.25;
					}
				}

			}
		}
#pragma omp parallel //firstprivate(dimY,dimX,midX,midY,hx,hy)
		{
	#pragma omp for
			for (size_t j = 1; j < dimY - 1; j++)
			{
				/*if (j == 1)
				{
				int tid1 = omp_get_num_threads();
				int tid = omp_get_num_threads();
				std::cout << "inside smooth for loop 2 " << tid1 << " " << tid << std::endl;
				}*/
				size_t l = (j & 0x1) + 1;
				for (size_t k = l; k < dimX - 1; k += 2)
				{
					if ((j == midY && k < midX) || j != midY)
						(*xgrd)(k, j) = (hx*hy*(*fgrd)(k, j) +  ((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) + ((*xgrd)(k, j + 1)
						+ (*xgrd)(k, j - 1))) * 0.25;


				}
			}
		}
    }

	//gettimeofday(&end, 0);
	//double elapsed = 0.000001 * ((double)((end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec));
	//std::cout << "***************:: Time spend for Smoothing= " << elapsed << '\n';
	/*cout << "====After smooth=== \n\n";
	for (size_t j = 0; j < dimX; j++)
	{
		for (size_t k = 0; k < dimY; k++)
		{
			cout << (*xgrd)(k, j) << " ";
		}
		cout << '\n';
	}*/
}

inline void restriction(const Grid * xgrd, const Grid * fgrd, Grid* rgrid)
{
    size_t xlen = (*xgrd).getXsize() - 1;
    size_t ylen = (*xgrd).getYsize() - 1;
    double hx = (*xgrd).getHx();
    double hy = (*xgrd).getHy();
    double	alpha = 1.0 / hx / hx;
    double	beta = 1.0 / hy / hy;
    double	center = (2.0 * alpha) + (2.0 * beta);
	size_t midY = ylen / 2;
	size_t midX = xlen / 2;

	//std::cout << "***************:: restriction= ";
	Grid tmpgrd(xlen + 1, ylen + 1, hx, hy, false);
	//size_t i = 1;
#pragma omp parallel //firstprivate(xlen,ylen,midX,midY,alpha,beta,center)
	{
#pragma omp for
		for (size_t i = 1; i < ylen; i++)
		{
			/*int tid1 = omp_get_num_threads();
			int tid = omp_get_num_threads();
			std::cout << "inside restriction residual for " << tid1 << " " << tid << std::endl;*/
			for (size_t j = 1; j < xlen; j++)
			{
				if ((i == midY && j < midX) || i != midY)
					tmpgrd(j, i) = (*fgrd)(j, i) + alpha*((*xgrd)(j + 1, i) + (*xgrd)(j - 1, i)) + beta * ((*xgrd)(j, i + 1)
					+ (*xgrd)(j, i - 1)) - (*xgrd)(j, i) * center;
			}


		}
	}

    size_t rxlen = (*rgrid).getXsize() - 1;
    size_t rylen = (*rgrid).getYsize() - 1;

	midY = rylen / 2;
	midX = rxlen / 2;
#pragma omp parallel //firstprivate(rxlen,rylen,midX,midY)
	{
#pragma omp for
		for (size_t i = 1; i < rylen; i++)
		{
			/*int tid1 = omp_get_num_threads();
			int tid = omp_get_num_threads();
			std::cout << "inside restriction for " << tid1 << " " << tid << std::endl;*/
			for (size_t j = 1; j < rxlen; j++)
			{
				if ((i == midY && j < midX) || i != midY)
					(*rgrid)(j, i) = (tmpgrd(2 * j - 1, 2 * i - 1) + tmpgrd(2 * j - 1, 2 * i + 1) +
					tmpgrd(2 * j + 1, 2 * i - 1) + tmpgrd(2 * j + 1, 2 * i + 1)) * 0.0625 +
					0.125 *(tmpgrd(2 * j, 2 * i - 1) + tmpgrd(2 * j, 2 * i + 1) +
					tmpgrd(2 * j - 1, 2 * i) + tmpgrd(2 * j + 1, 2 * i)) + 0.25 * tmpgrd(2 * j, 2 * i);
			}

		}
	}
   }

inline void interpolate(Grid * srcgrd, Grid * tgtgrd)
{
    //size_t len = (*srcgrd).getXsize() - 1;
    size_t txlen = (*tgtgrd).getXsize()-1;
	size_t tylen = (*tgtgrd).getYsize()-1;
	
	//ssize_t i;
	//std::cout << "***************:: interpolation= ";
#pragma omp parallel  //firstprivate(txlen,tylen) 
	{
#pragma omp for
		for (size_t i = 1; i < tylen; i += 2)
		{
			/*int tid1 = omp_get_num_threads();
			int tid = omp_get_num_threads();
			std::cout << "inside interpoalte for " << tid1 << " " << tid << std::endl;*/
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

//inline void resdualNorm(const Grid* xgrd, const Grid * fgrd, double* norm)
//{
//
//    size_t dimX = (*xgrd).getXsize() - 1;
//    size_t dimY = (*xgrd).getYsize() - 1;
//    double r = 0.0;
//    double hx = (*xgrd).getHx();
//    double hy = (*xgrd).getHy();
//   // double	alpha = 1.0;
//   // double	beta = 1.0;
//	//double	center = 4.0;
//
//
//	size_t midY = dimY / 2;
//	size_t midX = dimX / 2;
//	//size_t j;
//	double sum = 0.0;
//	//std::cout << "***************:: Residual= ";
//#pragma omp parallel reduction(+: sum) firstprivate(dimX,dimY,midX,midY,hx,hy,r)
//	{
//#pragma omp for
//		for (size_t j = 1; j < dimY; j++)
//		{
//			/*int tid1 = omp_get_num_threads();
//			int tid = omp_get_num_threads();
//			std::cout << "inside residual for " << tid1 << " " << tid << std::endl;*/
//			for (size_t k = 1; k < dimX; k++)
//			{
//				if ((j == midY && k < midX) || j != midY)
//				{
//					r = hx*hy*(*fgrd)(k, j) + ((*xgrd)(k + 1, j) + (*xgrd)(k - 1, j)) + ((*xgrd)(k, j + 1)
//						+ (*xgrd)(k, j - 1)) - (*xgrd)(k, j) * 4.0;
//
//					sum += r*r;
//				}
//			}
//		}
//	}
//
//        *norm = sqrt(sum / (dimX - 1) / (dimY - 1));
//}



inline double errorNorm(const Grid* xgrd, const Grid * sgrd)
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
				sum += r*r;
			}

		}
	}

    return sqrt(sum / dimX / dimY);
}

void mgsolve(size_t level)
{
    size_t gdim = pow(2, level) + 1;
    double newnorm = 1.0;
    double hsize = (XDOMHIGH - XDOMLOW) / (gdim - 1.0);
	size_t i = 0;
    init(hsize, level);
	sGrid = new Grid(gdim, gdim, hsize, hsize, true);
	size_t k;
#pragma omp parallel private(k) firstprivate(gdim,hsize)
	{
#pragma omp for
		for (k = 0; k < gdim; k++)
		{
			//#pragma omp parallel for
			for (size_t j = 0; j < gdim; j++)
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
		//std::cout << "Dirichlet:: Error L2 Norm for h as 1/" << gdim - 1 << " = " << newnorm << "\n\n";	

    }
   // vcycle = i;   
    
    std::cout << "Dirichlet:: Error L2 Norm for h as 1/" << gdim - 1 << " after " << i << " V-Cycle = " << newnorm << "\n\n";
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
	//std::cout << "333";
    gettimeofday(&start, 0);
	//std::cout << "22222";
    mgsolve(level);

    gettimeofday(&end, 0);
    double elapsed = 0.000001 * ((double)((end.tv_sec - start.tv_sec) * 1000000 + end.tv_usec - start.tv_usec));
    std::cout << "Dirichlet:: Time spend for Multigrid Solver = " << elapsed << '\n';

    double hsize = (*xGrids[0]).getHx();
    size_t gdim = (*xGrids[0]).getXsize();

    std::string fname1 = std::string("data/Dirichlet/solution_h_") + std::string(to_string(gdim - 1)) + std::string(".txt");
    std::ofstream	fOut1(fname1);
    std::string fnames1 = std::string("data/Dirichlet/exactsolution_h_") + std::string(to_string(gdim - 1)) + std::string(".txt");
    std::ofstream	fOutsolt1(fnames1);
	std::cout << "\n\nWriting solution to the file...\n\n";

//#pragma omp parallel for
		for (size_t y = 0; y < gdim; ++y) {
			for (size_t x = 0; x < gdim; ++x) {

				fOut1 << x*hsize - 1.0 << "\t" << y*hsize - 1.0 << "\t" << (*xGrids[0])(x, y) << std::endl;
				fOutsolt1 << x*hsize - 1.0 << "\t" << y*hsize - 1.0 << "\t" << (*sGrid)(x, y) << std::endl;
			}
			fOut1 << std::endl;
			fOutsolt1 << std::endl;
		}
	
		fOut1.close();
		fOutsolt1.close();
	
    std::cout << "\n\n =============== Dirichlet Boundary Value Problem 1 ends here ===================\n\n";

    return 0;
}
