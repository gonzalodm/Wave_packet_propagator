#ifndef WP_CUDA
#define WP_CUDA

#include "particle.h"

extern cuDoubleComplex *dev_psi;
extern cuDoubleComplex *dev_psin1;
extern cuDoubleComplex *dev_psin2;
extern double          *dev_x;
extern double          *dev_dens;
extern double          *dev_Vx;
extern UNINT            Ncores;
const UNINT      Nthreads = 512;

void init_cuda_subs(complex<double> *psi_grid, double *x_grid,
                    complex<double> *psin1_grid,
                    complex<double> *psin2_grid,
                    double *dens_grid,
                    double *Vx_grid, UNINT n_grid);

void propagate_cuda(UNINT n_grid, double mass, double dx, double dt);

void get_cuda_psi2(complex<double> *psi_grid, double *dens_grid, int n_grid);
void free_cuda_memory();

#endif
