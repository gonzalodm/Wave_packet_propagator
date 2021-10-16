#include "wp_cuda.h"

cuDoubleComplex *dev_psi;
cuDoubleComplex *dev_psin1;
cuDoubleComplex *dev_psin2;
double          *dev_x;
double          *dev_dens;
double          *dev_Vx;
UNINT Ncores;

__global__ void propagate_psi(cuDoubleComplex *psin_grid,
                              cuDoubleComplex *psi_grid, double *Vx_grid,
                              int n_grid, double u_term, double dt){

    int ind = threadIdx.x + blockIdx.x * blockDim.x;
    cuDoubleComplex idt = make_cuDoubleComplex(0.0e0,-dt);
    cuDoubleComplex u2VC;
    cuDoubleComplex uC;
    cuDoubleComplex s1;
    cuDoubleComplex s2;
    cuDoubleComplex aux1;

    if (ind > 0 && ind < n_grid-1){
        u2VC = make_cuDoubleComplex(2*u_term+Vx_grid[ind],0.0e0);
        uC   = make_cuDoubleComplex(-u_term,0.0e0);
        s1   = cuCadd(psi_grid[ind-1], psi_grid[ind+1]);
        s1   = cuCmul(uC, s1);
        s2   = cuCmul(u2VC, psi_grid[ind]);
        aux1 = cuCadd(s1,s2);
        aux1 = cuCmul(idt,aux1);
        psin_grid[ind] = aux1;
    }

    return;

}

__global__ void update_psi(cuDoubleComplex *psi_grid,
                           cuDoubleComplex *psin1_grid,
                           cuDoubleComplex *psin2_grid,
                           int n_grid){

    int ind = threadIdx.x + blockIdx.x * blockDim.x;
    cuDoubleComplex aux1;
    cuDoubleComplex aux2;
    cuDoubleComplex f_half= make_cuDoubleComplex(0.5e0,0.0e0);

    if (ind > 0 && ind < n_grid-1){
        aux1 = cuCmul(f_half, psin2_grid[ind]);
        aux2 = cuCadd(psin1_grid[ind],aux1);
        psi_grid[ind] = cuCadd(psi_grid[ind],aux2);
    }
    else if (ind == 0 || ind == n_grid-1){
        psi_grid[ind] = make_cuDoubleComplex(0.0e0,0.0e0);
    }

    return;
}

__global__ void calc_dens(cuDoubleComplex *psi_grid, double *dens_grid,
                          int n_grid){

    int    ind = threadIdx.x + blockIdx.x * blockDim.x;
    double vabs;

    if (ind > 0 && ind < n_grid){
        vabs            = cuCabs(psi_grid[ind]);
        dens_grid[ind] = pow(vabs,2.0);
    }

    return;
}

void init_cuda_subs(complex<double> *psi_grid, double *x_grid,
                    complex<double> *psin1_grid,
                    complex<double> *psin2_grid,
                    double *dens_grid,
                    double *Vx_grid, UNINT n_grid){

    double gaux = (double) n_grid;
    double taux = (double) Nthreads;

    Ncores = (UNINT) ceil(gaux/taux);

    cudaMalloc((void**) &dev_psi   , n_grid * sizeof(cuDoubleComplex));
    cudaMalloc((void**) &dev_psin1 , n_grid * sizeof(cuDoubleComplex));
    cudaMalloc((void**) &dev_psin2 , n_grid * sizeof(cuDoubleComplex));
    cudaMalloc((void**) &x_grid    , n_grid * sizeof(double));
    cudaMalloc((void**) &dev_dens  , n_grid * sizeof(double));
    cudaMalloc((void**) &dev_Vx    , n_grid * sizeof(double));

    cudaMemcpy(dev_psi, psi_grid, n_grid * sizeof(cuDoubleComplex),
               cudaMemcpyHostToDevice);
    cudaMemcpy(dev_psin1, psin1_grid, n_grid * sizeof(cuDoubleComplex),
               cudaMemcpyHostToDevice);
    cudaMemcpy(dev_psin2, psin2_grid, n_grid * sizeof(cuDoubleComplex),
               cudaMemcpyHostToDevice);
    cudaMemcpy(dev_x, x_grid, n_grid * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_dens, dens_grid, n_grid * sizeof(double),
               cudaMemcpyHostToDevice);
    cudaMemcpy(dev_Vx, Vx_grid, n_grid * sizeof(double),
               cudaMemcpyHostToDevice);

    return;

}

void propagate_cuda(UNINT n_grid, double mass, double dx, double dt){

    double u_term = 0.5*(1.0/(mass*pow(dx,2.0)));

    propagate_psi<<<Ncores, Nthreads>>>(dev_psin1, dev_psi, dev_Vx, n_grid,
                                        u_term, dt);
    propagate_psi<<<Ncores, Nthreads>>>(dev_psin2, dev_psin1, dev_Vx, n_grid,
                                        u_term, dt);
    update_psi<<<Ncores, Nthreads>>>(dev_psi, dev_psin1,dev_psin2, n_grid);

    return;
}

void get_cuda_psi2(complex<double> *psi_grid, double *dens_grid, int n_grid){

    calc_dens<<<Ncores, Nthreads>>>(dev_psi, dev_dens, n_grid);

    cudaMemcpy(dens_grid, dev_dens, n_grid*sizeof(double),
               cudaMemcpyDeviceToHost);
    cudaMemcpy(psi_grid, dev_psi, n_grid*sizeof(cuDoubleComplex),
               cudaMemcpyDeviceToHost);
    return;
}

void free_cuda_memory(){

    cudaFree(dev_psi);
    cudaFree(dev_psin1);
    cudaFree(dev_psin2);
    cudaFree(dev_x);
    cudaFree(dev_dens);
    cudaFree(dev_Vx);

    return;
}
