#include "particle.h"
#include "wp_subs.h"
#include "wp_cuda.h"

int main(){

   double                     x_min;
   double                     x_max;
   vector<double>             x_grid;
   UNINT                      n_grid;
   double                     dt;
   double                     dx;
   double                     x0;
   double                     p0;
   double                     aw;
   double                     mass;
   int                        tot_steps;
   int                        print_t;
   int                        print_x;
   vector <double>            dens_grid;
   vector <double>            Vx_grid;
   vector < complex<double> > psi_grid;
   vector < complex<double> > psin1_grid;
   vector < complex<double> > psin2_grid;
   const double               pi = 3.141592653589793;
   ifstream                   inputf;
   ofstream                   outfile[6];


   init_output(outfile);
   read_input(x_min, x_max, n_grid, dt, tot_steps, print_t, print_x, x0, p0,
              aw, mass, inputf);
   Particle part1(1, mass, p0, x0, aw);
   init_grids(x_grid, dens_grid, psi_grid, psin1_grid, psin2_grid, Vx_grid,
              x_min, x_max, dx, n_grid, mass);
   wave_init(x_grid, dens_grid, psi_grid, Vx_grid, n_grid, aw, x0, p0);

   write_output(x_grid, dens_grid, psi_grid, Vx_grid, n_grid, 0, dt,
                print_x, outfile);

   init_cuda_subs(& *psi_grid.begin(), & *x_grid.begin(),
                  & *psin1_grid.begin(), & *psin2_grid.begin(),
                  & *dens_grid.begin(),
                  & *Vx_grid.begin(), n_grid);

   for (int tt=1; tt<=tot_steps; tt++){
      propagate_cuda(n_grid, mass, dx, dt);
      if (tt%print_t==0){
         get_cuda_psi2(& *psi_grid.begin(), & *dens_grid.begin(), n_grid);
         write_output(x_grid, dens_grid, psi_grid, Vx_grid, n_grid, tt, dt,
                      print_x, outfile);
      }
   }

   free_cuda_memory();
   for (int ii=0; ii<5; ii++){outfile[ii].close();}
   return 0;
};
