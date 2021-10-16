#ifndef WP_SUBS
#define WP_SUBS

#include "particle.h"

void init_output(ofstream* outfile);

void read_input(double& x_min, double& x_max, UNINT& n_grid, double& dt,
                int& tot_steps, int& print_t, int& print_x, double& x0,
                double& p0, double& aw, double& mass, ifstream& inputf);

double V_potential_in_x(double xi,double mass);

void init_grids(vector<double>& x_grid, vector<double>& dens_grid,
                vector<complex<double> >& psi_grid,
                vector<complex<double> >& psin1_grid,
                vector<complex<double> >& psin2_grid,
                vector<double>& Vx_grid, double x_min, double x_max,
                double& dx, UNINT n_grid, double mass);

void wave_init(vector<double>& x_grid, vector<double>& dens_grid,
               vector<complex<double> >& psi_grid, vector<double>& Vx_grid,
               UNINT n_grid, double ax, double x0, double p0);

void write_output(vector<double>& x_grid, vector<double>& dens_grid,
                  vector<complex<double> >& psi_grid, vector<double>& Vx_grid,
                  UNINT n_grid, int t_step, double dt, int print_x, ofstream* outfile);


#endif
