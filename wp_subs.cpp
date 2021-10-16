#include "wp_subs.h"

//##############################################################################
void init_output(ofstream* outfile){
    outfile[0].open("xgrid.out");
    outfile[1].open("psi_re.out");
    outfile[2].open("psi_im.out");
    outfile[3].open("density.out");
    outfile[4].open("time.out");
    outfile[5].open("V_grid.out");

    return;
}
//##############################################################################
void read_input(double& x_min, double& x_max, UNINT& n_grid, double& dt,
                int& tot_steps, int& print_t, int& print_x, double& x0,
                double& p0, double& aw, double& mass, ifstream& inputf){

    string str;

    inputf.open("input.in");
    if (!inputf) {
        cout << "Unable to open file";
        exit(1); // terminate with error
    }

    string keys[] = {
        "x_min",
        "x_max",
        "n_grid",
        "x0",
        "p0",
        "aw",
        "mass",
        "dt",
        "tot_steps",
        "print_t",
        "print_x"
    };

    while (getline(inputf, str))
    {
    //cout << str << "\n";
    for(int jj=0; jj<12; jj++)
    {
      size_t found = str.find(keys[jj]);
      if (found != string::npos)
      {
        stringstream   linestream(str);
        string         data;
        getline(linestream, data, '=');
        if(jj==0)       linestream >> x_min;
        else if(jj==1)  linestream >> x_max;
        else if(jj==2)  linestream >> n_grid;
        else if(jj==3)  linestream >> x0;
        else if(jj==4)  linestream >> p0;
        else if(jj==5)  linestream >> aw;
        else if(jj==6)  linestream >> mass;
        else if(jj==7)  linestream >> dt;
        else if(jj==8)  linestream >> tot_steps;
        else if(jj==9)  linestream >> print_t;
        else if(jj==10) linestream >> print_x;
      }
    }

    }

    return;
}
//##############################################################################
double V_potential_in_x(double xi, double mass){
    double Vx;
    // double x0 = 0.40e0;
    double x0 = 0.5e0;
    double xf = 0.41e0;

    // if (xi >= x0 && xi <= xf){
    // Vx = 3.7e5;
    // }
    // else{
        // Vx = 0.0e0;
    // }

    double omega = 1000.0;
    Vx = 0.5*mass*pow(omega,2)*pow(xi-x0,2);

    return Vx;
}
//##############################################################################
void init_grids(vector<double>& x_grid, vector<double>& dens_grid,
                vector<complex<double> >& psi_grid,
                vector<complex<double> >& psin1_grid,
                vector<complex<double> >& psin2_grid,
                vector<double>& Vx_grid, double x_min, double x_max,
                double& dx, UNINT n_grid, double mass){

    double xi = x_min;
    dx = (x_max-x_min)/double(n_grid-1);

    for(int ii=0; ii<n_grid; ii++){
        x_grid.push_back(xi);
        dens_grid.push_back(0.0e0);
        psi_grid.push_back((0.0e0,0.0e0));
        psin1_grid.push_back((0.0e0,0.0e0));
        psin2_grid.push_back((0.0e0,0.0e0));
        Vx_grid.push_back(V_potential_in_x(xi,mass));
        xi += dx;
    }

    return;
}
//##############################################################################
void wave_init(vector<double>& x_grid, vector<double>& dens_grid,
               vector<complex<double> >& psi_grid, vector<double>& Vx_grid,
               UNINT n_grid, double ax, double x0, double p0){

    const double pi = 3.141592653589793;
    double       nf;

    nf     = pow(2.0e0*ax/pi, 0.25e0);
    for(int ii=1; ii<n_grid-1; ii++){
        double xi     = x_grid[ii];
        double e_aux  = exp(-ax*pow((xi-x0),2.0e0));
        double psi_re = nf * e_aux * cos(p0*(xi-x0));
        double psi_im = nf * e_aux * sin(p0*(xi-x0));
        psi_grid[ii]  = complex<double>(psi_re, psi_im);
        dens_grid[ii] = pow(psi_re, 2.0e0) + pow(psi_im, 2.0e0);
    }

    return;
}
//##############################################################################
void write_output(vector<double>& x_grid, vector<double>& dens_grid,
                  vector<complex<double> >& psi_grid, vector<double>& Vx_grid,
                  UNINT n_grid, int t_step, double dt, int print_x, ofstream* outfile){

    if (t_step == 0){
        for (int ii; ii<n_grid; ii++){
            if(ii%print_x==0){
                outfile[0]<< scientific << x_grid[ii] <<endl;
                outfile[5]<< scientific << Vx_grid[ii] <<endl;
            }
        }
    }

    for (int ii; ii<n_grid; ii++){
        if(ii%print_x==0){
            outfile[1]<<real(psi_grid[ii]) <<endl;
            outfile[2]<<imag(psi_grid[ii]) <<endl;
            outfile[3]<<dens_grid[ii]      <<endl;
        }
    }

    outfile[4]<< scientific << t_step * dt <<endl;

    return;
}
//##############################################################################
