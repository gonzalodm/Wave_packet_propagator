#ifndef CLASS_Particle
#define CLASS_Particle

#include <iostream>
#include <vector>
#include <string>
#include <complex>
#include <time.h>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <limits>
#include <cuda.h>
#include <cuComplex.h>
#include <ctime>
#include <cstdlib>

using namespace std;
typedef unsigned int UNINT;

class Particle
{
public:
   Particle(int ind, double mass, double p0, double x0, double aw);

   int    GetInd()   const;
   double GetMass()  const;
   double GetP0()  const;
   double GetX0()  const;
   double GetAw()  const;
   void   SetInd(int ind);
   void   SetMass(double mass);
   void   SetP0(double p0);
   void   SetX0(double x0);
   void   SetAw(double aw);

private:
   int    par_ind;
   double par_mass;
   double par_p0;
   double par_x0;
   double par_aw;

};

#endif
