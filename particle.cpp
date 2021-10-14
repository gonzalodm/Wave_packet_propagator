#include "particle.h"

Particle::Particle(int ind, double mass, double p0, double x0, double aw){

   par_ind  = ind;
   par_mass = mass;
   par_p0   = p0;
   par_x0   = x0;
   par_aw   = aw;

}

int Particle::GetInd()const{
   return par_ind;
}

double Particle::GetMass()const{
   return par_mass;
}

double Particle::GetP0()const{
   return par_p0;
}

double Particle::GetX0()const{
   return par_x0;
}

double Particle::GetAw()const{
   return par_aw;
}

void Particle::SetInd(int ind){
   par_ind = ind;
}
void Particle::SetMass(double mass){
   par_mass = mass;
}
void Particle::SetP0(double p0){
   par_p0 = p0;
}
void Particle::SetX0(double x0){
   par_x0 = x0;
}
void Particle::SetAw(double aw){
   par_aw = aw;
}
