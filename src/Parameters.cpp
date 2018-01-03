#include "Parameters.h"
#include <TMath.h>


const double mu_0 = 4*TMath::Pi()*TMath::Power(10,-7) ;                         //vacuum permeability in H/m
const double mu_r = 1.000022;                                                   //relative permeability for Aluminum
const double mu = mu_0*mu_r;
const double sigma = 3.5*TMath::Power(10,7);                                    //electrical conductivity in S/m

const int r_range = 250;                                                        //number of "bins" on the r direction (plot and integral calc)
const int z_range = 250;                                                        //number of "bins" on the z direction (plot and integral calc)

// Default values of non-const parameters (that can be changed using the GUI)
double a = 25.0/1000;                                                           //radius of the coil in m
double b = 20.0/1000;                                                           //distance between the coil and the slab in m
double d = 20.0/1000;                                                           //thickness of the slab in m
double omega = 50.0;                                                              //                                                               // 
double I = 1;                                                                  //
int order = 2000;                                                                 //order of the Gauss-Laguerre quadrature
