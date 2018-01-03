#include <Parameters.h>
#include <TMath.h>

const double a = 25.0/1000;                                                     //radius of the coil in m
const double b = 20.0/1000;                                                     //distance between the coil and the slab in m
const double d = 20.0/1000;                                                     //thickness of the slab in m
const double mu_0 = 4*TMath::Pi()*TMath::Power(10,-7) ;                         //vacuum permeability in H/m
const double mu_r = 1.000022;                                                   //relative permeability for Aluminum
const double mu = mu_0*mu_r;
const double sigma = 3.5*TMath::Power(10,7);                                    //electrical conductivity in S/m
const double omega = 50;                                                        //
const double j = 1;                                                             //
const double p2 = omega*mu*sigma*a*a;
const double I0 = 1;                                                            //
const int order = 20;                                                           //order of the Gauss-Laguerre quadrature
