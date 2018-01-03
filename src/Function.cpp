#include "Function.h"
#include <TMath.h>
#include <iostream>
#include "BesselFirstKind.h"
#include "Parameters.h"
#include <cmath>
#include <complex>

using namespace std;





/*
double F(double xi, double r, double z, double t) {
    double f, divider;
    double q, z_p;
    const double p2 = omega*mu*sigma*a*a;

    double k = xi/b;

    z_p = z-b;
    q = TMath::Sqrt(k*k+j*p2);

    f = (BESSJ1(k)*BESSJ1(k*r/(a)))*k;
    f *= ((k*mu_r + q)*TMath::Exp((d-z_p)/a) - (k*mu_r - q)*TMath::Exp(-q*(d-z_p)/a));
    divider = ((k*mu_r + q)*TMath::Exp(q*d/a) - (k*mu_r - q)*TMath::Exp(-q*d/a));

    if(isinf(divider)) f=0;
    else f /= divider;

    return f;
}

double F(double xi, double r, double z, double t, double a_i, double b_i, double d_i, double omega_i) {
    double f, divider;
    double q, p2_i, z_p;
    double k = xi/b_i;

    p2_i = omega_i*mu*sigma*a_i*a_i;
    z_p = z-b_i;
    q = TMath::Sqrt(k*k+j*p2_i);  // complex term

    f = (BESSJ1(k)*BESSJ1(k*r/(a_i))/b_i)*xi;
    f *= ((k*mu_r + q)*TMath::Exp((d_i-z_p)/a_i) - (k*mu_r - q)*TMath::Exp(-q*(d_i-z_p)/a_i));
    divider = ((k*mu_r + q)*TMath::Exp(q*d_i/a_i) - (k*mu_r - q)*TMath::Exp(-q*d_i/a_i));

    if(isinf(divider)) f=0;
    else f /= divider;

    return f;

}
*/
