#include "Function.h"
#include <TMath.h>
#include <iostream>
#include "BesselFirstKind.h"
#include "Parameters.h"
#include <cmath>
#include <complex>
#include <vector>

using namespace TMath;

double R1(double d_opt, double rho_q, double phi_q, double m) {
    double R1_output = m * mu_r * Cos( rho_q * d_opt * Sin( phi_q ) );

    R1_output += rho_q * Cos( phi_q + rho_q * d_opt * Sin( phi_q ) );

    R1_output *= Exp( rho_q * d_opt * Cos( phi_q ) );

    return R1_output;
}

double R2(double d_opt, double rho_q, double phi_q, double m) {
    double R2_output = m * mu_r * Cos( rho_q * d_opt * Sin( phi_q ) );

    R2_output -= rho_q * Cos( phi_q - rho_q * d_opt * Sin( phi_q ) );

    R2_output *= Exp( -rho_q * d_opt * Cos( phi_q ) );

    return R2_output;
}

double I1(double d_opt, double rho_q, double phi_q, double m) {
    double I1_output = m * mu_r * Sin( rho_q * d_opt * Sin( phi_q ) );

    I1_output += rho_q * Sin( phi_q + rho_q * d_opt * Sin( phi_q ) );

    I1_output *= Exp( rho_q * d_opt * Cos( phi_q ) );

    return I1_output;
}

double I2(double d_opt, double rho_q, double phi_q, double m) {
    double I2_output = m * mu_r * Sin( rho_q * d_opt * Sin( phi_q ) );

    I2_output += rho_q * Sin( phi_q - rho_q * d_opt * Sin( phi_q ) );

    I2_output *= - Exp( -rho_q * d_opt * Cos( phi_q ) );

    return I2_output;
}

std::vector< double > A2Z(double x, double y, double prm[], double* xi, double* wi) {
    std::vector< double > output;
    double rho_A, phi_A, rho_B, phi_B, rho_q, phi_q;
    double NR, NI, DR, DI;
    double GBR, GBI;
    double INTR = 0.0;
    double INTI = 0.0;
    double RA2, IA2;
    double wt = prm[3] * prm[5];

    // Parameters array (prm[]) is the list of all the parameters used to compute the Eddy EddyCurrents, resp.
    // index            parameter
    // 0                a
    // 1                b
    // 2                d
    // 3                omega
    // 4                I0
    // 5                t

    double p = Sqrt( prm[3] * mu * sigma );

    // TODO check import of xi and wi tables and their size calculation
    int GL_order = sizeof(xi)/sizeof(xi[0]);

    // TEMP TBR
    std::cout << "array size : " << GL_order << '\n';

    // Sum over all the elements wi*f(xi)
    for (size_t act_order = 0; act_order < GL_order; act_order++) {
        // Calculate rho_A, rho_B, phi_A, phi_B, rho_q and phi_q
        rho_A = Sqrt( Power(xi[act_order], 2) + Sqrt(2)*xi[act_order]*p +  Power(p, 2) );
        rho_B = Sqrt( Power(xi[act_order], 2) - Sqrt(2)*xi[act_order]*p +  Power(p, 2) );
        phi_A = ATan( -p / ( Sqrt(2)*xi[act_order] + p ) );
        phi_B = ATan( p / ( Sqrt(2)*xi[act_order] - p ) );
        rho_q = Sqrt( rho_A * rho_B );
        phi_q = ( phi_A + phi_B ) / 2;

        // Calculate NR, NI, DR, DI
        // args =  double d_opt, double rho_q, double phi_q, double m
        // d_opt is either 'd'  (prm[2]) or 'd+y' (prm[2] + y)

        NR = R1(prm[2] + y, rho_q, phi_q, xi[act_order]) + R2(prm[2] + y, rho_q, phi_q, xi[act_order]);
        NI = I1(prm[2] + y, rho_q, phi_q, xi[act_order]) + I2(prm[2] + y, rho_q, phi_q, xi[act_order]);

        DR = R1(prm[2], rho_q, phi_q, xi[act_order]) - R2(prm[2], rho_q, phi_q, xi[act_order]);
        DI = I1(prm[2], rho_q, phi_q, xi[act_order]) - I2(prm[2], rho_q, phi_q, xi[act_order]);

        GBR = ( NR * DR + NI * DI ) / ( Power(DR, 2) + Power(DI, 2) );
        GBI = ( NI * DR - NR * DI ) / ( Power(DR, 2) + Power(DI, 2) );

        INTR += wi[act_order] * Cos( xi[act_order] * x ) * GBR;
        INTI += wi[act_order] * Cos( xi[act_order] * x ) * GBI;
    }

    RA2 = Cos( wt ) * INTR - Sin( wt ) * INTI;
    RA2 *= mu * prm[4] / Pi();

    IA2 = Cos( wt ) * INTI + Sin( wt ) * INTR;
    IA2 *= mu * prm[4] / Pi();

    output.push_back(RA2);
    output.push_back(IA2);

    return output;
}


// TEMP TBR
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
