#ifndef FUNCTION_H
#define FUNCTION_H
#include <vector>

std::vector< double > A2Z(double x, double y, std::vector< double > prm, double* xi, double* wi);

double R1(double rho_q, double phi_q, double d_opt, double m);
double R2(double rho_q, double phi_q, double d_opt, double m);
double I1(double rho_q, double phi_q, double d_opt, double m);
double I2(double rho_q, double phi_q, double d_opt, double m);

// TEMP TBR
/*
double F(double xi, double r, double z, double t);
double F(double xi, double r, double z, double t, double a_i, double b_i, double d_i, double omega_i);
*/
#endif
