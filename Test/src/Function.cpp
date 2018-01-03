#include <Function.h>
#include <TMath.h>
#include <BesselFirstKind.h>
#include <Parameters.h>

double F(double xi, double r, double z, double t) {
    double f;
    double q, p2, z_p;

    z_p = z-b;
    q = TMath::Sqrt(xi*xi/(b*b)+j*p2);

    f = (BESSJ1(xi/b)*BESSJ1(xi*r/(a*b))/b);
    f *= ((xi*mu_r/b + q)*TMath::Exp((d-z_p)/a) - (xi*mu_r/b - q)*TMath::Exp(-q*(d-z_p)/a))*xi;
    f /= ((xi*mu_r/b + q)*TMath::Exp(q*d/a) - (xi*mu_r/b - q)*TMath::Exp(-q*d/a));

    return f;
}
