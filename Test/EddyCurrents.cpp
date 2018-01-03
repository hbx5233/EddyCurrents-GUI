#include <TMath.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <string>
#include <TH2.h>
#include <TCanvas.h>
#include <Parameters.h>
#include <GaussLaguerre.h>
#include <Function.h>

using namespace std;

int main() {
    // Definitions
    double          *w, *x;
    double          r, z, t;
    double          J_phi = 0.0;
    TH2F            *plot;
    TCanvas			*canvas;

    // Setup of the problem :
    // Induced current inside a plate having as an excitation a circular current loop above the plate
    // in a parallel position to the dividing surface alpha*alpha_prime

    // Get the Gauss-Laguerre roots and weights
    w = new double[order];                                                      // weights
    x = new double[order];                                                      // roots

    cgqf ( order, 5, 0.0, 0.0, 0.0, 1, x, w);                                   //5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a)). See GaussLaguerre.cpp.

    plot = new TH2F("TH2Plot","TH2Plot", 	20, 0.0, 2.0*a,
                                                20, b, b+d);

    t = 0.5;
    z = b;
    r = 0.0;

    for (size_t z_index = 0; z_index < 20; z_index++) {
        for (size_t r_index = 0; r_index < 20; r_index++) {
            for(size_t ActualOrder = 1; ActualOrder <= order; ActualOrder++ ) {
                J_phi = 0.0;

                for (size_t i = 0; i < ActualOrder; i++) {                                        // Sum (i = 1 to n) wi * f(xi)
                    J_phi += w[i]*F(x[i], r, z, t);
                }

                J_phi *= -j*p2*I0*TMath::Exp(j*omega*t)/(a*a);                              // Times the front factor
                plot->Fill(r, z, J_phi);
            }
            r += 2.0*a/20;
        }
        z += (b+d)/20;
    }

    canvas = new TCanvas("c", "c", 450, 300);
    canvas->cd();
	canvas->SetGridx(1);
	canvas->SetGridy(1);
	plot->Draw("colz");
	canvas->Update();
    canvas->SaveAs("plot.pdf");
}
