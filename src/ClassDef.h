// rootcint -f Dict.cxx -c src/ClassDef.h src/LinkDef.h
// g++ `root-config --cflags --glibs` -o EddyCurrents EddyCurrents.cxx Dict.cxx src/Function.cpp src/Parameters.cpp src/GaussLaguerre.cpp src/BesselFirstKind.cpp

#include <dirent.h> //<<<< problematic libs
#include "ClassDefRootList.h"
#include "Function.h"
#include "BesselFirstKind.h"
#include "GaussLaguerre.h"
#include "Parameters.h"

using namespace std;

class MyMainFrame : public TGMainFrame {
	//===========================================================================
	//=																			=
	//=							 Variables declaration							=
	//=																			=
	//===========================================================================
	public:
		TGTextView						*TextOutputFrame;
		TGTextButton					*extractCanvasButton, *Plot2DDistribButton;
		TGTextButton					*PlotIntVsVarButton, *PlotSkinDepthButton;
		TGNumberEntry					*aNumberEntry, *bNumberEntry, *dNumberEntry, *omegaNumberEntry;
		TGNumberEntry					*jNumberEntry, *INumberEntry, *orderNumberEntry;
		TGNumberEntry					*var_rangeMinEntry, *var_rangeMaxEntry, *var_rangeEntry;
		TGComboBox						*VarCombo;
		TRootEmbeddedCanvas				*embeddedCanvas;
		TCanvas							*canvas;
		TH1F 							*TH1Buffer;
		TGraph							*TGraphBuffer;
		TH2F 							*TH2Buffer, *dummy_his;
		TH3F 							*TH3Buffer;
		TGHProgressBar 					*ProgressBar;
		TGStatusBar 					*StatusBar;
		TFile 							*rootFile;
		TTree 							*tree;
		char 							StatusBarText[40], output[400];

		int 							wi_size;
		double          				*wi, *xi;
		double 							var_range, var_rangeMax, var_rangeMin;
		vector<double> 					varList;

		//===========================================================================
		//=																			=
		//=							 Functions declaration							=
		//=																			=
		//===========================================================================
		MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h);
		virtual ~MyMainFrame() {
			// Clean up main frame...
			Cleanup();
		}

		//===========================================================================
		//=						 StatusBarUpdate function							=
		//===========================================================================
		void StatusBarUpdate() {
		}

		//===========================================================================
		//=						 ParametersUpdate function							=
		//===========================================================================
		void ParametersUpdate() {
			varList[0] = aNumberEntry->GetNumber()/1000;						//a
			varList[1]  = bNumberEntry->GetNumber()/1000;						//b
			varList[2]  = dNumberEntry->GetNumber()/1000;						//d
			varList[3]  = omegaNumberEntry->GetNumber()*2*TMath::Pi();			//omega
			varList[5]  = INumberEntry->GetNumber();							//I
			varList[6]  = orderNumberEntry->GetNumber();						//order
		}

		//===========================================================================
		//=						 UpdateIntRange function							=
		//===========================================================================
		void UpdateIntRange() {
			var_rangeMin = var_rangeMinEntry->GetNumber();
			var_rangeMax = var_rangeMaxEntry->GetNumber();
			var_range = var_rangeEntry->GetNumber();
		}

		//===========================================================================
		//=				 ComputeGaussLaguerreQuadrature function					=
		//===========================================================================
		void ComputeGaussLaguerreQuadrature(int GLorder) {
			// Get the Gauss-Laguerre roots and weights
			delete xi;
			delete wi;

			wi = new double[GLorder];                                // weights
			xi = new double[GLorder];  								// roots

			wi_size = GLorder;

			cgqf (GLorder, 5, 0.0, 0.0, 0.0, 1, xi, wi);                     //5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a)). See GaussLaguerre.cpp.
		}

		//===========================================================================
		//=							 PlotSkinDepth function							=
		//===========================================================================
		void PlotSkinDepth() {
			double t = 0.5;
			//Get Js = Current at the surface

			this->ParametersUpdate();
			this->UpdateIntRange();
			// Definitions
			double          r, z, integral;
			double          J_phi = 0.0;
			double			z_min = varList[1];

			// BEWARE OF THE VARIABLE INTEGRATION RANGE IN R!!!!!!!!!!!!!!!!!!!!
			double			r_max = 4.0*varList[0];
			double 			p2_i = varList[3]*mu*sigma*varList[0]*varList[0];

			if (wi_size != (int) varList[6]) this->ComputeGaussLaguerreQuadrature((int) varList[6]);

			integral = 0.0;
			r = 4.0*varList[0]/((double) 2*r_range);

			for (size_t r_index = 0; r_index < r_range; r_index++) {
				J_phi = 0.0;

				for (size_t i = 0; i < ((int) varList[6]); i++) {                                        // Sum (i = 1 to n) wi * f(xi)
					J_phi += wi[i]*F(xi[i], r, z_min, t, varList[0], varList[1], varList[2], varList[4], varList[3]);
				}

				J_phi *= -varList[4]*p2_i*varList[5]*TMath::Exp(varList[4]*varList[3]*t)/(varList[0]*varList[0]);                              // Times the front factor

				integral += J_phi;

				r += 4.0*varList[0]/((double) r_range);
			}

			double Js = integral;

			std::cout << "Js = " << Js<< '\n';
			z = z_min;
			z += varList[2]/((double) z_range);
			int counter = 0;
			//Find z for which int J_phi = 0.37*Js => it's the skin depth
			for (size_t z_index = 0; z_index < z_range; z_index++) {
				integral = 0.0;
				r = 4.0*varList[0]/((double) 2*r_range);
				counter++;
				for (size_t r_index = 0; r_index < r_range; r_index++) {
					J_phi = 0.0;

					for (size_t i = 0; i < ((int) varList[6]); i++) {                                        // Sum (i = 1 to n) wi * f(xi)
						J_phi += wi[i]*F(xi[i], r, z, t, varList[0], varList[1], varList[2], varList[4], varList[3]);
					}

					J_phi *= -varList[4]*p2_i*varList[5]*TMath::Exp(varList[4]*varList[3]*t)/(varList[0]*varList[0]);                              // Times the front factor

					integral += J_phi;

					r += 4.0*varList[0]/((double) r_range);
				}
				///if(integral <= 0.37*Js) break;
				std::cout << "0.37*Js = " << Js*0.37 << "\t int " << integral << "\t z " << (z-varList[1])*1000 <<'\n';
				z += varList[2]/((double) z_range);
			}
			std::cout << "integral = " << integral << '\n';
			std::cout << "stopped at depth = " << (z-varList[1])*1000 << '\n';
			std::cout << "counter = " << counter << '\n';
		}

		//===========================================================================
		//=					 Plot2DDistribOfJPhi function							=
		//===========================================================================
		void Plot2DDistribOfJPhi() {
			// Definitions
			double          r, z, t;
			double          J_phi = 0.0;

			this->ParametersUpdate();

			// Setup of the problem :
			// Induced current inside a plate having as an excitation a circular current loop above the plate
			// in a parallel position to the dividing surface alpha*alpha_prime

			// Re-compute GaussLaguerre quadrature rule if necessary :
			if (wi_size != varList[6]) this->ComputeGaussLaguerreQuadrature((int) varList[6]);

			delete TH2Buffer;
			TH2Buffer = new TH2F("TH2Plot","TH2Plot", 	r_range, 0.0, 8.0*varList[0],
														z_range, varList[1], (varList[1]+varList[2]));

			t = 0.5;
			z = varList[1]+varList[2]/((double) z_range)/2;
			const double p2 = varList[3]*mu*sigma*varList[0]*varList[0];

			this->ProgressBarInit(z_range);

			for (size_t z_index = 0; z_index < z_range; z_index++) {
				r = 0.0+8.0*varList[0]/((double) r_range)/2;

				for (size_t r_index = 0; r_index < r_range; r_index++) {
					J_phi = 0.0;

					for (size_t i = 0; i < varList[6]; i++) {                    // Sum (i = 1 to n) wi * f(xi)
						J_phi += wi[i]*F(xi[i], r, z, t);
					}

					J_phi *= -varList[4]*p2*varList[5]*TMath::Exp(varList[4]*varList[3]*t)/(varList[0]*varList[0]);          // Times the front factor

					TH2Buffer->Fill(r, z, J_phi);

					r += 8.0*varList[0]/((double) r_range);
				}
				z += varList[2]/((double) z_range);
				this->ProgressBarUpdate(z_index);
			}

			canvas->cd();
			TH2Buffer->SetStats(kFALSE);
			canvas->SetGridx(1);
			canvas->SetGridy(1);
			TH2Buffer->Draw("colz");
			canvas->Update();
		}

		//===========================================================================
		//=						 J_Phi integral function							=
		//===========================================================================
		double JPhiIntegral(vector<double> var, double t_i) {
			// Definitions
			double          r, z, integral;
			double          J_phi = 0.0;
			double			z_min = var[1];

			// BEWARE OF THE VARIABLE INTEGRATION RANGE IN R!!!!!!!!!!!!!!!!!!!!
			double			r_max = 4.0*var[0];
			double 			p2_i = var[3]*mu*sigma*var[0]*var[0];

			if (wi_size != (int) var[6]) this->ComputeGaussLaguerreQuadrature((int) var[6]);

			z = z_min+var[2]/((double) z_range)/2;
			integral = 0.0;

			for (size_t z_index = 0; z_index < z_range; z_index++) {
				r = 4.0*var[0]/((double) 2*r_range);

				for (size_t r_index = 0; r_index < r_range; r_index++) {
					J_phi = 0.0;

					for (size_t i = 0; i < ((int) var[6]); i++) {                                        // Sum (i = 1 to n) wi * f(xi)
						J_phi += wi[i]*F(xi[i], r, z, t_i, var[0], var[1], var[2], var[4], var[3]);
					}

					J_phi *= -var[4]*p2_i*var[5]*TMath::Exp(var[4]*var[3]*t_i)/(var[0]*var[0]);                              // Times the front factor

					integral += J_phi;

					r += 4.0*var[0]/((double) r_range);
				}
				z += var[2]/((double) z_range);
			}

			return integral;
		}

		//===========================================================================
		//=								 PlotIntVSVar								=
		//===========================================================================
		void PlotIntVsVar() {
			int var_i = VarCombo->GetSelected()-1;
			int ScaleFactor = 1;

			if(var_i < 3) ScaleFactor = 1/1000;
			else if(var_i == 3) ScaleFactor = 2*TMath::Pi();

			this->ParametersUpdate();
			this->UpdateIntRange();

			double x[(size_t) var_range], y[(size_t) var_range];

			delete TGraphBuffer;

			double t = 1;
			double var_value = var_rangeMin;

			this->ProgressBarInit(var_range);

			for (size_t index = 0; index < var_range+1; index++) {
				varList[var_i] = var_value*ScaleFactor;

				x[index] = var_value;
				y[index] = JPhiIntegral(varList, t);

				var_value += ((double) (var_rangeMax-var_rangeMin))/var_range;
				this->ProgressBarUpdate(index);
			}

			TGraphBuffer = new TGraph(var_range, x, y);

			canvas->cd();
			canvas->SetGridx(1);
			canvas->SetGridy(1);
			TGraphBuffer->Draw("ACP");
			canvas->Update();
		}

		//===========================================================================
		//=						 Write text in TextOutputFrame						=
		//===========================================================================
		void WriteInTextFrame(string InputText, bool ClearFirst = true) {
			if (ClearFirst) TextOutputFrame->Clear();
			TextOutputFrame->AddLine(InputText.c_str());
		}

		//===========================================================================
		//=							 Extract Canvas function						=
		//===========================================================================
		void ExtractCanvas() {
		  canvas->DrawClone();
		}

		//===========================================================================
		//=								 ProgressBarInit							=
		//===========================================================================
		void ProgressBarInit(int Max) {
			ProgressBar->SetMax(Max);
			ProgressBar->SetPosition(0);
		}

		//===========================================================================
		//=								 ProgressBarUpdate							=
		//===========================================================================
		void ProgressBarUpdate(int Value) {
			if (Value%((int) (((double) ProgressBar->GetMax())/10 + 0.5)) == 0) {
				ProgressBar->Increment(ProgressBar->GetMax()/10.0);
				gSystem->ProcessEvents();
			}
			if (Value >= ProgressBar->GetMax()-1) {
				ProgressBar->SetPosition(Value+1);
				gSystem->ProcessEvents();
			}
		}

		//===========================================================================
		//=								 Exit Function								=
		//===========================================================================
		void DoExit() {
			delete TH1Buffer;
			gApplication->Terminate(0);
		}

	ClassDef(MyMainFrame, 0)
};
