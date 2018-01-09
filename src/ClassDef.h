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
		TGNumberEntry					*I0NumberEntry, *tNumberEntry, *orderNumberEntry;
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

		double          				*wi, *xi;
		double 							var_range, var_rangeMax, var_rangeMin;
		vector<double> 					prm_list;

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

			if (prm_list[1] != bNumberEntry->GetNumber()/1000 || prm_list[6] != orderNumberEntry->GetNumber()) {	// Update b and recalculate the GL parameters if b is updated
				this->ComputeGaussLaguerreQuadrature(orderNumberEntry->GetNumber());
			}

			prm_list[0] = aNumberEntry->GetNumber()/1000;						//a
			prm_list[1] = bNumberEntry->GetNumber()/1000;						//b
			prm_list[2] = dNumberEntry->GetNumber()/1000;						//d
			prm_list[3] = omegaNumberEntry->GetNumber()*2*TMath::Pi();			//omega
			prm_list[4] = I0NumberEntry->GetNumber();							//I0
			prm_list[5] = tNumberEntry->GetNumber();							//t
			prm_list[6] = orderNumberEntry->GetNumber();						//order
		}

		//===========================================================================
		//=					 ParametersOverload function							=
		//===========================================================================
		void ParametersOverload(int param_ID, double new_value) {
			prm_list[param_ID] = new_value;

			// Recalculate GL if required
			if( param_ID == 6 || param_ID == 1 ) this->ComputeGaussLaguerreQuadrature((int) prm_list[6]);
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
		void ComputeGaussLaguerreQuadrature(int GL_order) {
			// Get the Gauss-Laguerre roots and weights, if needed
			delete xi;
			delete wi;

			wi = new double[GL_order];                                // weights
			xi = new double[GL_order];  								// roots

			cgqf (GL_order, 5, 0.0, 0.0, 0.0, prm_list[1], xi, wi);                     //5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a)). See GaussLaguerre.cpp.
		}

		//===========================================================================
		//=							 PlotSkinDepth function							=
		//===========================================================================
		void PlotSkinDepth() {
		/*	double t = 0.5;
			//Get Js = Current at the surface

			this->ParametersUpdate();
			this->UpdateIntRange();
			// Definitions
			double          r, z, integral;
			double          J_phi = 0.0;
			double			z_min = prm_list[1];

			// BEWARE OF THE VARIABLE INTEGRATION RANGE IN R!!!!!!!!!!!!!!!!!!!!
			double			r_max = 4.0*prm_list[0];
			double 			p2_i = prm_list[3]*mu*sigma*prm_list[0]*prm_list[0];

			this->ComputeGaussLaguerreQuadrature((int) prm_list[6]);

			integral = 0.0;
			r = 4.0*prm_list[0]/((double) 2*r_range);

			for (size_t r_index = 0; r_index < r_range; r_index++) {
				J_phi = 0.0;

				for (size_t i = 0; i < ((int) prm_list[6]); i++) {                                        // Sum (i = 1 to n) wi * f(xi)
					J_phi += wi[i]*F(xi[i], r, z_min, t, prm_list[0], prm_list[1], prm_list[2], prm_list[4], prm_list[3]);
				}

				J_phi *= -prm_list[4]*p2_i*prm_list[4]*TMath::Exp(prm_list[4]*prm_list[3]*t)/(prm_list[0]*prm_list[0]);                              // Times the front factor

				integral += J_phi;

				r += 4.0*prm_list[0]/((double) r_range);
			}

			double Js = integral;

			std::cout << "Js = " << Js<< '\n';
			z = z_min;
			z += prm_list[2]/((double) z_range);
			int counter = 0;
			//Find z for which int J_phi = 0.37*Js => it's the skin depth
			for (size_t z_index = 0; z_index < z_range; z_index++) {
				integral = 0.0;
				r = 4.0*prm_list[0]/((double) 2*r_range);
				counter++;
				for (size_t r_index = 0; r_index < r_range; r_index++) {
					J_phi = 0.0;

					for (size_t i = 0; i < ((int) prm_list[6]); i++) {                                        // Sum (i = 1 to n) wi * f(xi)
						J_phi += wi[i]*F(xi[i], r, z, t, prm_list[0], prm_list[1], prm_list[2], prm_list[4], prm_list[3]);
					}

					J_phi *= -prm_list[4]*p2_i*prm_list[5]*TMath::Exp(prm_list[4]*prm_list[3]*t)/(prm_list[0]*prm_list[0]);                              // Times the front factor

					integral += J_phi;

					r += 4.0*prm_list[0]/((double) r_range);
				}
				///if(integral <= 0.37*Js) break;
				std::cout << "0.37*Js = " << Js*0.37 << "\t int " << integral << "\t z " << (z-prm_list[1])*1000 <<'\n';
				z += prm_list[2]/((double) z_range);
			}
			std::cout << "integral = " << integral << '\n';
			std::cout << "stopped at depth = " << (z-prm_list[1])*1000 << '\n';
			std::cout << "counter = " << counter << '\n';
			*/
		}

		//===========================================================================
		//=					 Plot2DDistribOfA2 function							=
		//===========================================================================
	void Plot2DDistribOfA2() {
				// Definitions
			double          x, y, A2;

			this->ParametersUpdate();

			// Setup of the problem :
			// Induced current inside a plate having as an excitation a wire above the plate
			// in a parallel position to the dividing surface alpha*alpha_prime

			// Re-compute GaussLaguerre quadrature rule if necessary :
			this->ComputeGaussLaguerreQuadrature((int) prm_list[6]);

			delete TH2Buffer;
			// TEMP check def of TH2
			TH2Buffer = new TH2F("TH2Plot","TH2Plot", 	x_range, 0.0, 8.0*prm_list[0],
														y_range, prm_list[1], (prm_list[1]+prm_list[2]));

			// TODO Clean the procedure to scan the z range
			y = prm_list[1]+prm_list[2]/((double) y_range)/2;

			// TODO change or load the var list into the prm[] table to call the A2 Function

			this->ProgressBarInit(y_range);

			for (size_t y_index = 0; y_index < y_range; y_index++) {
				x = 0.0+8.0*prm_list[0]/((double) x_range)/2;

				for (size_t x_index = 0; x_index < x_range; x_index++) {
					// TEMP test real part and then norm for A2
					// TODO define the PRM ARRAY INPUT
					A2 = A2Z(x, y, prm_list, xi, wi).at(0); // Real part

					TH2Buffer->Fill(x, y, A2);

					// TODO Clean the procedure to scan the r range
					x += 8.0*prm_list[0]/((double) x_range);
				}
				y += prm_list[2]/((double) y_range);
				this->ProgressBarUpdate(y_index);
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
/*
			// Definitions
			double          r, z, integral;
			double          J_phi = 0.0;
			double			z_min = var[1];

			// BEWARE OF THE VARIABLE INTEGRATION RANGE IN R!!!!!!!!!!!!!!!!!!!!
			double			r_max = 4.0*var[0];
			double 			p2_i = var[3]*mu*sigma*var[0]*var[0];

			this->ComputeGaussLaguerreQuadrature((int) var[6]);

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
			*/
		}

		//===========================================================================
		//=								 PlotIntVSVar								=
		//===========================================================================
		void PlotIntVsVar() {
			/*
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
				prm_list[var_i] = var_value*ScaleFactor;

				x[index] = var_value;
				y[index] = JPhiIntegral(prm_list, t);

				var_value += ((double) (var_rangeMax-var_rangeMin))/var_range;
				this->ProgressBarUpdate(index);
			}

			TGraphBuffer = new TGraph(var_range, x, y);

			canvas->cd();
			canvas->SetGridx(1);
			canvas->SetGridy(1);
			TGraphBuffer->Draw("ACP");
			canvas->Update();
			*/
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
