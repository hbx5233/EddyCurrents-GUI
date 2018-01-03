// rootcint -f Dict.cxx -c src/ClassDef.h src/LinkDef.h
// g++ `root-config --cflags --glibs` -o EddyCurrents EddyCurrents.cxx Dict.cxx

#include <dirent.h>
#include <TSystem.h>
#include <TArray.h>
#include <TMath.h>
#include <TGWindow.h>
#include <TGFrame.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TArray.h>
#include <TImage.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTimeStamp.h>
#include <TLegend.h>
#include <TF1.h>
#include <TGTextEntry.h>
#include <TGNumberEntry.h>
#include <TGTextView.h>
#include <TGShutter.h>
#include <TGStatusBar.h>
#include <TGButtonGroup.h>
#include <TROOT.h>
#include <TDatime.h>
#include <TGTab.h>
#include <TGLabel.h>
#include <TGComboBox.h>
#include <TApplication.h>
#include <TGProgressBar.h>
#include <TRootEmbeddedCanvas.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <stdlib.h>
#include <algorithm>

//#include <TGRadioButton.h>

#include "src/ClassDef.h"

using namespace std;

//===========================================================================
//=																			=
//=							 GUI Construction								=
//=																			=
//===========================================================================
MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h) : TGMainFrame(p, w, h) {

	//===========================================================================
	//=				 Variables initialization and constant definition			=
	//===========================================================================
	wi = new double[order];
	xi = new double[order];
	wi_size = -1;

	varList.push_back(a);
	varList.push_back(b);
	varList.push_back(d);
	varList.push_back(omega*2*TMath::Pi());
	varList.push_back(j);
	varList.push_back(I);
	varList.push_back(order);

	var_range = 50;
	var_rangeMin = 1;
	var_rangeMax = 5;

	//===========================================================================
	//=				 Initializing the histograms								=
	//===========================================================================
	TGraphBuffer = new TGraph(1);
	TH1Buffer = new TH1F("TH1Buffer","TH1Buffer",1,0.,1.);
	TH2Buffer = new TH2F("TH2Buffer","TH2Buffer",1,0.,1.,1,0.,1.);
	dummy_his = new TH2F("dummy_his","dummy_his",1,0.,1.,1,0.,1.);
	TH3Buffer = new TH3F("TH3Buffer","TH3Buffer",1,0.,1.,1,0.,1.,1,0.,1.);

	//===========================================================================
	//=				 Creatings buffers											=
	//===========================================================================

	//===========================================================================
	//=				File Name and number horz Top frame 						=
	//===========================================================================
	TGHorizontalFrame *horzFrameTopParam= new TGHorizontalFrame(this);
	AddFrame(horzFrameTopParam, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 2, 2, 2, 2));

	TGLabel *aParamLabel = new TGLabel(horzFrameTopParam, "a (mm) : ");
	horzFrameTopParam->AddFrame(aParamLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 0, 0, 0));

	aNumberEntry = new TGNumberEntry(horzFrameTopParam, 1,1);
	horzFrameTopParam->AddFrame(aNumberEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));
	aNumberEntry->Resize(75,20);
	aNumberEntry->SetNumber(a*1000);

	TGLabel *bParamLabel = new TGLabel(horzFrameTopParam, "b (mm) : ");
	horzFrameTopParam->AddFrame(bParamLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 0, 0, 0));

	bNumberEntry = new TGNumberEntry(horzFrameTopParam, 1,1);
	horzFrameTopParam->AddFrame(bNumberEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));
	bNumberEntry->Resize(75,20);
	bNumberEntry->SetNumber(b*1000);

	TGLabel *dParamLabel = new TGLabel(horzFrameTopParam, "d (mm) : ");
	horzFrameTopParam->AddFrame(dParamLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 0, 0, 0));

	dNumberEntry = new TGNumberEntry(horzFrameTopParam, 1,1);
	horzFrameTopParam->AddFrame(dNumberEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));
	dNumberEntry->Resize(75,20);
	dNumberEntry->SetNumber(d*1000);

	TGLabel *omegaParamLabel = new TGLabel(horzFrameTopParam, "f (Hz) : ");
	horzFrameTopParam->AddFrame(omegaParamLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 0, 0, 0));

	omegaNumberEntry = new TGNumberEntry(horzFrameTopParam, 1,1);
	horzFrameTopParam->AddFrame(omegaNumberEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));
	omegaNumberEntry->Resize(75,20);
	omegaNumberEntry->SetNumber(omega);

	TGHorizontalFrame *horzFrameTopParam2= new TGHorizontalFrame(this);
	AddFrame(horzFrameTopParam2, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 2, 2, 2, 2));

	TGLabel *jParamLabel = new TGLabel(horzFrameTopParam2, "j : ");
	horzFrameTopParam2->AddFrame(jParamLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 0, 0, 0));

	jNumberEntry = new TGNumberEntry(horzFrameTopParam2, 1,1);
	horzFrameTopParam2->AddFrame(jNumberEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));
	jNumberEntry->Resize(75,20);
	jNumberEntry->SetNumber(j);

	TGLabel *IParamLabel = new TGLabel(horzFrameTopParam2, "I : ");
	horzFrameTopParam2->AddFrame(IParamLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 0, 0, 0));

	INumberEntry = new TGNumberEntry(horzFrameTopParam2, 1,1);
	horzFrameTopParam2->AddFrame(INumberEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));
	INumberEntry->Resize(75,20);
	INumberEntry->SetNumber(I);

	TGLabel *orderParamLabel = new TGLabel(horzFrameTopParam2, "order : ");
	horzFrameTopParam2->AddFrame(orderParamLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 0, 0, 0));

	orderNumberEntry = new TGNumberEntry(horzFrameTopParam2, 1,1);
	horzFrameTopParam2->AddFrame(orderNumberEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));
	orderNumberEntry->Resize(75,20);
	orderNumberEntry->SetNumber(order);

	TGHorizontalFrame *horzFrameTop= new TGHorizontalFrame(this);
	AddFrame(horzFrameTop, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 2, 2, 2, 2));

	Plot2DDistribButton = new TGTextButton(horzFrameTop, " Plot 2D distrib ");
	Plot2DDistribButton->Connect("Pressed()", "MyMainFrame", this, "Plot2DDistribOfJPhi()");
	horzFrameTop->AddFrame(Plot2DDistribButton, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));

	PlotSkinDepthButton = new TGTextButton(horzFrameTop, " Skin Depth for Al ");
	PlotSkinDepthButton->Connect("Pressed()", "MyMainFrame", this, "PlotSkinDepth()");
	horzFrameTop->AddFrame(PlotSkinDepthButton, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));

	PlotIntVsVarButton = new TGTextButton(horzFrameTop, " Integral as fct of : ");
	PlotIntVsVarButton->Connect("Pressed()", "MyMainFrame", this, "PlotIntVsVar()");
	horzFrameTop->AddFrame(PlotIntVsVarButton, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));

	VarCombo = new TGComboBox(horzFrameTop);
	horzFrameTop->AddFrame(VarCombo, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));
	VarCombo->Resize(75, 20);
	VarCombo->AddEntry("a", 1);
	VarCombo->AddEntry("b", 2);
	VarCombo->AddEntry("d", 3);
	VarCombo->AddEntry("w", 4);
	VarCombo->AddEntry("j", 5);
	VarCombo->AddEntry("I", 6);
	VarCombo->AddEntry("order", 7);
	VarCombo->Select(1);

	extractCanvasButton = new TGTextButton(horzFrameTop, " Extract canvas ");
	extractCanvasButton->Connect("Pressed()", "MyMainFrame", this, "ExtractCanvas()");
	horzFrameTop->AddFrame(extractCanvasButton, new TGLayoutHints(kLHintsRight | kLHintsExpandY, 2, 0, 0, 0));

	TGHorizontalFrame *horzFrameTopIntRange= new TGHorizontalFrame(this);
	AddFrame(horzFrameTopIntRange, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 2, 2, 2, 2));

	TGLabel *IntRminLabel = new TGLabel(horzFrameTopIntRange, "Integral range :    Min var value : ");
	horzFrameTopIntRange->AddFrame(IntRminLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 0, 0, 0));

	var_rangeMinEntry = new TGNumberEntry(horzFrameTopIntRange, 1,1);
	horzFrameTopIntRange->AddFrame(var_rangeMinEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));
	var_rangeMinEntry->Resize(75,20);
	var_rangeMinEntry->SetNumber(var_rangeMin);
	var_rangeMinEntry->Connect("ValueSet(Double_t)", "MyMainFrame", this, "UpdateIntRange()");

	TGLabel *IntRmaxLabel = new TGLabel(horzFrameTopIntRange, "Max var value : ");
	horzFrameTopIntRange->AddFrame(IntRmaxLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 0, 0, 0));

	var_rangeMaxEntry = new TGNumberEntry(horzFrameTopIntRange, 1,1);
	horzFrameTopIntRange->AddFrame(var_rangeMaxEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));
	var_rangeMaxEntry->Resize(75,20);
	var_rangeMaxEntry->SetNumber(var_rangeMax);
	var_rangeMaxEntry->Connect("ValueSet(Double_t)", "MyMainFrame", this, "UpdateIntRange()");

	TGLabel *IntRLabel = new TGLabel(horzFrameTopIntRange, "Number of points : ");
	horzFrameTopIntRange->AddFrame(IntRLabel, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 2, 0, 0, 0));

	var_rangeEntry = new TGNumberEntry(horzFrameTopIntRange, 1,1);
	horzFrameTopIntRange->AddFrame(var_rangeEntry, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 0, 0));
	var_rangeEntry->Resize(75,20);
	var_rangeEntry->SetNumber(var_range);
	var_rangeEntry->Connect("ValueSet(Double_t)", "MyMainFrame", this, "UpdateIntRange()");

	//========================================= Tab and megaShutter field
	TGHorizontalFrame *shutterFrame= new TGHorizontalFrame(this);
	AddFrame(shutterFrame, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 2, 2, 2, 0));

	//========================================= Info Output
	TGHorizontalFrame *horzFrameCom= new TGHorizontalFrame(this, 10,150);
	AddFrame(horzFrameCom, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 2, 2, 5, 0));

	TextOutputFrame = new TGTextView(horzFrameCom, 50, 150);
	horzFrameCom->AddFrame(TextOutputFrame, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 2, 2, 0, 0));
	TextOutputFrame->SetBackground(0xF0F0F0);
	TextOutputFrame->AddLine("No file loaded");
	TextOutputFrame->Update();

	//========================================= Embedded Canvas
	embeddedCanvas = new TRootEmbeddedCanvas("embeddedCanvas",this,800,400);
	AddFrame(embeddedCanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 0));


	//========================================= Status Bar
	StatusBar = new TGStatusBar(this,20,50);
	AddFrame(StatusBar, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 2, 2, 0, 0));
	int parts [7] = { 20, 20, 6, 20, 6, 12, 15 };
	StatusBar->SetParts(parts,7);
	StatusBar->Draw3DCorner(kFALSE);


	//========================================= Progress Bar
	ProgressBar = new TGHProgressBar(this,TGProgressBar::kFancy | TGProgressBar::kSolidFill,20);
	AddFrame(ProgressBar, new TGLayoutHints(kLHintsExpandX | kLHintsTop, 2, 2, 0, 0));
	ProgressBar->SetBarColor("lightblue");
	ProgressBar->ShowPos(kTRUE);
	ProgressBar->SetMin(0);
	ProgressBar->SetMax(1);
	ProgressBar->SetPosition(1);


	//========================================= Exit Button
	TGTextButton *exitButton = new TGTextButton(this, " Exit ");
	exitButton->Connect("Pressed()", "MyMainFrame", this, "DoExit()");
	AddFrame(exitButton, new TGLayoutHints(kLHintsExpandX | kLHintsBottom, 2, 2, 2, 2));


	// Linking the canvas
	canvas = embeddedCanvas->GetCanvas();
	canvas->cd();

	// Loading the logo
	TImage *img = TImage::Open("nEDM.jpg");
	img->Draw("");

	// Set a name to the main frame
	SetWindowName("Eddy currents calculator");
	MapSubwindows();

	// Initialize the layout algorithm via Resize()
	Resize(GetDefaultSize());

	// Set minimum width, height
	this->SetWMSizeHints(400, 400, 1500, 1000, 0, 0);
	this->MapRaised();

	// Map main frame
	MapWindow();
}

void EddyCurrents() {
  // Popup the GUI...
  new MyMainFrame(gClient->GetRoot(), 200, 200);
}

int main(int argc, char **argv) {
	//gStyle->SetPalette(100);
	TApplication theApp("App",&argc,argv);
	EddyCurrents();
	theApp.Run();
	return 0;
}
