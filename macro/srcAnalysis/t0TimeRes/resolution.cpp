#include <iostream>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"

#include "BmnTrigWaveDigit.h"
#include "BmnTrigDigit.h"

using namespace std;


int doProj( TH2D * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE);
double wlk( double *x , double *p);


int main(int argc, char ** argv)
{

	if (argc !=2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tresolution /path/to/skimmed/res/file\n";
		return -1;
	}


	// Now with our quality parameters, we can go ahead and load the digi file, and cut on carbon-in/carbon-out:
	TString skimFile = argv[1];
	TFile * inFile = new TFile(skimFile);
	if( inFile->IsZombie() ){
		cerr << "Could not open file " << skimFile << "\n";
		return -3;
	}
	TTree * inTree = (TTree*) inFile->Get("events");

	double tdiff_t0_t03, tdiff_t0_mcp2, tdiff_t0_mcp3, tdiff_t0_t0Tqdc;
	double t03_width, mcp2_width, mcp3_width, t0Tqdc_width;

	inTree->SetBranchAddress("t0_t03"	,&tdiff_t0_t03		);
	inTree->SetBranchAddress("t0_mcp2"	,&tdiff_t0_mcp2		);
	inTree->SetBranchAddress("t0_mcp3"	,&tdiff_t0_mcp3		);
	inTree->SetBranchAddress("t0_t0Tqdc"	,&tdiff_t0_t0Tqdc	);

	inTree->SetBranchAddress("t03_wi"	,&t03_width		);
	inTree->SetBranchAddress("mcp2_wi"	,&mcp2_width		);
	inTree->SetBranchAddress("mcp3_wi"	,&mcp3_width		);
	inTree->SetBranchAddress("t0Tqdc_wi"	,&t0Tqdc_width		);


	TH2D * hT03	= new TH2D("hT03"	,"T0 Time - T03 Time vs T03 Width"	,1000,0,20,2500,-50,50);
	TH2D * hMCP3	= new TH2D("hMCP3"	,"T0 Time - MCP3 Time vs MCP3 Width"	,1000,0,20,2500,-50,50);
	TH2D * hMCP2	= new TH2D("hMCP2"	,"T0 Time - MCP2 Time vs MCP2 Width"	,1000,0,20,2500,-50,50);

	for( int i = 0; i < inTree->GetEntries() ; i++){
		if( i % 1000 == 0) cout << "\tPutting entry into histogram: " << i << "\n";
		inTree->GetEntry(i);
		hT03	->Fill( t03_width	,	tdiff_t0_t03	);
		hMCP3	->Fill(	mcp3_width	,	tdiff_t0_mcp3	);
		hMCP2	->Fill(	mcp2_width	,	tdiff_t0_mcp2	);
	}
	
	//TApplication theApp("App",&argc,argv);
	TFile * outFile = new TFile("projections.root","RECREATE");
	std::vector<double> t03_widths 		,mcp3_widths		,mcp2_widths	;
	std::vector<double> t03_times		,mcp3_times		,mcp2_times	;
	std::vector<double> t03_timesErr	,mcp3_timesErr		,mcp2_timesErr	;
	std::vector<double> t03_res		,mcp3_res		,mcp2_res	;
	std::vector<double> t03_resErr		,mcp3_resErr		,mcp2_resErr	;
	int currBin = 0;
	cout << "Working on T03 Histogram\n";
	while( currBin < 1000){
		double xPt, yPt, yEr, ySig, ySigEr;
		int step = doProj( hT03 , currBin , true , 0, xPt, yPt, yEr, ySig, ySigEr );
		currBin += step ;
		if( xPt < 12) continue;
		t03_widths	.push_back(xPt);
		t03_times	.push_back(yPt*-1);
		t03_timesErr	.push_back(yEr);
		t03_res		.push_back(ySig);
		t03_resErr	.push_back(ySigEr);
	}
	cout << "Working on MCP3 Histogram\n";
	currBin = 0;
	while( currBin < 1000){
		double xPt, yPt, yEr, ySig, ySigEr;
		int step = doProj( hMCP3 , currBin , true , 1, xPt, yPt, yEr, ySig, ySigEr );
		currBin += step ;
		if( xPt < 13) continue;
		mcp3_widths	.push_back(xPt);
		mcp3_times	.push_back(yPt*-1);
		mcp3_timesErr	.push_back(yEr);
		mcp3_res	.push_back(ySig);
		mcp3_resErr	.push_back(ySigEr);
	}
	cout << "Working on MCP2 Histogram\n";
	currBin = 0;
	while( currBin < 1000){
		double xPt, yPt, yEr, ySig, ySigEr;
		int step = doProj( hMCP2 , currBin , true , 2, xPt, yPt, yEr, ySig, ySigEr );
		currBin += step ;
		if( xPt < 12) continue;
		mcp2_widths	.push_back(xPt);
		mcp2_times	.push_back(yPt*-1);
		mcp2_timesErr	.push_back(yEr);
		mcp2_res	.push_back(ySig);
		mcp2_resErr	.push_back(ySigEr);
	}
	outFile->Close();
	
	TApplication theApp("App",&argc,argv);
	
	TCanvas *c1 = new TCanvas("c1");
	gStyle->SetOptFit(1);
        gStyle->SetStatX(0.8);
        gStyle->SetStatY(0.9);
        gStyle->SetStatW(0.2);
        gStyle->SetStatH(0.2);
	c1->Divide(3,1);
	
	c1->cd(1);
	int dim = mcp2_widths.size();
	TGraph *g = new TGraph(dim, &mcp2_widths[0], &mcp2_times[0]);
	TF1 * model1 = new TF1("timeWalk",wlk,0,20,2);
	TFitResultPtr one = g->Fit(model1,"QES");
	g->SetTitle("Timewalk of MCP2");
	g->GetXaxis()->SetTitle("MCP2 Width [ns]");
	g->GetYaxis()->SetTitleOffset(1.6);
	g->GetYaxis()->SetTitle("MCP2 Time - T0 Time [ns]");
	g->SetMarkerStyle(20);
	g->Draw("AP");
	c1->Update();


	c1->cd(2);
	dim = mcp3_widths.size();
	TGraph *m = new TGraph(dim, &mcp3_widths[0], &mcp3_times[0]);
	TF1 * model2 = new TF1("timeWalk",wlk,0,20,2);
	TFitResultPtr two = m->Fit(model2,"QES");
	m->SetTitle("Timewalk of MCP3");
	m->GetXaxis()->SetTitle("MCP3 Width [ns]");
	m->GetYaxis()->SetTitle("MCP3 Time - T0 Time [ns]");
	m->GetYaxis()->SetTitleOffset(1.6);
	m->SetMarkerStyle(20);
	m->Draw("AP");
	c1->Update();
	
	c1->cd(3);
	dim = t03_widths.size();
	TGraph *n = new TGraph(dim, &t03_widths[0], &t03_times[0]);
	TF1 * model3 = new TF1("timeWalk",wlk,0,20,2);
	TFitResultPtr three = n->Fit(model3,"QES");
	n->SetTitle("Timewalk of BMN-T0");
	n->GetXaxis()->SetTitle("BMN-T0 Width [ns]");
	n->GetYaxis()->SetTitle("BMN-T0 Time - T0 Time [ns]");
	n->SetMarkerStyle(20);
	n->Draw("AP");
	c1->Update();

	
	// Now with our fits, we can go back through and re-produce time-walk-corrected histograms:
	TH2D * hT03_fix	= new TH2D("hT03_fix"	,"T0 Time - T03 Time vs T03 Width"	,1000,0,20,2500,-50,50);
	for( int i = 0; i < inTree->GetEntries() ; i++){
		if( i % 1000 == 0) cout << "\tTime-walk correcting histograms: " << i << "\n";
		inTree->GetEntry(i);
		hT03_fix->Fill( t03_width	,	three->Parameter(0) + three->Parameter(1) / sqrt(t03_width)  + tdiff_t0_t03 	);
	}
	

	TCanvas * c2 = new TCanvas("c2");
	hT03_fix->Draw("colz");
	c2->Update();

	theApp.Run();
	////////////////////////////////////////////////
	// Plotting for visualization
	/*
	TApplication theApp("App",&argc,argv);

	TCanvas *c1 = new TCanvas("c1");
	hT03->Draw("colz");
	c1->Update();

	TCanvas *c2 = new TCanvas("c2");
	hMCP3->Draw("colz");
	c2->Update();

	TCanvas *c3 = new TCanvas("c3");
	hMCP2->Draw("colz");
	c3->Update();

	theApp.Run();
	*/
	////////////////////////////////////////////////


	inFile->Close();
	return 0;
}

int doProj( TH2D * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE){
	int cnt = 0;
	int step = 0;
	while ( cnt < 500 ){
		char temp[100];
		sprintf(temp,"slice_%d_%d",flag,bin);
		TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
		cnt = (int) pj->Integral();
		if( (bin+step) >= 1000) break;
		step+=1;
	}
	char temp[100];
	sprintf(temp,"slice_%d_%d",flag,bin);
	TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
	
	// Getting the mean x value in this range:
	hist->GetXaxis()->SetRange(bin,bin+step);
	x = hist->GetMean(1);
	hist->GetXaxis()->SetRange();
	TFitResultPtr f = pj->Fit("gaus","QES","",-50,50);
	y = f->Parameter(1);
	yE = f->ParError(1);
	sig = f->Parameter(2);
	sigE = f->ParError(2);
	
	if( write ){
		pj->Write();
	}
	return step;
}

double wlk( double *x , double *p){
	double var = *x;
	return p[0] + p[1] / sqrt(var);
}
