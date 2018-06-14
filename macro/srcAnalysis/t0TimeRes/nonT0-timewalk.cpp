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
#include "TGraphErrors.h"
#include "TStyle.h"

#include "BmnTrigWaveDigit.h"
#include "BmnTrigDigit.h"

using namespace std;


int doProj( TH2D * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE);
double wlk( double *x , double *p);
void walkCorr(	std::vector<double> *widths		,
		std::vector<double> *widthsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		double widthCut				,
		TH2D * hist				);
void fitGraph( 	std::vector<double> widths	,
		std::vector<double> times	,
		std::vector<double> widthsErr	,
		std::vector<double> timesErr	,
		double &par0			,
		double &par1			,
		TCanvas *c			,
		int cd				,
		TString name			);


int main(int argc, char ** argv)
{
	gStyle->SetOptFit(1);
        gStyle->SetStatX(0.8);
        gStyle->SetStatY(0.9);
        gStyle->SetStatW(0.2);
        gStyle->SetStatH(0.2);
	
	if (argc !=2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tnonT0-timewalk /path/to/skimmed/res/file\n";
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

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Load tree
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

	TH2D * hT03	= new TH2D("hT03"	,"T0 Time - T03 Time vs T03 Amplitude"	,1000,0,20,2500,-50,50);
	TH2D * hMCP3	= new TH2D("hMCP3"	,"T0 Time - MCP3 Time vs MCP3 Amplitude"	,1000,0,20,2500,-50,50);
	TH2D * hMCP2	= new TH2D("hMCP2"	,"T0 Time - MCP2 Time vs MCP2 Amplitude"	,1000,0,20,2500,-50,50);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Loop through input skimmed tree
	for( int i = 0; i < inTree->GetEntries() ; i++){
		if( i % 1000 == 0) cout << "\tPutting entry into histogram: " << i << "\n";
		inTree->GetEntry(i);
		hT03	->Fill( t03_width	,	tdiff_t0_t03	);
		hMCP3	->Fill(	mcp3_width	,	tdiff_t0_mcp3	);
		hMCP2	->Fill(	mcp2_width	,	tdiff_t0_mcp2	);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Go through and do time-walk correction once
	TFile * outFile = new TFile("projections.root","RECREATE");
	std::vector<double> t03_widths 		,mcp3_widths		,mcp2_widths	;
	std::vector<double> t03_widthsErr 	,mcp3_widthsErr		,mcp2_widthsErr	;
	std::vector<double> t03_times		,mcp3_times		,mcp2_times	;
	std::vector<double> t03_timesErr	,mcp3_timesErr		,mcp2_timesErr	;
	std::vector<double> t03_res		,mcp3_res		,mcp2_res	;
	std::vector<double> t03_resErr		,mcp3_resErr		,mcp2_resErr	;
	
	cout << "Working on T03 Histogram\n";
	walkCorr(	&t03_widths	,
			&t03_widthsErr	,
			&t03_times	,
			&t03_timesErr	,
			&t03_res	,
			&t03_resErr	,
			14		,
			hT03		);

	cout << "Working on MCP3 Histogram\n";
	walkCorr(	&mcp3_widths	,
			&mcp3_widthsErr	,
			&mcp3_times	,
			&mcp3_timesErr	,
			&mcp3_res	,
			&mcp3_resErr	,
			13		,
			hMCP3		);
	
	cout << "Working on MCP2 Histogram\n";
	walkCorr(	&mcp2_widths	,
			&mcp2_widthsErr	,
			&mcp2_times	,
			&mcp2_timesErr	,
			&mcp2_res	,
			&mcp2_resErr	,
			12		,
			hMCP2		);
	
	outFile->Close();
	
	
	TApplication theApp("App",&argc,argv);
	TCanvas *c1 = new TCanvas("First Iteration of Time-Walk");
	c1->Divide(3,1);
	double mcp2_par0	,mcp2_par1;
	double mcp3_par0	,mcp3_par1;
	double t03_par0		,t03_par1;
	fitGraph(mcp2_widths 	,mcp2_times	,mcp2_widthsErr	,mcp2_timesErr	,mcp2_par0	,mcp2_par1	,c1	,1, "MCP2");
	fitGraph(mcp3_widths 	,mcp3_times	,mcp3_widthsErr	,mcp3_timesErr	,mcp3_par0	,mcp3_par1	,c1	,2, "MCP3");
	fitGraph(t03_widths 	,t03_times	,t03_widthsErr	,t03_timesErr	,t03_par0	,t03_par1	,c1	,3, "T03");


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Now we can refill histograms and plot for visualization
	TH2D * hT03_fix	= new TH2D("hT03_fix"	,"Corrected: T0 Time - T03 Time vs T03 Amplitude"	,1000,0,20,250,-5,5);
	TH2D * hMCP3_fix= new TH2D("hMCP3_fix"	,"Corrected: T0 Time - MCP3 Time vs MCP3 Amplitude"	,1000,0,20,250,-5,5);
	TH2D * hMCP2_fix= new TH2D("hMCP2_fix"	,"Corrected: T0 Time - MCP2 Time vs MCP3 Amplitude"	,1000,0,20,250,-5,5);
	for( int i = 0; i < inTree->GetEntries() ; i++){
		if( i % 1000 == 0) cout << "\tTime-walk correcting histograms: " << i << "\n";
		inTree->GetEntry(i);
		hMCP2_fix->Fill(mcp2_width	,	mcp2_par0 + mcp2_par1 / sqrt(mcp2_width) 	+ tdiff_t0_mcp2);
		hMCP3_fix->Fill(mcp3_width	,	mcp3_par0 + mcp3_par1 / sqrt(mcp3_width) 	+ tdiff_t0_mcp3);
		hT03_fix->Fill( t03_width	,	t03_par0  + t03_par1  / sqrt(t03_width)   	+ tdiff_t0_t03);
	}

	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Plotting everything for visualization
	gStyle->SetOptFit(0);
	TCanvas *c2 = new TCanvas("Before Time-Walk Corrections");
	c2->Divide(3,1);
	c2->cd(1);
	hMCP2->Draw("colz");
	hMCP2->GetXaxis()->SetTitle("MCP2 Amplitude [a.u.]");
	hMCP2->GetYaxis()->SetTitle("MCP2 Time - T0 Time [ns]");
	hMCP2->GetYaxis()->SetTitleOffset(1.5);
	hMCP2->SetStats(0);
	c2->Update();

	c2->cd(2);
	hMCP3->Draw("colz");
	hMCP3->GetXaxis()->SetTitle("MCP3 Amplitude [a.u.]");
	hMCP3->GetYaxis()->SetTitle("MCP3 Time - T0 Time [ns]");
	hMCP3->GetYaxis()->SetTitleOffset(1.5);
	hMCP3->SetStats(0);
	c2->Update();

	c2->cd(3);
	hT03->Draw("colz");
	hT03->GetXaxis()->SetTitle("T03 Amplitude [a.u.]");
	hT03->GetYaxis()->SetTitle("T03 Time - T0 Time [ns]");
	hT03->GetYaxis()->SetTitleOffset(1.5);
	hT03->SetStats(0);
	c2->Update();

	TCanvas *c3 = new TCanvas("After Time-Walk Corrections");
	c3->Divide(3,1);
	c3->cd(1);
	hMCP2_fix->Draw("colz");
	hMCP2_fix->GetXaxis()->SetTitle("MCP2 Ampltiude [a.u.]");
	hMCP2_fix->GetYaxis()->SetTitle("Residual [ns]");
	hMCP2_fix->GetYaxis()->SetTitleOffset(1.5);
	hMCP2_fix->SetStats(0);
	c3->Update();;

	c3->cd(2);
	hMCP3_fix->Draw("colz");
	hMCP3_fix->GetXaxis()->SetTitle("MCP3 Amplitude [a.u.]");
	hMCP3_fix->GetYaxis()->SetTitle("Residual [ns]");
	hMCP3_fix->GetYaxis()->SetTitleOffset(1.5);
	hMCP3_fix->SetStats(0);
	c3->Update();
	
	c3->cd(3);
	hT03_fix->Draw("colz");
	hT03_fix->GetXaxis()->SetTitle("T03 Amplitude [a.u.]");
	hT03_fix->GetYaxis()->SetTitle("Residual [ns]");
	hT03_fix->GetYaxis()->SetTitleOffset(1.5);
	hT03_fix->SetStats(0);
	c3->Update();


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Save fits to file
	ofstream outfile("/home/segarrae/software/srcJINR/build/bin/timeWalk/fitFuncs/mcp2_mcp3_t03-timeWlkParam.txt");
	outfile << mcp2_par0 << "\t" << mcp2_par1 << "\n"
		<< mcp3_par0 << "\t" << mcp3_par1 << "\n"
		<<  t03_par0 << "\t" <<  t03_par1 << "\n";
	outfile.close();


	theApp.Run();
	inFile->Close();

	// Cleanup:
	delete c1, c2, c3;
	delete inFile, outFile;
	delete hT03, hMCP3, hMCP2;
	delete hT03_fix, hMCP3_fix, hMCP2_fix;
	delete inTree;

	return 0;
}

int doProj( TH2D * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE){
	int cnt = 0;
	int step = 0;
	TCanvas * trash = new TCanvas("trash");
	while ( cnt < 300 ){
		char temp[100];
		sprintf(temp,"slice_%d_%d",flag,bin);
		TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
		cnt = (int) pj->Integral();
		delete pj;
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

	delete trash;
	delete pj;

	return step;
}

double wlk( double *x , double *p){
	double var = *x;
	return p[0] + p[1] / sqrt(var);
}

void walkCorr(	std::vector<double> *widths		,
		std::vector<double> *widthsErr		,
		std::vector<double> *times		,
		std::vector<double> *timesErr		,
		std::vector<double> *res		,
		std::vector<double> *resErr		,
		double widthCut				,
		TH2D * hist				){
	int currBin = 0;
	while( currBin < 1000){
		double xPt, yPt, yEr, ySig, ySigEr;
		int step = doProj( hist , currBin , false , 0, xPt, yPt, yEr, ySig, ySigEr );
		currBin += step ;
		if( xPt < widthCut) continue;
		widths		->push_back(xPt);
		widthsErr	->push_back(0);
		times		->push_back(yPt*-1);
		timesErr	->push_back(yEr);
		res		->push_back(ySig);
		resErr		->push_back(ySigEr);
	}
}

void fitGraph( 	std::vector<double> widths	,
		std::vector<double> times	,
		std::vector<double> widthsErr	,
		std::vector<double> timesErr	,
		double &par0			,
		double &par1			,
		TCanvas *c			,
		int cd				,
		TString name			){
	
	int dim = widths.size();
	TGraphErrors *g = new TGraphErrors(dim, &widths[0], &times[0], &widthsErr[0], &timesErr[0]);
	TF1 * model = new TF1("timeWalk",wlk,0,20,2);
	TFitResultPtr ptr = g->Fit(model,"QES");
	par0 = ptr->Parameter(0);
	par1 = ptr->Parameter(1);

	c->cd(cd);
	g->SetMarkerStyle(20);
	g->Draw("AP");
	g->SetTitle("Timewalk of "+name);
	g->GetXaxis()->SetTitle(name+" Ampltiude [a.u.]");
	g->GetYaxis()->SetTitleOffset(1.6);
	g->GetYaxis()->SetTitle(name+" Time - T0 Time [ns]");
	c->Update();

}
