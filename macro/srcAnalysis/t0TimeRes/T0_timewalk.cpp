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
#include "Rtypes.h"

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
void timeRes(	std::vector<double> widths		,
		std::vector<double> times		,
		std::vector<double> timesErr		,
		std::vector<double> res			,
		std::vector<double> resErr		,
		TH2D * hist				,
		TCanvas * c				,
		int cd					,
		TString name				,
		TString color				);


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
			<< "\tT0_timewalk /path/to/skimmed/res/file\n";
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
	double tdiff_t0_t03, tdiff_t0_mcp2, tdiff_t0_mcp3;
	double t03_width, mcp2_width, mcp3_width, t0_width;

	inTree->SetBranchAddress("t0_t03"	,&tdiff_t0_t03		);
	inTree->SetBranchAddress("t0_mcp2"	,&tdiff_t0_mcp2		);
	inTree->SetBranchAddress("t0_mcp3"	,&tdiff_t0_mcp3		);
	inTree->SetBranchAddress("t03_wi"	,&t03_width		);
	inTree->SetBranchAddress("mcp2_wi"	,&mcp2_width		);
	inTree->SetBranchAddress("mcp3_wi"	,&mcp3_width		);
	inTree->SetBranchAddress("t0_wi"	,&t0_width		);

	TH2D * hT03	= new TH2D("hT03"	,"T0 Time - T03 Time vs T0 Amplitude"	,2000,0,40,2500,-50,50);
	TH2D * hMCP3	= new TH2D("hMCP3"	,"T0 Time - MCP3 Time vs T0 Amplitude"	,2000,0,40,2500,-50,50);
	TH2D * hMCP2	= new TH2D("hMCP2"	,"T0 Time - MCP2 Time vs T0 Amplitude"	,2000,0,40,2500,-50,50);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Loop through input skimmed tree
	for( int i = 0; i < inTree->GetEntries() ; i++){
		if( i % 1000 == 0) cout << "\tPutting entry into histogram: " << i << "\n";
		inTree->GetEntry(i);
		if( t03_width > 13 && t03_width < 19 	 )
			hT03	->Fill( t0_width	,	tdiff_t0_t03	);
		if( mcp3_width > 13 && mcp3_width < 16	 )
			hMCP3	->Fill(	t0_width	,	tdiff_t0_mcp3	);
		if( mcp2_width > 14 && mcp2_width < 16.75)
			hMCP2	->Fill(	t0_width	,	tdiff_t0_mcp2	);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Go through and do time-walk correction once
	TFile * outFile = new TFile("projections.root","RECREATE");
	std::vector<double> t03_widths 		,mcp3_widths		,mcp2_widths	;
	std::vector<double> t03_widthsErr	,mcp3_widthsErr		,mcp2_widthsErr	;
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
			24		,
			hT03		);

	cout << "Working on MCP3 Histogram\n";
	walkCorr(	&mcp3_widths	,
			&mcp3_widthsErr	,
			&mcp3_times	,
			&mcp3_timesErr	,
			&mcp3_res	,
			&mcp3_resErr	,
			24		,
			hMCP3		);
	
	cout << "Working on MCP2 Histogram\n";
	walkCorr(	&mcp2_widths	,
			&mcp2_widthsErr	,
			&mcp2_times	,
			&mcp2_timesErr	,
			&mcp2_res	,
			&mcp2_resErr	,
			24		,
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
	TH2D * hT03_fix	= new TH2D("hT03_fix"	,"Corrected: T0 Time - T03 Time vs T0 Amplitude"	,2000,0,40,1000,-5,5);
	TH2D * hMCP3_fix= new TH2D("hMCP3_fix"	,"Corrected: T0 Time - MCP3 Time vs T0 Amplitude"	,2000,0,40,1000,-5,5);
	TH2D * hMCP2_fix= new TH2D("hMCP2_fix"	,"Corrected: T0 Time - MCP2 Time vs T0 Amplitude"	,2000,0,40,1000,-5,5);
	for( int i = 0; i < inTree->GetEntries() ; i++){
		if( i % 1000 == 0) cout << "\tTime-walk correcting histograms: " << i << "\n";
		inTree->GetEntry(i);
		if( mcp2_width > 14 && mcp2_width < 16.75)
			hMCP2_fix->Fill(t0_width	,	tdiff_t0_mcp2 - (mcp2_par0 + mcp2_par1 / sqrt(t0_width) )	);
		if( mcp3_width > 13 && mcp3_width < 16	 )
			hMCP3_fix->Fill(t0_width	,	tdiff_t0_mcp3 - (mcp3_par0 + mcp3_par1 / sqrt(t0_width) )	);
		if( t03_width > 13 && t03_width < 19 	 )
			hT03_fix->Fill( t0_width	,	tdiff_t0_t03  - (t03_par0  + t03_par1  / sqrt(t0_width) )	);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Now I want to take the fixed histograms, do slices along x, and extract the time-resolution
	
	TCanvas *c10 = new TCanvas("c10");
	c10->Divide(4,1);
	
	cout << "Getting Time Resolution for T0 - MCP2\n";
	timeRes(	mcp2_widths	,
			mcp2_times	,
			mcp2_timesErr	,
			mcp2_res	,
			mcp2_resErr	,
			hMCP2_fix	,
			c10		,
			1		,
			"T0-MCP2"	,
			"blue"		);

	cout << "Getting Time Resolution for T0 - MCP3\n";
	timeRes(	mcp3_widths	,
			mcp3_times	,
			mcp3_timesErr	,
			mcp3_res	,
			mcp3_resErr	,
			hMCP3_fix	,
			c10		,
			2		,
			"T0-MCP3"	,
			"blue"		);


	cout << "Getting Time Resolution for T0 - T03\n";
	timeRes(	t03_widths	,
			t03_times	,
			t03_timesErr	,
			t03_res		,
			t03_resErr	,
			hT03_fix	,
			c10		,
			3		,
			"T0-T03"	,
			"blue"		);
	
	TH2D * last = (TH2D*)inFile->Get("mcp3_t03_wiT0");
	std::vector<double> last_widths, last_times, last_timesErr, last_res, last_resErr;
	cout << "Getting Time Resolution for MCP2 - T03\n";
	timeRes(	last_widths	,
			last_times	,
			last_timesErr	,
			last_res	,
			last_resErr	,
			last		,
			c10		,
			4		,
			"MCP3-T03"	,
			"blue");
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Read in the original file to show the improvement
	TFile * inDigi = new TFile("resolution-noCuts.root");
	TTree * inDigT = (TTree *)inDigi->Get("events");
	
	double t_t0_mcp2, 	w_mcp2;
	double t_t0_mcp3, 	w_mcp3;
	double t_t0_t03 , 	w_t03;
	double t_t0_t0Ana, 	w_t0Ana;
	double t_t0_t0Tqdc, 	w_t0Tqdc;	
	double w_t0;
	inDigT->SetBranchAddress("t0_mcp2"	,&t_t0_mcp2	);	
	inDigT->SetBranchAddress("mcp2_wi"	,&w_mcp2	);	
	inDigT->SetBranchAddress("t0_mcp3"	,&t_t0_mcp3	);	
	inDigT->SetBranchAddress("mcp3_wi"	,&w_mcp3	);	
	inDigT->SetBranchAddress("t0_t03"	,&t_t0_t03	);	
	inDigT->SetBranchAddress("t03_wi"	,&w_t03		);	
	inDigT->SetBranchAddress("t0_t0Ana"	,&t_t0_t0Ana	);	
	inDigT->SetBranchAddress("t0Ana_wi"	,&w_t0Ana	);	
	inDigT->SetBranchAddress("t0_t0Tqdc"	,&t_t0_t0Tqdc	);	
	inDigT->SetBranchAddress("t0Tqdc_wi"	,&w_t0Tqdc	);		
	inDigT->SetBranchAddress("t0_wi"	,&w_t0		);
	
	TH2D * hT03_raw	= new TH2D("hT03_raw"	,"Uncorrected: T0 Time - T03 Time vs T0 Amplitude"	,2000,0,40,10000,-50,50);
	TH2D * hMCP3_raw= new TH2D("hMCP3_raw"	,"Uncorrected: T0 Time - MCP3 Time vs T0 Amplitude"	,2000,0,40,10000,-50,50);
	TH2D * hMCP2_raw= new TH2D("hMCP2_raw"	,"Uncorrected: T0 Time - MCP2 Time vs T0 Amplitude"	,2000,0,40,10000,-50,50);
	TH2D * hLAST_raw= new TH2D("hLAST_raw"	,"Uncorrected: MCP3 Time - T03 Time vs T0 Amplitude"	,2000,0,40,10000,-50,50);
	for( int i = 0 ; i < inDigT->GetEntries() ; i ++){
		inDigT->GetEntry(i);
		if( w_mcp2 > 14 && w_mcp2 < 16.75)
			hMCP2_raw->Fill(w_t0,	t_t0_mcp2);
		if( w_mcp3 > 13 && w_mcp3 < 16	 )
			hMCP3_raw->Fill(w_t0,	t_t0_mcp3);
		if( w_t03 > 13 && w_t03 < 19 	 )
			hT03_raw ->Fill(w_t0,	t_t0_t03 );
		if( (w_t03 > 13) && (w_t03 < 19) && (w_mcp3 > 13) && (w_mcp3 < 16) )
			hLAST_raw->Fill(w_t0,   t_t0_t03 - t_t0_mcp3);
	}
	
	std::vector<double> t03_w 		,mcp3_w		,mcp2_w		,last_w;
	std::vector<double> t03_t		,mcp3_t		,mcp2_t		,last_t;
	std::vector<double> t03_tErr		,mcp3_tErr	,mcp2_tErr	,last_tErr;
	std::vector<double> t03_r		,mcp3_r		,mcp2_r		,last_r;
	std::vector<double> t03_rErr		,mcp3_rErr	,mcp2_rErr	,last_rErr;
	
	cout << "Working on raw for T0-MCP2\n";
	timeRes(	mcp2_w		,
			mcp2_t		,
			mcp2_tErr	,
			mcp2_r		,
			mcp2_rErr	,
			hMCP2_raw	,
			c10		,
			1		,
			"T0-MCP2"	,
			"red"		);
	cout << "Working on raw for T0-MCP3\n";
	timeRes(	mcp3_w		,
			mcp3_t		,
			mcp3_tErr	,
			mcp3_r		,
			mcp3_rErr	,
			hMCP3_raw	,
			c10		,
			2		,
			"T0-MCP3"	,
			"red"		);
	cout << "Working on raw for T0-T03\n";
	timeRes(	t03_w		,
			t03_t		,
			t03_tErr	,
			t03_r		,
			t03_rErr	,
			hT03_raw	,
			c10		,
			3		,
			"T0-T03"	,
			"red"		);
	cout << "Working on raw for MCP3-T03\n";
	timeRes(	last_w		,
			last_t		,
			last_tErr	,
			last_r		,
			last_rErr	,
			hLAST_raw	,
			c10		,
			4		,
			"MCP3-T03"	,
			"red"		);
	
	TCanvas *c4 = new TCanvas("Before Any Time-Walk Corrections");
	c4->Divide(3,1);
	c4->cd(1);
	hMCP2_raw->Draw("colz");
	hMCP2_raw->GetXaxis()->SetTitle("T0 Amplitude [a.u.]");
	hMCP2_raw->GetYaxis()->SetTitle("Residual [ns]");
	hMCP2_raw->GetYaxis()->SetTitleOffset(1.5);
	hMCP2_raw->SetStats(0);
	c4->Update();

	c4->cd(2);
	hMCP3_raw->Draw("colz");
	hMCP3_raw->GetXaxis()->SetTitle("T0 Amplitude [a.u.]");
	hMCP3_raw->GetYaxis()->SetTitle("Residual [ns]");
	hMCP3_raw->GetYaxis()->SetTitleOffset(1.5);
	hMCP3_raw->SetStats(0);
	c4->Update();
	
	c4->cd(3);
	hT03_raw->Draw("colz");
	hT03_raw->GetXaxis()->SetTitle("T0 Amplitude [a.u.]");
	hT03_raw->GetYaxis()->SetTitle("Residual [ns]");
	hT03_raw->GetYaxis()->SetTitleOffset(1.5);
	hT03_raw->SetStats(0);
	c4->Update();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Plotting everything for visualization
	TCanvas *c2 = new TCanvas("Before Time-Walk Corrections");
	c2->Divide(3,1);
	c2->cd(1);
	hMCP2->Draw("colz");
	hMCP2->GetXaxis()->SetTitle("T0 Amplitude [a.u.]");
	hMCP2->GetYaxis()->SetTitle("MCP2 Time - T0 Time [ns]");
	hMCP2->GetYaxis()->SetTitleOffset(1.5);
	hMCP2->SetStats(0);
	c2->Update();

	c2->cd(2);
	hMCP3->Draw("colz");
	hMCP3->GetXaxis()->SetTitle("T0 Amplitude [a.u.]");
	hMCP3->GetYaxis()->SetTitle("MCP3 Time - T0 Time [ns]");
	hMCP3->GetYaxis()->SetTitleOffset(1.5);
	hMCP3->SetStats(0);
	c2->Update();

	c2->cd(3);
	hT03->Draw("colz");
	hT03->GetXaxis()->SetTitle("T0 Amplitude [a.u.]");
	hT03->GetYaxis()->SetTitle("T03 Time - T0 Time [ns]");
	hT03->GetYaxis()->SetTitleOffset(1.5);
	hT03->SetStats(0);
	c2->Update();
	
	TCanvas *c3 = new TCanvas("After Time-Walk Corrections");
	c3->Divide(3,1);
	c3->cd(1);
	hMCP2_fix->Draw("colz");
	hMCP2_fix->GetXaxis()->SetTitle("T0 Amplitude [a.u.]");
	hMCP2_fix->GetYaxis()->SetTitle("Residual [ns]");
	hMCP2_fix->GetYaxis()->SetTitleOffset(1.5);
	hMCP2_fix->SetStats(0);
	c3->Update();;

	c3->cd(2);
	hMCP3_fix->Draw("colz");
	hMCP3_fix->GetXaxis()->SetTitle("T0 Amplitude [a.u.]");
	hMCP3_fix->GetYaxis()->SetTitle("Residual [ns]");
	hMCP3_fix->GetYaxis()->SetTitleOffset(1.5);
	hMCP3_fix->SetStats(0);
	c3->Update();
	
	c3->cd(3);
	hT03_fix->Draw("colz");
	hT03_fix->GetXaxis()->SetTitle("T0 Amplitude [a.u.]");
	hT03_fix->GetYaxis()->SetTitle("Residual [ns]");
	hT03_fix->GetYaxis()->SetTitleOffset(1.5);
	hT03_fix->SetStats(0);
	c3->Update();


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Save fits to file
	ofstream outfile("/home/segarrae/software/srcJINR/build/bin/mcp2_mcp3_t03-timeWlkParam.txt");
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
	int thres = hist->GetEntries() / 7;
	while ( cnt < thres ){
		char temp[100];
		sprintf(temp,"slice_%d_%d",flag,bin);
		TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
		cnt = (int) pj->Integral();
		delete pj;
		if( (bin+step) >= 2000) break;
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
	while( currBin < 2000){
		double xPt, yPt, yEr, ySig, ySigEr;
		int step = doProj( hist , currBin , true , 0, xPt, yPt, yEr, ySig, ySigEr );
		currBin += step ;
		if( xPt > widthCut) continue;
		widths		->push_back(xPt);
		widthsErr	->push_back(0);
		times		->push_back(yPt);
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
	TF1 * model = new TF1("timeWalk",wlk,0,40,2);
	TFitResultPtr ptr = g->Fit(model,"QES");
	par0 = ptr->Parameter(0);
	par1 = ptr->Parameter(1);

	c->cd(cd);
	g->SetMarkerStyle(20);
	g->GetXaxis()->SetLimits(17,24);
	g->Draw("AP");
	g->SetTitle("Timewalk of T0 Using "+name);
	g->GetXaxis()->SetTitle("T0 Amplitude [a.u.]");
	g->GetYaxis()->SetTitleOffset(1.6);
	g->GetYaxis()->SetTitle(name+" Time - T0 Time [ns]");
	c->Update();

}

void timeRes(	std::vector<double> widths		,
		std::vector<double> times		,
		std::vector<double> timesErr		,
		std::vector<double> res			,
		std::vector<double> resErr		,
		TH2D * hist				,
		TCanvas * c				,
		int cd					,
		TString name				,
		TString color				){
	std::vector<double> wErr;
	widths.clear();
	times.clear();
	timesErr.clear();
	res.clear();
	resErr.clear();
	walkCorr(	&widths		,
			&wErr		,
			&times		,
			&timesErr	,
			&res		,
			&resErr		,
			24		,
			hist		);

	int dim = widths.size();
	TGraphErrors *g = new TGraphErrors(dim, &widths[0], &res[0],&timesErr[0], &resErr[0]);
	c->cd(cd);
	g->SetMarkerStyle(20);
	g->GetYaxis()->SetRangeUser(0.0,0.250);
	g->GetXaxis()->SetLimits(17,23);
	if( color == "red"){
		g->Draw("P");
		g->SetMarkerColor(kRed);
	}
	else{
		g->Draw("AP");
		g->SetMarkerColor(kBlue);
	}
	g->SetTitle("Time Resolution of "+name);
	g->GetXaxis()->SetTitle("T0 Amplitude [a.u.]");
	g->GetYaxis()->SetTitleOffset(1.6);
	c->Update();
}
