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
#include "TVectorT.h"
#include "TRandom3.h"

using namespace std;

int doProj( TH2D * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE);
double wlk( double *x , double *p);
void fitGraph( 	std::vector<double> widths	,
		std::vector<double> times	,
		std::vector<double> widthsErr	,
		std::vector<double> timesErr	,
		double &par0			,
		double &par1			,
		double &par2			,
		double &par0Er			,
		double &par1Er			,
		double &par2Er			,
		TCanvas *c			,
		int cd				);

int main(int argc, char ** argv)
{
	if (argc < 2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\ttofTimeWalk /path/to/histFile\n";
		return -1;
	}

	ofstream f_call;
	TString NameCallFile = "TOF400_TimeWalk_RUN7_SRC.dat";
	f_call.open(NameCallFile.Data());
	f_call << "Plane\tStrip\tPts-For-Fit\tShift-in-Amp\tPar0\tPar1\tPar2\tPar0_Error\tPar1_Error\tPar2_Error" << endl << "==================================================================================================================" << endl;

	TFile * dataFile = new TFile(argv[1]);
	TString hName;

	TFile * outFile = new TFile("test.root","RECREATE");

	for( int plane = 0 ; plane < 20 ; plane++){

		//TApplication theApp("App",&argc,argv);
		hName = Form("Plane %i Fits",plane);
		TCanvas *c1 = new TCanvas("c1",hName);
		c1->Divide(4,2);
		TCanvas *c2 = new TCanvas("c2",hName);
		c2->Divide(4,2);
		TCanvas *c3 = new TCanvas("c3",hName);
		c3->Divide(4,2);
		TCanvas *c4 = new TCanvas("c4",hName);
		c4->Divide(4,2);
		TCanvas *c5 = new TCanvas("c5",hName);
		c5->Divide(4,2);
		TCanvas *c6 = new TCanvas("c6",hName);
		c6->Divide(4,2);

		cout << "Working on plane " << plane << "...\n";
		for( int strip = 0 ; strip < 48 ; strip++){
			cout << "\tWorking on strip " << strip << "...\n";
			hName = Form("hToF_All_2D_%i_%i",plane,strip);
			TH2D * hToF_v_Amp = (TH2D*)dataFile -> Get(hName);
			
			if( (hToF_v_Amp->GetEntries() / 1000.) < 4 ){
				cout << "Not enough entries, plane or strip likely broken. Quitting on " << plane << " " << strip << "...\n";
				f_call << std::setprecision(6);
				f_call << plane << "\t" << strip << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << -1 << "\t" << endl;
				delete hToF_v_Amp;
				continue;
			}
			
			// We need to sweep across the histogram and find bin centers where we can get 5 data points
			std::vector<double> amps;
			std::vector<double> amps_er;
			std::vector<double> meanT;
			std::vector<double> meanT_er;
			std::vector<double> resT;
			std::vector<double> resT_er;
			int currBin = 0;
			while( currBin < hToF_v_Amp->GetNbinsX()){
				cout << "\t\tOn bin " << currBin << "\n";
				double xPt, yPt, yEr, ySig, ySigEr;
				int step = doProj( hToF_v_Amp , currBin , true , strip + 48*plane, xPt, yPt, yEr, ySig, ySigEr );
				currBin += step;
				if( ySigEr != -1){
					amps.push_back( xPt );
					amps_er.push_back( 0 );
					meanT.push_back( yPt );
					meanT_er.push_back( yEr );
					resT.push_back( ySig );
					resT_er.push_back( ySigEr );
				}
				cout << "\t\t with results: " << xPt << " " << yPt << " " << yEr << " " << ySig << " " << ySigEr << "\n";
			}	
		
			double par0, par1, par2, par0Er, par1Er, par2Er;
			if( strip < 8)
				fitGraph( amps, meanT, amps_er, meanT_er, par0, par1, par2, par0Er, par1Er, par2Er, c1, strip );
			else if( strip < 16)
				fitGraph( amps, meanT, amps_er, meanT_er, par0, par1, par2, par0Er, par1Er, par2Er, c2, strip );
			else if( strip < 24)
			 	fitGraph( amps, meanT, amps_er, meanT_er, par0, par1, par2, par0Er, par1Er, par2Er, c3, strip );
			else if( strip < 32)
			 	fitGraph( amps, meanT, amps_er, meanT_er, par0, par1, par2, par0Er, par1Er, par2Er, c4, strip );
			else if( strip < 40)
			 	fitGraph( amps, meanT, amps_er, meanT_er, par0, par1, par2, par0Er, par1Er, par2Er, c5, strip );
			else if( strip < 48)
				fitGraph( amps, meanT, amps_er, meanT_er, par0, par1, par2, par0Er, par1Er, par2Er, c6, strip );

			
			f_call << plane << "\t" << strip << "\t" << amps.size() << "\t" << amps.at(0) << "\t" << par0 << "\t" << par1 << "\t" << par2 << "\t" << par0Er << "\t" << par1Er << "\t" << par2Er << "\t" << endl;
		
		}
		hName = Form("plane_%i_%i_fits.pdf",plane,0);
		c1->Print(hName,"pdf");	
		hName = Form("plane_%i_%i_fits.pdf",plane,1);
		c2->Print(hName,"pdf");	
		hName = Form("plane_%i_%i_fits.pdf",plane,2);
		c3->Print(hName,"pdf");	
		hName = Form("plane_%i_%i_fits.pdf",plane,3);
		c4->Print(hName,"pdf");	
		hName = Form("plane_%i_%i_fits.pdf",plane,4);
		c5->Print(hName,"pdf");	
		hName = Form("plane_%i_%i_fits.pdf",plane,5);
		c6->Print(hName,"pdf");	
	
		delete c1, c2, c3, c4, c5, c6;
	}

	f_call.close();
	dataFile->Close();
	outFile->Close();

}


int doProj( TH2D * hist , int bin , bool write , int flag , double &x, double &y, double &yE,double &sig, double &sigE){
	int cnt = 0;
	int step = 0;

	int thres = 1000;
	int dataPoints = hist->GetEntries() / thres;

	TCanvas * trash = new TCanvas("trash");
	while ( cnt < thres ){
		char temp[100];
		sprintf(temp,"slice_%d_%d",flag,bin);
		TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
		cnt = (int) pj->Integral();
		delete pj;
		if( (bin+step) >= hist->GetNbinsX() ) break;
		step+=1;
	}

	if ( cnt < 500){
		cout << "\t\t\t**skipping bin due to only " << cnt << " events**\n";
		x    = -1;
		y    = -1; 
		yE   = -1; 
		
		sig  = -1; 
		sigE = -1; 
	}
	else{

		char temp[100];
		sprintf(temp,"slice_%d_%d",flag,bin);
		TH1D * pj = hist->ProjectionY(temp,bin,bin+step);
		
		// Getting the mean x value in this range:
		hist->GetXaxis()->SetRange(bin,bin+step);
		x = hist->GetMean(1);
		hist->GetXaxis()->SetRange();

		// Do initial fit for projection
		double testMax = pj->GetXaxis()->GetBinCenter(pj->GetMaximumBin());
		TF1 * testFit = new TF1("testFit","gaus",-1,testMax+0.500);
		pj->Fit("testFit","QESRN");
		double init_par[3];
		testFit->GetParameters(&init_par[0]);
		const double * initError = testFit->GetParErrors();
		

		// Do second fit of projection
		TF1 * finFit = new TF1("finFit","gaus",-1,init_par[1]+0.5*init_par[2]);
		pj->Fit("finFit","QESRN");
		double fin_par[3], fin_parEr[3];
		finFit->GetParameters(&fin_par[0]);
		const double * finErrors = finFit->GetParErrors();
		fin_parEr[0] = finErrors[0];
		fin_parEr[1] = finErrors[1];
		fin_parEr[2] = finErrors[2];

		// Need to check that second fit is better
		int i = 1;
		while( fin_par[2] > init_par[2] ){
			cout << "\t\t\t*** trying to refit ... ***\n";
			finFit->SetRange( -1 , init_par[1] + (0.5+0.05*i)*init_par[2] );
			pj->Fit("finFit","QESRN");
			finFit->GetParameters(&fin_par[0]);
			const double * finEr = finFit->GetParErrors();
			fin_parEr[0] = finEr[0];
			fin_parEr[1] = finEr[1];
			fin_parEr[2] = finEr[2];
			if( (0.5+0.05*i) >= 2){ 
				cout << "\t\t\t***refit FAILED... ***\n";
				fin_par[2] = init_par[2];
				fin_par[1] = init_par[1];
				fin_par[0] = init_par[0];
				fin_parEr[0] = initError[0];
				fin_parEr[1] = initError[1];
				fin_parEr[2] = initError[2];
				break;
			}
			i++;
		}
		if ( fin_par[2] > 0.5 ){ // error with fit in this bin
			cout << "\t\t\t***issue with strip, fitting FAILED... ***\n";
			if( write) pj->Write();
			x    = -1;
			y    = -1; 
			yE   = -1; 
			
			sig  = -1; 
			sigE = -1; 

			delete pj;
			delete trash;
			return step;
		}
		
		y = fin_par[1];
		yE = fin_parEr[1];
		
		sig = fin_par[2];
		sigE = fin_parEr[2];

		if( write )
			pj->Write();
		
		delete pj;
	}

	delete trash;
	
	return step;
}
double wlk( double *x , double *p){
	double var = *x;
	//return p[0] + p[1] / sqrt(var);
	return p[0]+p[1]*exp(-var/p[2]);
}

void fitGraph( 	std::vector<double> widths	,
		std::vector<double> times	,
		std::vector<double> widthsErr	,
		std::vector<double> timesErr	,
		double &par0			,
		double &par1			,
		double &par2			,
		double &par0Er			,
		double &par1Er			,
		double &par2Er			,
		TCanvas *c			,
		int cd				){
	
	int dim = widths.size();
	std::vector<double> shiftedWidths;
	for(int i = 0 ; i < dim ; i++){
		shiftedWidths.push_back( widths.at(i) - widths.at(0) );
	}


	TGraphErrors *g = new TGraphErrors(dim, &shiftedWidths[0], &times[0], &widthsErr[0], &timesErr[0]);
	TF1 * model = new TF1("timeWalk",wlk,6,35,3);
	model->SetParameter(0,times.at(dim-1));
	model->SetParameter(1,times.at(0));
	model->SetParameter(2,3);

	TFitResultPtr ptr = g->Fit(model,"QES");
	par0 = ptr->Parameter(0);
	par0Er = ptr->Error(0);
	par1 = ptr->Parameter(1);
	par1Er = ptr->Error(1);
	par2 = ptr->Parameter(2);
	par2Er = ptr->Error(2);

	int pl = cd/8;
	int id = cd - pl*8;
	c->cd(id+1);
	g->SetMarkerStyle(20);
	gStyle->SetOptFit(1);
	//g->GetXaxis()->SetLimits(17,24);
	g->Draw("AP");
	TString name;
	name = Form("Timewalk of Strip %i" , cd);
	g->SetTitle(name);
	g->GetXaxis()->SetTitle("ToF Strip Amplitude (#sqrt{L*R}-Min) [a.u.]");
	g->GetYaxis()->SetTitleOffset(1.5);
	g->GetYaxis()->SetTitle("Mean Time Peak [ns]");
	c->Update();

}
