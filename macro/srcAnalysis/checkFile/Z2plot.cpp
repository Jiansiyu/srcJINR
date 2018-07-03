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
void findIdx( TClonesArray* data, int &index , double refT);

int main(int argc, char ** argv)
{

	if (argc !=3)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tskimForNonT0 /path/to/digi/file [fit? = 1]\n";
		return -1;
	}

	string line;
	string file = string(argv[1]);
	bool doFit = atoi(argv[2]);
	string run_number = file.substr( file.find(".")-9 , 4  );
	
	// Open up the checkedRun file to grab pedestals and where carbon peak is:
	string home = std::getenv("VMCWORKDIR");
	string openName = home+"/build/bin/qualityCheck/checked-" + run_number + ".txt";
	ifstream params (openName);
	double pedBC1, pedBC2, pedBC3, pedBC4;
	double C12in_mean, C12out_mean, C12in_sig, C12out_sig;
	if( params.is_open() ){
		params >> pedBC1 >> pedBC2 >> pedBC3 >> pedBC4;
		params >> C12in_mean >> C12in_sig >> C12out_mean >> C12out_sig;
		params.close();
	}
	else{
		cerr << "Could not find matching fileQuality.txt file.\n"
			<< "\tExpected: " << openName << "\n";
		return -2;
	}

	// Now with our quality parameters, we can go ahead and load the digi file, and cut on carbon-in/carbon-out:
	TString digiFile = argv[1];
	TFile * inFile = new TFile(digiFile);
	if( inFile->IsZombie() ){
		cerr << "Could not open file " << digiFile << "\n";
		return -3;
	}
	TTree * inTree = (TTree*) inFile->Get("cbmsim");

	// Load the tree
	TClonesArray * bc1Data = new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * bc2Data = new TClonesArray("BmnTrigWaveDigit");		// This is T0 analog signal to TQDC
	TClonesArray * bc3Data = new TClonesArray("BmnTrigWaveDigit");
	TClonesArray * bc4Data = new TClonesArray("BmnTrigWaveDigit");	
	
	TClonesArray * t0Data 		= new TClonesArray("BmnTrigDigit");	// This is the T0 LVDS signal
	
	inTree->SetBranchAddress("TQDC_BC1"	,&bc1Data	);
	inTree->SetBranchAddress("TQDC_T0"	,&bc2Data	);
	inTree->SetBranchAddress("TQDC_BC3"	,&bc3Data	);
	inTree->SetBranchAddress("TQDC_BC4"	,&bc4Data	);
	inTree->SetBranchAddress("T0"		,&t0Data	);



	TH2D * BC3BC4	= new TH2D("BC3BC4","",2000,0,2000,2000,0,2000);
	TH1D * sum	= new TH1D("sum"   ,"",2000,0,2000);
	TH1D * sumScale	= new TH1D("sumScale","",1000,0,80);


	const int nEvents = inTree->GetEntries();
	for( int event = 0 ; event < nEvents ; event++){
		if( event % 1000 == 0) cout << "\tWorking on event: " << event << "\n";
		bc1Data	->	Clear();
		bc2Data	->	Clear();
		bc3Data	->	Clear();
		bc4Data	->	Clear();
		t0Data	->	Clear();

		inTree->GetEvent(event);
		
		// First check that there is a T0 and if so, grab the data
		if( t0Data->GetEntriesFast() != 1) continue;
		BmnTrigDigit * t0Signal = (BmnTrigDigit*) t0Data->At(0);
		double t0Time = t0Signal->GetTime();
		double t0Width = t0Signal->GetAmp();

		// Now cut on Carbon-in-Carbon-out:
		if( bc1Data->GetEntriesFast() && bc2Data->GetEntriesFast() && bc3Data->GetEntriesFast() && bc4Data->GetEntriesFast() ){
			int bc1Idx, bc2Idx, bc3Idx, bc4Idx;
			findIdx(bc1Data, bc1Idx, t0Time);
			findIdx(bc2Data, bc2Idx, t0Time);
			findIdx(bc3Data, bc3Idx, t0Time);
			findIdx(bc4Data, bc4Idx, t0Time);
		
			BmnTrigWaveDigit * signalBC1 = (BmnTrigWaveDigit *) bc1Data->At(bc1Idx);			
			BmnTrigWaveDigit * signalBC2 = (BmnTrigWaveDigit *) bc2Data->At(bc2Idx);			
			BmnTrigWaveDigit * signalBC3 = (BmnTrigWaveDigit *) bc3Data->At(bc3Idx);			
			BmnTrigWaveDigit * signalBC4 = (BmnTrigWaveDigit *) bc4Data->At(bc4Idx);			
			double adcBC1 = signalBC1->GetPeak() - pedBC1;
			double adcBC2 = signalBC2->GetPeak() - pedBC2;
			double adcBC3 = signalBC3->GetPeak() - pedBC3;
			double adcBC4 = signalBC4->GetPeak() - pedBC4;
				
			if( fabs( sqrt( adcBC1 * adcBC2 ) - C12in_mean) > 2*C12in_sig ) continue;
			BC3BC4->Fill(adcBC3,adcBC4);
			sum->Fill( sqrt(adcBC3 * adcBC4) );
			
			// Scale for nominal trigger
			//sumScale->Fill( sqrt(adcBC3 * adcBC4) * 36./8.84208e+02 );
			
			// Scale for beam trigger
			sumScale->Fill( sqrt(adcBC3 * adcBC4) * 36./7.03804e+02 );

			// Scale for or trigger
			//sumScale->Fill( sqrt(adcBC3 * adcBC4) * 36./7.19375e+02 );
		}
	
	


	}
	TFile * outFile = new TFile("Z2Plot.root","RECREATE");
	outFile->cd();
	sumScale->Write();

	TFitResultPtr f1 = sumScale->Fit("gaus","QES","",1,6);
	TFitResultPtr f2 = sumScale->Fit("gaus","QES","",6,13);
	TFitResultPtr f3 = sumScale->Fit("gaus","QES","",13,20);
	TFitResultPtr f4 = sumScale->Fit("gaus","QES","",20,30);
	TFitResultPtr f5 = sumScale->Fit("gaus","QES","",30,43);
	TFitResultPtr f6 = sumScale->Fit("gaus","QES","",45,60);
	TFitResultPtr f7 = sumScale->Fit("gaus","QES","",58,72);
	
	double par[21];
	par[0] = f1->Parameter(0);
	par[1] = f1->Parameter(1);
	par[2] = f1->Parameter(2);
	par[3] = f2->Parameter(0);
	par[4] = f2->Parameter(1);
	par[5] = f2->Parameter(2);
	par[6] = f3->Parameter(0);
	par[7] = f3->Parameter(1);
	par[8] = f3->Parameter(2);
	par[9] = f4->Parameter(0);
	par[10] = f4->Parameter(1);
	par[11] = f4->Parameter(2);
	par[12] = f5->Parameter(0);
	par[13] = f5->Parameter(1);
	par[14] = f5->Parameter(2);
	par[15] = f6->Parameter(0);
	par[16] = f6->Parameter(1);
	par[17] = f6->Parameter(2);
	par[18] = f7->Parameter(0);
	par[19] = f7->Parameter(1);
	par[20] = f7->Parameter(2);

	cout << "Parameters:\n";
	for( int k = 0 ; k < 7 ; k++){
		cout << par[3*k+0] << " " << par[3*k+1] << " " << par[3*k+2] << "\n";
	}
	cout << "Entries: " << sumScale->GetEntries() << "\n";

	TF1 * total = new TF1("total","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)+gaus(18)",0,80);
	total->SetParameters(par);
	TFitResultPtr tot = sumScale->Fit(total,"QES");

	TF1 * g1	= new TF1("g1","gaus",0,80);
	g1->SetLineStyle(7);
	g1->SetLineColor(1);
	g1->SetParameter(0,tot->Parameter(0));
	g1->SetParameter(1,tot->Parameter(1));
	g1->SetParameter(2,tot->Parameter(2));
	TF1 * g2	= new TF1("g2","gaus",0,80);
	g2->SetParameter(0,tot->Parameter(3));
	g2->SetParameter(1,tot->Parameter(4));
	g2->SetParameter(2,tot->Parameter(5));
	g2->SetLineStyle(7);
	g2->SetLineColor(1);
	TF1 * g3	= new TF1("g3","gaus",0,80);
	g3->SetParameter(0,tot->Parameter(6));
	g3->SetParameter(1,tot->Parameter(7));
	g3->SetParameter(2,tot->Parameter(8));
	g3->SetLineStyle(7);
	g3->SetLineColor(1);
	TF1 * g4	= new TF1("g4","gaus",0,80);
	g4->SetParameter(0,tot->Parameter(9));
	g4->SetParameter(1,tot->Parameter(10));
	g4->SetParameter(2,tot->Parameter(11));
	g4->SetLineStyle(7);
	g4->SetLineColor(1);
	TF1 * g5	= new TF1("g5","gaus",0,80);
	g5->SetParameter(0,tot->Parameter(12));
	g5->SetParameter(1,tot->Parameter(13));
	g5->SetParameter(2,tot->Parameter(14));
	g5->SetLineStyle(7);
	g5->SetLineColor(1);
	TF1 * g6	= new TF1("g6","gaus",0,80);
	g6->SetParameter(0,tot->Parameter(15));
	g6->SetParameter(1,tot->Parameter(16));
	g6->SetParameter(2,tot->Parameter(17));
	g6->SetLineStyle(7);
	g6->SetLineColor(1);
	TF1 * g7	= new TF1("g7","gaus",0,80);
	g7->SetParameter(0,tot->Parameter(18));
	g7->SetParameter(1,tot->Parameter(19));
	g7->SetParameter(2,tot->Parameter(20));
	g7->SetLineStyle(7);
	g7->SetLineColor(1);




	TApplication theApp("App",&argc,argv);
	TCanvas * c = new TCanvas("c1");
	BC3BC4->Draw("colz");
	c->Update();

	TCanvas * d = new TCanvas("d1");
	sum->Draw();
	d->Update();


	TCanvas * c2 = new TCanvas("c2");
	sumScale->Draw();
	if( doFit){
		g1->Draw("same");
		g2->Draw("same");
		g3->Draw("same");
		g4->Draw("same");
		g5->Draw("same");
		g6->Draw("same");
		g7->Draw("same");
	}
	c2->Update();

	outFile->cd();
	g1->Write();
	g2->Write();
	g3->Write();
	g4->Write();
	g5->Write();
	g6->Write();
	g7->Write();
	total->Write();
	outFile->Close();

	//inFile->Close();
	
	theApp.Run();
	

	return 0;
}

void findIdx( TClonesArray* data, int &index , double refT){
	double minT = 1e4;
	for( int m = 0 ; m < data->GetEntriesFast() ; m++){
		BmnTrigWaveDigit * signal = (BmnTrigWaveDigit*)data->At(m);
		double time = fabs(signal->GetTime() - refT);
		if( time < minT){
			minT = time;
			index = m;
		}
	}
}
