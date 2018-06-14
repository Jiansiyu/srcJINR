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
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TCut.h"

#include "BmnTrigWaveDigit.h"
#include "BmnTrigDigit.h"

using namespace std;

double pedBC1, pedBC2, pedBC3, pedBC4;
double C12in_mean, C12out_mean, C12in_sig, C12out_sig;

void findIdx( TClonesArray* data, int &index , double refT);
void findIdx_tdc( TClonesArray* data, int &index , double refT);
void findPeak( TTree * inTree, double &peak, double &sig, TString branch);

int main(int argc, char ** argv)
{

	if (argc !=2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tskimForT0 /path/to/digi/file\n";
		return -1;
	}

	string line;
	string file = string(argv[1]);
	string run_number = file.substr( file.find(".")-9 , 4  );
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Open up the checkedRun file to grab pedestals and where carbon peak is:
	string openName = "/home/segarrae/software/srcJINR/build/bin/qualityCheck/checked-" + run_number + ".txt";
	ifstream params (openName);
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
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Open up the time-walk function parameters for MCP2, MCP3, T03
	double mcp2_par0, mcp2_par1, mcp3_par0, mcp3_par1, t03_par0, t03_par1;
	openName = "/home/segarrae/software/srcJINR/build/bin/timeWalk/fitFuncs/mcp2_mcp3_t03-timeWlkParam.txt";
	ifstream wlk (openName);
	if( wlk.is_open() ){
		string junk;
		wlk >> mcp2_par0 >> mcp2_par1;
		wlk >> mcp3_par0 >> mcp3_par1;
		wlk >> t03_par0  >> t03_par1 ;
		wlk.close();
	}
	else{
		cerr << "Could not find walk function file.\n"
			<< "\t Expected " << openName << "\n";
		return -4;
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Setup the input file, load the tree, and setup output file
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
	TClonesArray * t0AnaData	= new TClonesArray("BmnTrigDigit");	// This is the T0 analog signal to TDC
	TClonesArray * mcp2Data 	= new TClonesArray("BmnTrigDigit");	// This is MCP-2 signal
	TClonesArray * mcp3Data 	= new TClonesArray("BmnTrigDigit");	// This is MCP-3 signal (has more events than MCP-2
	TClonesArray * t03Data	 	= new TClonesArray("BmnTrigDigit");	// This is BMN-T0 signal
	
	inTree->SetBranchAddress("TQDC_BC1"	,&bc1Data	);
	inTree->SetBranchAddress("TQDC_T0"	,&bc2Data	);
	inTree->SetBranchAddress("TQDC_BC3"	,&bc3Data	);
	inTree->SetBranchAddress("TQDC_BC4"	,&bc4Data	);
	inTree->SetBranchAddress("T0"		,&t0Data	);
	inTree->SetBranchAddress("T0Analog"	,&t0AnaData	);
	inTree->SetBranchAddress("T03"		,&t03Data	);
	inTree->SetBranchAddress("MCP2"		,&mcp2Data	);
	inTree->SetBranchAddress("MCP3"		,&mcp3Data	);


	// Setup output file
	TFile * outFile = new TFile("sk_T0.root","RECREATE");
	TTree * outTree = new TTree("events","");
	double tdiff_t0_mcp2, mcp2_wi;
	double tdiff_t0_mcp3, mcp3_wi;
	double tdiff_t0_t03 , t03_wi;
	double width_t0;
	outTree->Branch("t0_wi"		,&width_t0		,"t0_wi/D");
	outTree->Branch("t03_wi"	,&t03_wi		,"t03_wi/D");
	outTree->Branch("mcp2_wi"	,&mcp2_wi		,"mcp2_wi/D");
	outTree->Branch("mcp3_wi"	,&mcp3_wi		,"mcp3_wi/D");

	outTree->Branch("t0_mcp2"	,&tdiff_t0_mcp2		,"t0_mcp2/D");
	outTree->Branch("t0_mcp3"	,&tdiff_t0_mcp3		,"t0_mcp3/D");
	outTree->Branch("t0_t03"	,&tdiff_t0_t03		,"t0_t03/D");

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Before hand, let's calculate the peak of MCP2, MCP3, T03 in order to sharply cut on them when
	// looking at T0 - (other det)
	double mcp2Wid_peak, mcp2Wid_sig;
	double mcp3Wid_peak, mcp3Wid_sig;
	double t03Wid_peak, t03Wid_sig;
	findPeak( inTree, mcp2Wid_peak, mcp2Wid_sig, "MCP2[0]");
	findPeak( inTree, mcp3Wid_peak, mcp3Wid_sig, "MCP3[0]");
	findPeak( inTree, t03Wid_peak,  t03Wid_sig,  "T03[0]");
	

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Now loop through the digi file. First we cut on 12-C in - out and then we look at T0 and other detectors.
	// Since we time-walk corrected MCP2, MCP3, T03, we are interested in how un-cut T0 looks like on these other
	// detectors. One question is do I cut on the other widths, or not?? I've already time-walk corrected them
	// but I can see how it looks if I do/dont...
	const int nEvents = inTree->GetEntries();
	double cnt = 0.;
	double cnt_t0andMCP2 = 0.;
	double cnt_t0andMCP3 = 0.;
	double cnt_t0andt03 = 0.;
	double cnt_t0andt0Ana = 0.;
	double cnt_t0andTqdc = 0.;
	TH2D * mcp2_t03_wiMCP = new TH2D("mcp2_t03_wiMCP","Corrected Time of (MCP2 - T03) vs MCP2 Width",2000,0,40,1000,-5,5);
	TH2D * mcp2_t03_wiT03 = new TH2D("mcp2_t03_wiT03","Corrected Time of (MCP2 - T03) vs T03 Width",2000,0,40,1000,-5,5);
	TH2D * mcp3_t03_wiMCP = new TH2D("mcp3_t03_wiMCP","Corrected Time of (MCP3 - T03) vs MCP3 Width",2000,0,40,1000,-5,5);
	TH2D * mcp3_t03_wiT03 = new TH2D("mcp3_t03_wiT03","Corrected Time of (MCP3 - T03) vs T03 Width",2000,0,40,1000,-5,5);
	TH2D * mcp3_t03_wiT0 = new TH2D("mcp3_t03_wiT0","Corrected Time of (MCP3 - T03) vs T0 Width",2000,0,40,1000,-5,5);
	for( int event = 0 ; event < nEvents ; event++){
		width_t0 = -666;
		tdiff_t0_mcp2	=-666;
		tdiff_t0_mcp3	=-666;
		tdiff_t0_t03	=-666;
		mcp3_wi		=-666;
		mcp2_wi		=-666;
		t03_wi		=-666;

		if( event % 1000 == 0) cout << "\tWorking on event: " << event << "\n";
		bc1Data	->	Clear();
		bc2Data	->	Clear();
		bc3Data	->	Clear();
		bc4Data	->	Clear();
		t0Data	->	Clear();
		t0AnaData->	Clear();
		t03Data	->	Clear();
		mcp2Data->	Clear();
		mcp3Data->	Clear();

		inTree->GetEvent(event);
		
		// First check that there is a T0 and if so, grab the data
		if( t0Data->GetEntriesFast() != 1) continue;
		BmnTrigDigit * t0Signal = (BmnTrigDigit*) t0Data->At(0);
		double t0Time = t0Signal->GetTime();
		width_t0 = t0Signal->GetAmp();

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
				
			//if( fabs( sqrt( adcBC1 * adcBC2 ) - C12in_mean) > 2*C12in_sig ) continue;
			if( fabs( sqrt( adcBC3 * adcBC4 ) - C12out_mean) > 2*C12out_sig ) continue;
		}
	
		//	Now look at MCP2 vs T03 -- no cuts on either widths
		if( mcp2Data->GetEntriesFast() && t03Data->GetEntriesFast() ){
			int mcp2Idx, t03Idx;
			
			findIdx_tdc( mcp2Data, mcp2Idx, t0Time);
			BmnTrigDigit * signal1 = (BmnTrigDigit *) mcp2Data->At(mcp2Idx);
			double width_mcp2 = signal1->GetAmp();
			
			findIdx_tdc( t03Data,  t03Idx,  t0Time);
			BmnTrigDigit * signal2 = (BmnTrigDigit *) t03Data->At(t03Idx);
			double width_t03 = signal2->GetAmp();

			if( width_t03 < 13) continue;
			if( width_t03 > 19) continue;
			if( width_mcp2 < 14) continue;
			if( width_mcp2 > 16.75) continue;
			
			double tdiff_mcp2_t03 = ( (t0Time - signal1->GetTime()) + (mcp2_par0 + mcp2_par1 / sqrt(width_mcp2)) ) - ( (t0Time - signal2->GetTime()) + (t03_par0  + t03_par1  / sqrt(width_t03))   );
			mcp2_t03_wiMCP->Fill( width_mcp2, tdiff_mcp2_t03);
			mcp2_t03_wiT03->Fill( width_t03 , tdiff_mcp2_t03);
		}
		//	Now look at MC3 vs T03 -- no cuts on either widths
		if( mcp3Data->GetEntriesFast() && t03Data->GetEntriesFast() ){
			int mcp3Idx, t03Idx;
	
			findIdx_tdc( mcp3Data, mcp3Idx, t0Time);
			BmnTrigDigit * signal1 = (BmnTrigDigit *) mcp3Data->At(mcp3Idx);
			double width_mcp3 = signal1->GetAmp();
		
			findIdx_tdc( t03Data,  t03Idx,  t0Time);
			BmnTrigDigit * signal2 = (BmnTrigDigit *) t03Data->At(t03Idx);
			double width_t03 = signal2->GetAmp();
			
			if( width_t03 < 13) continue;
			if( width_t03 > 19) continue;
			if( width_mcp3 < 13) continue;
			if( width_mcp3 > 16) continue;

			double tdiff_mcp3_t03 = ( (t0Time - signal1->GetTime()) + (mcp3_par0 + mcp3_par1 / sqrt(width_mcp3)) ) - ( (t0Time - signal2->GetTime()) + (t03_par0  + t03_par1  / sqrt(width_t03))   );
			mcp3_t03_wiMCP->Fill( width_mcp3, tdiff_mcp3_t03);
			mcp3_t03_wiT03->Fill( width_t03 , tdiff_mcp3_t03);
			mcp3_t03_wiT0->Fill( width_t0 , tdiff_mcp3_t03);
		}
	
		// Now I can look at T0 - MCP, T0 - T03 with no cuts on widths at first
		if( mcp2Data->GetEntriesFast()){
			int mcp2Idx;
				
			findIdx_tdc( mcp2Data, mcp2Idx, t0Time);
			BmnTrigDigit * signal = (BmnTrigDigit *) mcp2Data->At(mcp2Idx);
			mcp2_wi = signal->GetAmp();
		
			tdiff_t0_mcp2 = (t0Time - signal->GetTime()) + (mcp2_par0 + mcp2_par1 / sqrt(mcp2_wi));
		}
		if( mcp3Data->GetEntriesFast()){
			int mcp3Idx;
			
			findIdx_tdc( mcp3Data, mcp3Idx, t0Time);
			BmnTrigDigit * signal = (BmnTrigDigit *) mcp3Data->At(mcp3Idx);
			mcp3_wi = signal->GetAmp();
				
			tdiff_t0_mcp3 = (t0Time - signal->GetTime()) + (mcp3_par0 + mcp3_par1 / sqrt(mcp3_wi));
		}
		if( t03Data->GetEntriesFast()){
			int t03Idx;

			findIdx_tdc( t03Data, t03Idx, t0Time);
			BmnTrigDigit * signal = (BmnTrigDigit *) t03Data->At(t03Idx);
			t03_wi = signal->GetAmp();
			
			tdiff_t0_t03 = (t0Time - signal->GetTime()) + (t03_par0 + t03_par1 / sqrt(t03_wi) );
		}

	
		// This has 50/300k events, good!
		// as it should have nothing because these are read out independently
		/*
		if( mcp2Data->GetEntriesFast() && mcp3Data->GetEntriesFast() ){
			int mcp2Idx, mcp3Idx;
			findIdx_tdc( mcp2Data, mcp2Idx, t0Time);
			findIdx_tdc( mcp3Data, mcp3Idx, t0Time);
			BmnTrigDigit * signal1 = (BmnTrigDigit *) mcp2Data->At(mcp2Idx);
			BmnTrigDigit * signal2 = (BmnTrigDigit *) mcp3Data->At(mcp3Idx);
			width_mcp2 = signal1->GetAmp();
			width_mcp3 = signal2->GetAmp();
			tdiff_mcp2_mcp3 = (signal1->GetTime() - signal2->GetTime()) - (mcp2_par0 + mcp2_par1 / sqrt(width_mcp2)) + (mcp3_par0 + mcp3_par1 / sqrt(width_mcp3) );
		}
		*/


		outTree->Fill();	
	}



	outFile->cd();
	mcp2_t03_wiMCP->Write();
	mcp2_t03_wiT03->Write();
	mcp3_t03_wiMCP->Write();
	mcp3_t03_wiT03->Write();
	mcp3_t03_wiT0->Write();
	outTree->Write();
	outFile->Close();

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
void findIdx_tdc( TClonesArray* data, int &index , double refT){
	double minT = 1e4;
	for( int m = 0 ; m < data->GetEntriesFast() ; m++){
		BmnTrigDigit * signal = (BmnTrigDigit*)data->At(m);
		double time = fabs(signal->GetTime() - refT);
		if( time < minT){
			minT = time;
			index = m;
		}
	}
}

void findPeak( TTree * inTree, double &peak, double &sig, TString branch){
	TH1D * h = new TH1D("trash","trash",1000,0,100);
	TString w;
	//w.Form(" fabs (sqrt( (TQDC_BC1[0]->GetPeak() - %f) * (TQDC_T0[0]->GetPeak() - %f ) ) - %f ) < 2 * %f ",pedBC1,pedBC2,C12in_mean,C12in_sig);
	TString z;
	z.Form(" fabs (sqrt( (TQDC_BC3[0]->GetPeak() - %f) * (TQDC_BC4[0]->GetPeak() - %f ) ) - %f ) < 2 * %f ",pedBC3,pedBC4,C12out_mean,C12out_sig);
	
	inTree->Draw(branch+"->GetAmp() >> trash",z);
	TFitResultPtr init = h->Fit("gaus","QES","",0,100);
	double start = init->Parameter(1) - init->Parameter(2);
	TFitResultPtr fin  = h->Fit("gaus","QES","",start , 100.);

	peak = fin->Parameter(1);
	sig  = fin->Parameter(2);
	
	delete h;
}
