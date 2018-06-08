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
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TCut.h"

#include "BmnTrigWaveDigit.h"
#include "BmnTrigDigit.h"

using namespace std;
void findIdx( TClonesArray* data, int &index , double refT);
void findIdx_tdc( TClonesArray* data, int &index , double refT);

int main(int argc, char ** argv)
{

	if (argc !=2)
	{
		cerr << "Wrong number of arguments. Instead use\n"
			<< "\tprepResolution /path/to/digi/file\n";
		return -1;
	}

	string line;
	string file = string(argv[1]);
	string run_number = file.substr( file.find(".")-9 , 4  );
	
	// Open up the checkedRun file to grab pedestals and where carbon peak is:
	string openName = "/home/segarrae/software/srcJINR/build/bin/checked-" + run_number + ".txt";
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
	TFile * outFile = new TFile("resolution.root","RECREATE");
	TTree * outTree = new TTree("events","");
	double tdiff_t0_mcp2, width_mcp2;
	double tdiff_t0_mcp3, width_mcp3;
	double tdiff_t0_t03 , width_t03	;
	double tdiff_t0_t0Ana, width_t0Ana;
	double tdiff_t0_t0Tqdc, width_t0Tqdc;	

	outTree->Branch("t0_mcp2"	,&tdiff_t0_mcp2		,"t0_mcp2/D");
	outTree->Branch("mcp2_wi"	,&width_mcp2		,"mcp2_wi/D");
	
	outTree->Branch("t0_mcp3"	,&tdiff_t0_mcp3		,"t0_mcp3/D");
	outTree->Branch("mcp3_wi"	,&width_mcp3		,"mcp3_wi/D");
	
	outTree->Branch("t0_t03"	,&tdiff_t0_t03		,"t0_t03/D");
	outTree->Branch("t03_wi"	,&width_t03		,"t03_wi/D");

	outTree->Branch("t0_t0Ana"	,&tdiff_t0_t0Ana	,"t0_t0Ana/D");
	outTree->Branch("t0Ana_wi"	,&width_t0Ana		,"t0Ana_wi/D");

	outTree->Branch("t0_t0Tqdc"	,&tdiff_t0_t0Tqdc	,"t0_t0Tqdc/D");
	outTree->Branch("t0Tqdc_wi"	,&width_t0Tqdc		,"t0Tqdc_wi/D");


	// Before hand, let's calculate the peak of the T0 width
	TH1D * h = new TH1D("trash","trash",1000,0,100);
	TString w;
	//w.Form(" fabs (sqrt( (TQDC_BC1[0]->GetPeak() - %f) * (TQDC_T0[0]->GetPeak() - %f ) ) - %f ) < 2 * %f ",pedBC1,pedBC2,C12in_mean,C12in_sig);
	TString z;
	z.Form(" fabs (sqrt( (TQDC_BC3[0]->GetPeak() - %f) * (TQDC_BC4[0]->GetPeak() - %f ) ) - %f ) < 2 * %f ",pedBC3,pedBC4,C12out_mean,C12out_sig);
	//inTree->Draw("T0[0]->GetAmp() >> trash", w + "&&" + z);
	inTree->Draw("T0[0]->GetAmp() >> trash", z);
	TFitResultPtr init = h->Fit("gaus","QES","",0,100);
	double start = init->Parameter(1) - init->Parameter(2);
	TFitResultPtr fin  = h->Fit("gaus","QES","",start , 100.);
	double t0Wid_peak = fin->Parameter(1);
	double t0Wid_sig  = fin->Parameter(2);

	const int nEvents = inTree->GetEntries();
	double cnt = 0.;
	double cnt_t0andMCP2 = 0.;
	double cnt_t0andMCP3 = 0.;
	double cnt_t0andt03 = 0.;
	double cnt_t0andt0Ana = 0.;
	double cnt_t0andTqdc = 0.;
	for( int event = 0 ; event < nEvents ; event++){
		tdiff_t0_mcp2 	= -666;
		width_mcp2	= -666;
		tdiff_t0_mcp3 	= -666;
		width_mcp3	= -666;
		tdiff_t0_t03	= -666;
		width_t03	= -666;
		tdiff_t0_t0Ana	= -666;
		width_t0Ana	= -666;
		tdiff_t0_t0Tqdc	= -666;
		width_t0Tqdc	= -666;	
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
				
			//if( fabs( sqrt( adcBC1 * adcBC2 ) - C12in_mean) > 2*C12in_sig ) continue;
			if( fabs( sqrt( adcBC3 * adcBC4 ) - C12out_mean) > 2*C12out_sig ) continue;
		}
	
		// First let's cut on T0 width to take out any T0 time-walk dependence		
		if( fabs( t0Width - t0Wid_peak) > 0.1) continue;
		cnt++;

		// Now with a good quality-control cut, we can look at time resolution
		// 	Look at MCP
		if( mcp2Data->GetEntriesFast() ){
			int mcp2Idx;
			findIdx_tdc( mcp2Data, mcp2Idx, t0Time);
			BmnTrigDigit * signal = (BmnTrigDigit *) mcp2Data->At(mcp2Idx);
			tdiff_t0_mcp2 = t0Time - signal->GetTime();
			width_mcp2 = signal->GetAmp();
			cnt_t0andMCP2++;
		}
		if( mcp3Data->GetEntriesFast() ){
			int mcp3Idx;
			findIdx_tdc( mcp3Data, mcp3Idx, t0Time);
			BmnTrigDigit * signal = (BmnTrigDigit *) mcp3Data->At(mcp3Idx);
			tdiff_t0_mcp3 = t0Time - signal->GetTime();
			width_mcp3 = signal->GetAmp();
			cnt_t0andMCP3++;
		}
			// Look at BMN-T0
		if( t03Data->GetEntriesFast() ){
			int t03Idx;
			findIdx_tdc( t03Data, t03Idx, t0Time);
			BmnTrigDigit * signal = (BmnTrigDigit *) t03Data->At(t03Idx);
			tdiff_t0_t03 = t0Time - signal->GetTime();
			width_t03 = signal->GetAmp();
			cnt_t0andt03++;
		}
			// Look at T03Ana
		if( t0AnaData->GetEntriesFast() ){
			int t0AnaIdx;
			findIdx_tdc( t0AnaData, t0AnaIdx, t0Time);
			BmnTrigDigit * signal = (BmnTrigDigit *) t0AnaData->At(t0AnaIdx);
			tdiff_t0_t0Ana = t0Time - signal->GetTime();
			width_t0Ana = signal->GetAmp();
			cnt_t0andt0Ana++;
		}
			// Look at T0 in TQDC
		if( bc2Data->GetEntriesFast() ){
			int bc2Idx;
			findIdx( bc2Data, bc2Idx, t0Time);
			BmnTrigWaveDigit * signal = (BmnTrigWaveDigit *) bc2Data->At(bc2Idx);
			tdiff_t0_t0Tqdc = t0Time - signal->GetTime();
			width_t0Tqdc = signal->GetPeak();
			cnt_t0andTqdc++;
		}
	


		outTree->Fill();	
	}

	cout 	<< cnt << " " 
		<< cnt_t0andMCP2 << " " << cnt_t0andMCP3 << " " << cnt_t0andt03 << " " 
		<< cnt_t0andt0Ana<< " " << cnt_t0andTqdc << "\n";


	outFile->cd();
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
